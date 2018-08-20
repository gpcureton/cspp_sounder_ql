#!/usr/bin/env python
# encoding: utf-8
"""
sounder_image.py

Purpose: Create a plot of temperature, dewpoint or relative humidity,
         at a particular pressure level. Supports outputs from the following
         packages...

         * International ATOVS Processing Package (IAPP)
         * Microwave Integrated Retrieval System (MIRS)
         * CSPP Hyperspectral Retrieval (HSRTV) Package
         * NOAA Unique CrIS/ATMS Product System (NUCAPS)


Preconditions:
    * matplotlib (with basemap)
    * netCDF4 python module
    * h5py python module

Optional:
    *

Minimum commandline:

    python sounder_image.py  INPUTFILE DATATYPE

where...

    INPUTFILES: The fully qualified path to the input files. May be a
    file glob.

    DATATYPE: One of 'IAPP','MIRS', 'HSRTV' or 'NUCAPS'.


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2014-05-10.

Copyright 2014-2018, University of Wisconsin Regents.
Licensed under the GNU GPLv3.
"""

import os, sys, logging, traceback
from os import path,uname,environ
import string
import re
import uuid
from shutil import rmtree,copyfile
from glob import glob
from time import time
from datetime import datetime,timedelta
from collections import OrderedDict

import numpy as np
from numpy import ma
import copy

from scipy import vectorize

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

from mpl_toolkits.basemap import Basemap

import scipy.spatial as ss

import xarray
from netCDF4 import Dataset
from netCDF4 import num2date
import h5py

from thermo import dewhum

# every module should have a LOG object
LOG = logging.getLogger(__file__)

def uwyoming():
    sounding = SkewT.Sounding(filename="sounding.txt")

    sounding_keys = ['mixr',
                     'temp',
                     'thtv',
                     'SoundingDate',
                     'hght',
                     'StationNumber',
                     'dwpt',
                     'thte',
                     'sknt',
                     'relh',
                     'drct',
                     'thta',
                     'pres']

    sounding_inputs = {}

    for key in sounding_keys:
        sounding_inputs[key] = sounding[key]

    return sounding_inputs


def granuleFiles(prodGlob):
    '''
    Returns sorted lists of the geolocation and product files.
    '''

    prodDir = path.dirname(path.abspath(path.expanduser(prodGlob)))

    LOG.debug("Initial prodGlob = {}".format(prodGlob))

    prodGlob = path.basename(path.abspath(path.expanduser(prodGlob)))

    LOG.debug("prodDir = {}".format(prodDir))
    LOG.debug("prodGlob = {}".format(prodGlob))

    prodList = glob("%s/%s" % (prodDir,prodGlob))
    prodList.sort()

    return prodList

def generate_file_list(inputs, full=False):
    '''
    Trawl through the files and directories given at the command line, return a list of files.
    '''

    input_files = []

    for input in inputs:
        LOG.debug("bash glob input = {}".format(input))

    if full:
        input_files = list(set(inputs))
        input_files.sort()
        return input_files

    input_dirs = []
    input_files = []

    # Sort the command line inputs into directory and file inputs...
    for input in inputs:
        input = os.path.abspath(os.path.expanduser(input))
        if os.path.isdir(input):
            # Input file glob is of form "/path/to/files"
            LOG.debug("Input {} is a directory containing files...".format(input))
            input_dirs.append(input)
        elif os.path.isfile(input):
            # Input file glob is of form "/path/to/files/goes13_1_2015_143_1745.input"
            LOG.debug("Input {} is a file.".format(input))
            input_files.append(input)

    input_dirs = list(set(input_dirs))
    input_dirs.sort()
    input_files = list(set(input_files))
    input_files.sort()

    for dirs in input_dirs:
        LOG.debug("input dirs {}".format(dirs))
    for files in input_files:
        LOG.debug("input files {}".format(files))

    return input_files

def get_pressure_index(pressure,pres_0=850.,kd_dist=10.):
    '''
    Get the indices in the pressure level array closes to the desired
    pressure. If no pressure is specified, use the closest to 850mb.
    Uses kd-tree.
    '''

    try:
        pres_shape = pressure.shape

        nlevels = pres_shape[0]
        LOG.debug("nlevels= {}".format(nlevels))

        pres_mask = pressure.mask if ma.is_masked(pressure) \
                else np.zeros(pressure.shape,dtype='bool')


        pressure_masked = ma.array(pressure,mask=pres_mask)
        pressure_masked = pressure_masked.compressed()

        levels = np.arange(nlevels)
        levels_masked = ma.array(levels,mask=pres_mask)
        levels_masked = levels_masked.compressed()

        LOG.debug("There are {} masked pressure levels".\
                format(np.sum(pres_mask)))

        LOG.info("Searching for pressure {:4.2f}".format(pres_0))

        # Construct the kdtree
        points = zip(pressure_masked)
        tree = ss.KDTree(points)

        # get index of each point which distance from (lon_0,lat_0) is under 1
        a = tree.query_ball_point([pres_0], kd_dist)
        if a == []:
            LOG.warn("Specified pressure {:4.2f} not found...".\
                    format(pres_0))
            return None

        level = [levels_masked[i] for i in a][0]

        LOG.debug("Level index is {}".format(level))

        LOG.info("Retrieved pressure is {:4.2f} mb".format(pressure[level].squeeze()))

    except Exception, err :
        LOG.warn("{}".format(err))
        LOG.debug(traceback.format_exc())

    return level


def iapp_sounder(iapp_file_list,pres_0=850.):
    '''
    Pressure_Levels: Pressure for each level in mb (hPa)
    Temperature_Retrieval: Temperature profile in K
    WaterVapor_Retrieval: Water vapor profile (mixing ratio) in g/kg
    Dew_Point_Temp_Retrieval: Dew Point Temperature Profile in K
    '''

    '''
    float Latitude(Along_Track, Across_Track) ;
            Latitude:long_name = "Geodetic Latitude" ;
            Latitude:units = "degrees_north" ;
            Latitude:scale_factor = 1. ;
            Latitude:add_offset = 0. ;
            Latitude:Parameter_Type = "IAPP Output" ;
            Latitude:valid_range = -90.f, 90.f ;
            Latitude:_FillValue = -999.f ;

    float Longitude(Along_Track, Across_Track) ;
            Longitude:long_name = "Geodetic Longitude" ;
            Longitude:units = "degrees_east" ;

    float Pressure_Levels(Pres_Levels) ;
            Pressure_Levels:long_name = "Pressure Levels at which the retrieved
                                         values are calculated" ;
            Pressure_Levels:units = "hPa" ;

    float Temperature_Retrieval(Along_Track, Across_Track, Pres_Levels) ;
            Temperature_Retrieval:long_name = "Temperature Retrieval for
                                               the IAPP" ;
            Temperature_Retrieval:units = "degrees Kelvin" ;

    float WaterVapor_Retrieval(Along_Track, Across_Track, Pres_Levels) ;
            WaterVapor_Retrieval:long_name = "Water Vapor Retrieval for
                                             the IAPP" ;
            WaterVapor_Retrieval:units = "g/kg" ;

    float Dew_Point_Temp_Retrieval(Along_Track, Across_Track, Pres_Levels) ;
            Dew_Point_Temp_Retrieval:long_name = "Dew Point Temperature Retrieval
                                                  for the IAPP" ;
            Dew_Point_Temp_Retrieval:units = "degrees Kelvin" ;

    float Cloud_Top_Temperature_CO2(Along_Track, Across_Track, Num_of_FOVs) ;
        Cloud_Top_Temperature_CO2:long_name = "CO2 Slicing Cloud Top Temperature of the 3x3 box" ;
        Cloud_Top_Temperature_CO2:units = "K" ;

    float Cloud_Top_Pressure_CO2(Along_Track, Across_Track, Num_of_FOVs) ;
            Cloud_Top_Pressure_CO2:long_name = "CO2 Slicing Cloud Top Pressure of the 3x3 box" ;
            Cloud_Top_Pressure_CO2:units = "hPa" ;
    '''

    data_labels = [
                    'pres',
                    'temp',
                    'dwpt',
                    'wvap',
                    'relh',
                    'ctp',
                    'ctt',
                    'lat',
                    'lon'
                    ]

    sounding_inputs = {}
    for label in data_labels:
        sounding_inputs[label] = {}

    sounding_inputs['pres']['dset_name'] = 'Pressure_Levels'
    sounding_inputs['temp']['dset_name'] = 'Temperature_Retrieval'
    sounding_inputs['dwpt']['dset_name'] = 'Dew_Point_Temp_Retrieval'
    sounding_inputs['wvap']['dset_name'] = 'WaterVapor_Retrieval'
    sounding_inputs['ctp']['dset_name']  = 'Cloud_Top_Pressure_CO2'
    sounding_inputs['ctt']['dset_name']  = 'Cloud_Top_Temperature_CO2'
    sounding_inputs['lat']['dset_name']  = 'Latitude'
    sounding_inputs['lon']['dset_name']  = 'Longitude'
    sounding_inputs['relh'] = None

    LOG.debug("Input file list: {}".format(iapp_file_list))

    first_granule = True

    for grans in np.arange(len(iapp_file_list)):

        try:

            LOG.debug("Ingesting granule {} ...".format(grans))

            iapp_file = iapp_file_list[grans]
            LOG.debug("Ingesting granule {} ...".format(iapp_file))

            # Open the file object
            fileObj = Dataset(iapp_file)

            # Open the dataset nodes for each of the required datasets
            sounding_inputs['pres']['node'] = fileObj.variables['Pressure_Levels']
            sounding_inputs['temp']['node'] = fileObj.variables['Temperature_Retrieval']
            sounding_inputs['dwpt']['node'] = fileObj.variables['Dew_Point_Temp_Retrieval']
            sounding_inputs['wvap']['node'] = fileObj.variables['WaterVapor_Retrieval']
            sounding_inputs['ctp']['node'] = fileObj.variables['Cloud_Top_Pressure_CO2']
            sounding_inputs['ctt']['node'] = fileObj.variables['Cloud_Top_Temperature_CO2']
            sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
            sounding_inputs['lon']['node'] = fileObj.variables['Longitude']

            if first_granule:

                # Get the pressure scale
                pressure = sounding_inputs['pres']['node'][:]

                for label in data_labels:
                    if sounding_inputs[label] is not None:
                        LOG.info(">>> Processing {} ..."
                                .format(sounding_inputs[label]['dset_name']))
                        for attr_name in sounding_inputs[label]['node'].ncattrs():
                            attr = getattr(sounding_inputs[label]['node'],attr_name)
                            LOG.debug("{} = {}".format(attr_name,attr))
                            sounding_inputs[label][attr_name] = attr

                # Search for the pressure level closest to the input value
                pressure_scope = 10.
                LOG.debug("Scope of pressure level search is {} hPa"
                        .format(pressure_scope))
                level = get_pressure_index(pressure,pres_0=pres_0,
                        kd_dist=pressure_scope)
                attempts = 1
                while (level is None and attempts<10):
                    pressure_scope *= 1.1
                    LOG.debug("Scope of pressure level search is {} hPa"
                            .format(pressure_scope))
                    level = get_pressure_index(pressure,pres_0=pres_0,
                            kd_dist=pressure_scope)
                    attempts += 1

                if level is None:
                    raise Exception("No suitable pressure level found, aborting...")

                LOG.debug("Retrieved level = {}".format(level))
                sounding_inputs['pres_0'] = pressure[level]
                data_labels.remove('pres')

            # Add to the lat and lon and data arrays...
            for label in ['lat','lon']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:]
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            # Add to the ctp and   ctt and data arrays...
            cloud_top_qc = {'ctp':0., 'ctt':0.}
            for label in ['ctp','ctt']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:,5]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:,5]
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            # Add to the temp, dewpoint, water vapour and relative humidity data arrays...
            for label in ['temp','dwpt','wvap']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:,level]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:,level]
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            first_granule = False

            LOG.info("Closing {}".format(iapp_file))
            fileObj.close()

        except Exception, err :
            LOG.warn("There was a problem, closing {}".format(iapp_file))
            LOG.warn("{}".format(err))
            LOG.debug(traceback.format_exc())
            fileObj.close()
            LOG.info("Exiting...")
            sys.exit(1)

    # Contruct the pressure array
    pressure_array = pressure[level] * \
            np.ones(sounding_inputs['lat']['data'].shape,dtype='float')
    LOG.debug("pressure_array.shape =  {}".format(pressure_array.shape))

    # Construct the masks of the various datasets
    for label in ['temp','dwpt','wvap','ctp','ctt','lat','lon']:
        if sounding_inputs[label] is not None:
            LOG.debug(">>> Processing {} ..."
                    .format(sounding_inputs[label]['dset_name']))
            fill_value = sounding_inputs[label]['_FillValue']
            data = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            #data = ma.masked_equal(sounding_inputs[label]['data'],np.nan)
            LOG.debug("ma.is_masked({}) = {}".format(label,ma.is_masked(data)))

            if ma.is_masked(data):
                if data.mask.shape == ():
                    data_mask = np.ones(data.shape,dtype='bool')
                else:
                    data_mask = data.mask
            else:
                data_mask = np.zeros(data.shape,dtype='bool')

            LOG.debug("There are {}/{} masked values in {}".\
                    format(np.sum(data_mask),data.size,label))

            sounding_inputs[label]['mask'] = data_mask

    # Some extra sanity checking on the CTT and CTP dataets
    cloud_top_qc = {'ctp':0., 'ctt':150.}
    for label in ['ctp','ctt']:
        sounding_inputs[label]['mask'] += ma.masked_less(sounding_inputs[label]['data'], cloud_top_qc[label]).mask
    sounding_inputs['ctp']['mask'] += sounding_inputs['ctt']['mask']

    # Computing the relative humidity
    sounding_inputs['relh'] = {}
    LOG.debug("wvap.shape =  {}".format(sounding_inputs['wvap']['data'].shape))
    LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    relh_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']

    wvap = ma.array(sounding_inputs['wvap']['data'],mask=relh_mask)
    temp = ma.array(sounding_inputs['temp']['data'],mask=relh_mask)
    pressure_array = ma.array(pressure_array,mask=relh_mask)

    dewhum_vec = np.vectorize(dewhum)
    _,rh,_ = dewhum_vec(pressure_array,temp,wvap)

    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            else np.zeros(temp.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    # Implement dummy QC mask
    sounding_inputs['qc'] = {'mask':np.zeros(data_mask.shape,dtype='bool')}

    return sounding_inputs


def mirs_sounder(mirs_file_list,pres_0=850.):
    '''
    Player: Pressure for each layer in mb (hPa)
    PTemp: Temperature profile in K
    PVapor: Water vapor profile (mixing ratio) in g/kg
    '''

    '''
    float Latitude(Scanline, Field_of_view) ;
            Latitude:long_name = "Latitude of the view (-90,90)" ;

    float Longitude(Scanline,Field_of_view) ;
            Longitude:long_name = "Longitude of the view (-180,180)" ;

    float Player(P_Layer) ;
            Player:description = "Pressure for each layer in mb" ;

    float PTemp(Scanline, Field_of_view, P_Layer) ;
            PTemp:long_name = "Temperature profile in K" ;
            PTemp:scale = 1. ;

    float PVapor(Scanline, Field_of_view, P_Layer) ;
            PVapor:long_name = "Water vapor profile in g/kg" ;
            PVapor:scale = 1. ;

    short Qc(Scanline, Field_of_view, Qc_dim) ;
            Qc:description = "Qc(0): 0-good,1-usable with problem ,2-bad" ;
    '''

    data_labels = [
                    'pres',
                    'temp',
                    'dwpt',
                    'wvap',
                    'relh',
                    'lat',
                    'lon',
                    'qc'
                    ]

    sounding_inputs = {}
    for label in data_labels:
        sounding_inputs[label] = {}

    sounding_inputs['pres']['dset_name'] = 'Player'
    sounding_inputs['temp']['dset_name'] = 'PTemp'
    sounding_inputs['wvap']['dset_name'] = 'PVapor'
    sounding_inputs['lat']['dset_name']  = 'Latitude'
    sounding_inputs['lon']['dset_name']  = 'Longitude'
    sounding_inputs['qc']['dset_name']  = 'Qc'
    sounding_inputs['dwpt'] = None
    sounding_inputs['relh'] = None

    LOG.debug("Input file list: {}".format(mirs_file_list))

    first_granule = True

    for grans in np.arange(len(mirs_file_list)):

        try:

            LOG.debug("Ingesting granule {} ...".format(grans))

            mirs_file = mirs_file_list[grans]
            LOG.debug("Ingesting granule {} ...".format(mirs_file))

            # Open the file object
            fileObj = Dataset(mirs_file)

            sounding_inputs['pres']['node'] = fileObj.variables['Player']
            sounding_inputs['temp']['node'] = fileObj.variables['PTemp']
            sounding_inputs['wvap']['node'] = fileObj.variables['PVapor']
            sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
            sounding_inputs['lon']['node'] = fileObj.variables['Longitude']
            sounding_inputs['qc']['node'] = fileObj.variables['Qc']

            # Fill in some missing attributes
            fill_value = getattr(fileObj,'missing_value')
            sounding_inputs['pres']['units'] = 'hPa'
            sounding_inputs['temp']['units'] = 'K'
            sounding_inputs['wvap']['units'] = 'g/kg'

            if first_granule:

                # Get the pressure scale
                pressure = sounding_inputs['pres']['node'][:]

                for label in data_labels:
                    if sounding_inputs[label] is not None:
                        LOG.info(">>> Processing {} ..."
                                .format(sounding_inputs[label]['dset_name']))
                        sounding_inputs[label]['_FillValue'] = fill_value
                        for attr_name in sounding_inputs[label]['node'].ncattrs():
                            attr = getattr(sounding_inputs[label]['node'],attr_name)
                            LOG.debug("{} = {}".format(attr_name,attr))
                            sounding_inputs[label][attr_name] = attr

                # Search for the pressure level closest to the input value
                pressure_scope = 10.
                LOG.debug("Scope of pressure level search is {} hPa"
                        .format(pressure_scope))
                level = get_pressure_index(pressure,pres_0=pres_0,
                        kd_dist=pressure_scope)
                attempts = 1
                while (level is None and attempts<10):
                    pressure_scope *= 1.1
                    LOG.debug("Scope of pressure level search is {} hPa"
                            .format(pressure_scope))
                    level = get_pressure_index(pressure,pres_0=pres_0,
                            kd_dist=pressure_scope)
                    attempts += 1

                if level is None:
                    raise Exception("No suitable pressure level found, aborting...")

                LOG.debug("Retrieved level = {}".format(level))
                sounding_inputs['pres_0'] = pressure[level]
                data_labels.remove('pres')

            # Add to the lat and lon and data arrays...
            for label in ['lat','lon']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:]
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            # Add to the temp, dewpoint, water vapour and relative humidity data arrays...
            for label in ['temp','wvap']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:,level]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:,level]
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            for label in ['qc']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:,0]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:,0]
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            first_granule = False

            LOG.info("Closing {}".format(mirs_file))
            fileObj.close()

        except Exception, err :
            LOG.warn("There was a problem, closing {}".format(mirs_file))
            LOG.warn("{}".format(err))
            LOG.debug(traceback.format_exc())
            fileObj.close()
            LOG.info("Exiting...")
            sys.exit(1)

    # Contruct the pressure array
    pressure_array = pressure[level] * \
            np.ones(sounding_inputs['lat']['data'].shape,dtype='float')
    LOG.debug("pressure_array.shape =  {}".format(pressure_array.shape))

    # Construct the masks of the various datasets
    for label in ['temp','wvap','lat','lon']:
        if sounding_inputs[label] is not None:
            LOG.debug(">>> Processing {} ..."
                    .format(sounding_inputs[label]['dset_name']))
            data = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            LOG.debug("ma.is_masked({}) = {}".format(label,ma.is_masked(data)))

            if ma.is_masked(data):
                if data.mask.shape == ():
                    data_mask = np.ones(data.shape,dtype='bool')
                else:
                    data_mask = data.mask
            else:
                data_mask = np.zeros(data.shape,dtype='bool')

            LOG.debug("There are {}/{} masked values in {}".\
                    format(np.sum(data_mask),data.size,label))

            sounding_inputs[label]['mask'] = data_mask

    # Contruct the QC mask
    for label in ['qc']:
        if sounding_inputs[label] is not None:
            LOG.debug(">>> Processing {} ..."
                    .format(sounding_inputs[label]['dset_name']))
            fill_value = sounding_inputs[label]['_FillValue']
            data_fill = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            data_qc = ma.masked_greater(sounding_inputs[label]['data'],0)
            LOG.debug("{} has fill: {}".format(label,ma.is_masked(data_fill)))
            LOG.debug("{} has QC flag: {}".format(label,ma.is_masked(data_qc)))

            if ma.is_masked(data_fill):
                if data_fill.mask.shape == ():
                    data_fill_mask = np.ones(data_fill.shape,dtype='bool')
                else:
                    data_fill_mask = data_fill.mask
            else:
                data_fill_mask = np.zeros(data_fill.shape,dtype='bool')

            LOG.debug("There are {}/{} fill values in {}".\
                    format(np.sum(data_fill_mask),data_fill.size,label))

            if ma.is_masked(data_qc):
                if data_qc.mask.shape == ():
                    data_qc_mask = np.ones(data_qc.shape,dtype='bool')
                else:
                    data_qc_mask = data_qc.mask
            else:
                data_qc_mask = np.zeros(data_qc.shape,dtype='bool')

            LOG.debug("There are {}/{} qc values in {}".\
                    format(np.sum(data_qc_mask),data_qc.size,label))

            sounding_inputs[label]['mask'] = data_fill_mask + data_qc_mask

            LOG.debug("There are {}/{} masked values in total {}".format(
                np.sum(sounding_inputs[label]['mask']),sounding_inputs[label]['data'].size, label))

    # Computing the dewpoint temperature and relative humidity
    sounding_inputs['dwpt'] = {}
    sounding_inputs['relh'] = {}
    LOG.debug("wvap.shape =  {}".format(sounding_inputs['wvap']['data'].shape))
    LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    dwpt_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']
    relh_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']

    wvap = ma.array(sounding_inputs['wvap']['data'],mask=dwpt_mask)
    temp = ma.array(sounding_inputs['temp']['data'],mask=relh_mask)
    pressure_array = ma.array(pressure_array,mask=relh_mask)

    dewhum_vec = np.vectorize(dewhum)
    dwpt,rh,_ = dewhum_vec(pressure_array,temp,wvap)

    # Copy dewpoint temperature to output
    sounding_inputs['dwpt']['data'] = dwpt
    sounding_inputs['dwpt']['mask'] = dwpt.mask if  ma.is_masked(dwpt) \
            else np.zeros(dwpt.shape,dtype='bool')
    sounding_inputs['dwpt']['units'] = 'K'
    sounding_inputs['dwpt']['long_name'] = 'Dew-point Temperature in K'

    LOG.debug("ma.is_masked({}) = {}".format('dwpt',ma.is_masked(dwpt)))
    LOG.debug("There are {} masked values in dwpt".format(np.sum(dwpt.mask)))

    # Copy relative humidity to output
    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            else np.zeros(rh.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    return sounding_inputs


def hsrtv_sounder(dr_file_list,pres_0=850.):
    '''
    Plevs: Pressure for each level in mb (hPa)
    TAir: Temperature profile in K
    Dewpnt: Water vapor profile (mixing ratio) in g/kg

    /Latitude
        missing_value : -9999.
        units : degrees_north
        valid_range : -90.0,90.0
    /Longitude
        missing_value : -9999.
        units : degrees_east
        valid_range : -180.0,180.0
    /Plevs
        description: vertical coordinate axis (101 pressure levels)
        units : hPa
    /TAir
        description : retrieved temperature profile at 101 levels
        missing_value : -9999.
        units : K
    /Dewpnt
        description : dew point temperature profile
        missing_value : -9999.
        units : K
    /RelHum
        description : retrieved relative humidity at 101 levels
        missing_value : -9999.
        units : %
    /H2OMMR
        description : retrieved humidity (water vapor mixing ratio) profile at 101 levels
        missing_value : -9999.
        units : g/kg
    /CTP
        description : Cloud top pressure
        missing_value : -9999.
        units : hPa
    /CTT
        description : Cloud top temperature
        missing_value : -9999.
        units : K
    '''

    data_labels = [
                    'pres',
                    'temp',
                    'dwpt',
                    'wvap',
                    'relh',
                    'ctp',
                    'ctt',
                    'lat',
                    'lon'
                    ]

    sounding_inputs = {}
    for label in data_labels:
        sounding_inputs[label] = {}

    sounding_inputs['pres']['dset_name'] = 'Plevs'
    sounding_inputs['temp']['dset_name'] = 'TAir'
    sounding_inputs['dwpt']['dset_name'] = 'Dewpnt'
    sounding_inputs['relh']['dset_name'] = 'RelHum'
    sounding_inputs['wvap']['dset_name'] = 'H2OMMR'
    sounding_inputs['ctp']['dset_name']  = 'CTP'
    sounding_inputs['ctt']['dset_name']  = 'CTT'
    sounding_inputs['lat']['dset_name']  = 'Latitude'
    sounding_inputs['lon']['dset_name']  = 'Longitude'


    LOG.debug("Input file list: {}".format(dr_file_list))

    first_granule = True

    for grans in np.arange(len(dr_file_list)):

        try:

            LOG.debug("Ingesting granule {} ...".format(grans))

            dr_file = dr_file_list[grans]
            LOG.debug("Ingesting granule {} ...".format(dr_file))

            # Open the file object
            if h5py.is_hdf5(dr_file_list[grans]):
                fileObj = h5py.File(dr_file,'r')
            else:
                LOG.error("{} does not exist or is not a valid HDF5 file, aborting...".\
                        format(dr_file))
                sys.exit(1)

            # Open the dataset nodes for each of the required datasets
            sounding_inputs['pres']['node'] = fileObj['Plevs']
            sounding_inputs['temp']['node'] = fileObj['TAir']
            sounding_inputs['dwpt']['node'] = fileObj['Dewpnt']
            sounding_inputs['relh']['node'] = fileObj['RelHum']
            sounding_inputs['wvap']['node'] = fileObj['H2OMMR']
            sounding_inputs['ctp']['node'] = fileObj['CTP']
            sounding_inputs['ctt']['node'] = fileObj['CTT']
            sounding_inputs['lat']['node'] = fileObj['Latitude']
            sounding_inputs['lon']['node'] = fileObj['Longitude']

            if first_granule:

                # Get the pressure scale
                pressure = sounding_inputs['pres']['node'][:]

                #FIXME: All attributes for the HSRTV datasets appear to be strings,
                #       whether the attribute value is a string or not. So we can't
                #       determine the attribute type from the file, we just have
                #       *know* it. Why bother with self-descibing file formats then?
                for label in data_labels:
                    if sounding_inputs[label] is not None:
                        LOG.info(">>> Processing {} ...".\
                                format(sounding_inputs[label]['dset_name']))
                        for attr_name in sounding_inputs[label]['node'].attrs.keys():
                            attr = sounding_inputs[label]['node'].attrs[attr_name]
                            LOG.debug("{} = {}".format(attr_name,np.squeeze(attr)))
                            sounding_inputs[label][attr_name] = attr

                # Search for the pressure level closest to the input value
                pressure_scope = 10.
                LOG.debug("Scope of pressure level search is {} hPa"
                        .format(pressure_scope))
                level = get_pressure_index(pressure,pres_0=pres_0,
                        kd_dist=pressure_scope)
                attempts = 1
                while (level is None and attempts<10):
                    pressure_scope *= 1.1
                    LOG.debug("Scope of pressure level search is {} hPa"
                            .format(pressure_scope))
                    level = get_pressure_index(pressure,pres_0=pres_0,
                            kd_dist=pressure_scope)
                    attempts += 1

                if level is None:
                    raise Exception("No suitable pressure level found, aborting...")

                LOG.debug("Retrieved level = {}".format(level))
                sounding_inputs['pres_0'] = pressure[level]
                data_labels.remove('pres')

            # Add to the lat and lon and data arrays...
            for label in ['lat','lon','ctp','ctt']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:]
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            # Add to the temp, dewpoint, water vapour and relative humidity data arrays...
            for label in ['temp','dwpt','wvap','relh']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][level,:,:]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][level,:,:]
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            first_granule = False

            LOG.info("Closing {}".format(dr_file))
            fileObj.close()

        except Exception, err :
            LOG.warn("There was a problem, closing {}".format(dr_file))
            LOG.warn("{}".format(err))
            LOG.debug(traceback.format_exc())
            fileObj.close()
            LOG.info("Exiting...")
            sys.exit(1)

    # Contruct the pressure array
    pressure_array = pressure[level] * \
            np.ones(sounding_inputs['lat']['data'].shape,dtype='float')
    LOG.debug("pressure_array.shape =  {}".format(pressure_array.shape))

    # Construct the masks of the various datasets
    for label in ['temp','dwpt','wvap','relh','ctp','ctt','lat','lon']:
        if sounding_inputs[label] is not None:
            LOG.debug(">>> Processing {} ...".format(sounding_inputs[label]['dset_name']))
            fill_value = float(sounding_inputs[label]['missing_value'][0])
            data = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            LOG.debug("ma.is_masked({}) = {}".format(label,ma.is_masked(data)))

            if ma.is_masked(data):
                if data.mask.shape == ():
                    data_mask = np.ones(data.shape,dtype='bool')
                else:
                    data_mask = data.mask
            else:
                data_mask = np.zeros(data.shape,dtype='bool')

            LOG.debug("There are {}/{} masked values in {}".\
                    format(np.sum(data_mask),data.size,label))

            sounding_inputs[label]['mask'] = data_mask

    calc_relh = False
    #calc_relh = True

    if calc_relh:
        # Computing the relative humidity
        sounding_inputs['relh'] = {}
        LOG.debug("wvap.shape =  {}".format(sounding_inputs['wvap']['data'].shape))
        LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

        relh_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']

        wvap = ma.array(sounding_inputs['wvap']['data'],mask=relh_mask)
        temp = ma.array(sounding_inputs['temp']['data'],mask=relh_mask)
        pressure_array = ma.array(pressure_array,mask=relh_mask)

        dewhum_vec = np.vectorize(dewhum)
        _,rh,_ = dewhum_vec(pressure_array,temp,wvap)

        # Copy relative humidity to output
        sounding_inputs['relh']['data'] = rh
        sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
                else np.zeros(temp.shape,dtype='bool')
        sounding_inputs['relh']['units'] = '%'
        sounding_inputs['relh']['long_name'] = 'Relative Humidity'

        LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
        LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    # Implement dummy QC mask
    sounding_inputs['qc'] = {'mask':np.zeros(data_mask.shape,dtype='bool')}

    return sounding_inputs


def nucaps_sounder(nucaps_file_list,pres_0=850.):
    '''
    Pressure: Pressure for each layer in mb (hPa)
    Temperature: Temperature profile in K
    H2O_MR: Water vapor profile (mixing ratio) in g/g
    '''

    '''
    float Latitude(Number_of_CrIS_FORs) ;
            Latitude:long_name = "Retrieval latitude values for each CrIS FOR" ;
            Latitude:standard_name = "latitude" ;
            Latitude:units = "degrees_north" ;
            Latitude:parameter_type = "NUCAPS data" ;
            Latitude:valid_range = -90.f, 90.f ;
            Latitude:_FillValue = -9999.f ;

    float Longitude(Number_of_CrIS_FORs) ;
            Longitude:long_name = "Retrieval longitude values for each CrIS FOR" ;
            Longitude:standard_name = "longitude" ;
            Longitude:units = "degrees_east" ;
            Longitude:parameter_type = "NUCAPS data" ;
            Longitude:valid_range = -180.f, 180.f ;
            Longitude:_FillValue = -9999.f ;

    float Pressure(Number_of_CrIS_FORs, Number_of_P_Levels) ;
            Pressure:long_name = "Pressure" ;
            Pressure:units = "mb" ;
            Pressure:parameter_type = "NUCAPS data" ;
            Pressure:coordinates = "Time Latitude Longitude" ;
            Pressure:valid_range = 0.f, 2000.f ;
            Pressure:_FillValue = -9999.f ;

    float Temperature(Number_of_CrIS_FORs, Number_of_P_Levels) ;
            Temperature:long_name = "Temperature" ;
            Temperature:standard_name = "air_temperature" ;
            Temperature:units = "Kelvin" ;
            Temperature:parameter_type = "NUCAPS data" ;
            Temperature:coordinates = "Time Latitude Longitude Pressure" ;
            Temperature:valid_range = 0.f, 1000.f ;
            Temperature:_FillValue = -9999.f ;

    float H2O_MR(Number_of_CrIS_FORs, Number_of_P_Levels) ;
            H2O_MR:long_name = "water vapor mixing ratio" ;
            H2O_MR:units = "g/g" ;
            H2O_MR:parameter_type = "NUCAPS data" ;
            H2O_MR:coordinates = "Time Latitude Longitude Effective_Pressure" ;
            H2O_MR:valid_range = 0.f, 1.e+08f ;
            H2O_MR:_FillValue = -9999.f ;

    float Cloud_Top_Pressure(Number_of_CrIS_FORs, Number_of_Cloud_Layers) ;
            Cloud_Top_Pressure:long_name = "Cloud top pressure" ;
            Cloud_Top_Pressure:units = "mb" ;
            Cloud_Top_Pressure:parameter_type = "NUCAPS data" ;
            Cloud_Top_Pressure:coordinates = "Time Latitude Longitude" ;
            Cloud_Top_Pressure:valid_range = 0.f, 10000.f ;
            Cloud_Top_Pressure:_FillValue = -9999.f ;

    int Quality_Flag(Number_of_CrIS_FORs) ;
            Quality_Flag:long_name = "Quality flags for retrieval" ;
            Quality_Flag:parameter_type = "NUCAPS data" ;
            Quality_Flag:coordinates = "Time Latitude Longitude" ;
            Quality_Flag:valid_range = 0, 31 ;
            Quality_Flag:_FillValue = -9999 ;
            Quality_Flag:flag_values = 0, 1, 2, 4, 8, 9, 16, 17, 24, 25 ;
            Quality_Flag:flag_meanings = "accepted reject_physical reject_MIT reject_NOAA_reg
                                          reject_iMIT reject_phy_and_iMIT reject_iNOAA
                                          reject_phy_and_iNOAA reject_iMIT_and_iNOAA
                                          reject_phy_and_iMIT_and_iNOAA" ;
    '''

    data_labels = [
                    'pres',
                    'temp',
                    'dwpt',
                    'wvap',
                    'relh',
                    'lat',
                    'lon',
                    'qc'
                    ]

    sounding_inputs = {}
    for label in data_labels:
        sounding_inputs[label] = {}

    sounding_inputs['pres']['dset_name'] = 'Pressure'
    sounding_inputs['temp']['dset_name'] = 'Temperature'
    sounding_inputs['wvap']['dset_name'] = 'H2O_MR'
    sounding_inputs['lat']['dset_name']  = 'Latitude'
    sounding_inputs['lon']['dset_name']  = 'Longitude'
    sounding_inputs['qc']['dset_name']  = 'Quality_Flag'
    sounding_inputs['dwpt'] = None
    sounding_inputs['relh'] = None

    LOG.debug("Input file list: {}".format(nucaps_file_list))

    first_granule = True

    for grans in np.arange(len(nucaps_file_list)):

        try:

            LOG.debug("Ingesting granule {} ...".format(grans))

            nucaps_file = nucaps_file_list[grans]
            LOG.debug("Ingesting granule {} ...".format(nucaps_file))

            # Open the file object
            fileObj = xarray.open_dataset(nucaps_file, decode_times=False, mask_and_scale=False, group='/')

            # Open the dataset nodes for each of the required datasets
            sounding_inputs['pres']['node'] = fileObj['Pressure']
            sounding_inputs['temp']['node'] = fileObj['Temperature']
            sounding_inputs['wvap']['node'] = fileObj['H2O_MR']
            sounding_inputs['lat']['node'] = fileObj['Latitude']
            sounding_inputs['lon']['node'] = fileObj['Longitude']
            sounding_inputs['qc']['node'] = fileObj['Quality_Flag']


            # Fill in some missing attributes
            sounding_inputs['pres']['units'] = 'hPa'
            sounding_inputs['temp']['units'] = 'K'
            sounding_inputs['wvap']['units'] = 'g/kg'

            # Determine the shape of latitude dataset
            nucaps_shape = (int(sounding_inputs['lat']['node'].shape[0])/30, 30)

            if first_granule:

                # Get the pressure scale
                pressure = sounding_inputs['pres']['node'].data[0]

                for label in data_labels:
                    if sounding_inputs[label] is not None:
                        LOG.info(">>> Processing {} ...".format(sounding_inputs[label]['dset_name']))
                        LOG.debug("{}".format(sounding_inputs[label]['node'].attrs))
                        for attr_name in sounding_inputs[label]['node'].attrs.keys():
                            attr = sounding_inputs[label]['node'].attrs[attr_name]
                            LOG.debug("{} = {}".format(attr_name, attr))
                            sounding_inputs[label][attr_name] = attr

                # Search for the pressure level closest to the input value
                pressure_scope = 10.
                LOG.debug("Scope of pressure level search is {} hPa"
                        .format(pressure_scope))
                level = get_pressure_index(pressure,pres_0=pres_0,
                        kd_dist=pressure_scope)
                attempts = 1
                while (level is None and attempts<10):
                    pressure_scope *= 1.1
                    LOG.debug("Scope of pressure level search is {} hPa"
                            .format(pressure_scope))
                    level = get_pressure_index(pressure,pres_0=pres_0,
                            kd_dist=pressure_scope)
                    attempts += 1

                if level is None:
                    raise Exception("No suitable pressure level found, aborting...")

                LOG.debug("Retrieved level = {}".format(level))
                sounding_inputs['pres_0'] = pressure[level]
                data_labels.remove('pres')

            # Add to the lat and lon and data arrays...
            for label in ['lat','lon']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],
                        sounding_inputs[label]['node'].data.reshape(*nucaps_shape)))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'].data.reshape(*nucaps_shape)
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            # Add to the temp, dewpoint, water vapour and relative humidity data arrays...
            for label in ['temp','wvap']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'].data[:,level].reshape(*nucaps_shape)))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'].data[:,level].reshape(*nucaps_shape)
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            for label in ['qc']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'].data.reshape(*nucaps_shape)))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'].data.reshape(*nucaps_shape)
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            first_granule = False

            LOG.info("Closing {}".format(nucaps_file))
            fileObj.close()

        except Exception, err :
            LOG.warn("There was a problem, closing {}".format(nucaps_file))
            LOG.warn("{}".format(err))
            LOG.debug(traceback.format_exc())
            fileObj.close()
            LOG.info("Exiting...")
            sys.exit(1)

    # Convert water vapor from g/g to g/kg
    sounding_inputs['wvap']['data'] = 1000.*sounding_inputs['wvap']['data']

    # Contruct the pressure array
    pressure_array = pressure[level] * \
            np.ones(sounding_inputs['lat']['data'].shape,dtype='float')
    LOG.debug("pressure_array.shape =  {}".format(pressure_array.shape))

    # Construct the masks of the various datasets
    for label in ['temp','wvap','lat','lon']:
        if sounding_inputs[label] is not None:
            LOG.debug(">>> Processing {} ..."
                    .format(sounding_inputs[label]['dset_name']))
            fill_value = sounding_inputs[label]['_FillValue']
            data = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            LOG.debug("ma.is_masked({}) = {}".format(label,ma.is_masked(data)))

            if ma.is_masked(data):
                if data.mask.shape == ():
                    data_mask = np.ones(data.shape,dtype='bool')
                else:
                    data_mask = data.mask
            else:
                data_mask = np.zeros(data.shape,dtype='bool')

            LOG.debug("There are {}/{} masked values in {}".\
                    format(np.sum(data_mask),data.size,label))

            sounding_inputs[label]['mask'] = data_mask

    # Contruct the QC mask
    for label in ['qc']:
        if sounding_inputs[label] is not None:
            LOG.debug(">>> Processing {} ..."
                    .format(sounding_inputs[label]['dset_name']))
            fill_value = sounding_inputs[label]['_FillValue']
            data_fill = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            data_qc = ma.masked_greater(sounding_inputs[label]['data'],0)
            LOG.debug("{} has fill: {}".format(label,ma.is_masked(data_fill)))
            LOG.debug("{} has QC flag: {}".format(label,ma.is_masked(data_qc)))

            if ma.is_masked(data_fill):
                if data_fill.mask.shape == ():
                    data_fill_mask = np.ones(data_fill.shape,dtype='bool')
                else:
                    data_fill_mask = data_fill.mask
            else:
                data_fill_mask = np.zeros(data_fill.shape,dtype='bool')

            #data_fill_mask = np.zeros(data_fill.shape,dtype='bool') # DEBUG
            LOG.debug("There are {}/{} fill values in {}".\
                    format(np.sum(data_fill_mask),data_fill.size,label))

            if ma.is_masked(data_qc):
                if data_qc.mask.shape == ():
                    data_qc_mask = np.ones(data_qc.shape,dtype='bool')
                else:
                    data_qc_mask = data_qc.mask
                    # If the QC flags a total a bust, the QC flags might be crap...
                    if np.sum(data_qc_mask) == data_qc.size:
                        data_qc_mask = np.zeros(data_qc.shape,dtype='bool')
            else:
                data_qc_mask = np.zeros(data_qc.shape,dtype='bool')

            LOG.debug("There are {}/{} qc values in {}".\
                    format(np.sum(data_qc_mask),data_qc.size,label))

            sounding_inputs[label]['mask'] = data_fill_mask + data_qc_mask

            LOG.debug("There are {}/{} masked values in total {}".format(
                np.sum(sounding_inputs[label]['mask']),sounding_inputs[label]['data'].size, label))

    # Computing the dewpoint temperature and relative humidity
    sounding_inputs['dwpt'] = {}
    sounding_inputs['relh'] = {}
    LOG.debug("wvap.shape =  {}".format(sounding_inputs['wvap']['data'].shape))
    LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    dwpt_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']
    relh_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']

    wvap = ma.array(sounding_inputs['wvap']['data'],mask=relh_mask)
    temp = ma.array(sounding_inputs['temp']['data'],mask=relh_mask)
    pressure_array = ma.array(pressure_array,mask=relh_mask)

    dewhum_vec = np.vectorize(dewhum)
    dwpt,rh,_ = dewhum_vec(pressure_array,temp,wvap)

    # Copy dewpoint temperature to output
    sounding_inputs['dwpt']['data'] = dwpt
    sounding_inputs['dwpt']['mask'] = dwpt.mask if  ma.is_masked(dwpt) else np.zeros(dwpt.shape,dtype='bool')
    sounding_inputs['dwpt']['units'] = 'K'
    sounding_inputs['dwpt']['long_name'] = 'Dew-point Temperature in K'

    LOG.debug("ma.is_masked({}) = {}".format('dwpt',ma.is_masked(dwpt)))
    LOG.debug("There are {} masked values in dwpt".format(np.sum(dwpt.mask)))

    # Copy relative humidity to output
    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) else np.zeros(rh.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    return sounding_inputs


def plotData(data,data_mask,units,pngName,stride=1,dpi=200):
    '''
    Plot the input dataset in swath projection
    '''

    # General Setup
    dataShape = data.shape
    aspectRatio = float(dataShape[0])/float(dataShape[1])

    figWidth,figHeight = 4.,aspectRatio*4.
    ax_rect = [0.05, 0.15, 0.90, 0.80  ] # [left,bottom,width,height]

    fig = Figure(figsize=(figWidth,figHeight))
    canvas = FigureCanvas(fig)

    ax = fig.add_axes(ax_rect)

    data = ma.array(data,mask=data_mask)

    im = ax.imshow(data,interpolation='nearest')

    if units == '%':
        txt = ax.set_title('Reflectance')
    elif units == 'K':
        txt = ax.set_title('Brightness Temperature')
    else:
        txt = ax.set_title('scaled data')


    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    ax.set_aspect('equal')

    cax_rect = [0.05 , 0.05, 0.90 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')

    if units == '%':
        txt = cax.set_title('Reflectance (%)')
    elif units == 'K':
        txt = cax.set_title('Brightness Temperature (K)')
    else:
        txt = cax.set_title('scaled units')

    # Redraw the figure
    canvas.draw()
    canvas.print_figure(pngName,dpi=dpi)


def plotMapData(lat,lon,data,data_mask,pngName,**plot_options):

    # Copy the plot options to local variables
    title         = plot_options['title']
    cbar_title    = plot_options['cbar_title']
    units         = plot_options['units']
    stride        = plot_options['stride']
    lat_0         = plot_options['lat_0']
    lon_0         = plot_options['lon_0']
    pres_0        = plot_options['pres_0']
    latMin        = plot_options['latMin']
    lonMin        = plot_options['lonMin']
    latMax        = plot_options['latMax']
    lonMax        = plot_options['lonMax']
    plotMin       = plot_options['plotMin']
    plotMax       = plot_options['plotMax']
    bounding_lat  = plot_options['bounding_lat']
    scale         = plot_options['scale']
    map_res       = plot_options['map_res']
    proj          = plot_options['proj']
    #cmap          = plot_options['cmap']
    doScatterPlot = plot_options['scatterPlot']
    pointSize     = plot_options['pointSize']
    dpi           = plot_options['dpi']

    '''
    Plot the input dataset in mapped to particular projection
    '''

    # If our data is all missing, return
    if (np.sum(data_mask) == data.size):
        LOG.warn("Entire {} dataset is missing, aborting".\
                format(cbar_title))
        return -1

    # Determine if this is a global projection, so we can avoid having the
    # global "outset" plot.
    global_plot = False
    if (proj=='moll'
            or proj=='hammer'
            or proj=='robin'
            or proj=='eck4'
            or proj=='kav7'
            or proj=='mbtfpq'
            or proj=='sinu'
            or proj=='cyl'
            or proj=='merc'
            or proj=='mill'
            or proj=='gall'
            or proj=='cea'
            or proj=='vandg'):
        global_plot = True

    # Determine if this is a "rectangular" projection, so we can set sensible
    # latitude limits.
    rect_plot = False
    if (proj=='cyl'
            or proj=='merc'
            or proj=='mill'
            or proj=='gall'
            or proj=='cea'):
        rect_plot = True

    if (proj=='geos' or proj=='ortho' or global_plot):
        LOG.info("Setting lat_0 and lon_0 for {} projection.".format(proj))
        lat_0 = 0. if lat_0 is None else lat_0
        lon_0 = 0. if lon_0 is None else lon_0
        boundinglat = None

    if (proj=='npstere' and bounding_lat<-30.):
        LOG.warn("""North-Polar Stereographic bounding latitude (--bounding_lat)
         of {} degrees is too low, setting to -30 degrees latitude.""".format(bounding_lat))
        bounding_lat = -30.

    if (proj=='spstere' and bounding_lat>30.):
        LOG.warn("""South-Polar Stereographic bounding latitude (--bounding_lat)
         of {} degrees is too high, setting to 30 degrees latitude.""".format(bounding_lat))
        bounding_lat = 30.

    if (proj=='npaeqd' and bounding_lat<-30.):
        LOG.warn("""North-Polar Azimuthal Equidistant bounding latitude (--bounding_lat)
         of {} degrees is too low, setting to -30 degrees latitude.""".format(bounding_lat))
        bounding_lat = -30.

    if (proj=='spaeqd' and bounding_lat>30.):
        LOG.warn("""South-Polar Azimuthal Equidistant bounding latitude (--bounding_lat)
         of {} degrees is too high, setting to 30 degrees latitude.""".format(bounding_lat))
        bounding_lat = 30.

    if (proj=='nplaea' and bounding_lat<=0.):
        LOG.warn("""North-Polar Lambert Azimuthal bounding latitude (--bounding_lat)
         of {} degrees must be in northern hemisphere, setting to +1 degrees latitude.""".format(bounding_lat))
        bounding_lat = 1.

    if (proj=='splaea' and bounding_lat>=0.):
        LOG.warn("""South-Polar Lambert Azimuthal bounding latitude (--bounding_lat)
         of {} degrees must be in southern hemisphere, setting to -1 degrees latitude.""".format(bounding_lat))
        bounding_lat = -1.

    # Compute the central lat and lon if they are not specified
    LOG.debug("lat_0,lon_0= ({},{})".format(lat_0,lon_0))
    if (lat_0 is None) or (lon_0 is None):
        geo_shape = lat.shape
        nrows, ncols = geo_shape[0],geo_shape[1]
        LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
        # Non lat/lon pair given, use the central unmasked values.
        row_idx = int(nrows/2.)
        col_idx = int(ncols/2.)
        lat_0 = lat[row_idx,col_idx] if lat_0 is None else lat_0
        lon_0 = lon[row_idx,col_idx] if lon_0 is None else lon_0
        LOG.info("Incomplete lat/lon pair given, using ({:4.2f},{:4.2f})".
                format(lat_0,lon_0))

    LOG.info("Latitude extent of data:  ({:4.2f},{:4.2f})".format(np.min(lat),np.max(lat)))
    LOG.info("Longitude extent of data: ({:4.2f},{:4.2f})".format(np.min(lon),np.max(lon)))

    # General Setup
    if global_plot:
        figWidth,figHeight = 10.,6.
        ax_rect = [0.05, 0.17, 0.90, 0.72  ] # [left,bottom,width,height]
    else:
        figWidth,figHeight = 8.,8.
        ax_rect = [0.05, 0.15, 0.90, 0.75  ] # [left,bottom,width,height]

    fig = Figure(figsize=(figWidth,figHeight))
    canvas = FigureCanvas(fig)

    ax = fig.add_axes(ax_rect)

    # Common plotting options...
    plot_kw = {
        'ax'         : ax,
        'projection' : proj,
        'lon_0'      : lon_0,
        'lat_0'      : lat_0,
        'fix_aspect' : True,
        'resolution' : map_res
    }

    # Projection dependent plotting options
    if global_plot:
        if rect_plot:
            plot_kw['llcrnrlat'] =  -80. if latMin is None else latMin
            plot_kw['urcrnrlat'] =   80. if latMax is None else latMax
            plot_kw['llcrnrlon'] = -180. if lonMin is None else lonMin
            plot_kw['urcrnrlon'] =  180. if lonMax is None else lonMax
        else:
            pass
    else:
        LOG.info("scale = {}".format(scale))
        windowWidth = scale *(0.35*12000000.)
        windowHeight = scale *(0.50*9000000.)
        plot_kw['width'] = windowWidth
        plot_kw['height'] = windowHeight

    if (proj=='npstere' or proj=='spstere' or
            proj=='nplaea' or proj=='splaea' or
            proj=='npaeqd' or proj=='spaeqd'):
        plot_kw['boundinglat'] = bounding_lat

    if ((proj=='aea' or proj=='lcc' or proj=='eqdc') and (np.abs(lat_0)<1.e-6)):
        plot_kw['lat_1'] = 1.
        plot_kw['lat_2'] = 1.

    # Create the Basemap plotting object
    LOG.info("Plotting for {} projection...".format(proj))
    try:
        m = Basemap(**plot_kw)
    except ValueError, err :
            LOG.error("{} ({} projection), aborting.".format(err,proj))
            return 1
    except RuntimeError, err :
            LOG.error("{} ({} projection), aborting.".format(err,proj))
            return 1

    x,y=m(lon[::stride,::stride],lat[::stride,::stride])

    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    m.fillcontinents(color='0.85',zorder=0)
    m.drawparallels(np.arange( -90, 91,30), color = '0.25',
            linewidth = 0.5,labels=[1,0,1,0],fontsize=9,labelstyle="+/-")
    m.drawmeridians(np.arange(-180,180,30), color = '0.25',
            linewidth = 0.5,labels=[1,0,1,0],fontsize=9,labelstyle="+/-")

    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    plotMin = np.min(data) if plotMin is None else plotMin
    plotMax = np.max(data) if plotMax is None else plotMax
    LOG.debug("plotMin = {}".format(plotMin))
    LOG.debug("plotMax = {}".format(plotMax))

    if doScatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none',
                vmin=plotMin,vmax=plotMax, cmap=cm.jet)
    else:
        cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False,
                vmin=plotMin,vmax=plotMax, cmap=cm.jet)

    txt = ax.set_title(title,fontsize=11)

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    cax_rect = [0.05 , 0.05, 0.90 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')

    txt = cax.set_title(cbar_title)

    #
    # Add a small globe with the swath indicated on it #
    #
    if not global_plot:
        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        glax_rect = [0.81, 0.75, 0.18, 0.20 ] # [left,bottom,width,height]
        glax = fig.add_axes(glax_rect)

        m_globe = Basemap(lat_0=0.,lon_0=0.,\
            ax=glax,resolution='c',area_thresh=10000.,projection='robin')

        # If we previously had a zero size data array, increase the pointSize
        # so the data points are visible on the global plot
        if (np.shape(lon[::stride,::stride])[0]==2) :
            pointSize = 5.

        x,y=m_globe(lon[::stride,::stride],lat[::stride,::stride])
        swath = np.zeros(np.shape(x),dtype=int)

        m_globe.drawmapboundary(linewidth=0.1)
        m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
        m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

        p_globe = m_globe.scatter(x,y,s=0.5,c="red",axes=glax,edgecolors='none',zorder=2)
    else:
        pass


    # Redraw the figure
    canvas.draw()
    LOG.info("Writing image file {}".format(pngName))
    canvas.print_figure(pngName,dpi=dpi)

    return 0

def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    dataChoices = OrderedDict([
        ('iapp', 'International TOVS Processing Package'),
        ('mirs', 'Microwave Integrated Retrieval System'),
        ('hsrtv', 'CrIS, AIRS and IASI Hyperspectral Retrieval Software'),
        ('nucaps', 'NOAA Unique Combined Atmospheric Processing System')
        ])
    prodChoices=['temp','wvap','dwpt','relh','ctp','ctt']
    map_res_choice = ['c','l','i']
    map_proj_choice = OrderedDict([
                       ('lcc',     'Lambert Conformal'),
                       ('cea',     'Cylindrical Equal Area'),
                       ('ortho',   'Orthographic'),
                       ('geos',    'Geostationary'),
                       ('npstere', 'North-Polar Stereographic'),
                       ('spstere', 'South-Polar Stereographic'),
                       ('cyl',     'Cylindrical Equidistant'),
                       ('sinu',    'Sinusoidal'),
                       ('merc',    'Mercator'),
                       #('mbtfpq',  'McBryde-Thomas Flat-Polar Quartic'),
                       #('poly',    'Polyconic'),
                       #('moll',    'Mollweide'),
                       #('tmerc',   'Transverse Mercator'),
                       #('gall',    'Gall Stereographic Cylindrical'),
                       #('mill',    'Miller Cylindrical'),
                       #('stere',   'Stereographic'),
                       #('eqdc',    'Equidistant Conic'),
                       #(#'rotpole', 'Rotated Pole'),
                       #('hammer',  'Hammer'),
                       #('nsper',   'Near-Sided Perspective'),
                       #('eck4',    'Eckert IV'),
                       #('aea',     'Albers Equal Area'),
                       #('kav7',    'Kavrayskiy VII'),
                       #('nplaea',  'North-Polar Lambert Azimuthal'),
                       #('splaea',  'South-Polar Lambert Azimuthal'),
                       #('npaeqd',  'North-Polar Azimuthal Equidistant'),
                       #('spaeqd',  'South-Polar Azimuthal Equidistant'),
                       #('cass',    'Cassini-Soldner'),
                       #('laea',    'Lambert Azimuthal Equal Area'),
                       #('robin',   'Robinson')
                       ])


    defaults = {
                'dataset':'temp',
                'stride':1,
                'pressure':850.,
                'lat_0':None,
                'lon_0':None,
                'plotMin'  : None,
                'plotMax'  : None,
                'bounding_lat':0.,
                'scatter_plot':False,
                'pointSize':1,
                'scale':1,
                'map_res':'c',
                'proj':'lcc',
                'output_file':None,
                'outputFilePrefix' : None,
                'dpi':200
                }

    is_expert = False
    if '--expert' in sys.argv:
        expert_index = sys.argv.index('--expert')
        sys.argv[expert_index] = '--help'
        is_expert = True
    elif '-x' in sys.argv:
        expert_index = sys.argv.index('-x')
        sys.argv[expert_index] = '--help'
        is_expert = True
    else:
        pass

    version = 'CSPP Sounder QL v1.2'
    description = '''Create a plot of temperature, dewpoint or something else at a particular ''' \
                  '''pressure level. Supports IAPP, MIRS, HSRTV and NUCAPS files.'''
    epilog = ''
    parser = argparse.ArgumentParser(
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     add_help=False,
                                     epilog=epilog
                                     )

    # Mandatory/positional arguments

    parser.add_argument(
                      action="store",
                      dest="datatype",
                      type=str,
                      choices=dataChoices,
                      nargs='?',
                      help="""The type of the input sounder data file. Possible values """ \
                           """are:\n {}""".format(
                               "".join(["\t'{0:8s} ({1}),\n".format(tups[0]+"'",tups[1]) for
                               tups in zip(dataChoices.keys(),dataChoices.values())]))
                      )

    parser.add_argument(
                      action='store',
                      dest='input_file',
                      type=str,
                      nargs='*',
                      help = '''One or more input files.'''
    )

    # Optional arguments

    parser.add_argument('--dset',
                      action="store",
                      dest="dataset",
                      default=defaults["dataset"],
                      type=str,
                      choices=prodChoices,
                      help='''The sounder dataset to plot. Possible values are... {}.''' \
                           ''' [default: {}]'''.format(prodChoices.__str__()[1:-1],
                               defaults["dataset"])
                      )

    parser.add_argument('-S','--stride',
                      action="store",
                      dest="stride",
                      default=defaults["stride"],
                      type=int,
                      help='''Sample every STRIDE pixels in the band data.''' \
                      ''' [default: {}]'''.format(defaults["stride"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--pressure',
                      action="store",
                      dest="pressure",
                      default=defaults["pressure"],
                      type=float,
                      help='''The pressure level (in mbar) to plot.''' \
                      ''' [default: {}]'''.format(defaults["pressure"])
                      )

    parser.add_argument('--plot-min',
                      action="store",
                      dest="plotMin",
                      default=defaults["plotMin"],
                      type=float,
                      help="Minimum value to plot.".format(defaults["plotMin"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--plot-max',
                      action="store",
                      dest="plotMax",
                      default=defaults["plotMax"],
                      type=float,
                      help="Maximum value to plot.".format(defaults["plotMax"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--lat-0',
                      action="store",
                      dest="lat_0",
                      type=float,
                      help="Center latitude of plot." if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--lon-0',
                      action="store",
                      dest="lon_0",
                      type=float,
                      help="Center longitude of plot." if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--lat-min',
                      action="store",
                      dest="latMin",
                      type=float,
                      help="Minimum latitude to plot." if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--lat-max',
                      action="store",
                      dest="latMax",
                      type=float,
                      help="Maximum latitude to plot." if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--lon-min',
                      action="store",
                      dest="lonMin",
                      type=float,
                      help="Minimum longitude to plot." if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--lon-max',
                      action="store",
                      dest="lonMax",
                      type=float,
                      help="Maximum longitude to plot." if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--bounding-lat',
                      action="store",
                      dest="bounding_lat",
                      default=defaults["bounding_lat"],
                      type=float,
                      help='''The minimum/maximum latitude to plot for the polar projections.''' \
                      ''' [default: {}]'''.format(defaults["bounding_lat"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--scatter-plot',
                      action="store_true",
                      dest="doScatterPlot",
                      default=defaults["scatter_plot"],
                      help="Generate the plot using a scatterplot approach." if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-P','--point-size',
                      action="store",
                      dest="pointSize",
                      default=defaults["pointSize"],
                      type=float,
                      help='''Size of the plot point used to represent each pixel.''' \
                      ''' [default: {}]'''.format(defaults["pointSize"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-s','--scale',
                      action="store",
                      dest="scale",
                      default=defaults["scale"],
                      type=float,
                      help='''The scaling factor for the default viewport size of (w x h) = ''' \
                      '''(4200000.0 x 4500000.0) meters.
                      [default: {}]'''.format(defaults["pointSize"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-m','--map-res',
                      action="store",
                      dest="map_res",
                      default=defaults["map_res"],
                      type=str,
                      choices=map_res_choice,
                      help="""The map coastline resolution. Possible values are 'c' (coarse),'l'""" \
                           """ (low) and 'i' (intermediate). """ \
                           """[default: '{}']""".format(defaults["map_res"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('--proj',
                      action="store",
                      dest="proj",
                      default=defaults["proj"],
                      type=str,
                      choices=map_proj_choice,
                      help="""The map projection. Possible values are:\n {}[default: '{}']""".format(
                           "".join(["\t'{0:8s} ({1}),\n".format(tups[0]+"'",tups[1]) for
                           tups in zip(map_proj_choice.keys(),map_proj_choice.values())]),
                           defaults["proj"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-o','--output_file',
                      action="store",
                      dest="output_file",
                      default=defaults["output_file"],
                      type=str,
                      help='''The filename of the output Skew-T png file. [default: {}]'''.format(
                          defaults["output_file"])
                      )

    parser.add_argument('-O','--output-file-prefix',
                      action="store",
                      dest="outputFilePrefix",
                      default=defaults["outputFilePrefix"],
                      type=str,
                      help="""String to prefix to the automatically generated png name. """ \
                           """[default: {}]""".format(
                               defaults["outputFilePrefix"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default='200.',
                      type=float,
                      help='''The resolution in dots per inch of the output png file.''' \
                      ''' [default: {}]'''.format(defaults["dpi"]) if is_expert else argparse.SUPPRESS
                      )

    parser.add_argument('-V', '--version',
                      action='version',
                      version=version,
                      help='''Print the CSPP Sounder QL package version'''
                      )

    parser.add_argument("-v", "--verbosity",
                      dest='verbosity',
                      action="count",
                      default=2,
                      help='''Each occurrence increases verbosity one level from INFO. -v=DEBUG'''
                      )

    parser.add_argument("-q", "--quiet",
                      action="store_true",
                      dest='quiet',
                      default=False,
                      help='''Silence all console output'''
                      )

    parser.add_argument('-h', '--help',
                      action="help",
                      help='''Show this help message and exit'''
                      )

    parser.add_argument('-x', '--expert',
                      dest='is_expert',
                      action="store_true",
                      default=False,
                      help='''Display all help options, including the expert ones.'''
                      )
    args = parser.parse_args()

    # Do any parser checking for conditional arguments
    if args.datatype=='MIRS' or args.datatype=='NUCAPS':
        if args.dataset=='ctp' or args.dataset=='ctt':
            parser.error("""Datasets 'ctp' and 'ctt' can only be used with IAPP and HSRTV.""")

    # Set up the logging
    verbosity = 0 if args.quiet else args.verbosity
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    level = levels[min(verbosity,3)]

    if level == logging.DEBUG :
        console_logFormat = '%(asctime)s.%(msecs)03d (%(levelname)s) : %(filename)s : %(funcName)s : %(lineno)d:%(message)s'
        dateFormat='%Y-%m-%d %H:%M:%S'
    else:
        console_logFormat = '%(asctime)s.%(msecs)03d (%(levelname)s) : %(message)s'
        dateFormat='%Y-%m-%d %H:%M:%S'

    logging.basicConfig(level = level,
            format = console_logFormat,
            datefmt = dateFormat)

    return args

def main():
    '''
    The main method.
    '''

    # Read in the options
    options = _argparse()

    input_file = options.input_file
    datatype = options.datatype
    dataset = options.dataset
    stride = options.stride
    pressure = options.pressure
    lat_0  = options.lat_0
    lon_0  = options.lon_0
    latMin = options.latMin
    lonMin = options.lonMin
    latMax = options.latMax
    lonMax = options.lonMax
    bounding_lat = options.bounding_lat
    plotMin = options.plotMin
    plotMax = options.plotMax
    scale = options.scale
    doScatterPlot = options.doScatterPlot
    pointSize = options.pointSize
    proj = options.proj
    map_res = options.map_res
    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix
    dpi = options.dpi

    # Obtain list of the product files
    input_file_list = generate_file_list(input_file)

    if input_file_list==[]:
        LOG.warn('No files match the input file glob "{}", aborting.'
                .format(input_file))
        sys.exit(1)


    LOG.info("input file(s): {}".format(input_file_list))

    # Read in the input file, and return a dictionary containing the required
    # data

    if datatype == 'iapp' :
        sounding_inputs = iapp_sounder(input_file_list,pres_0=pressure)
    elif datatype == 'mirs' :
        sounding_inputs = mirs_sounder(input_file_list,pres_0=pressure)
    elif datatype == 'hsrtv' :
        sounding_inputs = hsrtv_sounder(input_file_list,pres_0=pressure)
    elif datatype == 'nucaps' :
        sounding_inputs = nucaps_sounder(input_file_list,pres_0=pressure)
    else:
        pass

    lats = sounding_inputs['lat']['data']
    lons = sounding_inputs['lon']['data']
    data = sounding_inputs[dataset]['data']
    data_mask = sounding_inputs[dataset]['mask'] + sounding_inputs['qc']['mask']
    units = sounding_inputs[dataset]['units']

    LOG.debug("There are {}/{} masked values in {}".format(np.sum(data_mask),data.size, dataset))

    pres_0 = int(sounding_inputs['pres_0'])

    input_file = path.basename(input_file_list[0])

    cbar_titles = {
            'temp':'temperature (K) @ {:4.2f} hPa'.format(pres_0),
            'wvap':'water vapor mixing ratio (g/kg) @ {:4.2f} hPa'.format(pres_0),
            'dwpt':'dewpoint temperature (K) @ {:4.2f} hPa'.format(pres_0),
            'relh':'relative humidity (%) @ {:4.2f} hPa'.format(pres_0),
            'ctp':'cloud top pressure (hPa)',
            'ctt':'cloud top temperature (K)'
            }

    # Determine the filename
    if dataset=='ctp' or dataset=='ctt':
        file_suffix = "{}_{}".format(datatype,dataset)
    else:
        file_suffix = "{}_{}_{}mb".format(datatype,dataset,pres_0)

    LOG.info('output_file = {}'.format(output_file))
    LOG.info('outputFilePrefix = {}'.format(outputFilePrefix))
    LOG.info('input_file = {}'.format(input_file))
    LOG.info('file_suffix = {}'.format(file_suffix))

    if output_file is None and outputFilePrefix is None :
        output_file = "{}.{}.png".format(input_file,file_suffix)
    if output_file is not None and outputFilePrefix is None :
        pass
    if output_file is None and outputFilePrefix is not None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)
    if output_file is not None and outputFilePrefix is not None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)

    # Default plot labels
    plot_options = {}
    plot_options['title'] = "{}\n\n".format(input_file)
    plot_options['cbar_title'] = cbar_titles[dataset]
    plot_options['units'] = sounding_inputs[dataset]['units']
    plot_options['stride'] = stride
    plot_options['lat_0'] = lat_0
    plot_options['lon_0'] = lon_0
    plot_options['pres_0'] = sounding_inputs['pres_0']
    plot_options['latMin'] = latMin
    plot_options['lonMin'] = lonMin
    plot_options['latMax'] = latMax
    plot_options['lonMax'] = lonMax
    plot_options['bounding_lat'] = bounding_lat
    plot_options['plotMin'] = plotMin
    plot_options['plotMax'] = plotMax
    plot_options['scale'] = scale
    plot_options['map_res'] = map_res
    plot_options['proj'] = proj
    #plot_options['cmap'] = cmap
    plot_options['scatterPlot'] = doScatterPlot
    plot_options['pointSize'] = pointSize
    plot_options['dpi'] = dpi

    # Create the plot
    retval = plotMapData(lats,lons,data,data_mask,output_file,**plot_options)

    return retval

if __name__=='__main__':
    sys.exit(main())
