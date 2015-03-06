#!/usr/bin/env python
# encoding: utf-8
"""
sounder_image.py

Purpose: Create a plot of temperature, dewpoint or relative humidity,
         at a particular pressure level. Supports outputs from the following 
         packages...
         
         * International ATOVS Processing Package (IAPP)
         * Microwave Integrated Retrieval System (MIRS)
         * CSPP Hyperspectral Retrieval (Dual Regression) Package
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
    directory or a file glob.

    DATATYPE: One of 'IAPP','MIRS', 'DR' or 'NUCAPS'.


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2014-05-10.
Copyright (c) 2014 University of Wisconsin Regents. All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

file_Date = '$Date: 2015-02-11 22:59:24 -0800 (Wed, 11 Feb 2015) $'
file_Revision = '$Revision: 2354 $'
file_Author = '$Author: geoffc $'
file_HeadURL = '$HeadURL: https://svn.ssec.wisc.edu/repos/jpss_adl/trunk/scripts/iapp/quicklooks/sounder_image.py $'
file_Id = '$Id: sounder_image.py 2354 2015-02-12 06:59:24Z geoffc $'

__author__ = 'Geoff Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id: sounder_image.py 2354 2015-02-12 06:59:24Z geoffc $'
__docformat__ = 'Epytext'

import os, sys, logging, traceback
from os import path,uname,environ
import string
import re
import uuid
from shutil import rmtree,copyfile
from glob import glob
from time import time
from datetime import datetime,timedelta

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

from netCDF4 import Dataset
from netCDF4 import num2date
import h5py

import skewt as skT
from thermo import mr_to_rh, rh_to_mr, dewpoint_magnus, dewpoint_AB
from thermo import MixR2VaporPress, DewPoint

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
        
        LOG.info("Retrieved pressure is {:4.2f} mb".\
                format(pressure[level].squeeze()))

    except Exception, err :
        LOG.warn("{}".format(err))
        LOG.debug(traceback.format_exc())

    return level


def iapp_sounder(iapp_file,pres_0=850.):
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
    '''

    data_labels = [
                    'pres',
                    'temp',
                    'dwpt',
                    'wvap',
                    'relh',
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
    sounding_inputs['lat']['dset_name']  = 'Latitude'
    sounding_inputs['lon']['dset_name']  = 'Longitude'
    sounding_inputs['relh'] = None


    # Open the file object
    fileObj = Dataset(iapp_file)

    try:

        sounding_inputs['pres']['node'] = fileObj.variables['Pressure_Levels']
        sounding_inputs['temp']['node'] = fileObj.variables['Temperature_Retrieval']
        sounding_inputs['dwpt']['node'] = fileObj.variables['Dew_Point_Temp_Retrieval']
        sounding_inputs['wvap']['node'] = fileObj.variables['WaterVapor_Retrieval']
        sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
        sounding_inputs['lon']['node'] = fileObj.variables['Longitude']

        pressure = sounding_inputs['pres']['node'][:]

        for label in data_labels:
            if sounding_inputs[label] != None:
                LOG.info(">>> Processing {} ...".format(sounding_inputs[label]['dset_name']))
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
        while (level==None and attempts<10):
            pressure_scope *= 1.1
            LOG.debug("Scope of pressure level search is {} hPa"
                    .format(pressure_scope))
            level = get_pressure_index(pressure,pres_0=pres_0,
                    kd_dist=pressure_scope)
            attempts += 1

        if level==None:
            raise Exception("No suitable pressure level found, aborting...")
             
        LOG.debug("Retrieved level = {}".format(level))
        sounding_inputs['pres_0'] = pressure[level]
        #row,col = 10,10

        sounding_inputs['lat']['data'] = sounding_inputs['lat']['node'][:,:]
        sounding_inputs['lon']['data'] = sounding_inputs['lon']['node'][:,:]
        sounding_inputs['temp']['data'] = sounding_inputs['temp']['node'][:,:,level]
        sounding_inputs['dwpt']['data'] = sounding_inputs['dwpt']['node'][:,:,level]
        sounding_inputs['wvap']['data'] = sounding_inputs['wvap']['node'][:,:,level]

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
    for label in ['temp','dwpt','wvap','lat','lon']:
        if sounding_inputs[label] != None:
            LOG.debug(">>> Processing {} ...".format(sounding_inputs[label]['dset_name']))
            #LOG.debug('{} --> sounding_inputs[label]['data'])
            fill_value = sounding_inputs[label]['_FillValue']
            #data = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            data = ma.masked_equal(sounding_inputs[label]['data'],np.nan)
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

    # Computing the relative humidity
    sounding_inputs['relh'] = {}
    LOG.debug("wvap.shape =  {}".format(sounding_inputs['wvap']['data'].shape))
    LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    relh_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']

    wvap = ma.array(sounding_inputs['wvap']['data'],mask=relh_mask)
    temp = ma.array(sounding_inputs['temp']['data'],mask=relh_mask)
    pressure_array = ma.array(pressure_array,mask=relh_mask)

    mr_to_rh_vec = np.vectorize(mr_to_rh)
    rh = mr_to_rh_vec(wvap,pressure_array,temp).real

    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            else np.zeros(temp.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    return sounding_inputs


def mirs_sounder(mirs_file,pres_0=850.):
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
    '''         

    data_labels = [
                    'pres',
                    'temp',
                    'dwpt',
                    'wvap',
                    'relh',
                    'lat',
                    'lon'
                    ]

    sounding_inputs = {}
    for label in data_labels:
        sounding_inputs[label] = {}

    sounding_inputs['pres']['dset_name'] = 'Player'
    sounding_inputs['temp']['dset_name'] = 'PTemp'
    sounding_inputs['wvap']['dset_name'] = 'PVapor'
    sounding_inputs['lat']['dset_name']  = 'Latitude'
    sounding_inputs['lon']['dset_name']  = 'Longitude'
    sounding_inputs['dwpt'] = None
    sounding_inputs['relh'] = None

    # Open the file object
    fileObj = Dataset(mirs_file)

    try:

        sounding_inputs['pres']['node'] = fileObj.variables['Player']
        sounding_inputs['temp']['node'] = fileObj.variables['PTemp']
        sounding_inputs['wvap']['node'] = fileObj.variables['PVapor']
        sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
        sounding_inputs['lon']['node'] = fileObj.variables['Longitude']

        pressure = sounding_inputs['pres']['node'][:]

        # Fill in some missing attributes
        fill_value = getattr(fileObj,'missing_value')
        sounding_inputs['pres']['units'] = 'hPa'
        sounding_inputs['temp']['units'] = 'K'
        sounding_inputs['wvap']['units'] = 'g/kg'

        for label in data_labels:
            if sounding_inputs[label] != None:
                LOG.info(">>> Processing {} ..."
                        .format(sounding_inputs[label]['dset_name']))
                sounding_inputs[label]['_FillValue'] = fill_value
                for attr_name in sounding_inputs[label]['node'].ncattrs(): 
                    attr = getattr(sounding_inputs[label]['node'],attr_name)
                    LOG.debug("{} = {}".format(attr_name,attr))
                    sounding_inputs[label][attr_name] = attr

        # Search for the pressure level closest to the input value
        LOG.debug("Shape of pressure is {}".format(pressure.shape))
        LOG.debug("pressure is {}".format(pressure))

        pressure_scope = 10.
        LOG.debug("Scope of pressure level search is {} hPa"
                .format(pressure_scope))
        level = get_pressure_index(pressure,pres_0=pres_0,
                kd_dist=pressure_scope)
        attempts = 1
        while (level==None and attempts<10):
            pressure_scope *= 1.1
            LOG.debug("Scope of pressure level search is {} hPa"
                    .format(pressure_scope))
            level = get_pressure_index(pressure,pres_0=pres_0,
                    kd_dist=pressure_scope)
            attempts += 1

        if level==None:
            raise Exception("No suitable pressure level found, aborting...")
             
        LOG.debug("Retrieved level = {}".format(level))
        sounding_inputs['pres_0'] = pressure[level]
        #row,col = 10,10

        sounding_inputs['lat']['data'] = sounding_inputs['lat']['node'][:,:]
        sounding_inputs['lon']['data'] = sounding_inputs['lon']['node'][:,:]
        sounding_inputs['temp']['data'] = sounding_inputs['temp']['node'][:,:,level]
        sounding_inputs['wvap']['data'] = sounding_inputs['wvap']['node'][:,:,level]

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
        if sounding_inputs[label] != None:
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

    # Computing the relative humidity
    sounding_inputs['relh'] = {}
    LOG.debug("wvap.shape =  {}".format(sounding_inputs['wvap']['data'].shape))
    LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    relh_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']

    wvap = ma.array(sounding_inputs['wvap']['data'],mask=relh_mask)
    temp = ma.array(sounding_inputs['temp']['data'],mask=relh_mask)
    pressure_array = ma.array(pressure_array,mask=relh_mask)

    mr_to_rh_vec = np.vectorize(mr_to_rh)
    rh = mr_to_rh_vec(wvap,pressure_array,temp).real

    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            else np.zeros(rh.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    # Computing the dewpoint temperature
    sounding_inputs['dwpt'] = {}
    LOG.debug("relh.shape =  {}".format(sounding_inputs['relh']['data'].shape))
    LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    dwpt_mask = sounding_inputs['relh']['mask']+sounding_inputs['temp']['mask']

    rh = ma.array(sounding_inputs['relh']['data'],mask=dwpt_mask)
    temp = ma.array(sounding_inputs['temp']['data'],mask=dwpt_mask)
    pressure_array = ma.array(pressure_array,mask=dwpt_mask)

    dewpoint_AB_vec = np.vectorize(dewpoint_AB)
    dwpt = dewpoint_AB_vec(temp,rh) + 273.15

    sounding_inputs['dwpt']['data'] = dwpt
    sounding_inputs['dwpt']['mask'] = dwpt.mask if  ma.is_masked(dwpt) \
            else np.zeros(dwpt.shape,dtype='bool')
    sounding_inputs['dwpt']['units'] = 'K'
    sounding_inputs['dwpt']['long_name'] = 'Dew-point Temperature in K'

    LOG.debug("ma.is_masked({}) = {}".format('dwpt',ma.is_masked(dwpt)))
    LOG.debug("There are {} masked values in dwpt".format(np.sum(dwpt.mask)))
    LOG.debug("dwpt.shape =  {}".format(sounding_inputs['dwpt']['data'].shape))
    LOG.debug("dwpt_mask.shape =  {}".format(sounding_inputs['dwpt']['mask'].shape))

    return sounding_inputs


def dual_regression_sounder(dr_file,pres_0=850.):
    '''
    Plevs: Pressure for each level in mb (hPa)
    TAir: Temperature profile in K
    Dewpnt: Water vapor profile (mixing ratio) in g/kg

    /Latitude
    /Longitude
    /Plevs
    /TAir
    /Dewpnt
    '''

    data_labels = [
                    'pres',
                    'temp',
                    'dwpt',
                    'wvap',
                    'relh',
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
    sounding_inputs['lat']['dset_name']  = 'Latitude'
    sounding_inputs['lon']['dset_name']  = 'Longitude'
    sounding_inputs['wvap'] = None

    # Open the file object
    if h5py.is_hdf5(dr_file):
        fileObj = h5py.File(dr_file,'r')
    else:
        LOG.error("{} does not exist or is not a valid HDF5 file, aborting...".\
                format(dr_file))
        sys.exit(1)

    try:

        sounding_inputs['pres']['node'] = fileObj['Plevs']
        sounding_inputs['temp']['node'] = fileObj['TAir']
        sounding_inputs['dwpt']['node'] = fileObj['Dewpnt']
        sounding_inputs['relh']['node'] = fileObj['RelHum']
        sounding_inputs['lat']['node'] = fileObj['Latitude']
        sounding_inputs['lon']['node'] = fileObj['Longitude']

        pressure = sounding_inputs['pres']['node'][:]

        #FIXME: All attributes for the DR datasets appear to be strings, 
        #       whether the attribute value is a string or not. So we can't 
        #       determine the attribute type from the file, we just have 
        #       *know* it. Why bother with self-descibing file formats then?
        for label in data_labels:
            if sounding_inputs[label] != None:
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
        while (level==None and attempts<10):
            pressure_scope *= 1.1
            LOG.debug("Scope of pressure level search is {} hPa"
                    .format(pressure_scope))
            level = get_pressure_index(pressure,pres_0=pres_0,
                    kd_dist=pressure_scope)
            attempts += 1

        if level==None:
            raise Exception("No suitable pressure level found, aborting...")
             
        LOG.debug("Retrieved level = {}".format(level))
        sounding_inputs['pres_0'] = pressure[level]

        sounding_inputs['lat']['data'] = sounding_inputs['lat']['node'][:,:]
        sounding_inputs['lon']['data'] = sounding_inputs['lon']['node'][:,:]
        sounding_inputs['temp']['data'] = sounding_inputs['temp']['node'][level,:,:]
        sounding_inputs['dwpt']['data'] = sounding_inputs['dwpt']['node'][level,:,:]
        sounding_inputs['relh']['data'] = sounding_inputs['relh']['node'][level,:,:]
        LOG.info("Copied the datasets from {}".format(dr_file))

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
    for label in ['temp','dwpt','relh','lat','lon']:
        if sounding_inputs[label] != None:
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


    return sounding_inputs


def nucaps_sounder(nucaps_file,pres_0=850.):
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
    '''         

    data_labels = [
                    'pres',
                    'temp',
                    'dwpt',
                    'wvap',
                    'relh',
                    'lat',
                    'lon'
                    ]

    sounding_inputs = {}
    for label in data_labels:
        sounding_inputs[label] = {}

    sounding_inputs['pres']['dset_name'] = 'Pressure'
    sounding_inputs['temp']['dset_name'] = 'Temperature'
    sounding_inputs['wvap']['dset_name'] = 'H2O_MR'
    sounding_inputs['lat']['dset_name']  = 'Latitude'
    sounding_inputs['lon']['dset_name']  = 'Longitude'
    sounding_inputs['dwpt'] = None
    sounding_inputs['relh'] = None

    # Open the file object
    fileObj = Dataset(nucaps_file)

    try:

        sounding_inputs['pres']['node'] = fileObj.variables['Pressure']
        sounding_inputs['temp']['node'] = fileObj.variables['Temperature']
        sounding_inputs['wvap']['node'] = fileObj.variables['H2O_MR']
        sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
        sounding_inputs['lon']['node'] = fileObj.variables['Longitude']

        pressure = sounding_inputs['pres']['node'][0,:]

        # Fill in some missing attributes
        sounding_inputs['pres']['units'] = 'hPa'
        sounding_inputs['temp']['units'] = 'K'
        sounding_inputs['wvap']['units'] = 'g/kg'

        for label in data_labels:
            if sounding_inputs[label] != None:
                LOG.info(">>> Processing {} ..."
                        .format(sounding_inputs[label]['dset_name']))
                for attr_name in sounding_inputs[label]['node'].ncattrs(): 
                    attr = getattr(sounding_inputs[label]['node'],attr_name)
                    LOG.debug("{} = {}".format(attr_name,attr))
                    sounding_inputs[label][attr_name] = attr

        # Search for the pressure level closest to the input value
        LOG.debug("Shape of pressure is {}".format(pressure.shape))
        LOG.debug("pressure is {}".format(pressure))

        pressure_scope = 10.
        LOG.debug("Scope of pressure level search is {} hPa"
                .format(pressure_scope))
        level = get_pressure_index(pressure,pres_0=pres_0,
                kd_dist=pressure_scope)
        attempts = 1
        while (level==None and attempts<10):
            pressure_scope *= 1.1
            LOG.debug("Scope of pressure level search is {} hPa"
                    .format(pressure_scope))
            level = get_pressure_index(pressure,pres_0=pres_0,
                    kd_dist=pressure_scope)
            attempts += 1

        if level==None:
            raise Exception("No suitable pressure level found, aborting...")
             
        LOG.debug("Retrieved level = {}".format(level))
        sounding_inputs['pres_0'] = pressure[level]
        #row,col = 10,10

        sounding_inputs['lat']['data'] = sounding_inputs['lat']['node'][:].reshape(4,30)
        sounding_inputs['lon']['data'] = sounding_inputs['lon']['node'][:].reshape(4,30)
        sounding_inputs['temp']['data'] = sounding_inputs['temp']['node'][:,level].reshape(4,30)
        # Convert water vapor from g/g to g/kg
        sounding_inputs['wvap']['data'] = 1000.*sounding_inputs['wvap']['node'][:,level].reshape(4,30)

        LOG.info("Closing {}".format(nucaps_file))
        fileObj.close()

    except Exception, err :
        LOG.warn("There was a problem, closing {}".format(nucaps_file))
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
        if sounding_inputs[label] != None:
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

    # Computing the relative humidity
    sounding_inputs['relh'] = {}
    LOG.debug("wvap.shape =  {}".format(sounding_inputs['wvap']['data'].shape))
    LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    relh_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']

    wvap = ma.array(sounding_inputs['wvap']['data'],mask=relh_mask)
    temp = ma.array(sounding_inputs['temp']['data'],mask=relh_mask)
    pressure_array = ma.array(pressure_array,mask=relh_mask)

    mr_to_rh_vec = np.vectorize(mr_to_rh)
    rh = mr_to_rh_vec(wvap,pressure_array,temp).real

    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            else np.zeros(rh.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    # Computing the dewpoint temperature
    sounding_inputs['dwpt'] = {}
    LOG.debug("relh.shape =  {}".format(sounding_inputs['relh']['data'].shape))
    LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    dwpt_mask = sounding_inputs['relh']['mask']+sounding_inputs['temp']['mask']

    rh = ma.array(sounding_inputs['relh']['data'],mask=dwpt_mask)
    temp = ma.array(sounding_inputs['temp']['data'],mask=dwpt_mask)
    pressure_array = ma.array(pressure_array,mask=dwpt_mask)

    dewpoint_AB_vec = np.vectorize(dewpoint_AB)
    dwpt = dewpoint_AB_vec(temp,rh) + 273.15

    sounding_inputs['dwpt']['data'] = dwpt
    sounding_inputs['dwpt']['mask'] = dwpt.mask if  ma.is_masked(dwpt) \
            else np.zeros(dwpt.shape,dtype='bool')
    sounding_inputs['dwpt']['units'] = 'K'
    sounding_inputs['dwpt']['long_name'] = 'Dew-point Temperature in K'

    LOG.debug("ma.is_masked({}) = {}".format('dwpt',ma.is_masked(dwpt)))
    LOG.debug("There are {} masked values in dwpt".format(np.sum(dwpt.mask)))
    LOG.debug("dwpt.shape =  {}".format(sounding_inputs['dwpt']['data'].shape))
    LOG.debug("dwpt_mask.shape =  {}".format(sounding_inputs['dwpt']['mask'].shape))

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
    scale         = plot_options['scale']
    map_res       = plot_options['map_res']
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

    # Compute the central lat and lon if they are not specified
    if (lat_0==None) and (lon_0==None):
        geo_shape = lat.shape
        nrows, ncols = geo_shape[0],geo_shape[1]
        LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
        # Non lat/lon pair given, use the central unmasked values.
        row_idx = int(nrows/2.)
        col_idx = int(ncols/2.)
        lat_0 = lat[row_idx,col_idx]
        lon_0 = lon[row_idx,col_idx]
        LOG.info("No lat/lon pair given, using ({:4.2f},{:4.2f})".
                format(lat_0,lon_0))

    if (latMin==None) and (latMax==None):
        LOG.info("Calculating lat extent...")
        latMin = np.min(lat)
        latMax = np.max(lat)
        LOG.info("Latitude extent: ({:4.2f},{:4.2f})".format(latMin,latMax))
    if (lonMin==None) and (lonMax==None):
        LOG.info("Calculating lon extent...")
        lonMin = np.min(lon)
        lonMax = np.max(lon)
        LOG.info("Longitude extent: ({:4.2f},{:4.2f})".format(lonMin,lonMax))

    # General Setup
    figWidth,figHeight = 8.,8.
    ax_rect = [0.05, 0.15, 0.90, 0.75  ] # [left,bottom,width,height]

    fig = Figure(figsize=(figWidth,figHeight))
    canvas = FigureCanvas(fig)

    ax = fig.add_axes(ax_rect)

    LOG.info("scale = {}".format(scale))
    windowWidth = scale *(0.35*12000000.)
    windowHeight = scale *(0.50*9000000.)
    m = Basemap(projection='lcc',lon_0=lon_0,lat_0=lat_0,
            width=windowWidth,height=windowHeight,
            #llcrnrlon=lonMin,
            #llcrnrlat=latMin,
            #urcrnrlon=lonMax,
            #urcrnrlat=latMax,
            ax=ax,fix_aspect=True,resolution=map_res)

    x,y=m(lon[::stride,::stride],lat[::stride,::stride])

    m.drawcoastlines()
    #m.drawstates()
    m.drawcountries()
    m.fillcontinents(color='0.85',zorder=0)
    m.drawparallels(np.arange( -90, 91,30), color = '0.25', 
            linewidth = 0.5,labels=[1,0,1,0])
    m.drawmeridians(np.arange(-180,180,30), color = '0.25', 
            linewidth = 0.5,labels=[1,0,1,0])

    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    plotMin = np.min(data) if plotMin==None else plotMin
    plotMax = np.max(data) if plotMax==None else plotMax
    LOG.debug("plotMin = {}".format(plotMin))
    LOG.debug("plotMax = {}".format(plotMax))

    if doScatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none',
                vmin=plotMin,vmax=plotMax)
    else:
        cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False,
                vmin=plotMin,vmax=plotMax)

    txt = ax.set_title(title,fontsize=11)

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    #ax.set_aspect('equal')

    cax_rect = [0.05 , 0.05, 0.90 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')

    txt = cax.set_title(cbar_title)

    #
    # Add a small globe with the swath indicated on it #
    #

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

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

    # Redraw the figure
    canvas.draw()
    LOG.info("Writing image file {}".format(pngName))
    canvas.print_figure(pngName,dpi=dpi)


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    dataChoices=['IAPP','MIRS','DR','NUCAPS']
    prodChoices=['temp','dwpt','relh']
    map_res_choice = ['c','l','i']

    defaults = {
                'input_file':None,
                'datatype':None,
                'dataset':'temp',
                'stride':1,
                'pressure':850.,
                'lat_0':None,
                'lat_0':None,
                'plotMin'  : None,
                'plotMax'  : None,
                'scatter_plot':False,
                'pointSize':1,
                'scale':1,
                'map_res':'c',
                'output_file':None,
                'outputFilePrefix' : None,
                'dpi':200
                }

    description = '''Create a contour plot of temperature, dewpoint or something 
                     else at a particular pressure level. Supports IAPP, MIRS, DR 
                     and NUCAPS files.'''

    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser(
                                     version=version,
                                     description=description
                                     )

    # Mandatory/positional arguments
    
    parser.add_argument(
                      action='store',
                      dest='input_file',
                      type=str,
                      help='''The fully qualified path to a single level-1d input file.'''
                      )

    parser.add_argument(
                      action="store",
                      dest="datatype",
                      default=defaults["datatype"],
                      type=str,
                      choices=dataChoices,
                      help='''The type of the input sounder data file.\n\n
                              Possible values are...
                              {}.
                           '''.format(dataChoices.__str__()[1:-1])
                      )

    # Optional arguments 

    parser.add_argument('--dset',
                      action="store",
                      dest="dataset",
                      default=defaults["dataset"],
                      type=str,
                      choices=prodChoices,
                      help='''The sounder dataset to plot.
                              Possible values are...
                              {}.
                              [default: {}]
                           '''.format(prodChoices.__str__()[1:-1],
                               defaults["dataset"])
                      )

    parser.add_argument('-S','--stride',
                      action="store",
                      dest="stride",
                      default=defaults["stride"],
                      type=int,
                      help='''Sample every STRIDE pixels in the band data. 
                      [default: {}]'''.format(defaults["stride"])
                      )

    parser.add_argument('--pressure',
                      action="store",
                      dest="pressure",
                      default=defaults["pressure"],
                      type=float,
                      help='''The pressure level (in mbar) to plot.. 
                      [default: {}]'''.format(defaults["pressure"])
                      )

    parser.add_argument('--plotMin',
                      action="store",
                      dest="plotMin",
                      default=defaults["plotMin"],
                      type=float,
                      help="Minimum value to plot.".format(defaults["plotMin"])
                      )

    parser.add_argument('--plotMax',
                      action="store",
                      dest="plotMax",
                      default=defaults["plotMax"],
                      type=float,
                      help="Maximum value to plot.".format(defaults["plotMax"])
                      )

    parser.add_argument('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default='200.',
                      type=float,
                      help='''The resolution in dots per inch of the output 
                      png file. 
                      [default: {}]'''.format(defaults["dpi"])
                      )

    parser.add_argument('--lat_0',
                      action="store",
                      dest="lat_0",
                      type=float,
                      help="Center latitude of plot."
                      )

    parser.add_argument('--lon_0',
                      action="store",
                      dest="lon_0",
                      type=float,
                      help="Center longitude of plot."
                      )

    parser.add_argument('--latMin',
                      action="store",
                      dest="latMin",
                      type=float,
                      help="Minimum latitude to plot."
                      )

    parser.add_argument('--latMax',
                      action="store",
                      dest="latMax",
                      type=float,
                      help="Maximum latitude to plot."
                      )

    parser.add_argument('--lonMin',
                      action="store",
                      dest="lonMin",
                      type=float,
                      help="Minimum longitude to plot."
                      )

    parser.add_argument('--lonMax',
                      action="store",
                      dest="lonMax",
                      type=float,
                      help="Maximum longitude to plot."
                      )

    parser.add_argument('--scatter_plot',
                      action="store_true",
                      dest="doScatterPlot",
                      default=defaults["scatter_plot"],
                      help="Generate the plot using a scatterplot approach."
                      )

    parser.add_argument('-P','--pointSize',
                      action="store",
                      dest="pointSize",
                      default=defaults["pointSize"],
                      type=float,
                      help='''Size of the plot point used to represent each pixel. 
                      [default: {}]'''.format(defaults["pointSize"])
                      )

    parser.add_argument('-s','--scale',
                      action="store",
                      dest="scale",
                      default=defaults["scale"],
                      type=float,
                      help='''The scaling factor for the default viewport size 
                      of (w x h) = (4200000.0 x 4500000.0) meters. 
                      [default: {}]'''.format(defaults["pointSize"])
                      )

    parser.add_argument('-m','--map_res',
                      action="store",
                      dest="map_res",
                      default=defaults["map_res"],
                      type=str,
                      choices=map_res_choice,
                      help="""The map coastline resolution. Possible values are 
                      'c' (coarse),'l' (low) and 'i' (intermediate). 
                      [default: '{}']""".format(defaults["map_res"])
                      )

    parser.add_argument('-o','--output_file',
                      action="store",
                      dest="output_file",
                      default=defaults["output_file"],
                      type=str,
                      help='''The filename of the output png file. 
                      '''
                      )

    parser.add_argument('-O','--output_file_prefix',
                      action="store",
                      dest="outputFilePrefix",
                      default=defaults["outputFilePrefix"],
                      type=str,
                      help="""String to prefix to the automatically generated 
                      png names, which are of the form 
                      <N_Collection_Short_Name>_<N_Granule_ID>.png. 
                      [default: {}]""".format(defaults["outputFilePrefix"])
                      )

    parser.add_argument("-z", "--verbose",
                      dest='verbosity',
                      action="count", 
                      default=0,
                      help='''each occurrence increases verbosity 1 level from 
                      ERROR: -z=WARNING -zz=INFO -zzz=DEBUG'''
                      )


    args = parser.parse_args()

    # Set up the logging
    console_logFormat = '%(asctime)s : (%(levelname)s):%(filename)s:%(funcName)s:%(lineno)d:  %(message)s'
    dateFormat = '%Y-%m-%d %H:%M:%S'
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[args.verbosity], 
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
    lonMax = options.lonMax
    plotMin = options.plotMin
    plotMax = options.plotMax
    scale = options.scale
    doScatterPlot = options.doScatterPlot
    pointSize = options.pointSize
    map_res = options.map_res
    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix
    dpi = options.dpi

    
    # Read in the input file, and return a dictionary containing the required
    # data

    dataChoices=['IAPP','MIRS','DR','NUCAPS']

    if datatype == 'IAPP' :
        sounding_inputs = iapp_sounder(input_file,pres_0=pressure)
    elif datatype == 'MIRS' :
        sounding_inputs = mirs_sounder(input_file,pres_0=pressure)
    elif datatype == 'DR' :
        sounding_inputs = dual_regression_sounder(input_file,pres_0=pressure)
    elif datatype == 'NUCAPS' :
        sounding_inputs = nucaps_sounder(input_file,pres_0=pressure)
    else:
        pass

    lats = sounding_inputs['lat']['data']
    lons = sounding_inputs['lon']['data']
    data = sounding_inputs[dataset]['data']
    data_mask = sounding_inputs[dataset]['mask']
    units = sounding_inputs[dataset]['units']

    pres_0 = int(sounding_inputs['pres_0'])

    input_file = path.basename(input_file)

    cbar_titles = {
            'temp':'temperature (K) @ {:4.2f} hPa'.format(pres_0),
            'dwpt':'dewpoint temperature (K) @ {:4.2f} hPa'.format(pres_0),
            'relh':'relative humidity (%) @ {:4.2f} hPa'.format(pres_0)
            }

    # Determine the filename
    file_suffix = "{}_{}_{}mb".format(datatype,dataset,pres_0)

    if output_file==None and outputFilePrefix==None :
        output_file = "{}.{}.png".format(input_file,file_suffix)
    if output_file!=None and outputFilePrefix==None :
        pass
    if output_file==None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)
    if output_file!=None and outputFilePrefix!=None :
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
    plot_options['plotMin'] = plotMin
    plot_options['plotMax'] = plotMax
    plot_options['scale'] = scale
    plot_options['map_res'] = map_res
    #plot_options['cmap'] = cmap
    plot_options['scatterPlot'] = doScatterPlot
    plot_options['pointSize'] = pointSize
    plot_options['dpi'] = dpi

    # Create the plot
    plotMapData(lats,lons,data,data_mask,output_file,**plot_options)

    return 0


if __name__=='__main__':
    sys.exit(main())  
