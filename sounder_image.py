#!/usr/bin/env python
# encoding: utf-8
"""
sounder_image.py

Purpose: Create a plot of temperature, dewpoint or relative humidity,
         at a particular pressure level. Supports HRPT, MIRS 
         and Dual Regression files.

Preconditions:
    * matplotlib (with basemap)
    * netCDF4 python module
    * h5py python module

Optional:
    * 

Minimum commandline:

    python sounder_image.py  --input_files=INPUTFILES --datatype=DATATYPE

where...

    INPUTFILES: The fully qualified path to the input files. May be a 
    directory or a file glob.

    DATATYPE: One of 'HRPT','MIRS' or 'DR'.


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

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'Geoff Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
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


def hrpt_sounder(hrpt_file,pres_0=850.):
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
    fileObj = Dataset(hrpt_file)

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

        LOG.info("Closing {}".format(hrpt_file))
        fileObj.close()

    except Exception, err :
        LOG.warn("There was a problem, closing {}".format(hrpt_file))
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
            LOG.debug(sounding_inputs[label]['data'])
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
    rh = mr_to_rh_vec(wvap,pressure_array,temp)

    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            else np.zeros(temp.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {} masked values in relh".format(np.sum(rh.mask)))

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
    rh = mr_to_rh_vec(wvap,pressure_array,temp)
    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            else np.zeros(rh.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {} masked values in relh".format(np.sum(rh.mask)))

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

    if doScatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none',
                vmin=plotMin,vmax=plotMax)
    else:
        cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False,
                vmin=plotMin,vmax=plotMax)

    txt = ax.set_title(title,fontsize=11)
    #if units == '%':
        #txt = ax.set_title('Reflectance')
    #elif units == 'K':
        #txt = ax.set_title('Brightness Temperature')
    #else:
        #txt = ax.set_title('scaled data')

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    #ax.set_aspect('equal')

    cax_rect = [0.05 , 0.05, 0.90 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')

    txt = cax.set_title(cbar_title)

    # Redraw the figure
    canvas.draw()
    LOG.info("Writing image file {}".format(pngName))
    canvas.print_figure(pngName,dpi=dpi)


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    dataChoices=['HRPT','MIRS','DR']
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
                'scatter_plot':False,
                'pointSize':1,
                'scale':1,
                'map_res':'c',
                'map_plot':False,
                'output_file':None,
                'dpi':200
                }

    description = '''Create a contour plot of temperature, dewpoint or something 
                     else at a particular pressure level. Supports HRPT and MIRS 
                     files.'''

    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser(
                                     version=version,
                                     description=description
                                     )

    # Mandatory arguments

    parser.add_argument('-i','--input_file',
                      action='store',
                      dest='input_file',
                      type=str,
                      required=True,
                      help="The fully qualified path to the input file."
                      )

    parser.add_argument('-t','--datatype',
                      action="store",
                      dest="datatype",
                      default=defaults["datatype"],
                      type=str,
                      required=True,
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
                           '''.format(prodChoices.__str__()[1:-1])
                      )

    parser.add_argument('--map_plot',
                      action="store_true",
                      dest="doMapPlots",
                      default=defaults["map_plot"],
                      help='''Plot the band arrays to png files, using a map 
                      projection (has no effect if --make_plots is not set). 
                      [default: {}]'''.format(defaults["map_plot"])
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
                      type=float,
                      help="Minimum value to plot."
                      )

    parser.add_argument('--plotMax',
                      action="store",
                      dest="plotMax",
                      type=float,
                      help="Maximum value to plot."
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
                      #default=defaults["output_file"],
                      type=str,
                      help='''The filename of the output png file. 
                      '''
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
    dpi = options.dpi

    
    # Read in the input file, and return a dictionary containing the required
    # data

    dataChoices=['HRPT','MIRS','DR']

    if datatype == 'HRPT' :
        sounding_inputs = hrpt_sounder(input_file,pres_0=pressure)
    elif datatype == 'MIRS' :
        sounding_inputs = mirs_sounder(input_file,pres_0=pressure)
    elif datatype == 'DR' :
        sounding_inputs = dual_regression_sounder(input_file,pres_0=pressure)
    else:
        pass

    #sys.exit(0)

    lats = sounding_inputs['lat']['data']
    lons = sounding_inputs['lon']['data']
    data = sounding_inputs[dataset]['data']
    data_mask = sounding_inputs[dataset]['mask']
    units = sounding_inputs[dataset]['units']

    if output_file == None:
        pres_0 = int(sounding_inputs['pres_0'])
        output_file = "{}_{}_{}mb.png".format(datatype,dataset,pres_0)
        LOG.info("output_file = {}".format(output_file))

    input_file = path.basename(input_file)

    cbar_titles = {
            'temp':'temperature (K) @ {:4.2f} hPa'.format(pres_0),
            'dwpt':'dewpoint temperature (K) @ {:4.2f} hPa'.format(pres_0),
            'relh':'relative humidity (%) @ {:4.2f} hPa'.format(pres_0)
            }

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
