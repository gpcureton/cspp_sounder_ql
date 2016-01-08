#!/usr/bin/env python
# encoding: utf-8
"""
IAPP.py

Purpose: A collection of routines for the reading and manipulation of data from
         the International ATOVS Processing Package (IAPP).

Preconditions:
    * netCDF4 python module
    * h5py python module


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2015-10-20.
Copyright (c) 2015 University of Wisconsin Regents. All rights reserved.

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

import scipy.spatial as ss

from netCDF4 import Dataset
from netCDF4 import num2date
import h5py

from thermo import dewhum

# every module should have a LOG object
LOG = logging.getLogger(__file__)

def sounder_image(iapp_file_list,pres_0=850.):
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
                    if sounding_inputs[label] != None:
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
            for label in ['ctp','ctt']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:,4]))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:,4]
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
        if sounding_inputs[label] != None:
            LOG.debug(">>> Processing {} ..."
                    .format(sounding_inputs[label]['dset_name']))
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

    dewhum_vec = np.vectorize(dewhum)
    _,rh,_ = dewhum_vec(pressure_array,temp,wvap)

    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            else np.zeros(temp.shape,dtype='bool')
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    return sounding_inputs


def sounder_skewT(iapp_file_list,lat_0=None,lon_0=None):
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

    LOG.debug("Input file list: {}".format(iapp_file_list))

    iapp_file = iapp_file_list[0]

    # Open the file object
    fileObj = Dataset(iapp_file)

    try:

        sounding_inputs['pres']['node'] = fileObj.variables['Pressure_Levels']
        sounding_inputs['temp']['node'] = fileObj.variables['Temperature_Retrieval']
        sounding_inputs['dwpt']['node'] = fileObj.variables['Dew_Point_Temp_Retrieval']
        sounding_inputs['wvap']['node'] = fileObj.variables['WaterVapor_Retrieval']
        sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
        sounding_inputs['lon']['node'] = fileObj.variables['Longitude']

        lat = sounding_inputs['lat']['node'][:,:]
        lon = sounding_inputs['lon']['node'][:,:]

        for label in data_labels:
            if sounding_inputs[label] != None:
                LOG.info(">>> Processing {} ...".format(sounding_inputs[label]['dset_name']))
                for attr_name in sounding_inputs[label]['node'].ncattrs(): 
                    attr = getattr(sounding_inputs[label]['node'],attr_name)
                    LOG.debug("{} = {}".format(attr_name,attr))
                    sounding_inputs[label][attr_name] = attr

        row,col = get_geo_indices(lat,lon,lat_0=lat_0,lon_0=lon_0)

        if row==None or col==None:
            raise Exception("No suitable lat/lon coordinates found, aborting...")
             
        LOG.debug("Retrieved row,col = {},{}".format(row,col))
        sounding_inputs['lat_0'] = lat[row,col]
        sounding_inputs['lon_0'] = lon[row,col]
        #row,col = 10,10

        sounding_inputs['pres']['data'] = sounding_inputs['pres']['node'][:]

        sounding_inputs['temp']['data'] = sounding_inputs['temp']['node'][row,col,:]
        sounding_inputs['dwpt']['data'] = sounding_inputs['dwpt']['node'][row,col,:]
        sounding_inputs['wvap']['data'] = sounding_inputs['wvap']['node'][row,col,:]

        LOG.info("Closing {}".format(iapp_file))
        fileObj.close()

    except Exception, err :
        LOG.warn("There was a problem, closing {}".format(iapp_file))
        LOG.warn("{}".format(err))
        LOG.debug(traceback.format_exc())
        fileObj.close()
        LOG.info("Exiting...")
        sys.exit(1)

    # Construct the masks of the various datasets, and mask the data accordingly
    for label in ['temp','dwpt','wvap']:
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
            sounding_inputs[label]['data'] = ma.array(data,mask=data_mask)

    # Computing the relative humidity
    sounding_inputs['relh'] = {}
    dewhum_vec = np.vectorize(dewhum)
    _,rh,_ = dewhum_vec(
            sounding_inputs['pres']['data'],
            sounding_inputs['temp']['data'],
            sounding_inputs['wvap']['data']
            )

    sounding_inputs['relh']['data'] = rh
    sounding_inputs['relh']['units'] = '%'
    sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    return sounding_inputs

