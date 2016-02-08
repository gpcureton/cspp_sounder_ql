
#!/usr/bin/env python
# encoding: utf-8
"""
NUCAPS.py

Purpose: Provide read methods for various datasets associated with the 
NOAA Unique CrIS/ATMS Product System (NUCAPS) Package.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2015-04-07.
Copyright (c) 2012-2016 University of Wisconsin Regents. All rights reserved.

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

from netCDF4 import Dataset
from netCDF4 import num2date
import h5py

# every module should have a LOG object
LOG = logging.getLogger(__file__)

from ql_common import get_pressure_index


def sounder(nucaps_file_list,pres_0=850.):
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

    LOG.debug("Input file list: {}".format(nucaps_file_list))

    first_granule = True

    for grans in np.arange(len(nucaps_file_list)):

        try:

            LOG.debug("Ingesting granule {} ...".format(grans))

            nucaps_file = nucaps_file_list[grans]
            LOG.debug("Ingesting granule {} ...".format(nucaps_file))

            # Open the file object
            fileObj = Dataset(nucaps_file)

            # Open the dataset nodes for each of the required datasets
            sounding_inputs['pres']['node'] = fileObj.variables['Pressure']
            sounding_inputs['temp']['node'] = fileObj.variables['Temperature']
            sounding_inputs['wvap']['node'] = fileObj.variables['H2O_MR']
            sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
            sounding_inputs['lon']['node'] = fileObj.variables['Longitude']

            # Fill in some missing attributes
            sounding_inputs['pres']['units'] = 'hPa'
            sounding_inputs['temp']['units'] = 'K'
            sounding_inputs['wvap']['units'] = 'g/kg'

            if first_granule:

                # Get the pressure scale
                pressure = sounding_inputs['pres']['node'][0,:]

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
                        sounding_inputs[label]['data'],
                        sounding_inputs[label]['node'][:].reshape(4,30)))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:].reshape(4,30)
                    LOG.debug("first arrays...")

                LOG.debug("Intermediate {} shape = {}".format(
                    label,sounding_inputs[label]['data'].shape))

            # Add to the temp, dewpoint, water vapour and relative humidity data arrays...
            for label in ['temp','wvap']:
                try :
                    sounding_inputs[label]['data']  = np.vstack((
                        sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,level].reshape(4,30)))
                    LOG.debug("subsequent arrays...")
                except KeyError :
                    sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,level].reshape(4,30)
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
