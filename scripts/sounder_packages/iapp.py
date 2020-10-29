#!/usr/bin/env python
# encoding: utf-8
"""
IAPP.py

Purpose: Provide read methods for various datasets associated with the
International ATOVS Processing Package (IAPP).

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

import os
import sys
import logging
import traceback
from os import path, uname, environ
import string
import re
import uuid
from shutil import rmtree, copyfile
from glob import glob
from time import time
from datetime import datetime, timedelta

import numpy as np
from numpy import ma
import copy

from ql_common import get_pressure_index, get_geo_indices, get_pressure_level
from ql_common import Datafile_NetCDF
from thermo import dewhum, rh_to_mr

from .sounder_packages import Sounder_Packages

# every module should have a LOG object
LOG = logging.getLogger(__file__)

# Dimension names
dimension_name = {}
dimension_name['nrows'] = ['Along_Track']
dimension_name['ncols'] = ['Across_Track']
dimension_name['nlevels'] = ['Pres_Levels']
dimension_name['nfov'] = ['Num_of_FOVs']

# Files dataset names for each dataset type
dset_name = {}
dset_name['pres'] = 'Pressure_Levels'
dset_name['pres_array'] = None
dset_name['lat']  = 'Latitude'
dset_name['lon']  = 'Longitude'
dset_name['ctp_fov']  = 'Cloud_Top_Pressure_CO2'
dset_name['ctp']  = None
dset_name['ctt_fov']  = 'Cloud_Top_Temperature_CO2'
dset_name['ctt']  = None
dset_name['temp'] = 'Temperature_Retrieval'
dset_name['2temp'] = None
dset_name['cold_air_aloft'] = None
dset_name['dwpt'] = 'Dew_Point_Temp_Retrieval'
dset_name['relh'] = None
dset_name['wvap'] = 'WaterVapor_Retrieval'

# Dataset types (surface/top; pressure level...)
dset_type = {}
dset_type['pres'] = 'profile'
dset_type['pres_array'] = 'level'
dset_type['lat']  = 'single'
dset_type['lon']  = 'single'
dset_type['ctp_fov']  = 'single'
dset_type['ctp']  = 'level'
dset_type['ctt_fov']  = 'single'
dset_type['ctt']  = 'level'
dset_type['temp'] = 'level'
dset_type['2temp'] = 'level'
dset_type['cold_air_aloft'] = 'level'
dset_type['dwpt'] = 'level'
dset_type['relh'] = 'level'
dset_type['wvap'] = 'level'

# Dataset dependencies for each dataset

dset_deps = {}
dset_deps['pres'] = []
dset_deps['pres_array'] = []
dset_deps['lat']  = []
dset_deps['lon']  = []
dset_deps['ctp_fov']  = []
dset_deps['ctp']  = ['ctp_fov']
dset_deps['ctt_fov']  = []
dset_deps['ctt']  = ['ctt_fov']
dset_deps['temp'] = []
dset_deps['2temp'] = ['temp']
dset_deps['cold_air_aloft'] = ['temp', 'wvap']
dset_deps['dwpt'] = []
dset_deps['wvap'] = []
dset_deps['relh'] = ['pres_array', 'temp','wvap']

# The class method used to read/calculate each dataset.
dset_method = {}
dset_method['pres'] = []
dset_method['lat']  = []
dset_method['lon']  = []
dset_method['ctp_fov']  = []
dset_method['ctp']  = []
dset_method['ctt_fov']  = []
dset_method['ctt']  = []
dset_method['temp'] = []
dset_method['2temp'] = []
#dset_method['cold_air_aloft'] = []
dset_method['dwpt'] = []
dset_method['relh'] = []
dset_method['wvap'] = []


class Iapp(Sounder_Packages):
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

    def __init__(self, *args, **kwargs):
        Sounder_Packages.__init__(self, *args, **kwargs)

        file_list, data_name, plot_type = args

        self.pkg_name = "International TOVS Processing Package (IAPP)"
        pres_0 = kwargs['pres_0']
        # elev_0 = kwargs['elev_0']
        lat_0 = kwargs['lat_0']
        lon_0 = kwargs['lon_0']
        # footprint = kwargs['footprint']

        self.file_list = file_list
        self.plot_type = plot_type
        self.have_geo = False
        self.dimensions = {}
        self.datasets = {}

        # Construct a list of the required datasets...
        dsets_to_read = [data_name] + self.get_dset_deps(data_name)
        LOG.info("Initial reverse list of datasets to read : {}".format(dsets_to_read))
        dsets_to_read.reverse()
        LOG.info("Dupe Datasets to read : {}".format(dsets_to_read))

        self.dsets_to_read = []
        for dset in dsets_to_read:
            if dset not in self.dsets_to_read:
                self.dsets_to_read.append(dset)
            else:
                pass

        self.dsets_to_read = ['lat','lon'] + self.dsets_to_read
        LOG.info("Final Datasets to read : {}".format(self.dsets_to_read))

        # Define any special methods for computing the required dataset
        dset_method['cold_air_aloft'] = self.cold_air_aloft
        dset_method['2temp'] = self.temp_times_two
        dset_method['ctt'] = self.cloud_top_temperature_fovavg
        dset_method['ctp'] = self.cloud_top_pressure_fovavg
        dset_method['relh'] = self.relative_humidity
        dset_method['pres_array'] = self.pressure_array

        # Read in the pressure dataset, and get the dimensions...
        LOG.info(">>> Reading in the pressure dataset...")
        dfile_obj = Datafile_NetCDF(file_list[0])
        data_obj = dfile_obj.Dataset(dfile_obj,dset_name['pres'])
        self.pressure = data_obj.dset[:]
        LOG.debug(f"self.pressure.shape = {self.pressure.shape}")

        for key in dimension_name:
            found_dim = False
            for key_var in dimension_name[key]:
                try:
                    LOG.debug(f"Looking for '{key_var}' as '{key}'...")
                    self.dimensions[key] = dfile_obj.dimensions[key_var]
                    LOG.debug(f"self.dimensions['{key}'] = {self.dimensions[key]}")
                    found_dim = True
                    break
                except Exception as err :
                    LOG.warn(f"\tThere was a problem reading dimensions '{key}'")
                    LOG.warn("\t{}".format(err))
                    LOG.debug(traceback.format_exc())
            if not found_dim:
                LOG.error(f"Failed to find dimension match for '{key}' for {os.path.basename(file_list[0])}, is this a valid {type(self).__name__.upper()} file?")
                dfile_obj.close()
                sys.exit(1)
            LOG.debug("dimension {} = {}".format(key,self.dimensions[key]))

        dfile_obj.close()

        # Determine the elevation dataset...
        # self.elev_0 = elev_0
        # if elev_0 != None:


            # # Create a cube of the temp_gdas and relh_gdas datasets...
            # LOG.info("Constructing a cube dataset...")
            # self.construct_cube_pass(file_list,['temp_gdas','relh_gdas'])

            # # Construct the pressure cube
            # #cube_shape = self.datasets['temp_gdas']['data'].shape
            # #nlevels,nrows,ncols = cube_shape[0],cube_shape[1],cube_shape[2]
            # (nlevels,nrows,ncols) = self.datasets['temp']['data'].shape
            # self.datasets['pressure'] = {}
            # self.datasets['pressure']['data'] = np.broadcast_to(np.broadcast_to(self.pressure,(nrows,nlevels)),(ncols,nrows,nlevels)).T
            # self.datasets['pressure']['attrs'] = {}
            # LOG.info("pressure cube shape = {}".format(self.datasets['pressure']['data'].shape))

            # # Compute the wvap_gdas dataset
            # self.datasets['wvap_gdas'] = {}
            # self.datasets['wvap_gdas']['attrs'] = {}
            # self.datasets['wvap_gdas']['data'] = np.array([1]) # placeholder
            # level = -4
            # press = self.datasets['pressure']['data']
            # relh = ma.array(self.datasets['relh_gdas']['data'],mask=self.datasets['relh_gdas']['data_mask'])
            # temp = ma.array(self.datasets['temp_gdas']['data'],mask=self.datasets['temp_gdas']['data_mask'])
            # # Compute the elevation cube
            # self.elevation = get_elevation(press,temp,relh)

            # # Get the cube indicies that correspond to the value of match_val...
            # match_val = elev_0
            # cube_idx = get_level_indices(self.elevation,match_val)
            # self.cube_idx = cube_idx

            # elev_level = -9999. * np.ones(self.elevation[0].shape,dtype='float')
            # elev_level[(cube_idx[1],cube_idx[2])] = self.elevation[cube_idx]
            # LOG.info("\nelev_level = \n{}".format(elev_level))


        if plot_type == 'image':

            LOG.info("Preparing a 'level' plot...")

            # if elev_0 != None:
                # # Contruct a pass of the required datasets at the desired pressure.
                # self.construct_level_pass(file_list,None,cube_idx,None,None)
            # else:
            if True:
                # Determine the level closest to the required pressure
                self.level,self.pres_0 = get_pressure_level(self.pressure,pres_0)
                # Contruct a pass of the required datasets at the desired pressure.
                self.construct_level_pass(file_list,self.level,None,None,None)

            LOG.debug("\n\n>>> Intermediate dataset manifest...\n")
            LOG.debug(self.print_dataset_manifest(self))

        # elif plot_type == 'slice':

            # if dset_type[data_name] == 'single':
                # LOG.warn("Cannot plot a slice of the single level dataset '{}', aborting...".format(data_name))
                # sys.exit(1)


            # LOG.info("Preparing a 'slice' plot...")
            # self.level,self.pres_0 = None,0

            # total_dsets_to_read = list(self.dsets_to_read)
            # self.dsets_to_read = ['lat','lon']

            # # Getting the full latitude and longitude
            # LOG.info(">>> Reading in the lat and lon arrays...")
            # self.construct_level_pass(file_list,None,None,None,None)

            # LOG.debug("\n\n>>> Intermediate dataset manifest...\n")
            # LOG.debug(self.print_dataset_manifest(self))

            # LOG.info("Computing the row/col and lon_0/lat_0 from the lon/lat pass...")
            # self.row,self.col,self.lat_0,self.lon_0 \
                    # = get_geo_indices(self.datasets['lat']['data'],
                                      # self.datasets['lon']['data'],
                                      # lat_0=None,lon_0=None)
            # if footprint != None:
                # self.col = footprint

            # for dsets in ['lat_row','lat_col','lon_row','lon_col']:
                # self.datasets[dsets] = {}


            # self.datasets['lat_row']['data'] = self.datasets['lat']['data'][self.row,:]
            # self.datasets['lat_col']['data'] = self.datasets['lat']['data'][:,self.col]
            # self.datasets['lat_row']['attrs'] =  dict(self.datasets['lat']['attrs'])
            # self.datasets['lat_col']['attrs'] =  dict(self.datasets['lat']['attrs'])

            # self.datasets['lon_row']['data'] = self.datasets['lon']['data'][self.row,:]
            # self.datasets['lon_col']['data'] = self.datasets['lon']['data'][:,self.col]
            # self.datasets['lon_row']['attrs'] =  dict(self.datasets['lon']['attrs'])
            # self.datasets['lon_col']['attrs'] =  dict(self.datasets['lon']['attrs'])

            # if elev_0 != None:
                # self.datasets['elev_slice'] = {}
                # self.datasets['elev_slice']['data'] = self.elevation[:,:,self.col]
                # self.datasets['elev_slice']['attrs'] = {}

            # self.dsets_to_read = list(total_dsets_to_read)
            # self.dsets_to_read.remove('lat')
            # self.dsets_to_read.remove('lon')

            # LOG.info(">>> Reading in the remaining dataset slices..")
            # slice_along_col = True
            # if slice_along_col:
                # LOG.debug("along column: (row,col,level)  = ({}, {}, {})".format(None,self.col,None))
                # # This is to slice along a column (across the rows)
                # self.construct_slice_pass(file_list, None, None, self.col)
            # else:
                # LOG.debug("along row: (row,col,level)  = ({}, {}, {})".format(self.row,None,None))
                # # This is to slice along a row (across the columns)
                # self.construct_slice_pass(file_list, None, self.row, None)

            # self.dsets_to_read = list(total_dsets_to_read)

        # Rationalise the satellite name...
        satellite_names = {'noaa15': u'NOAA-15',
                           'noaa16': u'NOAA-16',
                           'noaa18': u'NOAA-18',
                           'noaa19': u'NOAA-19',
                           'metopa': u'Metop-A',
                           'metopb': u'Metop-B'}
        for filename in self.datasets['file_attrs'].keys():
            satellite_name = self.datasets['file_attrs'][filename]['Satellite_Name']
            self.datasets['file_attrs'][filename]['Satellite_Name'] = satellite_names[satellite_name]

        LOG.debug("\n\n>>> Final dataset manifest...\n")
        LOG.debug(self.print_dataset_manifest(self))


    def construct_cube_pass(self,file_list,dsets_to_read=None):
        '''
        Read in each of the desired datasets. Each of the desired datasets may
        be either a "base" dataset, which is read directly from the file, or a
        "derived" dataset, which is constructed from previously read datasets.

        For each granule, all base and derived datasets should be completed, to
        ensure that any dataset interdependencies to not fail due to differing array
        shapes.
        '''

        LOG.info("Contructing a CUBE pass...")
        #LOG.info("Input file list: {}".format(file_list))

        first_granule = True
        self.datasets['file_attrs'] = {}

        if dsets_to_read == None:
            dsets_to_read = self.dsets_to_read

        for dset in dsets_to_read:
            self.datasets[dset] = {}

        this_granule_data = {}
        this_granule_mask = {}

        level = slice(None)
        row = slice(None)
        col = slice(None)

        # Loop through each of the granules...
        for grans in np.arange(len(file_list)):

            file_name = file_list[grans]

            try:

                LOG.info("")
                LOG.info("\tReading in granule {} : {}".format(grans,os.path.basename(file_name)))

                #LOG.info("\tOpening file {}".format(file_name))
                dfile_obj = Datafile_NetCDF(file_name)
                self.datasets['file_attrs'][os.path.basename(file_name)] \
                        = dfile_obj.attrs
                self.datasets['file_attrs'][os.path.basename(file_name)]['dt_obj'] \
                        = self.get_granule_dt(file_name)

                # Loop through each of the desired datasets
                for dset in dsets_to_read:
                    LOG.info("")
                    LOG.info("\t\tFor dataset {}".format(dset))

                    # Choose the correct "get" method for the dataset
                    if dset_method[dset] == []:
                        #LOG.info("\tUsing read method for {}".format(dset))
                        get_data = self.get_data
                    else:
                        #LOG.info("\tUsing compute method for {}".format(dset))
                        get_data = dset_method[dset]

                    #LOG.info("\t\tReading in granule {} of {}".format(grans,dset))

                    data = get_data(dfile_obj,dset,level,row,col,
                            this_granule_data,this_granule_mask)
                    this_granule_data[dset] = data

                    LOG.info("\t\tthis_granule_data['{}'].shape = {}".format(
                        dset,this_granule_data[dset].shape))

                    # Determine the missing/fill value
                    missing_value = None
                    for fill_val_key in dfile_obj.fill_val_keys:
                        if fill_val_key in  self.datasets[dset]['attrs'].keys():
                            missing_value = float(self.datasets[dset]['attrs'][fill_val_key])
                            LOG.debug("Setting missing_value to {} from {} dataset attributes"
                                    .format(missing_value, dset))
                    if missing_value == None:
                        missing_value = dfile_obj.fill_value
                    LOG.info("\t\tMissing value = {}".format(missing_value))

                    data_mask = ma.masked_equal(data,missing_value).mask
                    LOG.info("\t\tdata_mask.shape = {}".format(data_mask.shape))
                    if data_mask.shape == ():
                        data_mask = np.zeros(data.shape,dtype='bool')
                    this_granule_mask[dset] = data_mask

                    try :
                        self.datasets[dset]['data'] = \
                                np.hstack((self.datasets[dset]['data'],this_granule_data[dset]))
                        self.datasets[dset]['data_mask'] = \
                                np.hstack((self.datasets[dset]['data_mask'],this_granule_mask[dset]))

                        LOG.info("\t\tsubsequent arrays...")

                    except KeyError as err :

                        LOG.info("\t\tFirst arrays...")
                        LOG.info("\t\tCreating new data array for {}".format(dset))
                        #LOG.info(traceback.format_exc())

                        self.datasets[dset]['data'] = this_granule_data[dset]
                        self.datasets[dset]['data_mask'] = this_granule_mask[dset]


                    LOG.info("\t\tIntermediate {} shape = {}".format(
                        dset,self.datasets[dset]['data'].shape))
                    LOG.info("\tIntermediate {} mask shape = {}".format(
                        dset,self.datasets[dset]['data_mask'].shape))

                #LOG.info("\tClosing file {}".format(file_name))
                dfile_obj.close()

            except Exception as err :
                LOG.warn("\tThere was a problem, closing {}".format(file_name))
                LOG.warn("\t{}".format(err))
                LOG.debug(traceback.format_exc())
                LOG.info("\tClosing file {}".format(file_name))
                dfile_obj.close()
                LOG.info("\tExiting...")
                sys.exit(1)


    def construct_slice_pass(self,file_list,level,row,col):
        '''
        Read in each of the desired datasets. Each of the desired datasets may
        be either a "base" dataset, which is read directly from the file, or a
        "derived" dataset, which is constructed from previously read datasets.

        For each granule, all base and derived datasets should be completed, to
        ensure that any dataset interdependencies to not fail due to differing array
        shapes.
        '''

        LOG.info("\tContructing a SLICE pass...")
        #LOG.debug("\tInput file list: {}".format(file_list))
        LOG.debug("\t(row,col,level)  = ({}, {}, {})".format(row,col,level))

        first_granule = True
        self.datasets['file_attrs'] = {}
        for dset in self.dsets_to_read:
            self.datasets[dset] = {}

        this_granule_data = {}
        this_granule_mask = {}

        # Loop through each of the granules...
        for grans in np.arange(len(file_list)):

            file_name = file_list[grans]

            try:

                LOG.info("")
                LOG.info("\tReading in granule {} : {}".format(grans,os.path.basename(file_name)))

                #LOG.info("\tOpening file {}".format(file_name))
                dfile_obj = Datafile_NetCDF(file_name)
                self.datasets['file_attrs'][os.path.basename(file_name)] \
                        = dfile_obj.attrs
                self.datasets['file_attrs'][os.path.basename(file_name)]['dt_obj'] \
                        = self.get_granule_dt(file_name)

                # Loop through each of the desired datasets
                for dset in self.dsets_to_read:
                    LOG.info("")
                    LOG.info("\t\tFor dataset {}".format(dset))

                    # Choose the correct "get" method for the dataset
                    if dset_method[dset] == []:
                        LOG.info("\tUsing read method for {}".format(dset))
                        get_data = self.get_data
                    else:
                        LOG.info("\tUsing compute method for {}".format(dset))
                        get_data = dset_method[dset]

                    LOG.info("\tReading in granule {} of {}".format(grans,dset))

                    data = get_data(dfile_obj,dset,level,row,col,this_granule_data,this_granule_mask)

                    # Transpose to make the "level" dimension first, like HSRTV, for slice datasets
                    #this_granule_data[dset] = data.transpose()
                    this_granule_data[dset] = data[:]

                    LOG.info("\t\tthis_granule_data['{}'].shape = {}".format(
                        dset,this_granule_data[dset].shape))

                    # Determine the missing/fill value
                    missing_value = None
                    for fill_val_key in dfile_obj.fill_val_keys:
                        if fill_val_key in  self.datasets[dset]['attrs'].keys():
                            missing_value = float(self.datasets[dset]['attrs'][fill_val_key])
                            LOG.debug("Setting missing_value to {} from {} dataset attributes"
                                    .format(missing_value, dset))
                    if missing_value == None:
                        missing_value = dfile_obj.fill_value
                    LOG.info("\t\tMissing value = {}".format(missing_value))

                    data_mask = ma.masked_equal(data,missing_value).mask
                    LOG.info("\t\tdata_mask.shape = {}".format(data_mask.shape))
                    if data_mask.shape == ():
                        data_mask = np.zeros(data.shape,dtype='bool')
                    this_granule_mask[dset] = data_mask

                    try :
                        self.datasets[dset]['data'] = \
                                np.hstack((self.datasets[dset]['data'],this_granule_data[dset]))
                        self.datasets[dset]['data_mask'] = \
                                np.hstack((self.datasets[dset]['data_mask'],this_granule_mask[dset]))

                        LOG.info("\t\tsubsequent arrays...")

                    except KeyError as err :

                        LOG.info("\t\tFirst arrays...")
                        LOG.info("\t\tCreating new data array for {}".format(dset))
                        #LOG.info(traceback.format_exc())

                        self.datasets[dset]['data'] = this_granule_data[dset]
                        self.datasets[dset]['data_mask'] = this_granule_mask[dset]


                    LOG.info("\t\tIntermediate {} shape = {}".format(
                        dset,self.datasets[dset]['data'].shape))
                    #LOG.info("\tIntermediate {} mask shape = {}".format(
                        #dset,self.datasets[dset]['data_mask'].shape))

                #LOG.info("\tClosing file {}".format(file_name))
                dfile_obj.close()

            except Exception as err :
                LOG.warn("\tThere was a problem, closing {}".format(file_name))
                LOG.warn("\t{}".format(err))
                LOG.debug(traceback.format_exc())
                LOG.info("\tClosing file {}".format(file_name))
                dfile_obj.close()
                LOG.info("\tExiting...")
                sys.exit(1)

    def construct_level_pass(self,file_list,level,elev_idx,row,col):
        '''
        Read in each of the desired datasets. Each of the desired datasets may
        be either a "base" dataset, which is read directly from the file, or a
        "derived" dataset, which is constructed from previously read datasets.

        For each granule, all base and derived datasets should be completed, to
        ensure that any dataset interdependencies to not fail due to differing array
        shapes.
        '''

        LOG.info("Contructing a LEVEL pass...")
        #LOG.info("Input file list: {}".format(file_list))
        LOG.info("(level,row,col)  = ({}, {}, {})".format(level,row,col))

        first_granule = True
        self.datasets['file_attrs'] = {}
        for dset in self.dsets_to_read:
            self.datasets[dset] = {}

        this_granule_data = {}
        this_granule_mask = {}

        # Loop through each of the granules...
        for grans in np.arange(len(file_list)):

            file_name = file_list[grans]

            try:

                LOG.info("")
                LOG.info("\tReading in granule {} : {}".format(grans,os.path.basename(file_name)))

                #LOG.info("\tOpening file {}".format(file_name))
                dfile_obj = Datafile_NetCDF(file_name)
                self.datasets['file_attrs'][os.path.basename(file_name)] \
                        = dfile_obj.attrs
                self.datasets['file_attrs'][os.path.basename(file_name)]['dt_obj'] \
                        = self.get_granule_dt(file_name)

                # Loop through each of the desired datasets
                for dset in self.dsets_to_read:
                    LOG.info("")
                    LOG.info("\t\tFor dataset {}".format(dset))

                    # Choose the correct "get" method for the dataset
                    if dset_method[dset] == []:
                        LOG.info("\tUsing read method for {}".format(dset))
                        get_data = self.get_data
                    else:
                        LOG.info("\tUsing compute method for {}".format(dset))
                        get_data = dset_method[dset]

                    #LOG.info("\t\tReading in granule {} of {}".format(grans,dset))

                    data = get_data(dfile_obj,dset,level,row,col,
                            this_granule_data,this_granule_mask)
                    this_granule_data[dset] = data

                    LOG.info("\t\tthis_granule_data['{}'].shape = {}".format(
                        dset,this_granule_data[dset].shape))

                    # Determine the missing/fill value
                    missing_value = None
                    for fill_val_key in dfile_obj.fill_val_keys:
                        if fill_val_key in  self.datasets[dset]['attrs'].keys():
                            missing_value = float(self.datasets[dset]['attrs'][fill_val_key])
                            LOG.debug("Setting missing_value to {} from {} dataset attributes"
                                    .format(missing_value, dset))
                    if missing_value == None:
                        missing_value = dfile_obj.fill_value
                    LOG.info("\t\tMissing value = {}".format(missing_value))

                    data_mask = ma.masked_equal(data,missing_value).mask
                    #LOG.info("\t\tdata_mask.shape = {}".format(data_mask.shape))
                    if data_mask.shape == ():
                        data_mask = np.zeros(data.shape,dtype='bool')
                    this_granule_mask[dset] = data_mask

                    try :
                        self.datasets[dset]['data'] = \
                                np.vstack((self.datasets[dset]['data'],this_granule_data[dset]))
                        self.datasets[dset]['data_mask'] = \
                                np.vstack((self.datasets[dset]['data_mask'],this_granule_mask[dset]))

                        LOG.info("\t\tsubsequent arrays...")

                    except KeyError as err :

                        LOG.info("\t\tFirst arrays...")
                        LOG.info("\t\tCreating new data array for {}".format(dset))
                        #LOG.info(traceback.format_exc())

                        self.datasets[dset]['data'] = this_granule_data[dset]
                        self.datasets[dset]['data_mask'] = this_granule_mask[dset]


                    LOG.info("\t\tIntermediate {} shape = {}".format(
                        dset,self.datasets[dset]['data'].shape))
                    #LOG.info("\tIntermediate {} mask shape = {}".format(
                        #dset,self.datasets[dset]['data_mask'].shape))

                #LOG.info("\tClosing file {}".format(file_name))
                dfile_obj.close()

            except Exception as err :
                LOG.warn("\tThere was a problem, closing {}".format(file_name))
                LOG.warn("\t{}".format(err))
                LOG.debug(traceback.format_exc())
                LOG.info("\tClosing file {}".format(file_name))
                dfile_obj.close()
                LOG.info("\tExiting...")
                sys.exit(1)


    def get_data(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        This method reads a single granule of the required dataset, and returns
        a single array (which may be a surface/top, pressure level, or slice dataset.
        '''

        LOG.info("\t\tRetrieving {}".format(data_name))

        # The dataset dimensions from the file dimension block.
        nrows = self.dimensions['nrows']
        ncols = self.dimensions['ncols']
        nlevels = self.dimensions['nlevels']
        LOG.info("\t\t(nrows,ncols,nlevels)  = ({}, {}, {})".format(nrows,ncols,nlevels))

        LOG.info("\t\t(row,col,level)  = ({}, {}, {})".format(row,col,level))

        level = slice(level) if level == None else level
        row = slice(row) if row == None else row
        col = slice(col) if col == None else col

        LOG.info("\t\t(row,col,level) slice objects = ({}, {}, {})".format(row,col,level))

        LOG.info("\t\tplot type is {}".format(self.plot_type))

        # Get a reference to the dataset object...
        data_obj = dfile_obj.Dataset(dfile_obj,dset_name[data_name])
        LOG.info("\t\tdata_obj.dset.shape = {}".format(data_obj.dset.shape))

        # Determine whether this is a single pressure level (i.e.: ctp) or a profile dataset, and
        # get the correct array slice of the dataset array...
        if dset_type[data_name] == 'single':
            LOG.info("\t\tgetting a 'single' dataset...")
            dset = data_obj.dset[row,col].squeeze()
        elif dset_type[data_name] == 'level':
            LOG.info("\t\tgetting a 'level' dataset...")
            # Determine whether this is an level or slice plot, and get the correct array slice of
            # the data cube...
            if self.plot_type == 'image':
                LOG.info("\t\t> getting a pressure level...")
                dset = data_obj.dset[row,col,level].squeeze()
            else:
                LOG.info("\t\t> getting a pressure slice...")
                dset = data_obj.dset[row,col,level].squeeze().T

        LOG.info("\t\tDataset {} has shape {}".format(data_name,dset.shape))

        self.datasets[data_name]['attrs'] =  dict(data_obj.attrs)

        return dset

    def pressure_array(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to return an array of the pressure, whether it be at a single
        pressure level, or a vertical slice.
        '''

        LOG.info("\t\tComputing {}".format(data_name))

        # The dataset dimensions from the file dimension block.
        nrows = self.dimensions['nrows']
        ncols = self.dimensions['ncols']
        nlevels = self.dimensions['nlevels']
        LOG.info("\t\t(nrows,ncols,nlevels)  = ({}, {}, {})".format(nrows,ncols,nlevels))

        LOG.info("\t\t(row,col,level)  = ({}, {}, {})".format(row,col,level))

        level = slice(level) if level == None else level
        row = slice(row) if row == None else row
        col = slice(col) if col == None else col

        LOG.info("\t\t(row,col,level) slice objects = ({}, {}, {})".format(row,col,level))

        LOG.info("\t\ttemp has shape {}".format(this_granule_data['temp'].shape))

        LOG.info("\t\tplot type is {}".format(self.plot_type))

        # Create a pressure cube of dimensions (nrows,ncols,nlevels)
        self.pressure_array = np.broadcast_to(np.broadcast_to(self.pressure,(ncols,nlevels)),(nrows,ncols,nlevels))
        LOG.info("\t\tpressure_array.shape = {}".format(self.pressure_array.shape))

        # Determine whether this is an level or slice plot, and get the correct array slice of the
        # pressure cube...
        if self.plot_type == 'image':
            LOG.info("\t\t> getting a pressure level...")
            dset = self.pressure_array[row,col,level].squeeze()
        else:
            LOG.info("\t\t> getting a pressure slice...")
            dset = self.pressure_array[row,col,level].squeeze().T

        LOG.info("\t\tdset.shape = {}".format(dset.shape))

        dset_mask = np.zeros(dset.shape,dtype='bool')

        dset = ma.array(dset,mask=dset_mask)

        self.datasets[data_name]['attrs'] = self.datasets['lat']['attrs']

        return dset

    def cold_air_aloft(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to return the temperature, binned into three categories:

        temp < -65 degC             --> 0
        -65 degC < temp < -60 degC  --> 1
        temp > -60 degC             --> 2

        The resulting product is known as the "cold-air-aloft" and is of interest
        to aviation.
        '''

        LOG.info("\t\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        level = slice(level) if level == None else level
        row = slice(row) if row == None else row
        col = slice(col) if col == None else col

        LOG.info("\t\tComputing {}".format(data_name))
        dset_mask = this_granule_mask['temp']

        LOG.debug("\t\tTemp has shape {}".format(this_granule_data['temp'].shape))

        dset = this_granule_data['temp'][:,:].squeeze() - 273.16

        low_idx = np.where(dset < -65.)
        mid_idx = np.where((dset > -65.) * (dset < -60.))
        hi_idx = np.where(dset > -60.)
        dset[low_idx] = 0.
        dset[mid_idx] = 1.
        dset[hi_idx] = 2.

        dset = ma.array(dset,mask=dset_mask)

        self.datasets[data_name]['attrs'] = self.datasets['temp']['attrs']

        return dset

    def temp_times_two(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to compute two times the temperature. This is not used for anything, it's just
        to test making custom datasets.
        '''
        LOG.info("\t\tComputing {}".format(data_name))
        dset = 2. * self.datasets['temp']['data'] + 0.01*self.pres_0

        self.datasets[data_name]['attrs'] =  self.datasets['temp']['attrs']

        return dset

    def relative_humidity(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to compute the relative humidity.
        '''
        LOG.info("\t\tComputing {}".format(data_name))
        wvap = self.datasets['wvap']['data']
        temp = self.datasets['temp']['data']
        pres = self.datasets['pres_array']['data']

        dewhum_vec = np.vectorize(dewhum)
        _,rh,_ = dewhum_vec(pres,temp,wvap)

        dset = rh
        self.datasets[data_name]['attrs'] =  self.datasets['temp']['attrs']

        return dset

    def cloud_top_temperature_fovavg(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to compute the cloud top temperature by averaging over the valid FOV for each
        FOR.
        '''
        LOG.info("\t\tComputing {}".format(data_name))
        nrows,ncols,nfov = self.datasets['ctt_fov']['data'].shape
        dset = -999. * np.ones((nrows,ncols))
        for row in range(nrows):
            for col in range(ncols):
                fovarray = self.datasets['ctt_fov']['data'][row,col,:]
                fovarray = ma.masked_equal(fovarray,np.nan)
                fovarray = ma.masked_less(fovarray,153.)
                fovarray = ma.masked_equal(fovarray,-999.0)
                if fovarray.mask.sum() == len(fovarray.mask):
                    pass
                else:
                    dset[row,col] = np.mean(fovarray)

        dset = ma.masked_equal(dset,-999.0)
        LOG.info("\t\tCTT: {},{}".format(np.min(dset),np.max(dset)))

        self.datasets[data_name]['attrs'] =  self.datasets['ctt_fov']['attrs']

        return dset

    def cloud_top_pressure_fovavg(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to compute the cloud top pressure by averaging over the valid FOV for each
        FOR.
        '''
        LOG.info("\t\tComputing {}".format(data_name))
        nrows,ncols,nfov = self.datasets['ctp_fov']['data'].shape
        dset = -999. * np.ones((nrows,ncols))
        for row in range(nrows):
            for col in range(ncols):
                fovarray = self.datasets['ctp_fov']['data'][row,col,:]
                fovarray = ma.masked_equal(fovarray,np.nan)
                fovarray = ma.masked_equal(fovarray,0.)
                fovarray = ma.masked_equal(fovarray,-999.0)
                if fovarray.mask.sum() == len(fovarray.mask):
                    pass
                else:
                    dset[row,col] = np.mean(fovarray)

        dset = ma.masked_equal(dset,-999.0)
        LOG.info("\t\tCTP: {},{}".format(np.min(dset),np.max(dset)))

        self.datasets[data_name]['attrs'] =  self.datasets['ctp_fov']['attrs']

        return dset

    def print_dataset_manifest(self,pkg_obj):
        dataset_names = list(pkg_obj.datasets.keys())
        dataset_names.remove('file_attrs')
        LOG.debug("datasets: {}".format(dataset_names))

        for key in dataset_names:
            LOG.debug("For dataset {}:".format(key))
            LOG.debug("\tdatasets['{}']['attrs'] = {}".format(
                    key,pkg_obj.datasets[key]['attrs']))
            LOG.debug("\tdatasets['{}']['data'].shape = {}".format(
                    key,pkg_obj.datasets[key]['data'].shape))
        print("")

    def get_dset_deps(self,dset):
        '''
        For a particular dataset, returns the prerequisite datasets. Works
        recursively to find multiple levels of dependencies.
        Note: Currently there is nothing in place to detect circular dependencies,
        so be careful of how you construct the dependencies.
        '''
        deps = dset_deps[dset]
        for dset in deps:
            deps = deps + self.get_dset_deps(dset)
        return deps


    def get_granule_dt(self,file_name):
        '''
        Computes a datetime object based on the filename.
        This assumes a filename for HSRTV output like...
            "metopa_L2_d20150330_t1510413_e1516077_c20150330152548432117_iapp.nc"
        '''
        file_name = os.path.basename(file_name)
        image_date_time = "_".join(file_name.split("_")[2:4])
        image_date_time = image_date_time[:-1]
        dt_image_date = datetime.strptime(image_date_time,'d%Y%m%d_t%H%M%S')
        return dt_image_date


def gphite( press, temp, wvap, z_sfc, n_levels, i_dir):

    '''
    version of 18.05.00

    PURPOSE:

     Routine to compute geopotential height given the atmospheric state.
       Includes virtual temperature adjustment.

    CREATED:

     19-Sep-1996 Received from Hal Woolf, recoded by Paul van Delst
     18-May-2000 Logic error related to z_sfc corrected by Hal Woolf

     ARGUMENTS:

        Input
       --------
       press    - REAL*4 pressure array (mb)

       temp     - REAL*4 temperature profile array (K)

       wvap     - REAL*4 water vapour profile array (g/kg)

       z_sfc    - REAL*4 surface height (m).  0.0 if not known.

       n_levels - INT*4 number of elements used in passed arrays

        i_dir   - INT*4 direction of increasing layer number

                    i_dir = +1, Level[0] == p[top]         } satellite/AC
                                Level[n_levels-1] == p[sfc]  }    case

                    i_dir = -1, Level[0] == p[sfc]         } ground-based
                                Level[n_levels-1] == p[top]  }    case

        Output
       --------
          z     - REAL*4 pressure level height array (m)

    COMMENTS:

      Dimension of height array may not not be the same as that of the
        input profile data.

    ======================================================================
    python version  of gphite.f
    ======================================================================
    '''

    # -- Parameters

    rog = 29.2898
    fac = 0.5 * rog

    # Determine the level of the lowest valid temp and wvap values in the
    # profile...

    profile_mask = temp.mask + wvap.mask

    valid_idx = np.where(profile_mask!=True)[0]

    lowest_valid_idx = valid_idx[-1]
    lowest_valid_temp = temp[lowest_valid_idx]
    lowest_valid_wvap = wvap[lowest_valid_idx]
    lowest_valid_pressure = press[lowest_valid_idx]

    #LOG.debug("Lowest valid idx is {} at pressure {}".format(lowest_valid_idx,lowest_valid_pressure))
    #LOG.debug("Lowest valid temp is {} at pressure {}".format(lowest_valid_temp,lowest_valid_pressure))
    #LOG.debug("Lowest valid wvap is {} at pressure {}".format(lowest_valid_wvap,lowest_valid_pressure))

    # Reset number of levels, and truncate profile arrays to exclude missing
    # data at the bottom of the profile.
    n_levels = lowest_valid_idx+1
    t = temp[:lowest_valid_idx+1]
    w = wvap[:lowest_valid_idx+1]
    p = press[:lowest_valid_idx+1]

    #-----------------------------------------------------------------------
    #  -- Calculate virtual temperature adjustment and exponential       --
    #  -- pressure height for level above surface.  Also set integration --
    #  -- loop bounds                                                    --
    #-----------------------------------------------------------------------

    z = np.ones(n_levels)*(-999.)

    #LOG.debug("n_levels = {}".format(n_levels))
    #LOG.debug("t.shape = {}".format(t.shape))
    #LOG.debug("t[n_levels-1] = {}".format(t[-1]))

    if ( i_dir > 0 ) :

        # Data stored top down

        v_lower = t[-1] * ( 1.0 + ( 0.00061 * w[-1] ) )

        #LOG.debug("v_lower = {}".format(v_lower))

        algp_lower = np.log( p[-1] )

        #LOG.debug("algp_lower = {}".format(algp_lower))

        i_start = n_levels-1
        i_end   = 0

    else:

        # Data stored bottom up

        v_lower = t[0] * ( 1.0 + ( 0.00061 * w[0] ) )

        algp_lower = np.log( p[0] )

        i_start = 1
        i_end   = n_levels-1

    #-----------------------------------------------------------------------
    #                     -- Assign surface height --
    #-----------------------------------------------------------------------

    hgt = z_sfc

    # .. Following added 18 May 2000 ... previously, z(n_levels) for downward
    #       (usual) case was not defined!

    if(i_dir > 0):
        z[-1] = z_sfc
    else:
        z[0] = z_sfc;

    # .. End of addition

    #-----------------------------------------------------------------------
    #             -- Loop over layers always from sf% -> top --
    #-----------------------------------------------------------------------

    level_idx = np.arange(n_levels)

    # Looping from ground level to TOA.

    for l in level_idx[::-1*i_dir]:

        # Apply virtual temperature adjustment for upper level

        v_upper = t[l]

        if  ( p[l] >= 300.0 ):
            v_upper = v_upper * ( 1.0 + ( 0.00061 * w[l] ) )

        # Calculate exponential pressure height for upper layer

        algp_upper = np.log( p[l] )

        # Calculate height

        hgt = hgt + ( fac*(v_upper+v_lower)*(algp_lower-algp_upper) )


        # Overwrite values for next layer

        v_lower = v_upper
        algp_lower = algp_upper

        # Store heights in same direction as other data

        z[l] = hgt   # geopotential height

        #LOG.debug("p = {:6.1f}, t = {:6.1f}, w = {:6.1f}, z = {:6.1f}".format(
                #p[l],t[l],w[l],z[l]))

    return z


def pressure_to_height(p, t, w, z_sfc=0.):

    gc = 0.98

    z1 = gphite( p, t, w, z_sfc, 101,1) * gc
    z = z1 * 3.28 # meters to feet

    return z


def get_elevation(pressure,temp,relh):
    '''
    Given cubes of the pressure (mb), temperature (K) and relative humidity (%),
    computes the elevation for each element of the cube.
    '''


    np.set_printoptions(precision=3)
    cube_mask = np.zeros(pressure.shape,dtype='bool')

    pressure = ma.array(pressure,mask=cube_mask)
    temp = ma.array(temp,mask=cube_mask)
    relh = ma.array(relh,mask=cube_mask)

    rh_to_mr_vec = np.vectorize(rh_to_mr)

    wvap = rh_to_mr_vec(relh,pressure,temp)

    level = -4
    LOG.debug("\npressure[{}] = \n{}\n".format(level,pressure[level]))
    LOG.debug("\ntemp[{}] = \n{}\n".format(level,temp[level]))
    LOG.debug("\nrelh[{}] = \n{}\n".format(level,relh[level]))
    LOG.debug("\nwvap[{}] = \n{}\n".format(level,wvap[level]))

    elevation = np.zeros(temp.shape,dtype='float')

    (nlevels,nrows,ncols) = elevation.shape

    for row in np.arange(nrows):
        for col in np.arange(ncols):
            try:
                elevation[:,row,col] = pressure_to_height(
                        pressure[:,row,col],
                        temp[:,row,col],
                        wvap[:,row,col],
                        z_sfc=0.
                        )
            except Exception as err :
                #LOG.debug("([:,{},{}]): {}".format(row,col,err))
                #LOG.debug(traceback.format_exc())
                #LOG.warn("elevation[:,{},{}] failed".format(row,col))
                elevation[:,row,col] = -9999.

    return elevation


def get_level_indices(data,match_val):
    '''
    Given cubes of the pressure (mb), temperature (K) and relative humidity (%),
    computes the elevation for each element of the cube.
    '''

    (nlevels,nrows,ncols) = data.shape

    # Determine the level which corresponds with the required elevation...

    data_level_idx = -9999 * np.ones(data[0].shape,dtype='int')

    # Calculate the profile index corresponding to the smallest residual beteween
    # match_val and the profile values (for each row and col)
    for row in np.arange(nrows):
        for col in np.arange(ncols):
            try:
                # Get the data profile and indices for this row/col, masked with
                # the fill values.
                data_profile = ma.masked_equal(data[:,row,col],-9999.)
                data_profile_idx = ma.array(np.arange(len(data_profile)),mask=data_profile.mask)

                # Compress the profile and mask to remove fill values.
                data_profile = data_profile.compressed()
                data_profile_idx = data_profile_idx.compressed()

                if len(data_profile) == 0:
                    raise ValueError('Entire profile is masked.')

                # Compute the absolute residual of each element of this profile
                # from the match_val.
                data_profile = np.abs(data_profile - match_val)

                # Compute the minimum value of the residuals, and determine its
                # index in the original profile.
                data_profile_min = data_profile.min()
                data_profile_min_idx = np.where(data_profile == data_profile_min)[0][0]
                #LOG.debug("([:,{},{}]) data_profile_min , data_profile_min_idx = {:.3f} , {}".format(
                        #row,col,data_profile_min,data_profile_min_idx))

                # Assign the index of the minimum residual for this row/col
                data_level_idx[row,col] = data_profile_min_idx

            except Exception as err :
                #LOG.debug("([:,{},{}]): {}".format(row,col,err))
                #LOG.debug("([:,{},{}]): Setting as fill value.".format(row,col))
                data_level_idx[row,col] = -9999

    data_level_idx = ma.masked_equal(data_level_idx,-9999)
    LOG.debug("\ndata_level_idx = \n{}".format(data_level_idx))

    # Create a template index list for a single pressure level
    cube_idx = list(np.where(np.ones(data[0].shape,dtype=int)==1))
    cube_idx = [np.zeros(cube_idx[0].shape,dtype=int),cube_idx[0],cube_idx[1]]
    LOG.debug("\ncube_idx = {}".format(cube_idx))
    LOG.debug("cube_idx[0] = {} has shape {}".format(cube_idx[0], cube_idx[0].shape))
    LOG.debug("cube_idx[1] = {} has shape {}".format(cube_idx[1], cube_idx[1].shape))
    LOG.debug("cube_idx[2] = {} has shape {}".format(cube_idx[2], cube_idx[2].shape))

    # Assign the "level" index list to the
    cube_mask = ma.masked_equal(np.ravel(data_level_idx),-9999).mask
    cube_idx[0] = ma.array(data_level_idx, mask=cube_mask).compressed()
    cube_idx[1] = ma.array(cube_idx[1]   , mask=cube_mask).compressed()
    cube_idx[2] = ma.array(cube_idx[2]   , mask=cube_mask).compressed()
    cube_idx = tuple(cube_idx)

    LOG.debug("\ncompressed cube_idx = {}".format(cube_idx))
    LOG.debug("compressed cube_idx[0] = {} has shape {}".format(cube_idx[0], cube_idx[0].shape))
    LOG.debug("compressed cube_idx[1] = {} has shape {}".format(cube_idx[1], cube_idx[1].shape))
    LOG.debug("compressed cube_idx[2] = {} has shape {}".format(cube_idx[2], cube_idx[2].shape))

    return cube_idx



# def sounder(iapp_file_list, pres_0=850.):

    # data_labels = [
                    # 'pres',
                    # 'temp',
                    # 'dwpt',
                    # 'wvap',
                    # 'relh',
                    # 'ctp',
                    # 'ctt',
                    # 'lat',
                    # 'lon'
                    # ]

    # sounding_inputs = {}
    # for label in data_labels:
        # sounding_inputs[label] = {}

    # sounding_inputs['pres']['dset_name'] = 'Pressure_Levels'
    # sounding_inputs['temp']['dset_name'] = 'Temperature_Retrieval'
    # sounding_inputs['dwpt']['dset_name'] = 'Dew_Point_Temp_Retrieval'
    # sounding_inputs['wvap']['dset_name'] = 'WaterVapor_Retrieval'
    # sounding_inputs['ctp']['dset_name']  = 'Cloud_Top_Pressure_CO2'
    # sounding_inputs['ctt']['dset_name']  = 'Cloud_Top_Temperature_CO2'
    # sounding_inputs['lat']['dset_name']  = 'Latitude'
    # sounding_inputs['lon']['dset_name']  = 'Longitude'
    # sounding_inputs['relh'] = None

    # LOG.debug("Input file list: {}".format(iapp_file_list))

    # first_granule = True

    # for grans in np.arange(len(iapp_file_list)):

        # try:

            # LOG.debug("Ingesting granule {} ...".format(grans))

            # iapp_file = iapp_file_list[grans]
            # LOG.debug("Ingesting granule {} ...".format(iapp_file))

            # # Open the file object
            # fileObj = Dataset(iapp_file)

            # # Open the dataset nodes for each of the required datasets
            # sounding_inputs['pres']['node'] = fileObj.variables['Pressure_Levels']
            # sounding_inputs['temp']['node'] = fileObj.variables['Temperature_Retrieval']
            # sounding_inputs['dwpt']['node'] = fileObj.variables['Dew_Point_Temp_Retrieval']
            # sounding_inputs['wvap']['node'] = fileObj.variables['WaterVapor_Retrieval']
            # sounding_inputs['ctp']['node'] = fileObj.variables['Cloud_Top_Pressure_CO2']
            # sounding_inputs['ctt']['node'] = fileObj.variables['Cloud_Top_Temperature_CO2']
            # sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
            # sounding_inputs['lon']['node'] = fileObj.variables['Longitude']

            # if first_granule:

                # # Get the pressure scale
                # pressure = sounding_inputs['pres']['node'][:]

                # for label in data_labels:
                    # if sounding_inputs[label] != None:
                        # LOG.info(">>> Processing {} ..."
                                # .format(sounding_inputs[label]['dset_name']))
                        # for attr_name in sounding_inputs[label]['node'].ncattrs():
                            # attr = getattr(sounding_inputs[label]['node'],attr_name)
                            # LOG.debug("{} = {}".format(attr_name,attr))
                            # sounding_inputs[label][attr_name] = attr

                # # Search for the pressure level closest to the input value
                # pressure_scope = 10.
                # LOG.debug("Scope of pressure level search is {} hPa"
                        # .format(pressure_scope))
                # level = get_pressure_index(pressure,pres_0=pres_0,
                        # kd_dist=pressure_scope)
                # attempts = 1
                # while (level==None and attempts<10):
                    # pressure_scope *= 1.1
                    # LOG.debug("Scope of pressure level search is {} hPa"
                            # .format(pressure_scope))
                    # level = get_pressure_index(pressure,pres_0=pres_0,
                            # kd_dist=pressure_scope)
                    # attempts += 1

                # if level==None:
                    # raise Exception("No suitable pressure level found, aborting...")

                # LOG.debug("Retrieved level = {}".format(level))
                # sounding_inputs['pres_0'] = pressure[level]
                # data_labels.remove('pres')

            # # Add to the lat and lon and data arrays...
            # for label in ['lat','lon']:
                # try :
                    # sounding_inputs[label]['data']  = np.vstack((
                        # sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:]))
                    # LOG.debug("subsequent arrays...")
                # except KeyError :
                    # sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:]
                    # LOG.debug("first arrays...")

                # LOG.debug("Intermediate {} shape = {}".format(
                    # label,sounding_inputs[label]['data'].shape))

            # # Add to the ctp and   ctt and data arrays...
            # for label in ['ctp','ctt']:
                # try :
                    # sounding_inputs[label]['data']  = np.vstack((
                        # sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:,4]))
                    # LOG.debug("subsequent arrays...")
                # except KeyError :
                    # sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:,4]
                    # LOG.debug("first arrays...")

                # LOG.debug("Intermediate {} shape = {}".format(
                    # label,sounding_inputs[label]['data'].shape))

            # # Add to the temp, dewpoint, water vapour and relative humidity data arrays...
            # for label in ['temp','dwpt','wvap']:
                # try :
                    # sounding_inputs[label]['data']  = np.vstack((
                        # sounding_inputs[label]['data'],sounding_inputs[label]['node'][:,:,level]))
                    # LOG.debug("subsequent arrays...")
                # except KeyError :
                    # sounding_inputs[label]['data']  = sounding_inputs[label]['node'][:,:,level]
                    # LOG.debug("first arrays...")

                # LOG.debug("Intermediate {} shape = {}".format(
                    # label,sounding_inputs[label]['data'].shape))

            # first_granule = False

            # LOG.info("Closing {}".format(iapp_file))
            # fileObj.close()

        # except Exception as err :
            # LOG.warn("There was a problem, closing {}".format(iapp_file))
            # LOG.warn("{}".format(err))
            # LOG.debug(traceback.format_exc())
            # fileObj.close()
            # LOG.info("Exiting...")
            # sys.exit(1)

    # # Contruct the pressure array
    # pressure_array = pressure[level] * \
            # np.ones(sounding_inputs['lat']['data'].shape,dtype='float')
    # LOG.debug("pressure_array.shape =  {}".format(pressure_array.shape))

    # # Construct the masks of the various datasets
    # for label in ['temp','dwpt','wvap','ctp','ctt','lat','lon']:
        # if sounding_inputs[label] != None:
            # LOG.debug(">>> Processing {} ..."
                    # .format(sounding_inputs[label]['dset_name']))
            # fill_value = sounding_inputs[label]['_FillValue']
            # #data = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            # data = ma.masked_equal(sounding_inputs[label]['data'],np.nan)
            # LOG.debug("ma.is_masked({}) = {}".format(label,ma.is_masked(data)))

            # if ma.is_masked(data):
                # if data.mask.shape == ():
                    # data_mask = np.ones(data.shape,dtype='bool')
                # else:
                    # data_mask = data.mask
            # else:
                # data_mask = np.zeros(data.shape,dtype='bool')

            # LOG.debug("There are {}/{} masked values in {}".\
                    # format(np.sum(data_mask),data.size,label))

            # sounding_inputs[label]['mask'] = data_mask

    # # Computing the relative humidity
    # sounding_inputs['relh'] = {}
    # LOG.debug("wvap.shape =  {}".format(sounding_inputs['wvap']['data'].shape))
    # LOG.debug("temp.shape =  {}".format(sounding_inputs['temp']['data'].shape))

    # relh_mask = sounding_inputs['wvap']['mask']+sounding_inputs['temp']['mask']

    # wvap = ma.array(sounding_inputs['wvap']['data'],mask=relh_mask)
    # temp = ma.array(sounding_inputs['temp']['data'],mask=relh_mask)
    # pressure_array = ma.array(pressure_array,mask=relh_mask)

    # dewhum_vec = np.vectorize(dewhum)
    # _,rh,_ = dewhum_vec(pressure_array,temp,wvap)

    # sounding_inputs['relh']['data'] = rh
    # sounding_inputs['relh']['mask'] = rh.mask if  ma.is_masked(rh) \
            # else np.zeros(temp.shape,dtype='bool')
    # sounding_inputs['relh']['units'] = '%'
    # sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    # LOG.debug("ma.is_masked({}) = {}".format('relh',ma.is_masked(rh)))
    # LOG.debug("There are {}/{} masked values in relh".format(np.sum(rh.mask),rh.size))

    # return sounding_inputs


# def sounder_skewT(iapp_file_list,lat_0=None,lon_0=None):
    # '''
    # Pressure_Levels: Pressure for each level in mb (hPa)
    # Temperature_Retrieval: Temperature profile in K
    # WaterVapor_Retrieval: Water vapor profile (mixing ratio) in g/kg
    # Dew_Point_Temp_Retrieval: Dew Point Temperature Profile in K
    # '''

    # '''
    # float Latitude(Along_Track, Across_Track) ;
            # Latitude:long_name = "Geodetic Latitude" ;
            # Latitude:units = "degrees_north" ;
            # Latitude:scale_factor = 1. ;
            # Latitude:add_offset = 0. ;
            # Latitude:Parameter_Type = "IAPP Output" ;
            # Latitude:valid_range = -90.f, 90.f ;
            # Latitude:_FillValue = -999.f ;

    # float Longitude(Along_Track, Across_Track) ;
            # Longitude:long_name = "Geodetic Longitude" ;
            # Longitude:units = "degrees_east" ;

    # float Pressure_Levels(Pres_Levels) ;
            # Pressure_Levels:long_name = "Pressure Levels at which the retrieved
                                         # values are calculated" ;
            # Pressure_Levels:units = "hPa" ;

    # float Temperature_Retrieval(Along_Track, Across_Track, Pres_Levels) ;
            # Temperature_Retrieval:long_name = "Temperature Retrieval for
                                               # the IAPP" ;
            # Temperature_Retrieval:units = "degrees Kelvin" ;

    # float WaterVapor_Retrieval(Along_Track, Across_Track, Pres_Levels) ;
            # WaterVapor_Retrieval:long_name = "Water Vapor Retrieval for
                                             # the IAPP" ;
            # WaterVapor_Retrieval:units = "g/kg" ;

    # float Dew_Point_Temp_Retrieval(Along_Track, Across_Track, Pres_Levels) ;
            # Dew_Point_Temp_Retrieval:long_name = "Dew Point Temperature Retrieval
                                                  # for the IAPP" ;
            # Dew_Point_Temp_Retrieval:units = "degrees Kelvin" ;
    # '''

    # data_labels = [
                    # 'pres',
                    # 'temp',
                    # 'dwpt',
                    # 'wvap',
                    # 'relh',
                    # 'lat',
                    # 'lon'
                    # ]

    # sounding_inputs = {}
    # for label in data_labels:
        # sounding_inputs[label] = {}

    # sounding_inputs['pres']['dset_name'] = 'Pressure_Levels'
    # sounding_inputs['temp']['dset_name'] = 'Temperature_Retrieval'
    # sounding_inputs['dwpt']['dset_name'] = 'Dew_Point_Temp_Retrieval'
    # sounding_inputs['wvap']['dset_name'] = 'WaterVapor_Retrieval'
    # sounding_inputs['lat']['dset_name']  = 'Latitude'
    # sounding_inputs['lon']['dset_name']  = 'Longitude'
    # sounding_inputs['relh'] = None

    # LOG.debug("Input file list: {}".format(iapp_file_list))

    # iapp_file = iapp_file_list[0]

    # # Open the file object
    # fileObj = Dataset(iapp_file)

    # try:

        # sounding_inputs['pres']['node'] = fileObj.variables['Pressure_Levels']
        # sounding_inputs['temp']['node'] = fileObj.variables['Temperature_Retrieval']
        # sounding_inputs['dwpt']['node'] = fileObj.variables['Dew_Point_Temp_Retrieval']
        # sounding_inputs['wvap']['node'] = fileObj.variables['WaterVapor_Retrieval']
        # sounding_inputs['lat']['node'] = fileObj.variables['Latitude']
        # sounding_inputs['lon']['node'] = fileObj.variables['Longitude']

        # lat = sounding_inputs['lat']['node'][:,:]
        # lon = sounding_inputs['lon']['node'][:,:]

        # for label in data_labels:
            # if sounding_inputs[label] != None:
                # LOG.info(">>> Processing {} ...".format(sounding_inputs[label]['dset_name']))
                # for attr_name in sounding_inputs[label]['node'].ncattrs():
                    # attr = getattr(sounding_inputs[label]['node'],attr_name)
                    # LOG.debug("{} = {}".format(attr_name,attr))
                    # sounding_inputs[label][attr_name] = attr

        # row,col = get_geo_indices(lat,lon,lat_0=lat_0,lon_0=lon_0)

        # if row==None or col==None:
            # raise Exception("No suitable lat/lon coordinates found, aborting...")

        # LOG.debug("Retrieved row,col = {},{}".format(row,col))
        # sounding_inputs['lat_0'] = lat[row,col]
        # sounding_inputs['lon_0'] = lon[row,col]
        # #row,col = 10,10

        # sounding_inputs['pres']['data'] = sounding_inputs['pres']['node'][:]

        # sounding_inputs['temp']['data'] = sounding_inputs['temp']['node'][row,col,:]
        # sounding_inputs['dwpt']['data'] = sounding_inputs['dwpt']['node'][row,col,:]
        # sounding_inputs['wvap']['data'] = sounding_inputs['wvap']['node'][row,col,:]

        # LOG.info("Closing {}".format(iapp_file))
        # fileObj.close()

    # except Exception as err :
        # LOG.warn("There was a problem, closing {}".format(iapp_file))
        # LOG.warn("{}".format(err))
        # LOG.debug(traceback.format_exc())
        # fileObj.close()
        # LOG.info("Exiting...")
        # sys.exit(1)

    # # Construct the masks of the various datasets, and mask the data accordingly
    # for label in ['temp','dwpt','wvap']:
        # if sounding_inputs[label] != None:
            # LOG.debug(">>> Processing {} ..."
                    # .format(sounding_inputs[label]['dset_name']))
            # fill_value = sounding_inputs[label]['_FillValue']
            # data = ma.masked_equal(sounding_inputs[label]['data'],fill_value)
            # LOG.debug("ma.is_masked({}) = {}".format(label,ma.is_masked(data)))

            # if ma.is_masked(data):
                # if data.mask.shape == ():
                    # data_mask = np.ones(data.shape,dtype='bool')
                # else:
                    # data_mask = data.mask
            # else:
                # data_mask = np.zeros(data.shape,dtype='bool')

            # LOG.debug("There are {}/{} masked values in {}".\
                    # format(np.sum(data_mask),data.size,label))

            # sounding_inputs[label]['mask'] = data_mask
            # sounding_inputs[label]['data'] = ma.array(data,mask=data_mask)

    # # Computing the relative humidity
    # sounding_inputs['relh'] = {}
    # dewhum_vec = np.vectorize(dewhum)
    # _,rh,_ = dewhum_vec(
            # sounding_inputs['pres']['data'],
            # sounding_inputs['temp']['data'],
            # sounding_inputs['wvap']['data']
            # )

    # sounding_inputs['relh']['data'] = rh
    # sounding_inputs['relh']['units'] = '%'
    # sounding_inputs['relh']['long_name'] = 'Relative Humidity'

    # return sounding_inputs

