#!/usr/bin/env python
# encoding: utf-8
"""
HSRTV.py

Purpose: Provide read methods for various datasets associated with the 
Hyperspectral Retrieval (HSRTV) Package.

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

from ql_common import get_pressure_index,get_geo_indices,get_pressure_level
from ql_common import Datafile_HDF5
from thermo import rh_to_mr


# every module should have a LOG object
LOG = logging.getLogger(__file__)

# Files dataset names for each dataset type
dset_name = {}
dset_name['pres'] = '/Plevs'
dset_name['pres_array'] = None
dset_name['lat']  = '/Latitude'
dset_name['lon']  = '/Longitude'
dset_name['ctp']  = '/CTP'
dset_name['ctt']  = '/CTT'
dset_name['temp'] = '/TAir'
dset_name['temp_gdas'] = '/GDAS_TAir'
#dset_name['temp_gdas_cube'] = None
dset_name['2temp'] = None
dset_name['cold_air_aloft'] = None
dset_name['dwpt'] = '/Dewpnt'
dset_name['relh'] = '/RelHum'
dset_name['relh_gdas'] = '/GDAS_RelHum'
#dset_name['relh_gdas_cube'] = None
dset_name['wvap'] = '/H2OMMR'
dset_name['wvap_gdas'] = None
#dset_name['wvap_gdas_cube'] = None

# Dataset types (surface/top; pressure level...)
dset_type = {}
dset_type['pres'] = 'profile'
dset_type['pres_array'] = 'level'
dset_type['lat']  = 'single'
dset_type['lon']  = 'single'
dset_type['ctp']  = 'single'
dset_type['ctt']  = 'single'
dset_type['temp'] = 'level'
dset_type['temp_gdas'] = 'level'
#dset_type['temp_gdas_cube'] = 'level'
dset_type['2temp'] = 'level'
dset_type['cold_air_aloft'] = 'level'
dset_type['dwpt'] = 'level'
dset_type['relh'] = 'level'
dset_type['relh_gdas'] = 'level'
#dset_type['relh_gdas_cube'] = 'level'
dset_type['wvap'] = 'level'
dset_type['wvap_gdas'] = 'level'
#dset_type['wvap_gdas_cube'] = 'level'

# Dataset dependencies for each dataset
dset_deps = {}
dset_deps['pres'] = []
dset_deps['pres_array'] = []
dset_deps['lat']  = []
dset_deps['lon']  = []
dset_deps['ctp']  = []
dset_deps['ctt']  = []
dset_deps['temp'] = []
dset_deps['temp_gdas'] = []
#dset_deps['temp_gdas_cube'] = []
dset_deps['2temp'] = ['temp']
dset_deps['cold_air_aloft'] = ['temp','wvap','temp_gdas','relh_gdas']
dset_deps['dwpt'] = []
dset_deps['relh'] = []
dset_deps['relh_gdas'] = []
#dset_deps['relh_gdas_cube'] = []
dset_deps['wvap'] = []
dset_deps['wvap_gdas'] = ['temp_gdas','pres_array','relh_gdas']
#dset_deps['wvap_gdas_cube'] = ['temp_gdas','pres_array','relh_gdas']

# The class method used to read/calculate each dataset.
dset_method = {}
dset_method['pres'] = []
dset_method['pres_array'] = []
dset_method['lat']  = []
dset_method['lon']  = []
dset_method['ctp']  = []
dset_method['ctt']  = []
dset_method['temp'] = []
dset_method['temp_gdas'] = []
#dset_method['temp_gdas_cube'] = []
dset_method['2temp'] = []
dset_method['cold_air_aloft'] = []
dset_method['dwpt'] = []
dset_method['relh'] = []
dset_method['relh_gdas'] = []
#dset_method['relh_gdas_cube'] = []
dset_method['wvap'] = []
dset_method['wvap_gdas'] = []
#dset_method['wvap_gdas_cube'] = []


class HSRTV():
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

    def __init__(self,file_list,data_name,plot_type,pres_0=None,elev_0=None,
            lat_0=None,lon_0=None,footprint=None):

        self.file_list = file_list
        self.have_geo = False
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
        dset_method['2temp'] = self.temp_times_two
        dset_method['cold_air_aloft'] = self.cold_air_aloft
        #dset_method['temp_gdas'] = self.temp_gdas
        #dset_method['relh_gdas'] = self.relh_gdas
        dset_method['wvap_gdas'] = self.wvap_gdas
        dset_method['pres_array'] = self.pressure_array

        # Read in the pressure dataset
        LOG.info(">>> Reading in the pressure dataset...")
        dfile_obj = Datafile_HDF5(file_list[0])
        data_obj = dfile_obj.Dataset(dfile_obj,dset_name['pres'])
        self.pressure = data_obj.dset[:]
        dfile_obj.close()

        # Determine the elevation dataset...
        self.elev_0 = elev_0
        if elev_0 != None:


            # Create a cube of the temp_gdas and relh_gdas datasets...
            LOG.info("Constructing a cube dataset...")
            self.construct_cube_pass(file_list,['temp_gdas','relh_gdas'])

            # Construct the pressure cube
            #cube_shape = self.datasets['temp_gdas']['data'].shape
            #nlevels,nrows,ncols = cube_shape[0],cube_shape[1],cube_shape[2]
            (nlevels,nrows,ncols) = self.datasets['temp_gdas']['data'].shape
            self.datasets['pressure'] = {}
            self.datasets['pressure']['data'] = np.broadcast_to(np.broadcast_to(self.pressure,(nrows,nlevels)),(ncols,nrows,nlevels)).T
            self.datasets['pressure']['attrs'] = {}
            LOG.info("pressure cube shape = {}".format(self.datasets['pressure']['data'].shape))

            # Compute the wvap_gdas dataset
            self.datasets['wvap_gdas'] = {}
            self.datasets['wvap_gdas']['attrs'] = {}
            self.datasets['wvap_gdas']['data'] = np.array([1]) # placeholder
            level = -4
            press = self.datasets['pressure']['data']
            relh = ma.array(self.datasets['relh_gdas']['data'],mask=self.datasets['relh_gdas']['data_mask'])
            temp = ma.array(self.datasets['temp_gdas']['data'],mask=self.datasets['temp_gdas']['data_mask'])

            # Compute the elevation cube
            self.elevation = get_elevation(press,temp,relh)

            # Get the cube indicies that correspond to the value of match_val...
            match_val = elev_0
            cube_idx = get_level_indices(self.elevation,match_val)

            elev_level = -9999. * np.ones(self.elevation[0].shape,dtype='float')
            elev_level[(cube_idx[1],cube_idx[2])] = self.elevation[cube_idx]
            LOG.info("\nelev_level = \n{}".format(elev_level))

        #sys.exit(0)

        # 
        if plot_type == 'image':

            LOG.info("Preparing a 'level' plot...")
            # Determine the level closest to the required pressure
            self.level,self.pres_0 = get_pressure_level(self.pressure,pres_0)

            # Contruct a pass of the required datasets at the desired pressure.
            self.construct_level_pass(file_list,self.level,None,None)
            LOG.debug("\n\n>>> Intermediate dataset manifest...\n")
            LOG.debug(self.print_dataset_manifest(self))

        elif plot_type == 'slice':

            LOG.info("Preparing a 'slice' plot...")
            self.level,self.pres_0 = None,0

            total_dsets_to_read = list(self.dsets_to_read)
            self.dsets_to_read = ['lat','lon']

            # Getting the full latitude and longitude
            LOG.info(">>> Reading in the lat and lon arrays...")
            self.construct_level_pass(file_list,None,None,None)

            LOG.debug("\n\n>>> Intermediate dataset manifest...\n")
            LOG.debug(self.print_dataset_manifest(self))

            LOG.info("Computing the row/col and lon_0/lat_0 from the lon/lat pass...")
            self.row,self.col,self.lat_0,self.lon_0 \
                    = get_geo_indices(self.datasets['lat']['data'],
                                      self.datasets['lon']['data'],
                                      lat_0=None,lon_0=None)
            if footprint != None:
                self.col = footprint

            for dsets in ['lat_row','lat_col','lon_row','lon_col']:
                self.datasets[dsets] = {}

            #self.col = 64 # GPC: FIXME

            self.datasets['lat_row']['data'] = self.datasets['lat']['data'][self.row,:]
            self.datasets['lat_col']['data'] = self.datasets['lat']['data'][:,self.col]
            self.datasets['lat_row']['attrs'] =  dict(self.datasets['lat']['attrs'])
            self.datasets['lat_col']['attrs'] =  dict(self.datasets['lat']['attrs'])

            self.datasets['lon_row']['data'] = self.datasets['lon']['data'][self.row,:]
            self.datasets['lon_col']['data'] = self.datasets['lon']['data'][:,self.col]
            self.datasets['lon_row']['attrs'] =  dict(self.datasets['lon']['attrs'])
            self.datasets['lon_col']['attrs'] =  dict(self.datasets['lon']['attrs'])

            if elev_0 != None:
                self.datasets['elev_slice'] = {}
                self.datasets['elev_slice']['data'] = self.elevation[:,:,self.col]
                self.datasets['elev_slice']['attrs'] = {}

            self.dsets_to_read = list(total_dsets_to_read)
            self.dsets_to_read.remove('lat')
            self.dsets_to_read.remove('lon')

            LOG.info(">>> Reading in the remaining dataset slices..")
            LOG.debug("(level,row,col)  = ({}, {}, {})".format(None,self.row,self.col))
            # This is to slice along the rows
            self.construct_slice_pass(file_list,None,None,self.col)
            # This is to slice along the cols
            #self.construct_slice_pass(file_list,None,self.row,None)

            self.dsets_to_read = list(total_dsets_to_read)


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
                dfile_obj = Datafile_HDF5(file_name)
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

                    missing_value = float(self.datasets[dset]['attrs']['missing_value'])
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

                    except KeyError,err :

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

            except Exception, err :
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
        LOG.debug("\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

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
                dfile_obj = Datafile_HDF5(file_name)
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
                        #LOG.info("\tUsing read method for {}".format(dset))
                        get_data = self.get_data
                    else:
                        #LOG.info("\tUsing compute method for {}".format(dset))
                        get_data = dset_method[dset]

                    #LOG.info("\tReading in granule {} of {}".format(grans,dset))

                    data = get_data(dfile_obj,dset,level,row,col,
                            this_granule_data,this_granule_mask)
                    this_granule_data[dset] = data

                    LOG.info("\t\tthis_granule_data['{}'].shape = {}".format(
                        dset,this_granule_data[dset].shape))

                    missing_value = float(self.datasets[dset]['attrs']['missing_value'])
                    #LOG.info("\t\tMissing value = {}".format(missing_value))
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

                    except KeyError,err :

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

            except Exception, err :
                LOG.warn("\tThere was a problem, closing {}".format(file_name))
                LOG.warn("\t{}".format(err))
                LOG.debug(traceback.format_exc())
                LOG.info("\tClosing file {}".format(file_name))
                dfile_obj.close()
                LOG.info("\tExiting...")
                sys.exit(1)


    def construct_level_pass(self,file_list,level,row,col):
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
                dfile_obj = Datafile_HDF5(file_name)
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

                    missing_value = float(self.datasets[dset]['attrs']['missing_value'])
                    #LOG.info("\t\tMissing value = {}".format(missing_value))
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

                    except KeyError,err :

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

            except Exception, err :
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

        LOG.info("\t\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        level = slice(level) if level == None else level
        row = slice(row) if row == None else row
        col = slice(col) if col == None else col

        data_obj = dfile_obj.Dataset(dfile_obj,dset_name[data_name])

        if dset_type[data_name] == 'single':
            LOG.info("\t\tgetting a 'single' dataset...")
            dset = data_obj.dset[row,col].squeeze()
        elif dset_type[data_name] == 'level': 
            LOG.info("\t\tgetting a 'level' dataset...")
            dset = data_obj.dset[level,row,col].squeeze()


        LOG.info("\t\tDataset {} has shape {}".format(data_name,dset.shape))

        self.datasets[data_name]['attrs'] =  dict(data_obj.attrs)

        return dset


    #def latlon_array(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        #'''
        #Custom method to return an array of the latitude or longitude, whether it be at a single
        #pressure level, or a vertical slice.
        #'''

        #LOG.info("\t\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        #level = slice(level) if level == None else level
        #row = slice(row) if row == None else row
        #col = slice(col) if col == None else col

        #LOG.info("\t\tComputing {}".format(data_name))

        #LOG.info("\t\tlat has shape {}".format(this_granule_data['lat'].shape))

        ## Contruct the pressure array.
        #nrows = this_granule_data['lat'].shape[0]
        #ncols = this_granule_data['lat'].shape[1]
        #nlevels = len(self.pressure)
        #self.pressure_array = np.broadcast_to(np.broadcast_to(self.pressure,(nrows,nlevels)),(ncols,nrows,nlevels)).T

        #LOG.info("\t\tpressure_array.shape = {}".format(self.pressure_array.shape))

        #LOG.info("\t\tgetting a 'level' dataset...")
        #dset = self.pressure_array[level,row,col].squeeze()
        #LOG.info("\t\tdset.shape = {}".format(dset.shape))
        #dset_mask = np.zeros(dset.shape,dtype='bool')

        #dset = ma.array(dset,mask=dset_mask)

        #self.datasets[data_name]['attrs'] = self.datasets['lat']['attrs']

        #return dset


    def pressure_array(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to return an array of the pressure, whether it be at a single
        pressure level, or a vertical slice.
        '''

        LOG.info("\t\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        level = slice(level) if level == None else level
        row = slice(row) if row == None else row
        col = slice(col) if col == None else col

        LOG.info("\t\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        LOG.info("\t\tComputing {}".format(data_name))

        LOG.info("\t\trelh_gdas has shape {}".format(this_granule_data['relh_gdas'].shape))

        # Determine whether this is an level or slice, and determine the size
        # of the relh_gdas array...
        nlevels = len(self.pressure)
        if this_granule_data['relh_gdas'].shape[0] == nlevels:
            LOG.info("\t\t> getting a pressure slice...")
            nrows = this_granule_data['relh_gdas'].shape[1]
            ncols = None
            #self.pressure_array = np.broadcast_to(np.broadcast_to(self.pressure,(nrows,nlevels)),(ncols,nrows,nlevels)).T
        else:
            LOG.info("\t\t> getting a pressure level...")
            nrows = this_granule_data['relh_gdas'].shape[0]
            ncols = this_granule_data['relh_gdas'].shape[1]
            #self.pressure_array = np.broadcast_to(np.broadcast_to(self.pressure,(nrows,nlevels)),(ncols,nrows,nlevels)).T


        # Contruct the pressure array.
        nrows = this_granule_data['relh_gdas'].shape[0]
        ncols = this_granule_data['relh_gdas'].shape[1]
        LOG.info("\t\t(nlevels,nrows,ncols)  = ({}, {}, {})".format(nlevels,nrows,ncols))
        self.pressure_array = np.broadcast_to(np.broadcast_to(self.pressure,(nrows,nlevels)),(ncols,nrows,nlevels)).T

        LOG.info("\t\tpressure_array.shape = {}".format(self.pressure_array.shape))

        LOG.info("\t\tgetting a 'level' dataset...")
        dset = self.pressure_array[level,row,col].squeeze()
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


    def wvap_gdas(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to return the water vapor mixing ratio.
        '''

        LOG.info("\t\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        level = slice(level) if level == None else level
        row = slice(row) if row == None else row
        col = slice(col) if col == None else col

        LOG.info("\t\tComputing {}".format(data_name))
        dset_mask = this_granule_mask['temp_gdas'] + this_granule_mask['relh_gdas']

        pressure  = this_granule_data['pres_array'][:,:].squeeze()
        temp_gdas = this_granule_data['temp_gdas'][:,:].squeeze()
        relh_gdas = this_granule_data['relh_gdas'][:,:].squeeze()

        LOG.info("\t\tpressure  has shape {}".format(pressure.shape))
        LOG.info("\t\tGDAS Temp has shape {}".format(this_granule_data['temp_gdas'].shape))
        LOG.info("\t\tGDAS RH   has shape {}".format(this_granule_data['relh_gdas'].shape))

        # Get the pressure 
        #slice_length = temp_gdas.shape[1]
        #pressure_row = slice(0,slice_length)
        #pressure = self.pressure_array[level,pressure_row,col].squeeze()
        #LOG.info("\t\tNew pressure  has shape {}".format(pressure.shape))

        rh_to_mr_vec = np.vectorize(rh_to_mr)
        dset = rh_to_mr_vec( relh_gdas, pressure, temp_gdas)

        dset = ma.array(dset,mask=dset_mask)
        self.datasets[data_name]['attrs'] = self.datasets['relh_gdas']['attrs']

        return dset


    def temp_times_two(self,dfile_obj,data_name,level,row,col,this_granule_data,this_granule_mask):
        '''
        Custom method to compute two times the temperature.
        '''
        LOG.info("\t\tComputing {}".format(data_name))
        dset = 2. * self.datasets['temp']['data'] + 0.01*self.pres_0

        self.datasets[data_name]['attrs'] =  self.datasets['temp']['attrs']

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
        print ""

    def get_data_old(self,data_file,data_name,level):

        # Contruct the pressure array
        pressure_array = pressure[level] * \
                np.ones(sounding_inputs['lat']['data'].shape,dtype='float')
        LOG.debug("pressure_array.shape =  {}".format(pressure_array.shape))

        # Construct the masks of the various datasets
        for label in ['temp','dwpt','wvap','relh','ctp','ctt','lat','lon']:
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


        return sounding_inputs


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
            "IASI_d20150331_t153700_M01.atm_prof_rtv.h5"
        '''
        file_name = os.path.basename(file_name)
        image_date_time = "_".join(file_name.split("_")[1:3])
        image_date_time = image_date_time.split(".")[0]
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
    #gc = 1.

    z_sfc = 3161.
    z1 = gphite( p, t, w, z_sfc, 101,1) * gc
    z = z1 #* 3.28 # meters to feet

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
            except Exception, err :
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

            except Exception, err :
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


def get_pressure_levels():
    # Construct the pressure cube
    nlevs,nrows,ncols = 6,3,5
    little_pressure = np.arange(nlevs)*200.
    pressure_array = np.broadcast_to(np.broadcast_to(little_pressure,(nrows,nlevs)),(ncols,nrows,nlevs)).T

    # Get the indices of a single level
    pressure_idx = list(np.where(pressure_array==200.))

    # Construct a list of level indices
    lev_idx = np.arange(nrows*ncols)
    lev_idx[lev_idx>5] = 2

    # Replace the level indices in pressure_idx with the custom level indices in lev_idx
    pressure_idx[0] = lev_idx
    pressure_idx = tuple(pressure_idx)

    pressure_array[pressure_idx]
    pressure_levels = pressure_array[pressure_idx].reshape((nrows,ncols))


def press2alt(p,z):
    '''
    
    Routine to get altimeter setting from pressure and elevation.

    Inputs/Outputs:

       Variable     Var Type     I/O   Description
      ----------   ----------   ----- -------------
       p               RA         O    Pressure (X)
       z               RA         I    Elevation in meters.
       alt             RA         I    Altimeter setting (X)


    User Notes:

    1.  No quality control is performed in this routine.

    '''

    flg = 99998.0
    flag = 1e37
    T0 = 288.0
    gamma = 0.0065
    g_Rgamma = 5.2532

    alt = p / ((T0-gamma * z)/T0)**g_Rgamma

    return alt
