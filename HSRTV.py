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

#from netCDF4 import Dataset
#from netCDF4 import num2date
#import h5py

# every module should have a LOG object
LOG = logging.getLogger(__file__)

from ql_common import get_pressure_index,get_geo_indices,get_pressure_level
from ql_common import Datafile_HDF5

# Files dataset names for each dataset type
dset_name = {}
dset_name['pres'] = '/Plevs'
dset_name['lat']  = '/Latitude'
dset_name['lon']  = '/Longitude'
dset_name['ctp']  = '/CTP'
dset_name['ctt']  = '/CTT'
dset_name['temp'] = '/TAir'
dset_name['2temp'] = None
dset_name['cold_air_aloft'] = None
dset_name['dwpt'] = '/Dewpnt'
dset_name['relh'] = '/RelHum'
dset_name['wvap'] = '/H2OMMR'

# Dataset types (surface/top; pressure level...)
dset_type = {}
dset_type['pres'] = 'profile'
dset_type['lat']  = 'single'
dset_type['lon']  = 'single'
dset_type['ctp']  = 'single'
dset_type['ctt']  = 'single'
dset_type['temp'] = 'level'
dset_type['2temp'] = 'level'
dset_type['cold_air_aloft'] = 'level'
dset_type['dwpt'] = 'level'
dset_type['relh'] = 'level'
dset_type['wvap'] = 'level'

# Dataset dependencies for each dataset
dset_deps = {}
dset_deps['pres'] = []
dset_deps['lat']  = []
dset_deps['lon']  = []
dset_deps['ctp']  = []
dset_deps['ctt']  = []
dset_deps['temp'] = []
dset_deps['2temp'] = ['temp']
dset_deps['cold_air_aloft'] = ['temp']
dset_deps['dwpt'] = []
dset_deps['relh'] = []
dset_deps['wvap'] = []

#dset_deps = {}
#dset_deps['pres'] = []
#dset_deps['lat']  = []
#dset_deps['lon']  = []
#dset_deps['ctp']  = []
#dset_deps['ctt']  = []
#dset_deps['temp'] = ['pres']
#dset_deps['dwpt'] = ['temp']
#dset_deps['relh'] = ['temp','wvap']
#dset_deps['wvap'] = ['temp','dwpt']

# The class method used to read/calculate each dataset.
dset_method = {}
dset_method['pres'] = []
dset_method['lat']  = []
dset_method['lon']  = []
dset_method['ctp']  = []
dset_method['ctt']  = []
dset_method['temp'] = []
dset_method['2temp'] = []
dset_method['cold_air_aloft'] = []
dset_method['dwpt'] = []
dset_method['relh'] = []
dset_method['wvap'] = []


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

    def __init__(self,file_list,data_name,plot_type,pres_0=850.,lat_0=None,lon_0=None):

        self.file_list = file_list
        self.have_geo = False
        self.datasets = {}

        # Construct a list of the required datasets...
        dsets_to_read = [data_name] + self.get_dset_deps(data_name)
        LOG.info("Dupe Datasets to read : {}".format(dsets_to_read))
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

        # Read in the pressure dataset
        LOG.info(">>> Reading in the pressure dataset...")
        dfile_obj = Datafile_HDF5(file_list[0])
        data_obj = dfile_obj.Dataset(dfile_obj,dset_name['pres'])
        self.pressure = data_obj.dset[:]
        dfile_obj.close()

        LOG.debug("pressure = {}".format(self.pressure))

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
            print ""
            LOG.info(">>> Reading in the lat and lon arrays...")
            self.construct_level_pass(file_list,None,None,None)

            LOG.debug("\n\n>>> Intermediate dataset manifest...\n")
            LOG.debug(self.print_dataset_manifest(self))

            LOG.info("Computing the row/col and lon_0/lat_0 from the lon/lat pass...")
            self.row,self.col,self.lat_0,self.lon_0 \
                    = get_geo_indices(self.datasets['lat']['data'],
                                      self.datasets['lon']['data'],
                                      lat_0=None,lon_0=None)
            

            for dsets in ['lat_row','lat_col','lon_row','lon_col']:
                self.datasets[dsets] = {}

            self.datasets['lat_row']['data'] = self.datasets['lat']['data'][self.row,:]
            self.datasets['lat_col']['data'] = self.datasets['lat']['data'][:,self.col]
            self.datasets['lat_row']['attrs'] =  dict(self.datasets['lat']['attrs'])
            self.datasets['lat_col']['attrs'] =  dict(self.datasets['lat']['attrs'])

            self.datasets['lon_row']['data'] = self.datasets['lon']['data'][self.row,:]
            self.datasets['lon_col']['data'] = self.datasets['lon']['data'][:,self.col]
            self.datasets['lon_row']['attrs'] =  dict(self.datasets['lon']['attrs'])
            self.datasets['lon_col']['attrs'] =  dict(self.datasets['lon']['attrs'])

            self.dsets_to_read = list(total_dsets_to_read)
            self.dsets_to_read.remove('lat')
            self.dsets_to_read.remove('lon')

            print ""
            LOG.info(">>> Reading in the remaining dataset slices..")
            LOG.debug("(level,row,col)  = ({}, {}, {})".format(None,self.row,self.col))
            # This is to slice along the rows
            self.construct_slice_pass(file_list,None,None,self.col)
            # This is to slice along the cols
            #self.construct_slice_pass(file_list,None,self.row,None)

            self.dsets_to_read = list(total_dsets_to_read)


        LOG.debug("\n\n>>> Final dataset manifest...\n")
        LOG.debug(self.print_dataset_manifest(self))


    def construct_slice_pass(self,file_list,level,row,col):
        '''
        Read in each of the desired datasets. Each of the desired datasets may
        be either a "base" dataset, which is read directly from the file, or a 
        "derived" dataset, which is constructed from previously read datasets.
        '''

        LOG.info("\tContructing a SLICE pass...")
        LOG.debug("\tInput file list: {}".format(file_list))
        LOG.debug("\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        first_granule = True
        self.datasets['file_attrs'] = {}
        for dset in self.dsets_to_read:
            self.datasets[dset] = {}

        # Loop through each of the granules...
        for grans in np.arange(len(file_list)):

            file_name = file_list[grans]

            try:

                LOG.info("\tOpening file {}".format(file_name))
                dfile_obj = Datafile_HDF5(file_name)
                self.datasets['file_attrs'][os.path.basename(file_name)] \
                        = dfile_obj.attrs
                self.datasets['file_attrs'][os.path.basename(file_name)]['dt_obj'] \
                        = self.get_granule_dt(file_name)

                # Loop through each of the desired datasets
                for dset in self.dsets_to_read:
                    LOG.info("\tFor dataset {}".format(dset))

                    # Choose the correct "get" method for the dataset
                    if dset_method[dset] == []:
                        LOG.info("\tUsing read method for {}".format(dset))
                        get_data = self.get_data
                    else:
                        LOG.info("\tUsing compute method for {}".format(dset))
                        get_data = dset_method[dset]

                    LOG.info("\tReading in granule {} of {}".format(grans,dset))

                    try :

                        # This is to slice along the rows
                        self.datasets[dset]['data'] = np.vstack((
                            self.datasets[dset]['data'],
                            get_data(dfile_obj,level,row,col)
                            ))
                        LOG.debug("\tsubsequent arrays...")

                    except KeyError :

                        # This is to slice along the rows
                        self.datasets[dset]['data'] = get_data(dfile_obj,dset,level,row,col)
                        LOG.debug("\tfirst arrays...")

                    LOG.debug("\tIntermediate {} shape = {}".format(
                        dset,self.datasets[dset]['data'].shape))

                dfile_obj.close()

            except Exception, err :
                LOG.warn("\tThere was a problem, closing {}".format(file_name))
                LOG.warn("\t{}".format(err))
                LOG.debug(traceback.format_exc())
                dfile_obj.close()
                LOG.info("\tExiting...")
                sys.exit(1)


    def construct_level_pass(self,file_list,level,row,col):
        '''
        Read in each of the desired datasets. Each of the desired datasets may
        be either a "base" dataset, which is read directly from the file, or a 
        "derived" dataset, which is constructed from previously read datasets.
        '''

        LOG.info("\tContructing a LEVEL pass...")
        LOG.debug("\tInput file list: {}".format(file_list))
        LOG.debug("\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        first_granule = True
        self.datasets['file_attrs'] = {}
        for dset in self.dsets_to_read:
            self.datasets[dset] = {}

        # Loop through each of the granules...
        for grans in np.arange(len(file_list)):

            file_name = file_list[grans]

            try:

                LOG.info("\tOpening file {}".format(file_name))
                dfile_obj = Datafile_HDF5(file_name)
                self.datasets['file_attrs'][os.path.basename(file_name)] \
                        = dfile_obj.attrs
                self.datasets['file_attrs'][os.path.basename(file_name)]['dt_obj'] \
                        = self.get_granule_dt(file_name)

                # Loop through each of the desired datasets
                for dset in self.dsets_to_read:
                    LOG.info("\tFor dataset {}".format(dset))

                    # Choose the correct "get" method for the dataset
                    if dset_method[dset] == []:
                        LOG.info("\tUsing read method for {}".format(dset))
                        get_data = self.get_data
                    else:
                        LOG.info("\tUsing compute method for {}".format(dset))
                        get_data = dset_method[dset]

                    LOG.info("\tReading in granule {} of {}".format(grans,dset))

                    try :

                        self.datasets[dset]['data'] = np.vstack((
                            self.datasets[dset]['data'],
                            get_data(dfile_obj,dset,level,row,col)))
                        LOG.debug("\tsubsequent arrays...")

                    except KeyError :

                        self.datasets[dset]['data'] = get_data(dfile_obj,dset,level,row,col)
                        LOG.debug("\tfirst arrays...")

                    LOG.debug("\tIntermediate {} shape = {}".format(
                        dset,self.datasets[dset]['data'].shape))

                dfile_obj.close()

            except Exception, err :
                LOG.warn("\tThere was a problem, closing {}".format(file_name))
                LOG.warn("\t{}".format(err))
                LOG.debug(traceback.format_exc())
                dfile_obj.close()
                LOG.info("\tExiting...")
                sys.exit(1)


    def get_data(self,dfile_obj,data_name,level,row,col):
        '''
        This method reads a single granule of the required dataset, and returns
        a single array (which may be a surface/top, pressure level, or slice dataset.
        '''

        LOG.debug("\t\t(level,row,col)  = ({}, {}, {})".format(level,row,col))

        level = slice(level) if level == None else level
        row = slice(row) if row == None else row
        col = slice(col) if col == None else col

        data_obj = dfile_obj.Dataset(dfile_obj,dset_name[data_name])

        if dset_type[data_name] == 'single':
            LOG.debug("\tgetting a 'single' dataset...")
            dset = data_obj.dset[row,col].squeeze()
        elif dset_type[data_name] == 'level': 
            LOG.debug("\tgetting a 'level' dataset...")
            dset = data_obj.dset[level,row,col].squeeze()


        LOG.debug("\t\tDataset {} has shape {}".format(data_name,dset.shape))

        self.datasets[data_name]['attrs'] =  dict(data_obj.attrs)

        LOG.info("For dataset {}:".format(data_name))
        LOG.info("\t\tdatasets['{}']['attrs'] = {}".format(
                data_name,self.datasets[data_name]['attrs']))

        return dset


    def temp_times_two(self,dfile_obj,data_name,level):
        '''
        Custom method to compute two times the temperature.
        '''
        LOG.info("\t\tComputing {}".format(data_name))
        dset = 2. * self.datasets['temp']['data'] + 0.01*self.pres_0

        self.datasets[data_name]['attrs'] =  self.datasets['temp']['attrs']

        return dset


    def cold_air_aloft(self,dfile_obj,data_name,level):
        '''
        Custom method to returns the temperature.
        '''
        LOG.info("\t\tComputing {}".format(data_name))
        dset_mask = self.datasets['temp']['data'].mask

        dset = self.datasets['temp']['data'] - 273.16

        low_idx = np.where(dset < -65.)
        mid_idx = np.where((dset > -65.) * (dset < -60.))
        hi_idx = np.where(dset > -60.)
        dset[low_idx] = 0.
        dset[mid_idx] = 1.
        dset[hi_idx] = 2.

        dset = ma.array(dset,mask=dset_mask)

        self.datasets[data_name]['attrs'] = self.datasets['temp']['attrs']

        return dset


    def print_dataset_manifest(self,pkg_obj):
        dataset_names = list(pkg_obj.datasets.keys())
        dataset_names.remove('file_attrs')
        LOG.info("datasets: {}".format(dataset_names))

        for key in dataset_names:
            LOG.info("For dataset {}:".format(key))
            LOG.info("\tdatasets['{}']['attrs'] = {}".format(
                    key,pkg_obj.datasets[key]['attrs']))
            LOG.info("\tdatasets['{}']['data'].shape = {}".format(
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


