#!/usr/bin/env python
# encoding: utf-8
"""
sounder_image.py

Purpose: Create a plot of temperature, dewpoint or relative humidity,
         at a particular pressure level. Supports outputs from the following
         packages...

         * International ATOVS Processing Package (IAPP)
         * Microwave Integrated Retrieval System (MiRS)
         * CSPP Hyperspectral Retrieval (HSRTV) Package
         * NOAA Unique CrIS/ATMS Product System (NUCAPS)
         * Hyper-Spectral Enterprise Algorithm Package (HEAP NUCAPS)


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

    DATATYPE: One of 'iapp','mirs', 'hsrtv', 'nucaps' OR 'heap'.


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2014-05-10.
Copyright (c) 2018 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

import os
import sys
import logging
import traceback
from cffi import FFI

#import string
#import re
#import uuid
#from shutil import rmtree, copyfile
#from glob import glob
#from time import time
#from datetime import datetime, timedelta

from args import argument_parser_image

import numpy as np
from numpy import ma
import copy

from scipy import vectorize

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

from netCDF4 import Dataset
from netCDF4 import num2date
import h5py

import sounder_image_data
from thermo import dewhum
from ql_common import granuleFiles
from ql_common import plotMapDataContinuous, plotMapDataContinuous_cartopy, plotMapDataDiscrete
from ql_common import set_plot_styles

import sounder_packages

# every module should have a LOG object
LOG = logging.getLogger(__file__)

def main():
    '''
    The main method.
    '''

    # Read in the command line options
    args, work_dir = argument_parser_image()


    input_files = args.inputs
    datatype = args.datatype
    dataset = args.dataset
    stride = args.stride
    pressure = args.pressure
    # elevation = args.elevation
    lat_0  = args.lat_0
    lon_0  = args.lon_0
    # latMin = args.latMin
    # lonMin = args.lonMin
    # latMax = args.latMax
    # lonMax = args.lonMax
    # yMin = args.yMin
    # yMax = args.yMax
    # footprint = args.footprint
    bounding_lat = args.bounding_lat
    extent = args.extent
    plotMin = args.plotMin
    plotMax = args.plotMax
    plot_type = args.plot_type
    # scale = args.scale
    doScatterPlot = args.doScatterPlot
    pointSize = args.pointSize
    proj = args.proj
    map_res = args.map_res
    output_file  = args.output_file
    outputFilePrefix  = args.outputFilePrefix
    dpi = args.dpi

    # Obtain matched lists of the geolocation and product files
    LOG.info("Input file(s):\n\t{}".format('\n\t'.join(input_files)))
    input_file_list = granuleFiles(input_files)

    if input_file_list==[]:
        LOG.warn('No files match the input file glob "{}", aborting.'
                .format(input_files))
        sys.exit(1)


    LOG.debug("Input file(s):\n\t{}".format('\n\t'.join(input_file_list)))

    # Read in the input file, and return a dictionary containing the required
    # data

    LOG.info("Input pressure = {}".format(pressure))
    # LOG.info("Input elevation = {}".format(elevation))

    # Get all of the required command line options for the required souder package...
    LOG.info("datatype = {}".format(datatype))
    LOG.info("dataset = {}".format(dataset))
    sounder_args = (input_file_list,dataset,plot_type)
    sounder_kwargs = {'pres_0':pressure,
                      # 'elev_0':elevation,
                      'lat_0':None,
                      'lon_0':None,
                      # 'footprint':footprint
                      }

    # Get a reference to the desired sounder package class, and instantiate
    sounder_package_ref = sounder_packages.sounder_package_ref[datatype]
    sounder_obj = sounder_package_ref(*sounder_args, **sounder_kwargs)

    LOG.debug("sounder_obj.pres_0 = {}".format(sounder_obj.pres_0))
    # LOG.debug("sounder_obj.elev_0 = {}".format(sounder_obj.elev_0))

    pres_0 = sounder_obj.pres_0
    # elev_0 = sounder_obj.elev_0
    lats = sounder_obj.datasets['lat']['data']
    lons = sounder_obj.datasets['lon']['data']
    data = sounder_obj.datasets[dataset]['data']
    data_mask = sounder_obj.datasets[dataset]['data_mask']

    # if plot_type == 'slice':
        # lat_col = sounder_obj.datasets['lat_col']['data']
        # lon_col = sounder_obj.datasets['lon_col']['data']
        # pressure = sounder_obj.pressure
        # args.footprint = sounder_obj.col
        # if elevation != None:
            # elevation = sounder_obj.datasets['elev_slice']['data']

    input_file = os.path.basename(input_file_list[0])

    # Get the dataset options
    try:
        dataset_options = sounder_image_data.Dataset_Options.data[dataset]
    except KeyError:
        dataset_options = sounder_image_data.Dataset_Options.data['unknown']
        dataset_options['name'] = dataset

    for key in dataset_options.keys():
        LOG.debug("dataset_options['{}'] = {}".format(key,dataset_options[key]))

    # Determine the filename
    if dataset=='ctp' or dataset=='ctt' or plot_type=='slice':
        LOG.info("pres_0 = {}".format(pres_0))
        LOG.info("elev_0 = {}".format(elev_0))
        vertical_type = "elev" if elev_0 != None else "pres"
        LOG.info("vertical_type = {}".format(vertical_type))
        file_suffix = "{}_{}_{}".format(datatype,dataset,vertical_type)
    else:
        if pres_0 != None:
            file_suffix = "{}_{}_{}mb".format(datatype,dataset,int(pres_0))
        else:
            file_suffix = "{}_{}_{}ft".format(datatype,dataset,int(elev_0))

    if output_file==None and outputFilePrefix==None :
        output_file = "{}.{}.png".format(input_file,file_suffix)
    if output_file!=None and outputFilePrefix==None :
        pass
    if output_file==None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)
    if output_file!=None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)

    dataset = args.dataset

    # Set the navigation
    # plot_nav_options = set_plot_navigation(lats,lons,sounder_obj, args)
    # for key in plot_nav_options.keys():
        # LOG.info("plot_nav_options['{}'] = {}".format(key,plot_nav_options[key]))

    # Set the plot styles
    plot_style_options = set_plot_styles(sounder_obj, dataset, dataset_options, args)

    # Get pointers to the desired plotting routines
    plot_image = plot_style_options['plot_image']
    plot_map = plot_style_options['plot_map']
    plot_slice = plot_style_options['plot_slice']

    #plot_style_options['version'] = cspp_geo_version
    # plot_style_options['yMin'] = yMin
    # plot_style_options['yMax'] = yMax

    plot_options = {}
    #plot_options['title'] = "{}\n\n".format(input_file)
    #plot_options['cbar_title'] = cbar_titles[dataset]
    #plot_options['units'] = sounding_inputs[dataset]['units']
    #plot_options['stride'] = stride
    plot_options['lat_0'] = lat_0
    plot_options['lon_0'] = lon_0
    #plot_options['pres_0'] = sounding_inputs['pres_0']
    # plot_options['latMin'] = latMin
    # plot_options['lonMin'] = lonMin
    # plot_options['latMax'] = latMax
    # plot_options['lonMax'] = lonMax
    plot_options['bounding_lat'] = bounding_lat
    plot_options['extent'] = extent
    # plot_options['yMin'] = yMin
    # plot_options['yMax'] = yMax
    #plot_options['plotMin'] = plotMin
    #plot_options['plotMax'] = plotMax
    # plot_options['scale'] = scale
    #plot_options['map_res'] = map_res
    plot_options['proj'] = proj
    plot_options['pkg_name'] = sounder_obj.pkg_name
    #plot_options['cmap'] = cmap
    #plot_options['scatterPlot'] = doScatterPlot
    #plot_options['pointSize'] = pointSize
    #plot_options['dpi'] = dpi

    # Create the plot
    if plot_type == 'image':
        retval = plot_map(lats, lons, data, data_mask,
                output_file, dataset_options, plot_style_options, plot_options)
    # if plot_type == 'slice':
        # retval = plot_slice(lat_col, lon_col, lats, lons, pressure, elevation,
                # data, data_mask, output_file, dataset_options, plot_style_options, plot_options)

    print("")
    return retval


if __name__=='__main__':
    sys.exit(main())
