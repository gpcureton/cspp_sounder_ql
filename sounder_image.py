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

#import scipy.spatial as ss

from netCDF4 import Dataset
from netCDF4 import num2date
import h5py

import sounder_image_data
from thermo import dewhum
from ql_common import granuleFiles
#from ql_common import get_pressure_index
from ql_common import plotMapDataContinuous,plotMapDataDiscrete
from ql_common import set_plot_styles
from ql_common import set_plot_navigation_bm as set_plot_navigation

import HSRTV

# every module should have a LOG object
LOG = logging.getLogger(__file__)


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    dataChoices=['IAPP','MIRS','HSRTV','NUCAPS']
    prodChoices=['temp','temp_gdas','wvap','wvap_gdas','dwpt','relh','relh_gdas',
                 'ctp','ctt','2temp','cold_air_aloft']
    map_res_choice = ['c','l','i']
    plot_type_choice = ['image','slice']
    map_proj_choice = {
                       'cea':     'Cylindrical Equal Area',
                       'lcc':     'Lambert Conformal',
                       'ortho':   'Orthographic',
                       'geos':    'Geostationary',
                       'npstere': 'North-Polar Stereographic',
                       'spstere': 'South-Polar Stereographic',
                       'cyl':     'Cylindrical Equidistant',
                       'sinu':    'Sinusoidal',
                       'merc':    'Mercator'
                       #'mbtfpq':  'McBryde-Thomas Flat-Polar Quartic',
                       #'poly':    'Polyconic',
                       #'moll':    'Mollweide',
                       #'tmerc':   'Transverse Mercator',
                       #'gall':    'Gall Stereographic Cylindrical',
                       #'mill':    'Miller Cylindrical',
                       #'stere':   'Stereographic',
                       #'eqdc':    'Equidistant Conic',
                       ##'rotpole': 'Rotated Pole',
                       #'hammer':  'Hammer',
                       #'nsper':   'Near-Sided Perspective',
                       #'eck4':    'Eckert IV',
                       #'aea':     'Albers Equal Area',
                       #'kav7':    'Kavrayskiy VII',
                       #'nplaea':  'North-Polar Lambert Azimuthal',
                       #'splaea':  'South-Polar Lambert Azimuthal',
                       #'npaeqd':  'North-Polar Azimuthal Equidistant',
                       #'spaeqd':  'South-Polar Azimuthal Equidistant',
                       #'cass':    'Cassini-Soldner',
                       #'laea':    'Lambert Azimuthal Equal Area',
                       #'robin':   'Robinson'
                       }


    #defaults = {
                #'input_file':None,
                #'datatype':None,
                #'dataset':'temp',
                #'stride':1,
                #'pressure':850.,
                #'lat_0':None,
                #'lat_0':None,
                #'plotMin'  : None,
                #'plotMax'  : None,
                #'bounding_lat':0.,
                #'scatter_plot':False,
                #'pointSize':1,
                #'scale':1,
                #'map_res':'c',
                #'proj':'lcc',
                #'output_file':None,
                #'outputFilePrefix' : None,
                #'dpi':200
                #}

    defaults = {
                'image_size' : [7.5, 7.5],
                'cbar_axis' : [0.10 , 0.05, 0.8 , 0.05],
                'cmap':None,
                'dataset':'temp',
                'datatype':None,
                'dpi':200,
                'font_scale':1.,
                'input_file':None,
                'bounding_lat':0.,
                'lat_0':None,
                'lon_0':None,
                'llcrnrx':None,
                'llcrnry':None,
                'urcrnrx':None,
                'urcrnry':None,
                'map_axis' : [0.10, 0.15, 0.80, 0.8],
                'map_res':'c',
                'output_file':None,
                'outputFilePrefix' : None,
                'plotMax'  : None,
                'plotMin'  : None,
                'plot_type'  : 'image',
                'pointSize':1,
                'pressure':850.,
                'proj':'lcc',
                'scale':1,
                'scatter_plot':False,
                'stride':1,
                'unnavigated':False
                }



    description = '''Create a plot of temperature, dewpoint or something 
                     else at a particular pressure level. Supports IAPP, MIRS, HSRTV 
                     and NUCAPS files.'''

    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser(
                                     #version=version,
                                     description=description
                                     )

    # Mandatory/positional arguments
    
    parser.add_argument(
                      action='store',
                      dest='input_file',
                      type=str,
                      help='''The fully qualified path to the input file(s). May 
                              be a file glob (which should be in quotes).'''
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

    parser.add_argument('--bounding_lat',
                      action="store",
                      dest="bounding_lat",
                      default=defaults["bounding_lat"],
                      type=float,
                      help='''The minimum/maximum latitude to plot for the polar projections. 
                      [default: {}]'''.format(defaults["bounding_lat"])
                      )

    parser.add_argument('--viewport',
                      action="store",
                      dest="viewport",
                      type=float,
                      nargs=4,
                      metavar=('LLCRNRX', 'LLCRNRY', 'URCRNRX', 'URCRNRY'),
                      help="""Lower-left and upper-right coordinates 
                      [*llcrnrx*, *llcrnry*, *urcrnrx*, *urcrnry*] of the projection 
                      viewport, in the range [-0.5,+0.5] (for navigated plots only)"""
                      )

    parser.add_argument('--image_size',
                      action="store",
                      dest="image_size",
                      default=defaults["image_size"],
                      type=float,
                      nargs=2,
                      metavar=('WIDTH', 'HEIGHT'),
                      help="""The size of the output image [*width*, *height*]
                      in inches. [default: '{}']""".format(defaults["image_size"])
                      )

    parser.add_argument('--map_axis',
                      action="store",
                      dest="map_axis",
                      default=defaults["map_axis"],
                      type=float,
                      nargs=4,
                      metavar=('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
                      help="""Set the map axes at position [*left*, *bottom*, *width*, *height*] 
                      where all quantities are in fractions of figure width and height. 
                      [default: '{}']""".format(defaults["map_axis"])
                      )

    parser.add_argument('--cbar_axis',
                      action="store",
                      dest="cbar_axis",
                      default=defaults["cbar_axis"],
                      type=float,
                      nargs=4,
                      metavar=('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
                      help="""Set the colorbar axes at position [*left*, *bottom*, *width*, *height*] 
                      where all quantities are in fractions of figure width and height. 
                      [default: '{}']""".format(defaults["cbar_axis"])
                      )

    parser.add_argument('--scatter_plot',
                      action="store_true",
                      dest="doScatterPlot",
                      default=defaults["scatter_plot"],
                      help="Generate the plot using a scatterplot approach."
                      )

    parser.add_argument('--logscale',
                      action="store_true",
                      dest="logscale",
                      help="""Plot the dataset using a logarithmic scale."""
                      )

    parser.add_argument('--no_logscale',
                      action="store_true",
                      dest="no_logscale",
                      help="""Plot the dataset using a linear scale."""
                      )

    parser.add_argument('-P','--pointSize',
                      action="store",
                      dest="pointSize",
                      default=defaults["pointSize"],
                      type=float,
                      help='''Size of the plot point used to represent each pixel. 
                      [default: {}]'''.format(defaults["pointSize"])
                      )

    parser.add_argument('--font_scale',
                      action="store",
                      dest="font_scale",
                      default=defaults["font_scale"],
                      type=float,
                      help='''The scale factor to apply to the default font size
                      for the plot labels. [default: {}]'''.format(defaults["font_scale"])
                      )

    parser.add_argument('--cmap',
                      action="store",
                      dest="cmap",
                      default=defaults["cmap"],
                      type=str,
                      help="""The matplotlib colormap to use. 
                      [default: '{}']""".format(defaults["cmap"])
                      )

    parser.add_argument('--plot_type',
                      action="store",
                      dest="plot_type",
                      default=defaults["plot_type"],
                      type=str,
                      choices=plot_type_choice,
                      help='''The type of plot. Possible values are {}.'''.format(
                          plot_type_choice.__str__()[1:-1])
                      )

    parser.add_argument('--plot_title',
                      action="store",
                      dest="plot_title",
                      type=str,
                      help='''The plot title. Must be placed in double quotes.
                      '''
                      )

    parser.add_argument('--cbar_title',
                      action="store",
                      dest="cbar_title",
                      type=str,
                      help='''The colourbar title. Must be placed in double quotes.
                      '''
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

    parser.add_argument('--proj',
                      action="store",
                      dest="proj",
                      default=defaults["proj"],
                      type=str,
                      choices=map_proj_choice,
                      help='''The map projection. Possible values are 
                      {{{}}}. [default: '{}']'''.format(
                          "".join(["'{0:8s} ({1}), ".format(tups[0]+"'",tups[1]) for tups in zip(map_proj_choice.keys(),map_proj_choice.values())]),
                          defaults["proj"]
                          )
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
                      png name. 
                      [default: {}]""".format(defaults["outputFilePrefix"])
                      )

    parser.add_argument("-v", "--verbose",
                      dest='verbosity',
                      action="count", 
                      default=2,
                      help='''each occurrence increases verbosity 1 level from 
                      INFO. -v=DEBUG'''
                      )

    parser.add_argument("-q", "--quiet",
                      action="store_true",
                      dest='quiet',
                      default=False,
                      help='''Silence all console output'''
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
    plot_type = options.plot_type
    scale = options.scale
    doScatterPlot = options.doScatterPlot
    pointSize = options.pointSize
    proj = options.proj
    map_res = options.map_res
    output_file  = options.output_file
    outputFilePrefix  = options.outputFilePrefix
    dpi = options.dpi

    # Obtain matched lists of the geolocation and product files
    input_file_list = granuleFiles(input_file)

    if input_file_list==[]:
        LOG.warn('No files match the input file glob "{}", aborting.'
                .format(input_file))
        sys.exit(1)


    LOG.info("input file(s): {}".format(input_file_list))

    # Read in the input file, and return a dictionary containing the required
    # data

    dataChoices=['IAPP','MIRS','HSRTV','NUCAPS']

    LOG.info("Input pressure = {}".format(pressure))
    
    hsrtv_obj = HSRTV.HSRTV(input_file_list,dataset,plot_type,pres_0=pressure)

    pres_0 = hsrtv_obj.pres_0
    lats = hsrtv_obj.datasets['lat']['data']
    lons = hsrtv_obj.datasets['lon']['data']
    data = hsrtv_obj.datasets[dataset]['data']
    data_mask = hsrtv_obj.datasets[dataset]['data_mask']

    if plot_type == 'slice':
        lat_col = hsrtv_obj.datasets['lat_col']['data']
        lon_col = hsrtv_obj.datasets['lon_col']['data']

    LOG.debug(hsrtv_obj.datasets['file_attrs'].items())

    input_file = path.basename(input_file_list[0])

    #print 'temp = np.array({})'.format(hsrtv_obj.datasets['temp']['data'][:,0])
    #print 'wvap = np.array({})'.format(hsrtv_obj.datasets['wvap']['data'][:,0])
    #print 'temp_gdas = np.array({})'.format(hsrtv_obj.datasets['temp_gdas']['data'][:,0])
    #print 'relh_gdas = np.array({})'.format(hsrtv_obj.datasets['relh_gdas']['data'][:,0])
    #print "Hello"

    # Get the dataset options
    try:
        dataset_options = sounder_image_data.Dataset_Options.data[dataset]
    except KeyError:
        dataset_options = sounder_image_data.Dataset_Options.data['unknown']
        dataset_options['name'] = dataset

    for key in dataset_options.keys():
        LOG.debug("dataset_options['{}'] = {}".format(key,dataset_options[key]))


    # Determine the filename
    if dataset=='ctp' or dataset=='ctt':
        file_suffix = "{}_{}".format(datatype,dataset)
    else:
        file_suffix = "{}_{}_{}mb".format(datatype,dataset,int(pres_0))

    if output_file==None and outputFilePrefix==None :
        output_file = "{}.{}.png".format(input_file,file_suffix)
    if output_file!=None and outputFilePrefix==None :
        pass
    if output_file==None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)
    if output_file!=None and outputFilePrefix!=None :
        output_file = "{}_{}.png".format(outputFilePrefix,file_suffix)

    dataset = options.dataset

    # Set the navigation 
    plot_nav_options = set_plot_navigation(lats,lons,hsrtv_obj,options)
    for key in plot_nav_options.keys():
        LOG.info("plot_nav_options['{}'] = {}".format(key,plot_nav_options[key]))


    # Set the plot styles 
    plot_style_options = set_plot_styles(hsrtv_obj,dataset,dataset_options,options)

    # Get pointers to the desired plotting routines
    plot_image = plot_style_options['plot_image']
    plot_map = plot_style_options['plot_map']
    plot_slice = plot_style_options['plot_slice']

    #plot_style_options['version'] = cspp_geo_version

    plot_options = {}
    #plot_options['title'] = "{}\n\n".format(input_file)
    #plot_options['cbar_title'] = cbar_titles[dataset]
    #plot_options['units'] = sounding_inputs[dataset]['units']
    #plot_options['stride'] = stride
    plot_options['lat_0'] = lat_0
    plot_options['lon_0'] = lon_0
    #plot_options['pres_0'] = sounding_inputs['pres_0']
    plot_options['latMin'] = latMin
    plot_options['lonMin'] = lonMin
    plot_options['latMax'] = latMax
    plot_options['lonMax'] = lonMax
    plot_options['bounding_lat'] = bounding_lat
    #plot_options['plotMin'] = plotMin
    #plot_options['plotMax'] = plotMax
    plot_options['scale'] = scale
    #plot_options['map_res'] = map_res
    plot_options['proj'] = proj
    #plot_options['cmap'] = cmap
    #plot_options['scatterPlot'] = doScatterPlot
    #plot_options['pointSize'] = pointSize
    #plot_options['dpi'] = dpi

    # Create the plot
    if plot_type == 'image':
        retval = plot_map(lats, lons, data, data_mask, 
                output_file, dataset_options, plot_style_options,plot_options)
    if plot_type == 'slice':
        retval = plot_slice(lat_col, lon_col, lats, lons, data, data_mask,
                output_file, dataset_options, plot_style_options,plot_options)

    print ""
    return retval


if __name__=='__main__':
    sys.exit(main())  
