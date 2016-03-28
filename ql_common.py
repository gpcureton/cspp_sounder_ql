#!/usr/bin/env python
# encoding: utf-8
"""
ql_common.py

Purpose: Common methods for datafile I/O, plots etc...

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2015-03-04.
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

from scipy import interpolate
import scipy.spatial as ss

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.ticker import LogLocator,LogFormatter
from matplotlib.mlab import griddata
from scipy.interpolate import griddata as griddata_sp

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

from mpl_toolkits.basemap import Basemap

from netCDF4 import Dataset
from netCDF4 import num2date
import h5py

#from basemap_utils import drawparallels,drawmeridians

# every module should have a LOG object
LOG = logging.getLogger(__file__)


def _tuple2args(parms):
    s = ' '.join( '+%s=%s' % (k,v) for (k,v) in parms )
    return s.encode('ascii')


class Datafile_NetCDF():

    def __init__(self,input_file):


        self.input_file = input_file

        if os.path.exists(self.input_file):
            LOG.debug('Opening {} with NetCDF...'.format(self.input_file))
            self.file_obj = Dataset(self.input_file)
        else:
            LOG.error('Input file {} does not exist, aborting...'.format(self.input_file))
            sys.exit(1)

        # Dictionary of file object attributes
        self.attrs = {}
        for attr_key in self.file_obj.ncattrs():
            self.attrs[attr_key] = getattr(self.file_obj,attr_key)
        
        # Ordered dictionary of dataset objects
        self.data_dict = self.file_obj.variables

        # List of dataset names
        self.datanames = self.data_dict.keys()
        self.datanames.sort()


    class Dataset():

        def __init__(selfd,L1_obj,dataname):

            selfd.dataname = dataname
            LOG.debug("selfd.dataname = {}".format(selfd.dataname))

            selfd.dset_obj = L1_obj.file_obj.variables[dataname]

            selfd.attrs = {}
            for attr_key in selfd.dset_obj.ncattrs():
                selfd.attrs[attr_key] = getattr(selfd.dset_obj,attr_key)
            
            selfd.dset = ma.masked_equal(selfd.dset_obj[:],selfd.attrs['_FillValue'])
            #selfd.dset = selfd.dset * selfd.attrs['scale_factor'] + selfd.attrs['add_offset']

    def close(self):
        LOG.debug('Closing {}...'.format(self.input_file))
        self.file_obj.close()


class Datafile_HDF5():

    def __init__(self,input_file):


        self.input_file = input_file

        if h5py.is_hdf5(self.input_file):
            LOG.debug('Opening {} with HDF5...'.format(self.input_file))
            self.file_obj = h5py.File(self.input_file,'r')
        else:
            LOG.error('Input file {} does not exist, aborting...'.format(self.input_file))
            sys.exit(1)

        # Dictionary of file object attributes
        self.attrs = {}
        for attr_key in self.file_obj.attrs.keys():
            self.attrs[attr_key] = self.file_obj.attrs[attr_key][0]
        
        # Dictionary of dataset objects
        self.data_dict = {}
        for dset_obj in self.file_obj.values():
            key = dset_obj.name
            self.data_dict[key] = dset_obj

        # List of dataset names
        self.datanames = self.file_obj.keys()
        self.datanames.sort()


    class Dataset():

        def __init__(selfd,dfile_obj,dataname):

            selfd.dataname = dataname
            LOG.debug("selfd.dataname = {}".format(selfd.dataname))

            selfd.dset_obj = dfile_obj.data_dict[dataname]

            selfd.attrs = {}
            for attr_key in selfd.dset_obj.attrs.keys():
                selfd.attrs[attr_key] = selfd.dset_obj.attrs[attr_key][0]
            
            if 'missing_value' in selfd.attrs.keys():
                selfd.dset = ma.masked_equal(selfd.dset_obj[:],
                        float(selfd.attrs['missing_value']))
            else:
                selfd.dset = selfd.dset_obj[:]


    def close(self):
        LOG.debug('Closing {}...'.format(self.input_file))
        self.file_obj.close()


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


def get_pressure_index(pressure,pres_0,kd_dist=10.):
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


def get_pressure_level(pressure,pres_0):
    '''
    Calculate the closest pressure level to the desired pressure. 
    '''
    pressure_scope = 10.
    LOG.debug("Scope of pressure level search is {} hPa".format(pressure_scope))

    level = get_pressure_index(pressure, pres_0, kd_dist=pressure_scope)

    attempts = 1
    while (level==None and attempts<10):
        pressure_scope *= 1.1
        LOG.debug("Scope of pressure level search is {} hPa".format(pressure_scope))
        level = get_pressure_index(pressure, pres_0, kd_dist=pressure_scope)
        attempts += 1

    if level==None:
        raise Exception("No suitable pressure level found, aborting...")
         
    LOG.debug("Retrieved level = {}".format(level))

    LOG.debug("Retrieved pressure = {}".format(pressure[level]))

    return level,pressure[level]


def get_geo_indices(lat,lon,lat_0=None,lon_0=None):
    '''
    Get the indices in the lat and lon arrays closes to the desired
    latitude and longitude. If no lat and lon are specified, try to select 
    the central indices. Uses kd-tree.
    '''

    try:
        geo_shape = lat.shape

        nrows, ncols = geo_shape[0],geo_shape[1]
        LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
        LOG.debug("lat_0,lon_0= ({},{})".format(lat_0,lon_0))

        if (lat_0==None) or (lon_0==None):
            # Non lat/lon pair given, use the central unmasked values.
            col_offset = 0
            row_idx = int(np.floor(nrows/2.))
            col_idx = int(np.floor(ncols/2.)) - col_offset
            lat_0 = lat[row_idx,col_idx] if lat_0==None else lat_0
            lon_0 = lon[row_idx,col_idx] if lon_0==None else lon_0
            LOG.info("No lat/lon pair given, using ({:4.2f},{:4.2f})".
                    format(lat_0,lon_0))

        lat_mask = lat.mask if ma.is_masked(lat) else np.zeros(lat.shape,dtype='bool')
        lon_mask = lon.mask if ma.is_masked(lon) else np.zeros(lon.shape,dtype='bool')
        geo_mask = lat_mask + lon_mask

        lat_masked = ma.array(lat,mask=geo_mask).ravel()
        lon_masked = ma.array(lon,mask=geo_mask).ravel()
        lat_masked = lat_masked.compressed()
        lon_masked = lon_masked.compressed()

        rows, cols = np.mgrid[0:nrows, 0:ncols]
        cols_masked = ma.array(cols,mask=geo_mask).ravel()
        rows_masked = ma.array(rows,mask=geo_mask).ravel()
        cols_masked = cols_masked.compressed()
        rows_masked = rows_masked.compressed()

        LOG.debug("There are {} masked latitude values".format(np.sum(lat_mask)))
        LOG.debug("There are {} masked longitude values".format(np.sum(lon_mask)))

        LOG.info("Searching for lat,lon coordinate ({:4.2f},{:4.2f})".\
                format(lat_0,lon_0))

        # Construct the kdtree
        points = zip(lon_masked, lat_masked)
        tree = ss.KDTree(points)

        # get index of each point which distance from (lon_0,lat_0) is under 1
        a = tree.query_ball_point([lon_0, lat_0], 1.0)
        if a == []:
            LOG.error("Specified coordinates ({:4.2f},{:4.2f}) not found...".\
                    format(lat_0,lon_0))
            return None,None

        row = [rows_masked[i] for i in a][0]
        col = [cols_masked[i] for i in a][0]

        LOG.debug("Row index is {}".format(row))
        LOG.debug("Col index is {}".format(col))
        
        lat_0 = lat[row,col].squeeze()
        lon_0 = lon[row,col].squeeze()

        LOG.debug("Retrieved latitude is {:4.2f}".format(lat_0))
        LOG.debug("Retrieved longitude is {:4.2f}".format(lon_0))

    except Exception, err :
        LOG.warn("{}".format(err))
        LOG.debug(traceback.format_exc())

    return row,col,lat_0,lon_0


def set_plot_navigation_proj(lats,lons,goes_l1_obj, options):
    """
    Collects the various navigation options and does any required tweaking
    before passing to the plotting method.

    Uses Proj to generate the map projection.
    """
    # Here are some attributes from the geocat Level-1 files...
    from mpl_toolkits.basemap.pyproj import Proj

    nrows,ncols = lats.shape[0],lats.shape[1]
    nrows_div2,ncols_div2 = np.floor_divide(nrows,2),np.floor_divide(ncols,2)
    LOG.debug('lats.shape = {}'.format(lats.shape))
    LOG.debug('nrows,ncols = ({},{})'.format(nrows,ncols))
    LOG.debug('nrows_div2,ncols_div2 = ({},{})'.format(nrows_div2,ncols_div2))

    corners = [(0,0), (0,-1), (-1,-1), (-1,0)]

    for crnr,crnr_idx in zip(range(4),corners):
        LOG.debug("({}) (lon,lat):({},{})".format(crnr,lons[crnr_idx],lats[crnr_idx]))


    subsatellite_longitude = goes_l1_obj.attrs['Subsatellite_Longitude']
    LOG.debug(' File subsatellite_Longitude = {}'.format(subsatellite_longitude))

    plot_nav_options = {}

    p1 = Proj("+proj=geos +h=35774290 +a= 6378137 +b= 6378137 +lon_0=-75 +units=meters +no_defs")

    # Western edge
    x,y = p1(lons[:,0],lats[:,0])
    west_x = np.average(x.compressed())
    # Eastern edge
    x,y = p1(lons[:,-1],lats[:,-1])
    east_x = np.average(x.compressed())
    # Northern edge
    x,y = p1(lons[0,:],lats[0,:])
    north_y = np.average(y.compressed())
    # Southern edge
    x,y = p1(lons[-1,:],lats[-1,:])
    south_y = np.average(y.compressed())

    corner_x_fallback = [west_x,east_x,east_x,west_x]
    corner_y_fallback = [north_y,north_y,south_y,south_y]

    LOG.debug("(west_x ,east_x ):({:10.1f},{:10.1f})".format(west_x,east_x))
    LOG.debug("(north_y,south_y):({:10.1f},{:10.1f})".format(north_y,south_y))

    crnr_x_names_proj = ['ulcrnrx_proj','urcrnrx_proj','lrcrnrx_proj','llcrnrx_proj']
    crnr_y_names_proj = ['ulcrnry_proj','urcrnry_proj','lrcrnry_proj','llcrnry_proj']

    # The maximum extent of the full disk in the x and y directions
    plot_nav_options['extent_x'] = p1(subsatellite_longitude+81.,0)[0]*2.
    plot_nav_options['extent_y'] = p1(subsatellite_longitude,81.)[1]*2.
    LOG.debug("(extent_x,extent_y):({},{})".format(plot_nav_options['extent_x'],plot_nav_options['extent_y']))

    # Generate the corner-origin coordinates from Basemap object...
    for crnr in range(4):
        crnr_idx = corners[crnr]
        lon,lat = lons[crnr_idx], lats[crnr_idx]
        if ma.is_masked(lon) and ma.is_masked(lat):
            x,y = corner_x_fallback[crnr],corner_y_fallback[crnr]
        else:
            x,y = p1(lon, lat)
        #if crnr_x_names_proj[crnr] is not None:
        plot_nav_options[crnr_x_names_proj[crnr]] = x
        #if crnr_y_names_proj[crnr] is not None:
        plot_nav_options[crnr_y_names_proj[crnr]] = y
        LOG.debug(
                "({}) (lon,lat):({},{}), (x,y): ({:10.1f},{:10.1f})".format(crnr,lon,lat,x,y))

    # Default plot options
    plot_nav_options['llcrnrx'] = options.llcrnrx
    plot_nav_options['llcrnry'] = options.llcrnry
    plot_nav_options['urcrnrx'] = options.urcrnrx
    plot_nav_options['urcrnry'] = options.urcrnry


    return plot_nav_options


def set_plot_navigation_bm(lats,lons,dfile_obj, options):
    """
    Collects the various navigation options and does any required tweaking
    before passing to the plotting method.

    Uses Basemap to generate the map projection.
    """

    nrows,ncols = lats.shape[0],lats.shape[1]
    nrows_div2,ncols_div2 = np.floor_divide(nrows,2),np.floor_divide(ncols,2)
    LOG.debug('lats.shape = {}'.format(lats.shape))
    LOG.debug('nrows,ncols = ({},{})'.format(nrows,ncols))
    LOG.debug('nrows_div2,ncols_div2 = ({},{})'.format(nrows_div2,ncols_div2))

    corners = [(0,0), (0,-1), (-1,-1), (-1,0)]

    for crnr,crnr_idx in zip(range(4),corners):
        LOG.debug("({}) (lon,lat):({},{})".format(crnr,lons[crnr_idx],lats[crnr_idx]))

    return {}

    #subsatellite_longitude = goes_l1_obj.attrs['Subsatellite_Longitude']
    #LOG.info('File subsatellite_Longitude = {:6.1f}'.format(subsatellite_longitude))

    m1 = Basemap(projection='geos',lon_0=subsatellite_longitude,resolution=None)

    plot_nav_options = {}

    # The maximum extent of the full disk in the x and y directions
    plot_nav_options['extent_x'] = m1.urcrnrx
    plot_nav_options['extent_y'] = m1.urcrnry
    LOG.debug("(extent_x,extent_y):({},{})".format(plot_nav_options['extent_x'],plot_nav_options['extent_y']))

    # Compute coordinates of sub-satellite point
    x_subsat,y_subsat = m1(subsatellite_longitude,0)
    LOG.debug("(x_subsat,y_subsat):({},{})".format(x_subsat,y_subsat))

    is_full_disk = False
    if ma.is_masked(lats[corners[0]]) and ma.is_masked(lats[corners[1]]) and \
       ma.is_masked(lats[corners[2]]) and ma.is_masked(lats[corners[3]]):
        is_full_disk = True

    if options.region == "FD":
        is_full_disk = True

    if is_full_disk:

        LOG.info("This image is full disk")

        plot_nav_options['urcrnrx_map'] = plot_nav_options['extent_x'] - x_subsat
        plot_nav_options['urcrnry_map'] = plot_nav_options['extent_y'] - y_subsat
        plot_nav_options['llcrnrx_map'] = 0. - x_subsat
        plot_nav_options['llcrnry_map'] = 0  - y_subsat

    else:
        LOG.info("This image is NOT full disk")

        # Get the north, south, east and west edges of the plot region (in map
        # coordinates). Check that the edge isnn't all missing, and step towards 
        # center as necessary.

        # Western edge
        we_col = 0
        while True:
            LOG.debug("Checking edge for we_col = {}".format(we_col))
            x,y = m1(lons[:,we_col],lats[:,we_col])
            x_compressed = x.compressed()
            #LOG.debug("x_compressed:{}".format(x_compressed))
            if (len(x_compressed) == 0):
                LOG.debug("No Disk! ({})".format(we_col))
                we_col = we_col + 1
            else : 
                LOG.debug("Checking complete: we_col = {}".format(we_col))
                break

        west_x = np.average(x_compressed)

        # Eastern edge

        ee_col = -1
        while True:
            LOG.debug("Checking edge for ee_col = {}".format(ee_col))
            x,y = m1(lons[:,ee_col],lats[:,ee_col])
            x_compressed = x.compressed()
            #LOG.debug("x_compressed:{}".format(x_compressed))
            if (len(x_compressed) == 0):
                LOG.debug("No Disk! ({})".format(ee_col))
                ee_col = ee_col - 1
            else : 
                LOG.debug("Checking complete: ee_col = {}".format(ee_col))
                break

        east_x = np.average(x_compressed)

        # Northern edge
        ne_row = 0
        while True:
            LOG.debug("Checking edge for ne_row = {}".format(ne_row))
            x,y = m1(lons[ne_row,:],lats[ne_row,:])
            y_compressed = y.compressed()
            #LOG.debug("y_compressed:{}".format(y_compressed))
            if (len(y_compressed) == 0):
                LOG.debug("No Disk! ({})".format(ne_row))
                ne_row = ne_row + 1
            else : 
                LOG.debug("Checking complete: ne_row = {}".format(ne_row))
                break

        north_y = np.average(y_compressed)

        # Southern edge
        se_row = -1
        while True:
            LOG.debug("Checking edge for se_row = {}".format(se_row))
            x,y = m1(lons[se_row,:],lats[se_row,:])
            y_compressed = y.compressed()
            #LOG.debug("y_compressed:{}".format(y_compressed))
            if (len(y_compressed) == 0):
                LOG.debug("No Disk! ({})".format(se_row))
                se_row = se_row - 1
            else : 
                LOG.debug("Checking complete: se_row = {}".format(se_row))
                break

        south_y = np.average(y_compressed)


        corner_x_fallback = [west_x,east_x,east_x,west_x]
        corner_y_fallback = [north_y,north_y,south_y,south_y]

        LOG.debug("(west_x ,east_x ):({:10.1f},{:10.1f})".format(west_x,east_x))
        LOG.debug("(north_y,south_y):({:10.1f},{:10.1f})".format(north_y,south_y))

        crnr_x_names = ['ulcrnrx','urcrnrx','lrcrnrx','llcrnrx']
        crnr_y_names = ['ulcrnry','urcrnry','lrcrnry','llcrnry']

        # Generate the center-origin coordinates from Proj object...
        for crnr in range(4):
            crnr_idx = corners[crnr]
            lon,lat = lons[crnr_idx], lats[crnr_idx]
            if ma.is_masked(lon) and ma.is_masked(lat):
                x,y = corner_x_fallback[crnr],corner_y_fallback[crnr]
            else:
                x,y = m1(lon, lat)
            plot_nav_options[crnr_x_names[crnr]] = x
            plot_nav_options[crnr_y_names[crnr]] = y
            LOG.debug(
                    "({}) (lon,lat):({},{}), (x,y): ({:10.1f},{:10.1f})".format(crnr,lon,lat,x,y))

        plot_nav_options['ulcrnrx_map'] = plot_nav_options['ulcrnrx'] - x_subsat
        plot_nav_options['ulcrnry_map'] = plot_nav_options['ulcrnry'] - y_subsat
        plot_nav_options['urcrnrx_map'] = plot_nav_options['urcrnrx'] - x_subsat
        plot_nav_options['urcrnry_map'] = plot_nav_options['urcrnry'] - y_subsat
        plot_nav_options['lrcrnrx_map'] = plot_nav_options['lrcrnrx'] - x_subsat
        plot_nav_options['lrcrnry_map'] = plot_nav_options['lrcrnry'] - y_subsat
        plot_nav_options['llcrnrx_map'] = plot_nav_options['llcrnrx'] - x_subsat
        plot_nav_options['llcrnry_map'] = plot_nav_options['llcrnry'] - y_subsat


    plot_nav_options['llcrnrx'] = plot_nav_options['llcrnrx_map'] if options.viewport is None else options.viewport[0] * plot_nav_options['extent_x']
    plot_nav_options['llcrnry'] = plot_nav_options['llcrnry_map'] if options.viewport is None else options.viewport[1] * plot_nav_options['extent_y'] 
    plot_nav_options['urcrnrx'] = plot_nav_options['urcrnrx_map'] if options.viewport is None else options.viewport[2] * plot_nav_options['extent_x'] 
    plot_nav_options['urcrnry'] = plot_nav_options['urcrnry_map'] if options.viewport is None else options.viewport[3] * plot_nav_options['extent_y'] 

    plot_nav_options['lon_0'] = subsatellite_longitude
    plot_nav_options['is_full_disk'] = is_full_disk

    return plot_nav_options


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


def plotMapDataContinuous(lat, lon, data, data_mask, pngName,
        dataset_options, plot_style_options, plot_options):
        
    # Copy the plot options to local variables
    #title         = plot_options['title']
    #cbar_title    = plot_options['cbar_title']
    #units         = plot_options['units']
    #stride        = plot_options['stride']
    lat_0         = plot_options['lat_0']
    lon_0         = plot_options['lon_0']
    #pres_0        = plot_options['pres_0']
    latMin        = plot_options['latMin']
    lonMin        = plot_options['lonMin']
    latMax        = plot_options['latMax']
    lonMax        = plot_options['lonMax']
    #plotMin       = plot_options['plotMin']
    #plotMax       = plot_options['plotMax']
    bounding_lat  = plot_options['bounding_lat']
    scale         = plot_options['scale']
    #map_res       = plot_options['map_res']
    proj          = plot_options['proj']
    #cmap          = plot_options['cmap']
    #doScatterPlot = plot_options['scatterPlot']
    #pointSize     = plot_options['pointSize']
    #dpi           = plot_options['dpi']

    # Copy the plot options to local variables
    #llcrnrx        = plot_nav_options['llcrnrx']
    #llcrnry        = plot_nav_options['llcrnry']
    #urcrnrx        = plot_nav_options['urcrnrx']
    #urcrnry        = plot_nav_options['urcrnry']

    #lon_0         = plot_nav_options['lon_0']

    title         = plot_style_options['title']
    cbar_title    = plot_style_options['cbar_title']
    stride        = plot_style_options['stride']
    plotMin       = plot_style_options['plotMin']
    plotMax       = plot_style_options['plotMax']
    plotLims      = plot_style_options['plotLims']

    map_axis      = plot_style_options['map_axis']
    cbar_axis     = plot_style_options['cbar_axis']
    image_size    = plot_style_options['image_size']

    map_res       = plot_style_options['map_res']
    cmap          = plot_style_options['cmap']
    doScatterPlot = plot_style_options['scatterPlot']
    pointSize     = plot_style_options['pointSize']
    font_scale    = plot_style_options['font_scale']
    log_plot      = plot_style_options['log_plot']
    dpi           = plot_style_options['dpi']

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
        lat_0 = 0. if lat_0==None else lat_0
        lon_0 = 0. if lon_0==None else lon_0
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
    if (lat_0==None) or (lon_0==None):
        geo_shape = lat.shape
        nrows, ncols = geo_shape[0],geo_shape[1]
        LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
        # Non lat/lon pair given, use the central unmasked values.
        row_idx = int(nrows/2.)
        col_idx = int(ncols/2.)
        lat_0 = lat[row_idx,col_idx] if lat_0==None else lat_0
        lon_0 = lon[row_idx,col_idx] if lon_0==None else lon_0
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
            plot_kw['llcrnrlat'] =  -80. if latMin==None else latMin
            plot_kw['urcrnrlat'] =   80. if latMax==None else latMax
            plot_kw['llcrnrlon'] = -180. if lonMin==None else lonMin
            plot_kw['urcrnrlon'] =  180. if lonMax==None else lonMax
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

    LOG.info("data.shape = {}".format(data.shape))
    LOG.info("data_mask.shape = {}".format(data_mask.shape))

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

    txt = ax.set_title(title,fontsize=11,y=1.04)

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


def plotMapDataDiscrete(lat, lon, data, data_mask, pngName,
        dataset_options, plot_style_options, plot_options):
        
    # Copy the plot options to local variables
    #title         = plot_options['title']
    #cbar_title    = plot_options['cbar_title']
    #units         = plot_options['units']
    #stride        = plot_options['stride']
    lat_0         = plot_options['lat_0']
    lon_0         = plot_options['lon_0']
    #pres_0        = plot_options['pres_0']
    latMin        = plot_options['latMin']
    lonMin        = plot_options['lonMin']
    latMax        = plot_options['latMax']
    lonMax        = plot_options['lonMax']
    #plotMin       = plot_options['plotMin']
    #plotMax       = plot_options['plotMax']
    bounding_lat  = plot_options['bounding_lat']
    scale         = plot_options['scale']
    #map_res       = plot_options['map_res']
    proj          = plot_options['proj']
    #cmap          = plot_options['cmap']
    #doScatterPlot = plot_options['scatterPlot']
    #pointSize     = plot_options['pointSize']
    #dpi           = plot_options['dpi']

    # Copy the plot options to local variables
    #llcrnrx        = plot_nav_options['llcrnrx']
    #llcrnry        = plot_nav_options['llcrnry']
    #urcrnrx        = plot_nav_options['urcrnrx']
    #urcrnry        = plot_nav_options['urcrnry']

    #lon_0         = plot_nav_options['lon_0']

    title         = plot_style_options['title']
    cbar_title    = plot_style_options['cbar_title']
    stride        = plot_style_options['stride']
    plotMin       = plot_style_options['plotMin']
    plotMax       = plot_style_options['plotMax']
    plotLims      = plot_style_options['plotLims']

    map_axis      = plot_style_options['map_axis']
    cbar_axis     = plot_style_options['cbar_axis']
    image_size    = plot_style_options['image_size']

    map_res       = plot_style_options['map_res']
    cmap          = plot_style_options['cmap']
    doScatterPlot = plot_style_options['scatterPlot']
    pointSize     = plot_style_options['pointSize']
    font_scale    = plot_style_options['font_scale']
    log_plot      = plot_style_options['log_plot']
    dpi           = plot_style_options['dpi']

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
        lat_0 = 0. if lat_0==None else lat_0
        lon_0 = 0. if lon_0==None else lon_0
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
    if (lat_0==None) or (lon_0==None):
        geo_shape = lat.shape
        nrows, ncols = geo_shape[0],geo_shape[1]
        LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
        # Non lat/lon pair given, use the central unmasked values.
        row_idx = int(nrows/2.)
        col_idx = int(ncols/2.)
        lat_0 = lat[row_idx,col_idx] if lat_0==None else lat_0
        lon_0 = lon[row_idx,col_idx] if lon_0==None else lon_0
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

    # Define the discrete colorbar tick locations
    fill_colours = dataset_options['fill_colours']
    LOG.info("fill_colours = {}".format(fill_colours))
    cmap = ListedColormap(fill_colours)

    numCats = np.array(fill_colours).size
    numBounds = numCats + 1

    tickPos = np.arange(float(numBounds))/float(numCats)
    tickPos = tickPos[0 :-1] + tickPos[1]/2.

    LOG.debug('plotLims = {},{}'.format(plotLims[0],plotLims[1]))
    vmin,vmax = plotLims[0],plotLims[-1]

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
            plot_kw['llcrnrlat'] =  -80. if latMin==None else latMin
            plot_kw['urcrnrlat'] =   80. if latMax==None else latMax
            plot_kw['llcrnrlon'] = -180. if lonMin==None else lonMin
            plot_kw['urcrnrlon'] =  180. if lonMax==None else lonMax
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

    plotMin = np.min(data) if plotMin==None else plotMin
    plotMax = np.max(data) if plotMax==None else plotMax
    LOG.debug("plotMin = {}".format(plotMin))
    LOG.debug("plotMax = {}".format(plotMax))

    if doScatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none',
                vmin=plotMin,vmax=plotMax,cmap=cmap)
    else:
        cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False,
                vmin=plotMin,vmax=plotMax,cmap=cmap)

    txt = ax.set_title(title,fontsize=11,y=1.04)

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.90 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=font_scale*9)
    ppl.setp(cax.get_xticklines(),visible=False)

    # Set the colourbar tick locations and ticklabels
    #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
    tickpos_data_coords = vmax*tickPos
    cb.set_ticks(tickpos_data_coords) # In data coords (0..3)
    tick_names = dataset_options['tick_names']
    ppl.setp(cb.ax,xticklabels=tick_names)

    # Colourbar title
    cax_title = ppl.setp(cax,title=cbar_title)
    ppl.setp(cax_title,fontsize=font_scale*10)

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


def plotSliceContinuous(lat, lon, lat_arr, lon_arr, pressure, elevation, data, data_mask, 
        pngName, dataset_options, plot_style_options, plot_options):
        
    # Copy the plot options to local variables
    #title         = plot_options['title']
    #cbar_title    = plot_options['cbar_title']
    #units         = plot_options['units']
    #stride        = plot_options['stride']
    lat_0         = plot_options['lat_0']
    lon_0         = plot_options['lon_0']
    #pres_0        = plot_options['pres_0']
    latMin        = plot_options['latMin']
    lonMin        = plot_options['lonMin']
    latMax        = plot_options['latMax']
    lonMax        = plot_options['lonMax']
    #plotMin       = plot_options['plotMin']
    #plotMax       = plot_options['plotMax']
    yMin          = plot_options['yMin']
    yMax          = plot_options['yMax']
    bounding_lat  = plot_options['bounding_lat']
    scale         = plot_options['scale']
    #map_res       = plot_options['map_res']
    proj          = plot_options['proj']
    #cmap          = plot_options['cmap']
    #doScatterPlot = plot_options['scatterPlot']
    #pointSize     = plot_options['pointSize']
    #dpi           = plot_options['dpi']

    # Copy the plot options to local variables
    #llcrnrx        = plot_nav_options['llcrnrx']
    #llcrnry        = plot_nav_options['llcrnry']
    #urcrnrx        = plot_nav_options['urcrnrx']
    #urcrnry        = plot_nav_options['urcrnry']

    #lon_0         = plot_nav_options['lon_0']

    title         = plot_style_options['title']
    cbar_title    = plot_style_options['cbar_title']
    stride        = plot_style_options['stride']
    plotMin       = plot_style_options['plotMin']
    plotMax       = plot_style_options['plotMax']
    plotLims      = plot_style_options['plotLims']

    map_axis      = plot_style_options['map_axis']
    cbar_axis     = plot_style_options['cbar_axis']
    image_size    = plot_style_options['image_size']

    map_res       = plot_style_options['map_res']
    cmap          = plot_style_options['cmap']
    doScatterPlot = plot_style_options['scatterPlot']
    pointSize     = plot_style_options['pointSize']
    font_scale    = plot_style_options['font_scale']
    log_plot      = plot_style_options['log_plot']
    dpi           = plot_style_options['dpi']

    '''
    Plot the input dataset in mapped to particular projection
    '''

    LOG.info("lat.shape = {}".format(lat.shape))
    LOG.info("lon.shape = {}".format(lon.shape))
    LOG.info("lat_arr.shape = {}".format(lat_arr.shape))
    LOG.info("lon_arr.shape = {}".format(lon_arr.shape))
    LOG.info("data.shape = {}".format(data.shape))

    (nlevels,nrows) = data.shape
    lat_slice = np.broadcast_to(lat,(nlevels,nrows))
    lat_slice = ma.masked_less(lat_slice,-90.)
    LOG.info("lat_slice.shape = {}".format(lat_slice.shape))

    LOG.info("pressure.shape = {}".format(pressure.shape))

    pressure_slice = np.broadcast_to(pressure,(nrows,nlevels)).T
    LOG.info("pressure_slice = {}".format(pressure_slice))
    LOG.info("pressure_slice.shape = {}".format(pressure_slice.shape))

    LOG.info("elevation = {}".format(elevation))
    LOG.info("elevation.shape = {}".format(elevation.shape if elevation != None else None))

    # If our data is all missing, return
    #data_mask = data.mask
    if (np.sum(data_mask) == data.size):
        LOG.warn("Entire {} dataset is missing, aborting".\
                format(cbar_title))
        return -1

    # General Setup
    figWidth,figHeight = 10.,6.
    ax_rect = [0.10, 0.22, 0.78, 0.67  ] # [left,bottom,width,height]

    fig = Figure(figsize=(figWidth,figHeight))
    canvas = FigureCanvas(fig)

    ax = fig.add_axes(ax_rect)

    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    plotMin = np.min(data) if plotMin==None else plotMin
    plotMax = np.max(data) if plotMax==None else plotMax
    LOG.debug("plotMin = {}".format(plotMin))
    LOG.debug("plotMax = {}".format(plotMax))

    ###########################################################################
    # Building interpolated dataset
    ###########################################################################
    
    num_x_indicies = 101
    num_y_indicies = 101

    z = ma.array(data.ravel(),mask=data_mask.ravel())
    x = lat_slice.ravel()
    LOG.info("z.shape = {}".format(z.shape))
    LOG.info("x.shape = {}".format(x.shape))

    LOG.info("zMin,zMax = {},{}".format(np.min(z),np.max(z)))

    if elevation != None:

        yMin = 0. if yMin == None else yMin
        yMax = 50000. if yMax == None else yMax
        LOG.info("yMin,yMax = {},{}".format(yMin,yMax))

        y = elevation.real.ravel()
        LOG.info("y = {}".format(y))
        LOG.info("y.shape = {}".format(y.shape))
        LOG.info("yMin,yMax = {},{}".format(np.min(y),np.max(y)))

        if doScatterPlot:
            im = ax.scatter(x,y,c=z,
                    axes=ax, vmin=plotMin,vmax=plotMax, edgecolors='none',cmap=cmap)
        else:
            points = np.vstack((x,y)).T

            #xi = np.linspace(np.min(x), np.max(x),num_x_indicies)
            xi = lat_slice[0]
            xi.sort()
            yi = np.linspace(np.min(y), np.max(y),num_y_indicies)

            zi = griddata(x, y, z, xi, yi, interp='linear') 

            #grid_x, grid_y = np.meshgrid(xi, yi)
            #zi = griddata_sp(points, z, (grid_x, grid_y), method='linear')

            zi = ma.array(zi,mask=np.isnan(zi),fill_value=-9999.)
            zi = ma.masked_less(zi,0.)

            mask_range = dataset_options['mask_ranges'] if dataset_options['mask_ranges']!=[] else [None,None]
            zi = ma.masked_inside(zi,mask_range[0],mask_range[1])

            LOG.info("zMin,zMax = {},{}".format(np.min(zi),np.max(zi)))

            #im = ax.contourf(xi, yi, zi, cmap=cmap)
            im = ax.pcolor(xi, yi, zi, vmin=plotMin,vmax=plotMax,cmap=cmap)

        y_label = ppl.setp(ax,ylabel='elevation [ft]')

    else:

        yMin = np.max(pressure_slice) if yMin == None else yMin
        yMax = 0. if yMax == None else yMax
        LOG.info("yMin,yMax = {},{}".format(yMin,yMax))

        y = pressure_slice.ravel()
        LOG.info("y = {}".format(y))
        LOG.info("y.shape = {}".format(y.shape))

        LOG.info("zMin,zMax = {},{}".format(np.min(z),np.max(z)))

        if doScatterPlot:
            im = ax.scatter(x,y,c=z,
                    axes=ax, vmin=plotMin,vmax=plotMax, edgecolors='none',cmap=cmap)
            LOG.info("zMin,zMax = {},{}".format(np.min(z),np.max(z)))
        else:
            points = np.vstack((x,y)).T

            #xi = np.linspace(np.min(x), np.max(x),num_x_indicies)
            xi = lat_slice[0]
            xi.sort()
            yi = np.linspace(np.min(y), np.max(y),num_y_indicies)

            grid_x, grid_y = np.meshgrid(xi, yi)
            zi = griddata_sp(points, z, (grid_x, grid_y), method='linear')

            zi = ma.array(zi,mask=np.isnan(zi),fill_value=-9999.)
            zi = ma.masked_less(zi,0.)

            mask_range = dataset_options['mask_ranges'] if dataset_options['mask_ranges']!=[] else [None,None]
            zi = ma.masked_inside(zi,mask_range[0],mask_range[1])

            LOG.info("zMin,zMax = {},{}".format(np.min(zi),np.max(zi)))
            LOG.info("plotMin,plotMax = {},{}".format(plotMin,plotMax))

            #im = ax.contourf(xi, yi, zi, cmap=cmap)
            im = ax.pcolor(xi, yi, zi, vmin=plotMin,vmax=plotMax, cmap=cmap)

        y_label = ppl.setp(ax,ylabel='pressure [mb]')

    ax.set_ylim(yMin,yMax)
    ax.set_xlim(np.min(lat_slice.ravel()),np.max(lat_slice.ravel()))

    x_label = ppl.setp(ax,xlabel='latitude [$^{\circ}$]')
    ppl.setp(x_label,fontsize=font_scale*10)
    ppl.setp(y_label,fontsize=font_scale*10)

    ax.grid(True)

    txt = ax.set_title(title,fontsize=11,y=1.04)

    #ppl.setp(ax.get_xticklines(),visible=False)
    #ppl.setp(ax.get_yticklines(),visible=False)

    #ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_xticklabels(),fontsize=font_scale*10)

    #ppl.setp(ax.get_yticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(),fontsize=font_scale*10)

    # add a colorbar axis
    cax_rect = [0.10 , 0.05, 0.78 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')

    ppl.setp(cax.get_xticklabels(),fontsize=font_scale*9)
    ppl.setp(cax.get_xticklines(),visible=True)

    # Colourbar title
    cax_title = ppl.setp(cax,title=cbar_title)
    ppl.setp(cax_title,fontsize=font_scale*10)

    #
    # Add a small globe with the swath indicated on it #
    #
    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance

    globe_proj = 'lcc'

    # Compute the central lat and lon if they are not specified
    if (lat_0==None) or (lon_0==None):
        geo_shape = lat_arr.shape
        nrows, ncols = geo_shape[0],geo_shape[1]
        LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
        # Non lat/lon pair given, use the central unmasked values.
        row_idx = int(nrows/2.)
        col_idx = int(ncols/2.)
        lat_0 = lat_arr[row_idx,col_idx] if lat_0==None else lat_0
        lon_0 = lon_arr[row_idx,col_idx] if lon_0==None else lon_0
        LOG.info("Incomplete lat/lon pair given, using ({:4.2f},{:4.2f})".
                format(lat_0,lon_0))

    LOG.info("Latitude extent of data:  ({:4.2f},{:4.2f})".format(np.min(lat_arr),np.max(lat_arr)))
    LOG.info("Longitude extent of data: ({:4.2f},{:4.2f})".format(np.min(lon_arr),np.max(lon_arr)))

    glax_rect = [0.85, 0.71, 0.18, 0.18 ] # [left,bottom,width,height]
    glax = fig.add_axes(glax_rect)

    #m_globe = Basemap(lat_0=0.,lon_0=0.,\
        #ax=glax,resolution='c',area_thresh=10000.,projection='lcc')

    LOG.info("scale = {}".format(scale))
    windowWidth = scale *(0.35*12000000.)
    windowHeight = scale *(0.50*9000000.)

    # Common plotting options...
    plot_kw = {
        'ax'         : glax,
        'projection' : globe_proj,
        'lon_0'      : lon_0,
        'lat_0'      : lat_0,
        'width'      : windowWidth,
        'height'     : windowHeight,
        'fix_aspect' : True,
        'resolution' : 'c'
    }


    if ((globe_proj=='aea' or globe_proj=='lcc' or globe_proj=='eqdc') and (np.abs(lat_0)<1.e-6)):
        plot_kw['lat_1'] = 1.
        plot_kw['lat_2'] = 1.

    # Create the Basemap plotting object
    LOG.info("Plotting for {} projection...".format(proj))
    try:
        m_globe = Basemap(**plot_kw)
    except ValueError, err :
            LOG.error("{} ({} projection), aborting.".format(err,proj))
            return 1
    except RuntimeError, err :
            LOG.error("{} ({} projection), aborting.".format(err,proj))
            return 1

    # If we previously had a zero size data array, increase the pointSize
    # so the data points are visible on the global plot
    if (np.shape(lon[::stride])[0]==2) :
        pointSize = 5.

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    #swath = np.zeros(np.shape(x),dtype=int)

    x,y=m_globe(lon_arr[::stride,::stride],lat_arr[::stride,::stride])
    p_globe = m_globe.scatter(x,y,s=0.5,c="#B1B1B1",axes=glax,edgecolors='none',zorder=2)

    x,y=m_globe(lon[::stride],lat[::stride])
    p_globe = m_globe.scatter(x,y,s=0.5,c="red",axes=glax,edgecolors='none',zorder=2)

    #p_globe = m_globe.scatter(x,y,s=0.5,c="red",axes=glax,edgecolors='none',zorder=2)
    #p_globe = m_globe.scatter(x,y,s=0.5,c="red",axes=glax,edgecolors='none',zorder=2)

    # Redraw the figure
    canvas.draw()
    LOG.info("Writing image file {}".format(pngName))
    canvas.print_figure(pngName,dpi=dpi)

    return 0


def plotSliceDiscrete(lat, lon, lat_arr, lon_arr, pressure, elevation, data, data_mask,
        pngName, dataset_options, plot_style_options, plot_options):
        
    # Copy the plot options to local variables
    #title         = plot_options['title']
    #cbar_title    = plot_options['cbar_title']
    #units         = plot_options['units']
    #stride        = plot_options['stride']
    lat_0         = plot_options['lat_0']
    lon_0         = plot_options['lon_0']
    #pres_0        = plot_options['pres_0']
    latMin        = plot_options['latMin']
    lonMin        = plot_options['lonMin']
    latMax        = plot_options['latMax']
    lonMax        = plot_options['lonMax']
    #plotMin       = plot_options['plotMin']
    #plotMax       = plot_options['plotMax']
    yMin          = plot_options['yMin']
    yMax          = plot_options['yMax']
    bounding_lat  = plot_options['bounding_lat']
    scale         = plot_options['scale']
    #map_res       = plot_options['map_res']
    proj          = plot_options['proj']
    #cmap          = plot_options['cmap']
    #doScatterPlot = plot_options['scatterPlot']
    #pointSize     = plot_options['pointSize']
    #dpi           = plot_options['dpi']

    # Copy the plot options to local variables
    #llcrnrx        = plot_nav_options['llcrnrx']
    #llcrnry        = plot_nav_options['llcrnry']
    #urcrnrx        = plot_nav_options['urcrnrx']
    #urcrnry        = plot_nav_options['urcrnry']

    #lon_0         = plot_nav_options['lon_0']

    title         = plot_style_options['title']
    cbar_title    = plot_style_options['cbar_title']
    stride        = plot_style_options['stride']
    plotMin       = plot_style_options['plotMin']
    plotMax       = plot_style_options['plotMax']
    plotLims      = plot_style_options['plotLims']

    map_axis      = plot_style_options['map_axis']
    cbar_axis     = plot_style_options['cbar_axis']
    image_size    = plot_style_options['image_size']

    map_res       = plot_style_options['map_res']
    cmap          = plot_style_options['cmap']
    doScatterPlot = plot_style_options['scatterPlot']
    pointSize     = plot_style_options['pointSize']
    font_scale    = plot_style_options['font_scale']
    log_plot      = plot_style_options['log_plot']
    dpi           = plot_style_options['dpi']

    '''
    Plot the input dataset in mapped to particular projection
    '''

    LOG.info("lat.shape = {}".format(lat.shape))
    LOG.info("lon.shape = {}".format(lon.shape))
    LOG.info("lat_arr.shape = {}".format(lat_arr.shape))
    LOG.info("lon_arr.shape = {}".format(lon_arr.shape))
    LOG.info("data.shape = {}".format(data.shape))

    (nlevels,nrows) = data.shape
    lat_slice = np.broadcast_to(lat,(nlevels,nrows))
    lat_slice = ma.masked_less(lat_slice,-90.)
    LOG.info("lat_slice.shape = {}".format(lat_slice.shape))

    LOG.info("pressure.shape = {}".format(pressure.shape))
    pressure_slice = np.broadcast_to(pressure,(nrows,nlevels)).T
    LOG.info("pressure_slice.shape = {}".format(pressure_slice.shape))

    LOG.info("elevation = {}".format(elevation))
    LOG.info("elevation.shape = {}".format(elevation.shape if elevation != None else None))

    # If our data is all missing, return
    #data_mask = data.mask
    if (np.sum(data_mask) == data.size):
        LOG.warn("Entire {} dataset is missing, aborting".\
                format(cbar_title))
        return -1

    # General Setup
    figWidth,figHeight = 10.,6.
    ax_rect = [0.10, 0.22, 0.78, 0.67  ] # [left,bottom,width,height]

    fig = Figure(figsize=(figWidth,figHeight))
    canvas = FigureCanvas(fig)

    ax = fig.add_axes(ax_rect)

    # Define the discrete colorbar tick locations
    fill_colours = dataset_options['fill_colours']
    LOG.info("fill_colours = {}".format(fill_colours))
    cmap = ListedColormap(fill_colours)

    numCats = np.array(fill_colours).size
    numBounds = numCats + 1

    tickPos = np.arange(float(numBounds))/float(numCats)
    tickPos = tickPos[0 :-1] + tickPos[1]/2.

    LOG.debug('plotLims = {},{}'.format(plotLims[0],plotLims[1]))
    vmin,vmax = plotLims[0],plotLims[-1]

    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    plotMin = np.min(data) if plotMin==None else plotMin
    plotMax = np.max(data) if plotMax==None else plotMax
    LOG.debug("plotMin = {}".format(plotMin))
    LOG.debug("plotMax = {}".format(plotMax))

    ###########################################################################
    # Building interpolated dataset
    ###########################################################################
    LOG.info("elevation = {}".format(elevation))
    
    num_x_indicies = 101
    num_y_indicies = 101

    z = ma.array(data.ravel(),mask=data_mask.ravel())
    x = lat_slice.ravel()
    LOG.info("z.shape = {}".format(z.shape))
    LOG.info("x.shape = {}".format(x.shape))

    LOG.info("zMin,zMax = {},{}".format(np.min(z),np.max(z)))

    if elevation != None:

        yMin = 0. if yMin == None else yMin
        yMax = 50000. if yMax == None else yMax
        LOG.info("yMin,yMax = {},{}".format(yMin,yMax))

        y = elevation.ravel()
        LOG.info("y = {}".format(y))
        LOG.info("y.shape = {}".format(y.shape))
        LOG.info("yMin,yMax = {},{}".format(np.min(y),np.max(y)))

        if doScatterPlot:
            im = ax.scatter(x,y,c=z,
                    axes=ax, vmin=plotMin,vmax=plotMax, edgecolors='none',cmap=cmap)
        else:
            points = np.vstack((x,y)).T

            #xi = np.linspace(np.min(x), np.max(x),num_x_indicies)
            xi = lat_slice[0]
            xi.sort()
            yi = np.linspace(np.min(y), np.max(y),num_y_indicies)

            zi = griddata(x, y, z, xi, yi, interp='nn') 

            #grid_x, grid_y = np.meshgrid(xi, yi)
            #zi = griddata_sp(points, z, (grid_x, grid_y), method='nearest')

            zi = ma.array(zi,mask=np.isnan(zi),fill_value=-9999.)
            zi = ma.masked_less(zi,0.)

            mask_range = dataset_options['mask_values'] if dataset_options['mask_values']!=[] else [None,None]
            zi = ma.masked_inside(zi,mask_range[0],mask_range[1])

            LOG.info("zMin,zMax = {},{}".format(np.min(zi),np.max(zi)))

            #im = ax.contourf(xi, yi, zi, cmap=cmap)
            im = ax.pcolor(xi, yi, zi, cmap=cmap)

        y_label = ppl.setp(ax,ylabel='elevation [ft]')

    else:

        yMin = np.max(pressure_slice) if yMin == None else yMin
        yMax = 0. if yMax == None else yMax
        LOG.info("yMin,yMax = {},{}".format(yMin,yMax))

        z = ma.array(data.ravel(),mask=data_mask.ravel())
        x = lat_slice.ravel()
        y = pressure_slice.ravel()

        if doScatterPlot:
            im = ax.scatter(x,y,c=z,
                    axes=ax, vmin=plotMin,vmax=plotMax, edgecolors='none',cmap=cmap)
            LOG.info("zMin,zMax = {},{}".format(np.min(z),np.max(z)))
        else:

            points = np.vstack((x,y)).T

            #xi = np.linspace(np.min(x), np.max(x),num_x_indicies)
            xi = lat_slice[0]
            xi.sort()
            yi = np.linspace(np.min(y), np.max(y),num_y_indicies)

            grid_x, grid_y = np.meshgrid(xi, yi)
            zi = griddata_sp(points, z, (grid_x, grid_y), method='nearest')

            zi = ma.array(zi,mask=np.isnan(zi),fill_value=-9999.)
            zi = ma.masked_less(zi,0.)

            LOG.info("zMin,zMax = {},{}".format(np.min(zi),np.max(zi)))

            #im = ax.contourf(xi, yi, zi, cmap=cmap)
            im = ax.pcolor(xi, yi, zi, cmap=cmap)

        y_label = ppl.setp(ax,ylabel='pressure [mb]')

    ax.set_ylim(yMin,yMax)
    ax.set_xlim(np.min(lat_slice.ravel()),np.max(lat_slice.ravel()))

    x_label = ppl.setp(ax,xlabel='latitude [$^{\circ}$]')
    ppl.setp(x_label,fontsize=font_scale*10)
    ppl.setp(y_label,fontsize=font_scale*10)

    ax.grid(True)

    txt = ax.set_title(title,fontsize=11,y=1.04)

    #ppl.setp(ax.get_xticklines(),visible=False)
    #ppl.setp(ax.get_yticklines(),visible=False)

    #ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_xticklabels(),fontsize=font_scale*10)
    
    #ppl.setp(ax.get_yticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(),fontsize=font_scale*10)

    # add a colorbar axis
    cax_rect = [0.10 , 0.04, 0.78 , 0.05 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')

    ppl.setp(cax.get_xticklabels(),fontsize=font_scale*9)
    ppl.setp(cax.get_xticklines(),visible=False)

    # Set the colourbar tick locations and ticklabels
    #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
    tickpos_data_coords = vmax*tickPos
    cb.set_ticks(tickpos_data_coords) # In data coords (0..3)
    tick_names = dataset_options['tick_names']
    ppl.setp(cb.ax,xticklabels=tick_names)

    # Colourbar title
    cax_title = ppl.setp(cax,title=cbar_title)
    ppl.setp(cax_title,fontsize=font_scale*10)

    #
    # Add a small globe with the swath indicated on it #
    #
    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance

    globe_proj = 'lcc'

    # Compute the central lat and lon if they are not specified
    if (lat_0==None) or (lon_0==None):
        geo_shape = lat_arr.shape
        nrows, ncols = geo_shape[0],geo_shape[1]
        LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
        # Non lat/lon pair given, use the central unmasked values.
        row_idx = int(nrows/2.)
        col_idx = int(ncols/2.)
        lat_0 = lat_arr[row_idx,col_idx] if lat_0==None else lat_0
        lon_0 = lon_arr[row_idx,col_idx] if lon_0==None else lon_0
        LOG.info("Incomplete lat/lon pair given, using ({:4.2f},{:4.2f})".
                format(lat_0,lon_0))

    LOG.info("Latitude extent of data:  ({:4.2f},{:4.2f})".format(np.min(lat_arr),np.max(lat_arr)))
    LOG.info("Longitude extent of data: ({:4.2f},{:4.2f})".format(np.min(lon_arr),np.max(lon_arr)))

    glax_rect = [0.85, 0.71, 0.18, 0.18 ] # [left,bottom,width,height]
    glax = fig.add_axes(glax_rect)

    #m_globe = Basemap(lat_0=0.,lon_0=0.,\
        #ax=glax,resolution='c',area_thresh=10000.,projection='lcc')

    LOG.info("scale = {}".format(scale))
    windowWidth = scale *(0.35*12000000.)
    windowHeight = scale *(0.50*9000000.)

    # Common plotting options...
    plot_kw = {
        'ax'         : glax,
        'projection' : globe_proj,
        'lon_0'      : lon_0,
        'lat_0'      : lat_0,
        'width'      : windowWidth,
        'height'     : windowHeight,
        'fix_aspect' : True,
        'resolution' : 'c'
    }


    if ((globe_proj=='aea' or globe_proj=='lcc' or globe_proj=='eqdc') and (np.abs(lat_0)<1.e-6)):
        plot_kw['lat_1'] = 1.
        plot_kw['lat_2'] = 1.

    # Create the Basemap plotting object
    LOG.info("Plotting for {} projection...".format(proj))
    try:
        m_globe = Basemap(**plot_kw)
    except ValueError, err :
            LOG.error("{} ({} projection), aborting.".format(err,proj))
            return 1
    except RuntimeError, err :
            LOG.error("{} ({} projection), aborting.".format(err,proj))
            return 1

    # If we previously had a zero size data array, increase the pointSize
    # so the data points are visible on the global plot
    if (np.shape(lon[::stride])[0]==2) :
        pointSize = 5.

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    #swath = np.zeros(np.shape(x),dtype=int)

    x,y=m_globe(lon_arr[::stride,::stride],lat_arr[::stride,::stride])
    p_globe = m_globe.scatter(x,y,s=0.5,c="#B1B1B1",axes=glax,edgecolors='none',zorder=2)

    x,y=m_globe(lon[::stride],lat[::stride])
    p_globe = m_globe.scatter(x,y,s=0.5,c="red",axes=glax,edgecolors='none',zorder=2)

    #p_globe = m_globe.scatter(x,y,s=0.5,c="red",axes=glax,edgecolors='none',zorder=2)
    #p_globe = m_globe.scatter(x,y,s=0.5,c="red",axes=glax,edgecolors='none',zorder=2)

    # Redraw the figure
    canvas.draw()
    LOG.info("Writing image file {}".format(pngName))
    canvas.print_figure(pngName,dpi=dpi)

    return 0


#def set_plot_styles(data_obj, dataset_options, options, plot_nav_options):
def set_plot_styles(dfile_obj, dset, dataset_options, options):
    """
    Collects the various plot formatting options and does any required tweaking
    before passing to the plotting method.
    """

    plot_style_options = {}
    plot_style_options['stride'] = options.stride
    plot_style_options['plotMin'] = dataset_options['values'][0] if options.plotMin==None else options.plotMin
    plot_style_options['plotMax'] = dataset_options['values'][-1] if options.plotMax==None else options.plotMax
    plot_style_options['map_res'] = options.map_res
    plot_style_options['map_axis'] = options.map_axis
    plot_style_options['cbar_axis'] = options.cbar_axis
    plot_style_options['image_size'] = options.image_size
    plot_style_options['scatterPlot'] = options.doScatterPlot
    plot_style_options['pointSize'] = options.pointSize
    plot_style_options['font_scale'] = options.font_scale
    plot_style_options['proj'] = options.proj
    plot_style_options['dpi'] = options.dpi

    filenames = dfile_obj.datasets['file_attrs'].keys()
    filenames.sort()
    dt_image_date = dfile_obj.datasets['file_attrs'][filenames[0]]['dt_obj']

    # Set the plot title
    if options.plot_title==None:
        if options.plot_type == 'image':
            if options.pressure != None:
                vert_lev_str = "" if not dataset_options['level'] else " @ {:4.2f} hPa".format(dfile_obj.pres_0)
            elif options.elevation != None:
                vert_lev_str = "" if not dataset_options['level'] else " @ {:4.2f} feet".format(dfile_obj.elev_0)
            else:
                pass
            
        elif options.plot_type == 'slice':
            vert_lev_str = "(footprint {})".format(options.footprint)

        plot_style_options['title'] = "{} , {} {}\n{}z".format(
                dfile_obj.datasets['file_attrs'][filenames[0]]['Mission_Name'],
                dataset_options['name'],
                vert_lev_str,
                dt_image_date.strftime('%Y-%m-%d %H:%M')
                )
    else:
        plot_style_options['title'] = options.plot_title

    LOG.info("plot_style_options['title'] = {}".format(plot_style_options['title']))

    # Set the colorbar label
    if 'units' in dfile_obj.datasets['file_attrs'][filenames[0]].keys():
        plot_style_options['units'] = dfile_obj.datasets['file_attrs'][filenames[0]]['units']
    else:
        plot_style_options['units'] = dataset_options['units']

    if options.cbar_title==None:
        quantity = '' if dataset_options['quantity']==None else dataset_options['quantity']
        units = '' if dataset_options['units']==None else '[{}]'.format(dataset_options['units'])
        plot_style_options['cbar_title'] = "{} {}".format(quantity,units)
    else:
        plot_style_options['cbar_title'] = options.cbar_title

    # Set the colormap
    if options.cmap == None:
        plot_style_options['cmap'] = dataset_options['cmap']
    else :
        try:
            plot_style_options['cmap'] = getattr(cm,options.cmap)
        except AttributeError:
            warn_str = """See http://matplotlib.org/users/colormaps.html for more options."""
            LOG.warning('Colormap {} does not exist, falling back to cubehelix_r\n{}'
                    .format(options.cmap,warn_str))
            plot_style_options['cmap'] = getattr(cm,'cubehelix_r')

    # Determine whether to plot on a log scale
    if 'logscale' in dataset_options.keys():
        if options.logscale and options.no_logscale:
            LOG.warning("Only one of the options --logscale and --no_logscale can be used, defaulting to linear scale.")
            plot_style_options['log_plot'] = False
        elif options.logscale:
            plot_style_options['log_plot'] = True
        elif options.no_logscale:
            plot_style_options['log_plot'] = False
        else:
            plot_style_options['log_plot'] = dataset_options['logscale']

        if plot_style_options['log_plot']:
            if plot_style_options['plotMin'] <= 0.:
                plot_style_options['plotMin'] = 0.1
            if plot_style_options['plotMax'] <= 0.:
                plot_style_options['plotMax'] = 1.0

    else:
        if options.logscale:
            LOG.warning('The dataset {} does not support log scaling, plotting on a linear scale'
                    .format(options.dataset))
        plot_style_options['log_plot'] = False


    if plot_style_options['plotMin'] > plot_style_options['plotMax']:
        LOG.warning('Plot limit --plotMin={} > --plotMax={}, reversing the limits'
                .format(plot_style_options['plotMin'],plot_style_options['plotMax']))
        plot_style_options['plotMin'],plot_style_options['plotMax'] = plot_style_options['plotMax'],plot_style_options['plotMin']

    plot_style_options['plotLims'] = [plot_style_options['plotMin'],plot_style_options['plotMax']]

    # If this is a navigated plot, set which axes parallels and meridians get 
    # labeled at...
    #if plot_nav_options != {}:
        #if plot_nav_options['is_full_disk']:
            #plot_style_options['parallel_axes'] = [0,0,0,0]
            #plot_style_options['meridian_axes'] = [0,0,0,0]
        #else:
            #plot_style_options['parallel_axes'] = [1,0,0,0]
            #plot_style_options['meridian_axes'] = [0,0,0,1]
    #else:
        #pass

    plot_style_options['plot_map'] = plotMapDataDiscrete if dataset_options['discrete'] else plotMapDataContinuous
    plot_style_options['plot_image'] = plot_image_discrete if dataset_options['discrete'] else plot_image_continuous
    plot_style_options['plot_slice'] = plotSliceDiscrete if dataset_options['discrete'] else plotSliceContinuous

    return plot_style_options


def plot_image_continuous(data,data_mask,png_file,
        dataset_options,plot_style_options):

    # Copy the plot options to local variables
    title         = plot_style_options['title']
    cbar_title    = plot_style_options['cbar_title']
    stride        = plot_style_options['stride']
    plotMin       = plot_style_options['plotMin']
    plotMax       = plot_style_options['plotMax']
    plotLims      = plot_style_options['plotLims']

    map_axis      = plot_style_options['map_axis']
    cbar_axis     = plot_style_options['cbar_axis']
    image_size    = plot_style_options['image_size']

    cmap          = plot_style_options['cmap']
    font_scale    = plot_style_options['font_scale']
    log_plot      = plot_style_options['log_plot']
    dpi           = plot_style_options['dpi']

    '''
    Plot the input dataset in in native data coordinates
    '''

    # If our data is all missing, return
    if (np.sum(data_mask) == data.size):
        LOG.warn("Entire {} dataset is missing, aborting".\
                format(cbar_title))
        return -1

    # Construct particular data mask
    particular_mask = np.zeros(data.shape,dtype=np.bool)
    for ranges in dataset_options['mask_ranges']:
        print "Mask range is {}".format(ranges)
        particular_mask += ma.masked_inside(data,ranges[0],ranges[1]).mask

    LOG.debug("data array has shape {}".format(data.shape))
    data_aspect = float(data.shape[1])/float(data.shape[0])
    LOG.debug("data array has aspect {}".format(data_aspect))

    # Create figure with default size, and create canvas to draw on
    fig = Figure(figsize=(image_size[0],image_size[1]))
    canvas = FigureCanvas(fig)

    fig.text(0.98, 0.01, plot_style_options['version'],fontsize=font_scale*5, color='gray',ha='right',va='bottom',alpha=0.9)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = map_axis # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect,axis_bgcolor='black')

    # Granule axis title
    ax_title = ppl.setp(ax,title=title)
    ppl.setp(ax_title,fontsize=font_scale*11)
    ppl.setp(ax_title,family="sans-serif")

    LOG.debug('data.shape = {}'.format(data.shape))
    LOG.debug('data_mask.shape = {}'.format(data_mask.shape))

    data_mask = data_mask + particular_mask
    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    LOG.debug('plotLims = {},{}'.format(plotLims[0],plotLims[1]))
    vmin,vmax = plotLims[0],plotLims[-1]

    cs_kwargs = {'interpolation':'nearest',
                 'axes':ax,
                 'vmin':vmin,
                 'vmax':vmax,
                 'cmap':cmap}
    if log_plot: cs_kwargs['norm'] = LogNorm(vmin=vmin,vmax=vmax)

    im = ax.imshow(data,**cs_kwargs)

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    # add a colorbar axis
    cax_rect = cbar_axis
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb_kwargs = {'cax':cax,'orientation':'horizontal'}
    if log_plot: cb_kwargs['ticks'] = LogLocator(subs=range(10))
    cb = fig.colorbar(im, **cb_kwargs)
    ppl.setp(cax.get_xticklabels(),fontsize=font_scale*9)

    # Colourbar title
    cax_title = ppl.setp(cax,title=cbar_title)
    ppl.setp(cax_title,fontsize=font_scale*10)

    # Redraw the figure
    canvas.draw()

    # Write the figure to file
    canvas.print_figure(png_file,dpi=dpi)
    LOG.info("Writing image file {}".format(png_file))


def plot_map_continuous(lat,lon,data,png_file,
        dataset_options,plot_style_options):
        #dataset_options,plot_nav_options,plot_style_options):

    # Copy the plot options to local variables
    #llcrnrx        = plot_nav_options['llcrnrx']
    #llcrnry        = plot_nav_options['llcrnry']
    #urcrnrx        = plot_nav_options['urcrnrx']
    #urcrnry        = plot_nav_options['urcrnry']

    #lon_0         = plot_nav_options['lon_0']

    title         = plot_style_options['title']
    cbar_title    = plot_style_options['cbar_title']
    stride        = plot_style_options['stride']
    plotMin       = plot_style_options['plotMin']
    plotMax       = plot_style_options['plotMax']
    plotLims      = plot_style_options['plotLims']

    map_axis      = plot_style_options['map_axis']
    cbar_axis     = plot_style_options['cbar_axis']
    image_size    = plot_style_options['image_size']

    map_res       = plot_style_options['map_res']
    cmap          = plot_style_options['cmap']
    doScatterPlot = plot_style_options['scatterPlot']
    pointSize     = plot_style_options['pointSize']
    font_scale    = plot_style_options['font_scale']
    log_plot      = plot_style_options['log_plot']
    proj          = plot_style_options['proj']
    dpi           = plot_style_options['dpi']


    '''
    Plot the input dataset in mapped to particular projection
    '''

    if ma.is_masked(data):
        if data.mask.shape == ():
            data_mask = np.ones(data.shape,dtype='bool')
        else:
            data_mask = data.mask
    else: 
        data_mask = np.zeros(data.shape,dtype='bool')

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
        lat_0 = 0. if lat_0==None else lat_0
        lon_0 = 0. if lon_0==None else lon_0
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
    if (lat_0==None) or (lon_0==None):
        geo_shape = lat.shape
        nrows, ncols = geo_shape[0],geo_shape[1]
        LOG.debug("nrows,ncols= ({},{})".format(nrows,ncols))
        # Non lat/lon pair given, use the central unmasked values.
        row_idx = int(nrows/2.)
        col_idx = int(ncols/2.)
        lat_0 = lat[row_idx,col_idx] if lat_0==None else lat_0
        lon_0 = lon[row_idx,col_idx] if lon_0==None else lon_0
        LOG.info("Incomplete lat/lon pair given, using ({:4.2f},{:4.2f})".
                format(lat_0,lon_0))

    LOG.info("Latitude extent of data:  ({:4.2f},{:4.2f})".format(np.min(lat),np.max(lat)))
    LOG.info("Longitude extent of data: ({:4.2f},{:4.2f})".format(np.min(lon),np.max(lon)))


    # Construct particular data mask
    particular_mask = np.zeros(data.shape,dtype=np.bool)
    for ranges in dataset_options['mask_ranges']:
        print "Mask range is {}".format(ranges)
        particular_mask += ma.masked_inside(data,ranges[0],ranges[1]).mask

    LOG.debug("data array has shape {}".format(data.shape))
    data_aspect = float(data.shape[1])/float(data.shape[0])
    LOG.debug("data array has aspect {}".format(data_aspect))

    # Create figure with default size, and create canvas to draw on
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
            plot_kw['llcrnrlat'] =  -80. if latMin==None else latMin
            plot_kw['urcrnrlat'] =   80. if latMax==None else latMax
            plot_kw['llcrnrlon'] = -180. if lonMin==None else lonMin
            plot_kw['urcrnrlon'] =  180. if lonMax==None else lonMax
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

    fig.text(0.98, 0.01, plot_style_options['version'],fontsize=font_scale*5, color='gray',ha='right',va='bottom',alpha=0.9)


    # Granule axis title
    ax_title = ppl.setp(ax,title=title)
    ppl.setp(ax_title,fontsize=font_scale*11)
    ppl.setp(ax_title,family="sans-serif")

    coastline_color = dataset_options['coastline_color']
    country_color = dataset_options['country_color']
    meridian_color = dataset_options['meridian_color']
    parallel_axes = plot_style_options['parallel_axes']
    meridian_axes = plot_style_options['meridian_axes']

    m.drawcoastlines(ax=ax,color=coastline_color,linewidth = 0.3)
    m.drawcountries(ax=ax,color=country_color,linewidth = 0.2)
    m.drawstates(ax=ax,color=country_color,linewidth = 0.2)
    m.fillcontinents(color='0.',zorder=0)

    drawparallels(m,np.arange( -90, 91,30), color = meridian_color, 
            linewidth = 0.5,fontsize=font_scale*6,labels=parallel_axes) # left, right, top or bottom
    drawmeridians(m,np.arange(-180,180,30), color = meridian_color, 
            linewidth = 0.5,fontsize=font_scale*6,labels=meridian_axes) # left, right, top or bottom

    LOG.debug('data.shape = {}'.format(data.shape))
    LOG.debug('data_mask.shape = {}'.format(data_mask.shape))

    data_mask = data_mask + particular_mask
    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    LOG.debug('plotLims = {},{}'.format(plotLims[0],plotLims[1]))
    vmin,vmax = plotLims[0],plotLims[-1]

    if doScatterPlot:
        cs_kwargs = {'s':pointSize,'c':data,'axes':ax,'edgecolors':'none',
                'vmin':vmin,'vmax':vmax,'cmap':cmap}
        if log_plot: cs_kwargs['norm'] = LogNorm(vmin=vmin,vmax=vmax)
        cs = m.scatter(x,y,**cs_kwargs)
    else:
        cs_kwargs = {'axes':ax,'edgecolors':'none','antialiased':False,
                'vmin':vmin,'vmax':vmax,'cmap':cmap}
        if log_plot: cs_kwargs['norm'] = LogNorm(vmin=vmin,vmax=vmax)
        cs = m.pcolor(x,y,data,**cs_kwargs)

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    # add a colorbar axis
    cax_rect = cbar_axis
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb_kwargs = {'cax':cax,'orientation':'horizontal'}
    if log_plot: cb_kwargs['ticks'] = LogLocator(subs=range(10))
    cb = fig.colorbar(cs, **cb_kwargs)
    ppl.setp(cax.get_xticklabels(),fontsize=font_scale*9)

    # Colourbar title
    cax_title = ppl.setp(cax,title=cbar_title)
    ppl.setp(cax_title,fontsize=font_scale*10)

    # Redraw the figure
    canvas.draw()

    # Write the figure to file
    canvas.print_figure(png_file,dpi=dpi)
    LOG.info("Writing image file {}".format(png_file))


def plot_image_discrete(data,data_mask,png_file,
        dataset_options,plot_style_options):

    # Copy the plot options to local variables
    title         = plot_style_options['title']
    cbar_title    = plot_style_options['cbar_title']
    stride        = plot_style_options['stride']
    plotMin       = plot_style_options['plotMin']
    plotMax       = plot_style_options['plotMax']
    plotLims      = plot_style_options['plotLims']

    map_axis      = plot_style_options['map_axis']
    cbar_axis     = plot_style_options['cbar_axis']
    image_size    = plot_style_options['image_size']

    cmap          = plot_style_options['cmap']
    font_scale    = plot_style_options['font_scale']
    log_plot      = plot_style_options['log_plot']
    dpi           = plot_style_options['dpi']

    '''
    Plot the input dataset in in native data coordinates
    '''

    # If our data is all missing, return
    if (np.sum(data_mask) == data.size):
        LOG.warn("Entire {} dataset is missing, aborting".\
                format(dataset))
        return -1

    # Construct particular data mask
    particular_mask = np.zeros(data.shape,dtype=np.bool)
    for vals in dataset_options['mask_values']:
        print "Mask value is {}".format(vals)
        particular_mask += ma.masked_equal(data,vals).mask

    LOG.debug("data array has shape {}".format(data.shape))
    data_aspect = float(data.shape[1])/float(data.shape[0])
    LOG.debug("data array has aspect {}".format(data_aspect))

    # Create figure with default size, and create canvas to draw on
    fig = Figure(figsize=(image_size[0],image_size[1]))
    canvas = FigureCanvas(fig)

    fig.text(0.98, 0.01, plot_style_options['version'],fontsize=font_scale*5, color='gray',ha='right',va='bottom',alpha=0.9)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = map_axis # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect,axis_bgcolor='black')

    # Granule axis title
    ax_title = ppl.setp(ax,title=title)
    ppl.setp(ax_title,fontsize=font_scale*11)
    ppl.setp(ax_title,family="sans-serif")

    # Define the discrete colorbar tick locations
    fill_colours = dataset_options['fill_colours']
    cmap = ListedColormap(fill_colours)

    numCats = np.array(fill_colours).size
    numBounds = numCats + 1

    tickPos = np.arange(float(numBounds))/float(numCats)
    tickPos = tickPos[0 :-1] + tickPos[1]/2.

    LOG.debug('data.shape = {}'.format(data.shape))
    LOG.debug('data_mask.shape = {}'.format(data_mask.shape))

    data_mask = data_mask + particular_mask
    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    LOG.debug('plotLims = {},{}'.format(plotLims[0],plotLims[1]))
    vmin,vmax = plotLims[0],plotLims[-1]

    cs_kwargs = {'interpolation':'nearest',
                 'axes':ax,
                 'vmin':vmin,
                 'vmax':vmax,
                 'cmap':cmap}
    if log_plot: cs_kwargs['norm'] = LogNorm(vmin=vmin,vmax=vmax)

    im = ax.imshow(data,**cs_kwargs)

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    # add a colorbar axis
    cax_rect = cbar_axis
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb_kwargs = {'cax':cax,'orientation':'horizontal'}
    cb = fig.colorbar(im, **cb_kwargs)
    ppl.setp(cax.get_xticklabels(),fontsize=font_scale*9)
    ppl.setp(cax.get_xticklines(),visible=False)

    # Set the colourbar tick locations and ticklabels
    #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
    tickpos_data_coords = vmax*tickPos
    cb.set_ticks(tickpos_data_coords) # In data coords (0..3)
    tick_names = dataset_options['tick_names']
    ppl.setp(cb.ax,xticklabels=tick_names)

    # Colourbar title
    cax_title = ppl.setp(cax,title=cbar_title)
    ppl.setp(cax_title,fontsize=font_scale*10)

    # Redraw the figure
    canvas.draw()

    # Write the figure to file
    canvas.print_figure(png_file,dpi=dpi)
    LOG.info("Writing image file {}".format(png_file))


def plot_map_discrete(lat,lon,data,data_mask,png_file,
        dataset_options,plot_nav_options,plot_style_options):
        
    # Copy the plot options to local variables
    llcrnrx        = plot_nav_options['llcrnrx']
    llcrnry        = plot_nav_options['llcrnry']
    urcrnrx        = plot_nav_options['urcrnrx']
    urcrnry        = plot_nav_options['urcrnry']

    lon_0         = plot_nav_options['lon_0']

    title         = plot_style_options['title']
    cbar_title    = plot_style_options['cbar_title']
    stride        = plot_style_options['stride']
    plotMin       = plot_style_options['plotMin']
    plotMax       = plot_style_options['plotMax']
    plotLims      = plot_style_options['plotLims']

    map_axis      = plot_style_options['map_axis']
    cbar_axis     = plot_style_options['cbar_axis']
    image_size    = plot_style_options['image_size']

    map_res       = plot_style_options['map_res']
    cmap          = plot_style_options['cmap']
    doScatterPlot = plot_style_options['scatterPlot']
    pointSize     = plot_style_options['pointSize']
    font_scale    = plot_style_options['font_scale']
    dpi           = plot_style_options['dpi']

    '''
    Plot the input dataset in mapped to particular projection
    '''

    # If our data is all missing, return
    if (np.sum(data_mask) == data.size):
        LOG.warn("Entire {} dataset is missing, aborting".\
                format(dataset))
        return -1

    # Construct particular data mask
    particular_mask = np.zeros(data.shape,dtype=np.bool)
    for vals in dataset_options['mask_values']:
        print "Mask value is {}".format(vals)
        particular_mask += ma.masked_equal(data,vals).mask

    LOG.debug("data array has shape {}".format(data.shape))
    data_aspect = float(data.shape[1])/float(data.shape[0])
    LOG.debug("data array has aspect {}".format(data_aspect))

    # Create figure with default size, and create canvas to draw on
    fig = Figure(figsize=(image_size[0],image_size[1]))
    canvas = FigureCanvas(fig)

    fig.text(0.98, 0.01, plot_style_options['version'],fontsize=font_scale*5, color='gray',ha='right',va='bottom',alpha=0.9)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = map_axis # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect,axis_bgcolor='black')

    # Granule axis title
    ax_title = ppl.setp(ax,title=title)
    ppl.setp(ax_title,fontsize=font_scale*11)
    ppl.setp(ax_title,family="sans-serif")

    # Define the discrete colorbar tick locations
    fill_colours = dataset_options['fill_colours']
    cmap = ListedColormap(fill_colours)

    numCats = np.array(fill_colours).size
    numBounds = numCats + 1

    tickPos = np.arange(float(numBounds))/float(numCats)
    tickPos = tickPos[0 :-1] + tickPos[1]/2.

    # Setup the map
    m = Basemap(projection='geos',lon_0=lon_0,ax=ax,fix_aspect=True,resolution=map_res,
            llcrnrx=llcrnrx,
            llcrnry=llcrnry,
            urcrnrx=urcrnrx,
            urcrnry=urcrnry
            )

    x,y=m(lon[::stride,::stride],lat[::stride,::stride])

    coastline_color = dataset_options['coastline_color']
    country_color = dataset_options['country_color']
    meridian_color = dataset_options['meridian_color']
    parallel_axes = plot_style_options['parallel_axes']
    meridian_axes = plot_style_options['meridian_axes']

    m.drawcoastlines(ax=ax,color=coastline_color,linewidth = 0.3)
    m.drawcountries(ax=ax,color=country_color,linewidth = 0.2)
    m.drawstates(ax=ax,color=country_color,linewidth = 0.2)
    m.fillcontinents(color='0.',zorder=0)

    drawparallels(m,np.arange( -90, 91,30), color = meridian_color, 
            linewidth = 0.5,fontsize=font_scale*6,labels=parallel_axes) # left, right, top or bottom
    drawmeridians(m,np.arange(-180,180,30), color = meridian_color, 
            linewidth = 0.5,fontsize=font_scale*6,labels=meridian_axes) # left, right, top or bottom

    LOG.debug('data.shape = {}'.format(data.shape))
    LOG.debug('data_mask.shape = {}'.format(data_mask.shape))

    data_mask = data_mask + particular_mask
    data = ma.array(data[::stride,::stride],mask=data_mask[::stride,::stride])

    LOG.debug('plotLims = {},{}'.format(plotLims[0],plotLims[1]))
    vmin,vmax = plotLims[0],plotLims[-1]

    if doScatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=data,axes=ax,edgecolors='none',
                vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        cs = m.pcolor(x,y,data,axes=ax,edgecolors='none',antialiased=False,
                vmin=vmin,vmax=vmax,cmap=cmap)

    ppl.setp(ax.get_xticklines(),visible=False)
    ppl.setp(ax.get_yticklines(),visible=False)
    ppl.setp(ax.get_xticklabels(), visible=False)
    ppl.setp(ax.get_yticklabels(), visible=False)

    # add a colorbar axis
    cax_rect = cbar_axis
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=font_scale*9)
    ppl.setp(cax.get_xticklines(),visible=False)

    # Set the colourbar tick locations and ticklabels
    #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
    tickpos_data_coords = vmax*tickPos
    cb.set_ticks(tickpos_data_coords) # In data coords (0..3)
    tick_names = dataset_options['tick_names']
    ppl.setp(cb.ax,xticklabels=tick_names)

    # Colourbar title
    cax_title = ppl.setp(cax,title=cbar_title)
    ppl.setp(cax_title,fontsize=font_scale*10)

    # Redraw the figure
    canvas.draw()

    # Write the figure to file
    canvas.print_figure(png_file,dpi=dpi)
    LOG.info("Writing image file {}".format(png_file))


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
