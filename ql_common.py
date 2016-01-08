#!/usr/bin/env python
# encoding: utf-8
"""
ql_geocat_common.py

Purpose: Common methods for level 1 and level 2 plots from geocat output.

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

from scipy import vectorize

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.ticker import LogLocator,LogFormatter

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

from mpl_toolkits.basemap import Basemap

from pyhdf.SD import SD
from netCDF4 import Dataset
from netCDF4 import num2date

from basemap_utils import drawparallels,drawmeridians

# every module should have a LOG object
LOG = logging.getLogger(__file__)


def _tuple2args(parms):
    s = ' '.join( '+%s=%s' % (k,v) for (k,v) in parms )
    return s.encode('ascii')

class GOES_HDF4():

    def __init__(self,input_file):


        self.input_file = input_file

        LOG.debug('Opening {} with GOES_HDF4...'.format(self.input_file))
        self.file_obj = SD(self.input_file)
        
        # Dictionary of file object attributes
        self.attrs = self.file_obj.attributes()

        # Dictionary of dataset shape and type attributes
        self.data_dict = self.file_obj.datasets()

        # List of dataset names
        self.datanames = self.data_dict.keys()
        self.datanames.sort()


    class Dataset():

        def __init__(selfd,L1_obj,dataname):

            selfd.dataname = dataname

            selfd.dset_obj = L1_obj.file_obj.select(dataname)

            selfd.attrs = selfd.dset_obj.attributes()
            selfd.dset = ma.masked_equal(selfd.dset_obj.get(),selfd.attrs['_FillValue'])
            
            selfd.dset = selfd.dset * selfd.attrs['scale_factor'] + selfd.attrs['add_offset']

    def close(self):
        LOG.debug('Closing {}...'.format(self.input_file))
        self.file_obj.end()


class GOES_NetCDF():

    def __init__(self,input_file):


        self.input_file = input_file

        if os.path.exists(self.input_file):
            LOG.debug('Opening {} with GOES_NetCDF...'.format(self.input_file))
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


def set_plot_navigation_bm(lats,lons,goes_l1_obj, options):
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


    subsatellite_longitude = goes_l1_obj.attrs['Subsatellite_Longitude']
    LOG.info('File subsatellite_Longitude = {:6.1f}'.format(subsatellite_longitude))

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


def set_plot_styles(goes_obj, data_obj, dataset_options, options, plot_nav_options):
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
    plot_style_options['dpi'] = options.dpi

    image_date_time = goes_obj.attrs['Image_Date_Time']
    dt_image_date = datetime.strptime(image_date_time,'%Y-%m-%dT%H:%M:%SZ')

    # Set the plot title
    if options.plot_title==None:
        plot_style_options['title'] = "{} Imager, {}\n{}z".format(
                goes_obj.attrs['Spacecraft_Name'],
                dataset_options['name'],
                dt_image_date.strftime('%Y-%m-%d %H:%M')
                )
    else:
        plot_style_options['title'] = options.plot_title

    # Set the colorbar label
    plot_style_options['units'] = dataset_options['units'] if data_obj.attrs['units']=="none" else data_obj.attrs['units']
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
    if plot_nav_options != {}:
        if plot_nav_options['is_full_disk']:
            plot_style_options['parallel_axes'] = [0,0,0,0]
            plot_style_options['meridian_axes'] = [0,0,0,0]
        else:
            plot_style_options['parallel_axes'] = [1,0,0,0]
            plot_style_options['meridian_axes'] = [0,0,0,1]
    else:
        pass

    plot_style_options['plot_map'] = plot_map_discrete if dataset_options['discrete'] else plot_map_continuous
    plot_style_options['plot_image'] = plot_image_discrete if dataset_options['discrete'] else plot_image_continuous

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


def plot_map_continuous(lat,lon,data,data_mask,png_file,
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
