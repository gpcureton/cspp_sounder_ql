#!/usr/bin/env python
# encoding: utf-8
"""
sounder_image_data.py

Purpose: Provide required data for sounder image products.

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

from matplotlib import cm as cm


class Dataset_Options:
    """
    This class contains static data for the interpretation of the various  
    sounder discrete and continuous datasets.
    """

    data = {}

    data['cold_air_aloft'] = {
                'name':'Cold Air Aloft',
                'quantity': 'cold air aloft',
                'level':True,
                'discrete':True,
                'values':(0,1,2),
                'mask_values':[],
                'units': "$^{\circ}$C",
                'fill_boundaries':[-0.5,0.5,1.5,2.5],
                'fill_colours':['#9800CB','#66CCFF','#B1B1B1'],
                'tick_names':['< -65','-65 to -60','> -60'],
                'cmap':None,
                'meridian_color':'black',
                'coastline_color':'black',
                'country_color':'black'
                }

    data['temp'] = {
                'name':'temperature',
                'quantity':'temperature',
                'level':True,
                'discrete':False,
                'values':[200., 280.],
                'mask_ranges':[0.,173.16],
                'logscale':False,
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.jet,
                'n_levels':256
                }
    data['temp_gdas'] = dict(data['temp'])
    data['temp_gdas']['name'] = 'GDAS temperature'

    data['2temp'] = {
                'name':'temperature',
                'quantity':'temperature',
                'level':True,
                'discrete':False,
                'values':[440.,600.],
                'mask_ranges':[],
                'logscale':False,
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.jet,
                'n_levels':256
                }


    data['wvap'] = {
                'name':'water vapor mixing ratio',
                'quantity':'water vapor mixing ratio',
                'level':True,
                'discrete':False,
                'values':[0.,2.],
                'mask_ranges':[],
                'logscale':False,
                'units': 'g/kg',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.jet,
                'n_levels':256
                }
    data['wvap_gdas'] = dict(data['wvap'])
    data['wvap_gdas']['name'] = 'GDAS water vapor mixing ratio'

    data['dwpt'] = {
                'name':'dewpoint temperature',
                'quantity':'dewpoint temperature',
                'level':True,
                'discrete':False,
                'values':[170.,280.],
                'mask_ranges':[],
                'logscale':False,
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.jet,
                'n_levels':256
                }

    data['relh'] = {
                'name':'relative humidity',
                'quantity':'relative humidity',
                'level':True,
                'discrete':False,
                'values':[0.,100.],
                'mask_ranges':[],
                'logscale':False,
                'units': '%',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.jet,
                'n_levels':256
                }
    data['relh_gdas'] = dict(data['relh'])
    data['relh_gdas']['name'] = 'GDAS relative humidity'

    data['ctp'] = {
                'name':'cloud top pressure',
                'quantity':'cloud top pressure',
                'level':False,
                'discrete':False,
                'values':[0.,1000.],
                'mask_ranges':[],
                'logscale':False,
                'units': 'hPa',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.jet,
                'n_levels':256
                }

    data['ctt'] = {
                'name':'cloud top temperature',
                'quantity':'cloud top temperature',
                'level':False,
                'discrete':False,
                'values':[190.,270.],
                'mask_ranges':[],
                'logscale':False,
                'units': 'K',
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.jet,
                'n_levels':256
                }

    data['unknown'] = {
                'name':None,
                'quantity':None,
                'level':False,
                'discrete':False,
                'values':[None,None],
                'mask_ranges':[],
                'logscale':False,
                'units': None,
                'coastline_color':'cyan',
                'country_color':'magenta',
                'meridian_color':'yellow',
                'cmap':cm.cubehelix,
                'n_levels':256
                }
