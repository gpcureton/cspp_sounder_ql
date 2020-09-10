#!/usr/bin/env python
# encoding: utf-8
"""
This is the base class defining the common characteristics of the sounder 
software packages.

Created by Geoff Cureton on 2016-11-17.
Copyright (c) 2016 University of Wisconsin SSEC. All rights reserved.
"""

import os
import sys
import re, string
import logging
import time
from datetime import datetime


LOG = logging.getLogger('sounder_packages')


class Sounder_Packages(object):
    '''
    The various software packages treat some datasets differently. Below is a summary...


    # Latitude and Longitude 

    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:,:]                   # HSRTV
    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:,:]                   # MIRS
    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:,:]                   # IAPP
    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:].reshape(4,30)       # NUCAPS

    # Profile datasets (temp, dwpt,...)
    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][level,:,:]             # HSRTV
    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:,:,level]             # MIRS
    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:,:,level]             # IAPP
    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:,level].reshape(4,30) # NUCAPS

    # CTP, CTT

    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:,:]                   # HSRTV
    sounding_inputs[label]['data'] = sounding_inputs[label]['node'][:,:,4]                 # IAPP

    '''



    def __init__(self, *args, **kwargs):
        LOG.info("Inside Sounder_Packages.__init__()")
        pass

    def get_pressure_levels(self):
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


    def press2alt(self,p,z):
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
