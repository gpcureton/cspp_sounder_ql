#!/usr/bin/env python
# encoding: utf-8
"""
This module provides classes to handle the output products for several sounder
software packages.

Created by Geoff Cureton on 2016-11-17.
Copyright (c) 2016 University of Wisconsin SSEC. All rights reserved.
"""

import logging

#from sounder_packages import Sounder_Packages
from hsrtv import HSRTV
from iapp import IAPP
from mirs import MIRS
#from nucaps import NUCAPS

LOG = logging.getLogger('__init__')

# Dictionary containing the names of the sounder package classes.
sounder_package_cfg = {
        'hsrtv': HSRTV,
        'iapp': IAPP,
        'mirs': MIRS,
        #'nucaps': NUCAPS
        }
