#!/usr/bin/env python
# encoding: utf-8
"""
exceptions.py

 * DESCRIPTION: This file contains a few custom exceptions

Created by Geoff Cureton on 2020-09-05.
Copyright (c) 2020 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

class UnsupportedInputType(RuntimeError):
    '''To catch unsupported input files'''
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.code)
