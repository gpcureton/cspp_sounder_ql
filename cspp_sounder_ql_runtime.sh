#!/bin/bash
# Sets the PATH and python interpreter locations for the CSPP Sounder QL package
#
# Copyright 2018, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

# python interpreter including numpy, h5py, pytables, scipy; used by CSPP scripts
export PY=${CSPP_SOUNDER_QL_HOME}/vendor/ShellB3/bin/python

# common modules location used by CSPP scripts
export PYTHONPATH=${CSPP_SOUNDER_QL_HOME}/scripts

export PATH=${PYTHONPATH}:${PATH}

# insurance
unset LD_LIBRARY_PATH
