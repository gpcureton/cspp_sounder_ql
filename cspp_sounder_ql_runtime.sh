#!/bin/bash
# Sets the PATH and python interpreter locations for the CSPP Sounder QL package
#
# Copyright 2014-2018, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

export DPE_VER=CSPP_SOUNDER_QL_1_1

# python interpreter including numpy, h5py, pytables, scipy; used by CSPP scripts
export PY=${CSPP_SOUNDER_QL_HOME}/vendor/env/bin/python

# common modules location used by CSPP scripts
export PYTHONPATH=${CSPP_SOUNDER_QL_HOME}/scripts

export PATH=${PYTHONPATH}:${PATH}

# insurance
export LD_LIBRARY_PATH=${CSPP_SOUNDER_QL_HOME}/vendor/env/lib
export LD_LIBRARY_PATH=${CSPP_SOUNDER_QL_HOME}/vendor:${LD_LIBRARY_PATH}
