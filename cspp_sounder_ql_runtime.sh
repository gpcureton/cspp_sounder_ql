#!/bin/bash
# $Id$
# Contains some package specific setup.
#
# Copyright 2011-2013, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

test -n "$CSPP_SOUNDER_QL_HOME" || echo "CSPP_IAPP_HOME is not set. Please set this environment variable to the install location of CSPP software packages."

export DPE_VER=CSPP_SOUNDER_QL_1_0

#
# derived CSPP default locations (site installs may move these under some circumstances)
#
#

#
# scripting environment settings
#

# python interpreter including numpy, h5py, pytables, scipy; used by CSPP scripts
export PY=${CSPP_SOUNDER_QL_HOME}/common/ShellB3/bin/python

# common modules location used by CSPP scripts
export PYTHONPATH=$CSPP_SOUNDER_QL_HOME/common:${CSPP_SOUNDER_QL_HOME}/scripts

#environment cleanups
unset LD_PRELOAD

test -x "$PY" || echo "Python interpreter not available; please source cspp_sounder_ql_env.sh"

# Linux execution configuration
export OSTYPE=`uname`

# make the stack size unlimited
ulimit -s unlimited

# Make the core file size unlimited, so that if the algorithm does have a
# segmentation fault, it'll write a core file that can be investigated.
ulimit -c unlimited

# Make the data size unlimited
ulimit -d unlimited
