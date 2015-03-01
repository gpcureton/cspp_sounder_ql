#!/bin/bash
# $Id: ql_iapp_level2_skewt.sh 2354 2015-02-12 06:59:24Z geoffc $
#
# Wrapper script for sounder_skewt.py, which generates 
# skew-T plots from IAPP level-2 files.
#
# Environment settings:
# CSPP_RT_HOME : the location of the CSPP_RT directory
#
# Copyright 2014-2014, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.
#
# file_Date = '$Date: 2015-02-11 22:59:24 -0800 (Wed, 11 Feb 2015) $'
# file_Revision = '$Revision: 2354 $'
# file_Author = '$Author: geoffc $'
# file_HeadURL = '$HeadURL: https://svn.ssec.wisc.edu/repos/jpss_adl/trunk/scripts/iapp/quicklooks/ql_iapp_level2_skewt.sh $'
# file_Id = '$Id: ql_iapp_level2_skewt.sh 2354 2015-02-12 06:59:24Z geoffc $'

if [ -z "$CSPP_IAPP_HOME" ]; then
    echo "CSPP_IAPP_HOME is not set, but is required for this script to operate."
    exit 9
fi

. ${CSPP_IAPP_HOME}/cspp_iapp_runtime.sh

usage() {
    $PY $CSPP_IAPP_HOME/iapp/quicklooks/sounder_skewt.py --help
}

if [ -z "$1" ]; then
    usage
    exit 3
fi

$PY $CSPP_IAPP_HOME/iapp/quicklooks/sounder_skewt.py  -zz "$@"

