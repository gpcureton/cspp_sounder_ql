#!/bin/bash
# Bash front-end script for CSPP Sounder QL skew-T python script
#
# Copyright 2014-2018, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_SOUNDER_QL_HOME" ]; then
    echo "CSPP_SOUNDER_QL_HOME must be set to the path where the CSPP software was installed."
    echo "i.e.: export CSPP_SOUNDER_QL_HOME=/home/me/Sounder_QL"
    exit 1
fi

. ${CSPP_SOUNDER_QL_HOME}/cspp_sounder_ql_runtime.sh

usage() {
    $PY -W ignore $CSPP_SOUNDER_QL_HOME/scripts/sounder_skewt.py --help
}

if [ -z "$1" ]; then
    usage
    exit 3
fi

$PY -W ignore ${CSPP_SOUNDER_QL_HOME}/scripts/sounder_skewt.py "$@"
