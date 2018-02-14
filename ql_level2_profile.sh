#!/bin/bash
# Bash front-end script for CSPP Sounder QL python script
#
# Copyright 2014-2018, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_SOUNDER_QL_HOME" ]; then
    echo "CSPP_SOUNDER_QL_HOME is not set, but is required for this script to operate."
    exit 9
fi

. ${CSPP_SOUNDER_QL_HOME}/cspp_sounder_ql_runtime.sh

usage() {
    $PY $CSPP_SOUNDER_QL_HOME/scripts/sounder_profile.py --help
}

if [ -z "$1" ]; then
    usage
    exit 3
fi

$PY $CSPP_SOUNDER_QL_HOME/scripts/sounder_profile.py  -zz "$@"
