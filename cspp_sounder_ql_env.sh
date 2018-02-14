#!/bin/bash
# Environment script for the CSPP Sounder QL package
#
# Copyright 2014-2018, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "${CSPP_SOUNDER_QL_HOME}" ]; then
    echo "CSPP_SOUNDER_QL_HOME must be set to the path where the CSPP software was installed."
    echo "i.e.: export CSPP_SOUNDER_QL_HOME=/home/me/Sounder_QL"
    exit 1
fi
export PATH=${PATH}:${CSPP_SOUNDER_QL_HOME}/vendor/env/bin
export PATH=${PATH}:${CSPP_SOUNDER_QL_HOME}/bin:${CSPP_SOUNDER_QL_HOME}/scripts
