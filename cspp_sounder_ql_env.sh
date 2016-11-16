#!/bin/bash
# $Id$
# Environment script for CSPP / SOUNDER_QL


test -n "$CSPP_SOUNDER_QL_HOME" || echo "CSPP_SOUNDER_QL_HOME is not set. Please set this environment variable to the install location of CSPP software packages."

#
# user path environment settings, making it easy to invoke wrapper scripts
#


export PATH=${CSPP_SOUNDER_QL_HOME}/common:$PATH
#export PATH=${CSPP_SOUNDER_QL_HOME}/common/ShellB3/bin:$PATH
export PATH=${CSPP_SOUNDER_QL_HOME}/common/shellb3-v1.1-2.7/ShellB3/bin:$PATH
export PATH=${CSPP_SOUNDER_QL_HOME}/scripts:$PATH


