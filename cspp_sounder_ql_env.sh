#!/bin/bash
# $Id$
# Environment script for CSPP / IAPP


test -n "$CSPP_SOUNDER_QL_HOME" || echo "CSPP_SOUNDER_QL_HOME is not set. Please set this environment variable to the install location of CSPP software packages."

#
# user path environment settings, making it easy to invoke wrapper scripts
#


export PATH=${CSPP_SOUNDER_QL_HOME}/common:$PATH
export PATH=${CSPP_SOUNDER_QL_HOME}/common/ShellB3/bin:$PATH
export PATH=${CSPP_SOUNDER_QL_HOME}/sounder_ql:$PATH


