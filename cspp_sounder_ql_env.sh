#!/bin/bash
# $Id$
# Environment script for CSPP / IAPP


test -n "$CSPP_IAPP_HOME" || echo "CSPP_IAPP_HOME is not set. Please set this environment variable to the install location of CSPP software packages. (When installed, \$CSPP_IAPP_HOME/ADL is a directory.)"

test -d "$CSPP_IAPP_HOME/common/IAPP_VENDOR" || echo "CSPP_IAPP_HOME does not appear to be set properly. See cspp_iapp_env.sh"

# revision string for this CSPP release, which we set if we have reasonable expectation that environment is correct
test -d "$CSPP_IAPP_HOME/common/IAPP_VENDOR" && export CSPP_REV="20120215"


#
# user path environment settings, making it easy to invoke wrapper scripts
#


export PATH=${CSPP_IAPP_HOME}/common:$PATH
export PATH=${CSPP_IAPP_HOME}/common/ShellB3/bin:$PATH
export PATH=${CSPP_IAPP_HOME}/iapp:$PATH


