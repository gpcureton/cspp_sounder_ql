#!/bin/bash
#
# Bumps the repo version number and commits, then tags.
#
# Usage:
#
#   tag_version.sh <version_number>
#
# where:
#
#   version_number: desired version number (symantic numbering preferred).
#
# Copyright 2019, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

rc=0

if [ $# != 1 ] ;
then

    if [ $# == 1 ] && [ $1 == "-h" ] ;
    then
        sleep 0
    else
        echo $'\n'"ERROR: Incorrect number of arguments ($#), must be 1."$'\n'
        rc=1
    fi

    echo "usage: tag_version.sh <version_number>"\
    $'\n\n'"where:"\
    $'\n\t'"version_number: desired version number (symantic numbering preferred)."
    exit $rc
fi

VERSION_TAG=$1

# Absolute path to the directory this script is in
#scriptdir=$(realpath $(dirname $0))
scriptdir=$(readlink -f $(dirname $0))

echo "The proposed tag is: "$VERSION_TAG
echo "The current directory is "$scriptdir

CSPP_SOUNDER_QL_VERSION_FILE=$scriptdir/VERSION

# Edit the VERSION file in the repo root
if [ -f $CSPP_SOUNDER_QL_VERSION_FILE ] && [ -w $CSPP_SOUNDER_QL_VERSION_FILE ];
then
    echo $CSPP_SOUNDER_QL_VERSION_FILE" exists!"
    VERSION=$(cat $CSPP_SOUNDER_QL_VERSION_FILE)
    echo "VERSION = "$VERSION
    echo "Deleting old file "$CSPP_SOUNDER_QL_VERSION_FILE
    rm -f $CSPP_SOUNDER_QL_VERSION_FILE
    echo "Creating new VERSION file..."
    echo $VERSION_TAG >> $CSPP_SOUNDER_QL_VERSION_FILE
else
    echo $CSPP_SOUNDER_QL_VERSION_FILE" does not exist or is not writable!"
    exit 1
fi

git commit -a -m "Bumped version to '$VERSION_TAG'"

if [ $? -eq 0 ] ;
then
    echo "Commit command succeeded"
else
    echo "Commit command failed, aborting."
    exit 1
fi

git tag -a 'cspp-sounder-ql-'$VERSION_TAG -m "Release $VERSION_TAG"

if [ $? -eq 0 ] ;
then
    echo "Tag command succeeded"
else
    echo "Tag command failed, aborting."
    exit 1
fi

exit 0
