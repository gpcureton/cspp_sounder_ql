#!/bin/bash
#
# Creates the CSPP Sounder Quicklook Package tarballs
#
# Usage:
#
#   cspp-sounder-ql_packaging.sh <pkg_tree_dir> <version> <tarball_dir>
#
# where:
#
#   pkg_tree_dir : Directory containing the package tree
#   version : The desired package number
#   tarball_dir : The desired directory to copy the tarballs to
#
# Copyright 2020, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

set -e

pkg_dir=$1
version=$2
tarball_location=$3

echo $pkg_dir
echo $version
echo $tarball_location

cur_dir=$PWD
echo "Current dir: "$cur_dir

pkg_dir_link=$pkg_dir-$version
echo "package dir link: "$pkg_dir_link

if [ -h $pkg_dir_link ];
then
    echo "link "$pkg_dir_link" exists!";
else
    echo "Creating link "$pkg_dir_link
    ln -s $pkg_dir $pkg_dir_link
fi

echo "Creating the package tarball (sans ShellB3)"
rm -f cspp-sounder-ql-$version.tar*
tar cvhf \
    cspp-sounder-ql-$version.tar \
    --exclude=__pycache__ \
    --exclude=*.pyc \
    --exclude=.git* \
    --exclude=oldstuff \
    --exclude=ShellB3 \
    --exclude=env \
    --exclude=tag_version.sh \
    --exclude=dummy_file \
    --exclude=.vim \
    $pkg_dir_link 

echo "Including ShellB3 into cspp-sounder-ql-"$version".tar"
tar rf cspp-sounder-ql-$version.tar $pkg_dir_link/vendor/ShellB3

echo "Zipping cspp-sounder-ql-"$version".tar"
gzip cspp-sounder-ql-$version.tar

echo "Creating checksum cspp-sounder-ql-"$version".tar.gz.sha1"
sha1sum cspp-sounder-ql-$version.tar.gz >> cspp-sounder-ql-$version.tar.gz.sha1

tarball_dir=$tarball_location/cspp-sounder-ql-$version
if [ -d $tarball_dir ];
then
    echo "Tarball directory "$tarball_dir" exists!";
else
    echo "Creating tarball dir "$tarball_dir
    mkdir $tarball_dir
fi

echo "Moving cspp-sounder-ql-"$version".tar.gz to "$tarball_dir
mv -f cspp-sounder-ql-$version.tar.gz* $tarball_dir

echo "Exiting..."
exit 0
