#!/bin/sh

# Copy relevant files to GRASS GIS addons directory

# fail fast
set -e

if [ $# -lt 1 ] ; then
    echo "Provide directory of the addon in a SVN copy as a parameter"
    echo "Usage: $0 /path/to/svn/addons/grass7/raster/r.futures"
    echo "Run in the directory which is a clone of the Git repository"
    exit 1
fi

cp -r r.futures/* $1
rm $1/r.futures.pga/Doxyfile
rm $1/r.futures.pga/incentive.svg
rm $1/r.futures.pga/futures_indent.sh
