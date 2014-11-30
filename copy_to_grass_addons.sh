#!/bin/sh

# Copy relevant files to GRASS GIS addons directory


if [ $# -lt 1 ] ; then
    echo "Provide directory of the addon in a SVN copy as a parameter"
    echo "Usage: $0 /path/to/svn/addons/grass7/raster/r.futures"
    echo "Run in the directory which is a clone of Git repository"
    exit 1
fi

cp *.cpp Makefile r.futures.html r_futures_detail.png r_futures.png $1
