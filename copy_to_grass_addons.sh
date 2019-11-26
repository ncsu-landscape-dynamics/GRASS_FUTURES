#!/bin/sh

# Copy relevant files to GRASS GIS addons directory

# fail fast
set -e

if [ $# -lt 1 ] ; then
    echo "Provide directory of the addon in a Git clone of grass-addons repo as a parameter"
    echo "Usage: $0 /path/to/git/addons/grass7/raster/r.futures"
    echo "Run in the directory which is a clone of the FUTURES Git repository"
    exit 1
fi

git archive HEAD r.futures/ | tar -x -C /tmp
rm /tmp/r.futures/r.futures.pga/Doxyfile
rm /tmp/r.futures/r.futures.pga/incentive.svg
rm /tmp/r.futures/r.futures.pga/futures_indent.sh
cp -r /tmp/r.futures/* $1
rm -rf /tmp/r.futures/
