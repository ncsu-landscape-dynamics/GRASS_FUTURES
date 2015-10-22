#!/bin/bash

# Copy relevant files to GRASS GIS addons directory

# fail fast
set -e

if [ $# -lt 1 ] ; then
    echo "Provide a GRASS GIS Subversion revision number"
    echo "Usage: $0 revision-number"
    echo "Example: $0 59542"
    echo "Run in the directory which is a clone of the Git repository"
    exit 1
fi

if [[ $1 == r* ]] ; then
    echo "Please do not include 'r' in the revision number"
    exit 1
fi

if [[ ! $1 =~ ^-?[0-9]+$ ]] ; then
    echo "Please provide a GRASS GIS Subversion revision number"
    exit 1
fi

git tag -a grass-svn-$1 -m \
    "r.futures synced with GRASS GIS Addons in r$1 (http://trac.osgeo.org/grass/changeset/$1)"
