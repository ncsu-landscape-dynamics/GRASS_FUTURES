#!/usr/bin/env bash

# The make step requires something like:
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PREFIX/lib"
# further steps additionally require:
# export PATH="$PATH:$PREFIX/bin"

# fail on non-zero return code from a subprocess
set -e

if [ -z "$1" ]
then
    echo "Usage: $0 PREFIX"
    exit 1
fi

export INSTALL_PREFIX=$1

# GRASS GIS

git clone https://github.com/OSGeo/grass.git --branch main --depth=1

cd grass

./configure \
    --prefix="$INSTALL_PREFIX/" \
    --enable-largefile \
    --with-cxx \
    --with-zstd \
    --with-bzlib \
    --with-blas \
    --with-lapack \
    --with-readline \
    --with-openmp \
    --with-pthread \
    --with-tiff \
    --with-freetype \
    --with-freetype-includes="/usr/include/freetype2/" \
    --with-proj-share=/usr/share/proj \
    --with-geos \
    --with-sqlite \
    --with-fftw \
    --with-netcdf \
    --without-pdal

make
make install
