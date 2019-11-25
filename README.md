# FUTURES in GRASS GIS

This is a repository for FUTURES model as GRASS GIS r.futures module.

<img align="right"  height="250" src="readme_grass_r_futures.png">

* [Installation](README.md#installation)
  * [Ubuntu Linux](README.md#ubuntu-linux)
  * [Windows](README.md#windows)
  * [Mac OS](README.md#mac-os)
* [About this repository](README.md#about-this-repository)
* [How to cite](README.md#how-to-cite)
* [Authors](README.md#authors)
* [License](README.md#license)
* [FAQs](FAQs.md)




## Installation

To get the officially released version
first install GRASS GIS (http://grass.osgeo.org/) and then install *r.futures* modules
from GRASS Addons using GUI (Settings -> Addons extensions -> Install extension from addons) or the following command:

    g.extension r.futures

You will also need r.sample.category addon:

    g.extension r.sample.category

To run r.futures.potential, you also need R (>= 3.5.0) and R packages MuMIn, lme4, optparse.
Once you install R, run in R:

    install.packages(c("MuMIn", "lme4", "optparse"))

Optionally, install SciPy to be used for certain parameters in r.futures.demand:

    pip install scipy
 
See notes for specific platforms.

### Ubuntu Linux
Install GRASS GIS from packages:

    sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
    sudo apt-get update
    sudo apt-get install grass

Install R from packages, minimum version is 3.5.0. Install SciPy using pip.

### Windows
It is recommended to first install R and then GRASS GIS, so that GRASS knows where to find R executables, because
typically they are not on your PATH.

Use [OSGeo4W](https://trac.osgeo.org/osgeo4w/) package manager and
install the latest stable GRASS GIS 7 version and the SciPy Python package using the *Advanced install* option.
A guide with step by step screen shot instructions is available [here](https://docs.google.com/presentation/d/1yEGpriBne7RvjB35HI6GecNwO1P1y2nRGqDThqoNnCE/present?usp=sharing).

### Mac OS
Install latest stable GRASS GIS from [downloads page](http://grassmac.wikidot.com/downloads)
and follow the instruction.

## About this repository

This repository is the primary repository for development of GRASS GIS
version of FUTURES model implementation.

If you want to use latest tested code, use the version in GRASS GIS Addons
Subversion repository which is installable using *g.extension* in GRASS GIS
and the code is available at:

 * https://github.com/OSGeo/grass-addons/tree/master/grass7/raster/r.futures

The corresponding documentation is available at:

 * https://grass.osgeo.org/grass7/manuals/addons/r.futures

If you want the latest code with ongoing development (which might be possibly
broken) use this repository.

The code which goes to GRASS GIS Addons is in the `r.futures` directory.

If you are a FUTURES developer and you are updating the version in GRASS GIS
Addons, use the script *copy_to_grass_addons.sh*.

If you are a FUTURES developer and you completely synced the version in
this repository and in Addons, create a Git tag using *create_tag_synced.sh*
script.


## How to cite
If you cite FUTURES as a model concept, please use [1]. If you are referring to the actual software, also use [2].
For example: *... we projected urban growth using FUTURES model (Meentemeyer et al, 2013) implemented as GRASS GIS r.futures addon (Petrasova et al, 2016) version 1.0.0 ...*


[1] Meentemeyer, R. K., Tang, W., Dorning, M. A., Vogler, J. B., Cunniffe, N. J., & Shoemaker, D. A. (2013). FUTURES: Multilevel Simulations of Emerging Urban–Rural Landscape Structure Using a Stochastic Patch-Growing Algorithm. Annals of the Association of American Geographers, 103(4), 785–807. http://doi.org/10.1080/00045608.2012.707591

[2] Petrasova, A., Petras, V., Van Berkel, D., Harmon, B. A., Mitasova, H., & Meentemeyer, R. K. (2016). Open Source Approach To Urban Growth Simulation. Int. Arch. Photogramm. Remote Sens. Spatial Inf. Sci., XLI-B7(July), 953–959. http://doi.org/10.5194/isprsarchives-XLI-B7-953-2016

See the [FUTURES website](https://cnr.ncsu.edu/geospatial/research/landscape-forecasting/futures/) for list of publications.

## Authors

 * Ross K. Meentemeyer
 * Wenwu Tang
 * Monica A. Dorning
 * John B. Vogler
 * Nik J. Cunniffe
 * Douglas A. Shoemaker
 * Jennifer A. Koch
 * Vaclav Petras
 * Anna Petrasova
 
Maintainers: Anna Petrasova, Vaclav Petras


## License

Copyright (C) 2013-2019 Meentemeyer et al.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

See the LICENSE file for details.

## FAQs
See [FAQs page](FAQs.md)
