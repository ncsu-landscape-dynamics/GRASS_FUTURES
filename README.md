FUTURES in GRASS GIS
====================

This is a repository for FUTURES model port to GRASS GIS,
namely for r.futures module.

![r.futures module, results and manual](readme_grass_r_futures.png)


Getting r.futures
-----------------

To get the officially released version
install GRASS GIS (http://grass.osgeo.org/) and then install an *r.futures*
from GRASS Addons using GUI or the following command:

    g.extension r.futures


About this repository
---------------------

This repository is the primary repository for development of GRASS GIS
version of FUTURES model implementation.

If you want to use latest tested code, use the version in GRASS GIS Addons
Subversion repository which is installable using *g.extension* in GRASS GIS
and the code is available at:

 * https://trac.osgeo.org/grass/browser/grass-addons/grass7/raster/r.futures

The correspong documentation is available at:

 * https://grass.osgeo.org/grass70/manuals/addons/r.futures

If you want the latest code with ongoing development (which might be possibly
broken) use this repository.

The code which goes to GRASS GIS Addons is in the `r.futures` directory.

If you are a FUTURES developer and you are updating the version in GRASS GIS
Addons, use the script *copy_to_grass_addons.sh*.

If you are a FUTURES developer and you completely synced the version in
this repository and in Addons, create a Git tag using *create_tag_synced.sh*
script.


Authors
-------

 * Ross K. Meentemeyer
 * Wenwu Tang
 * Monica A. Dorning
 * John B. Vogler
 * Nik J. Cunniffe
 * Douglas A. Shoemaker
 * Jennifer A. Koch
 * Vaclav Petras
 * Anna Petrasova

See the GRASS module manual page for details and references.


License
-------

Copyright (C) 2013-2015 Meentemeyer et al.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

See the LICENSE file for details.
