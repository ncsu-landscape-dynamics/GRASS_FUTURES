#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
##############################################################################
#
# MODULE:       r.futures.calib
#
# AUTHOR(S):    Anna Petrasova (kratochanna gmail.com)
#
# PURPOSE:      FUTURES patches calibration tool
#
# COPYRIGHT:    (C) 2015 by the GRASS Development Team
#
#		This program is free software under the GNU General Public
#		License (version 2). Read the file COPYING that comes with GRASS
#		for details.
#
##############################################################################

#%module
#% description: Script for calibrating patch characteristics used as input to r.futures.pga
#% keyword: raster
#% keyword: patch
#%end
#%option G_OPT_R_INPUT
#% key: development_start
#% description: Name of input binary raster map representing development in the beginning
#%end
#%option G_OPT_R_INPUT
#% key: development_end
#% description: Name of input binary raster map representing development in the end
#%end
#%option
#% type: integer
#% key: repeat
#% description: How many times is the simulation repeated
#% required: yes
#% answer: 10
#%end
#%option
#% key: compactness_mean
#% type: double
#% description: Patch compactness mean to be tested
#% required: yes
#%end
#%option
#% type: double
#% key: compactness_range
#% description: Patch compactness range to be tested
#% required: yes
#%end
#%option
#% type: double
#% key: discount_factor
#% description: Patch size discount factor
#% required: yes
#%end
#%option G_OPT_F_OUTPUT
#% key: patch_sizes
#% description: File with patch sizes
#% required: yes
#%end
#%option
#% type: double
#% key: patch_threshold
#% description: Minimum size of a patch in meters squared
#% required: yes
#% answer: 0
#%end
#%option G_OPT_R_INPUT
#% key: development_pressure
#% required: yes
#% description: Files containing the information to read in
#% guisection: FUTURES
#%end
#%option G_OPT_R_INPUT
#% key: cons_weight
#% required: no
#% label: Name of raster map representing development potential constraint weight for scenarios
#% description: Values must be between 0 and 1, 1 means no constraint
#% guisection: FUTURES
#%end
#%option G_OPT_R_INPUTS
#% key: predictors
#% required: yes
#% multiple: yes
#% description: Names of predictor variable raster maps
#% guisection: FUTURES
#%end
#%option
#% key: n_dev_neighbourhood
#% type: integer
#% description: Size of square used to recalculate development pressure
#% required: yes
#% guisection: FUTURES
#%end
#%option G_OPT_F_INPUT
#% key: devpot_params
#% required: yes
#% multiple: yes
#% label: Development potential parameters for each region
#% description: Each line should contain region ID followed by parameters. Values are separated by whitespace (spaces or tabs). First line is ignored, so it can be used for header
#% guisection: FUTURES
#%end
#%option G_OPT_F_INPUT
#% key: incentive_table
#% required: yes
#% description: File containing incentive lookup table (infill vs. sprawl)
#% guisection: FUTURES
#%end
#%option
#% key: num_neighbors
#% type: integer
#% required: yes
#% multiple: no
#% options: 4,8
#% description: The number of neighbors to be used for patch generation (4 or 8)
#% guisection: FUTURES
#%end
#%option
#% key: seed_search
#% type: integer
#% required: yes
#% multiple: no
#% options: 1,2
#% description: The way that the location of a seed is determined
#% guisection: FUTURES
#%end
#%option
#% key: development_pressure_approach
#% type: string
#% required: yes
#% multiple: no
#% options: occurrence,gravity,kernel
#% answer: gravity
#% description: Approaches to derive development pressure
#% guisection: FUTURES
#%end
#%option
#% key: gamma
#% type: double
#% required: yes
#% multiple: no
#% description: Required for development_pressure_approach 1 and 2
#% guisection: FUTURES
#%end
#%option
#% key: scaling_factor
#% type: double
#% required: yes
#% multiple: no
#% description: Required for development_pressure_approach 2 and 3
#% guisection: FUTURES
#%end
#%option
#% key: num_regions
#% type: integer
#% required: yes
#% multiple: no
#% description: Number of sub-regions (e.g., counties) to be simulated
#% guisection: FUTURES
#%end
#%option G_OPT_R_INPUT
#% key: subregions
#% required: yes
#% description: Raster map of subregions
#% guisection: FUTURES
#%end
#%option G_OPT_F_INPUT
#% key: demand
#% required: yes
#% description: Control file with number of cells to convert
#% guisection: FUTURES
#%end


import sys
import os
import atexit
import numpy as np

from grass.exceptions import CalledModuleError
import grass.script.core as gcore
import grass.script.raster as grast
import grass.script.utils as gutils


TMP = []
TMPFILE = None
CLEANUP = True


def cleanup():
    if CLEANUP:
        gcore.run_command('g.remove', flags='f', type=['raster', 'vector'], name=TMP)
        print TMPFILE
        gutils.try_remove(TMPFILE)


def main():
    dev_start = options['development_start']
    dev_end = options['development_end']
    repeat = int(options['repeat'])
    compactness_mean = float(options['compactness_mean'])
    compactness_range = float(options['compactness_range'])
    discount_factor = float(options['discount_factor'])
    patches_file = options['patch_sizes']
    threshold = float(options['patch_threshold'])
    # v.clean removes size <= threshold, we want to keep size == threshold
    threshold -= 1e-6

    # compute cell size
    region = gcore.region()
    res = (region['nsres'] + region['ewres'])/2.
    coeff = float(gcore.parse_command('g.proj', flags='g')['meters'])
    cell_size = res * res * coeff * coeff

    tmp_name = 'tmp_futures_calib_' + str(os.getpid()) + '_'
    global TMP, TMPFILE

    orig_patch_diff = tmp_name + 'orig_patch_diff'
    TMP.append(orig_patch_diff)
    tmp_patch_vect = tmp_name + 'tmp_patch_vect'
    TMP.append(tmp_patch_vect)
    temp_file = TMPFILE = gcore.tempfile(create=False)
    simulation_dev_end = tmp_name + 'simulation_dev_end'
    simulation_dev_diff = tmp_name + 'simulation_dev_diff'
    TMP.append(simulation_dev_end)
    TMP.append(simulation_dev_diff)

    gcore.message(_("Analyzing original patches..."))
    diff_development(dev_start, dev_end, options['subregions'], orig_patch_diff)
    patch_analysis(orig_patch_diff, tmp_patch_vect, threshold, temp_file)
    area, perimeter = np.loadtxt(fname=temp_file, unpack=True)
    compact = compactness(area, perimeter)
    write_patches_file(tmp_patch_vect, cell_size, patches_file)

    # area histogram
    area = area / cell_size
    bin_width = 1.  # automatic ways to determine bin width do not perform well in this case
    hist_bins_area_orig = int(np.ptp(area) / bin_width)
    hist_range_area_orig = (np.min(area), np.max(area))
    histogram_area_orig, _edges = np.histogram(area, bins=hist_bins_area_orig,
                                               range=hist_range_area_orig, density=True)
    histogram_area_orig = histogram_area_orig * 100  # to get percentage for readability

    # compactness histogram
    bin_width = 0.1
    hist_bins_compactness_orig = int(np.ptp(compact) / bin_width)
    hist_range_compactness_orig = (np.min(compact), np.max(compact))
    histogram_compactness_orig, _edges = np.histogram(compact, bins=hist_bins_compactness_orig,
                                                      range=hist_range_compactness_orig, density=True)
    histogram_compactness_orig = histogram_compactness_orig * 100  # to get percentage for readability
#    import matplotlib.pyplot as plt
#    width = 0.7 * (_edges[1] - _edges[0])
#    center = (_edges[:-1] + _edges[1:]) / 2
#    plt.bar(center, histogram_compactness_orig, align='center', width=width)
#    plt.show()

    sum_dist_area = 0
    sum_dist_compactness = 0
    for i in range(repeat):
        gcore.message(_("Running FUTURES simulation {i}/{repeat}...".format(i=i + 1, repeat=repeat)))
        run_simulation(development_start=dev_start, development_end=simulation_dev_end,
                       compactness_mean=compactness_mean, compactness_range=compactness_range,
                       discount_factor=discount_factor, patches_file=patches_file, fut_options=options)
        new_development(simulation_dev_end, simulation_dev_diff)
        patch_analysis(simulation_dev_diff, tmp_patch_vect, threshold, temp_file)
        sim_hist_area, sim_hist_compactness = create_histograms(temp_file, hist_bins_area_orig, hist_range_area_orig,
                                                                hist_bins_compactness_orig, hist_range_compactness_orig, cell_size)
        sum_dist_area += compare_histograms(histogram_area_orig, sim_hist_area)
        sum_dist_compactness += compare_histograms(histogram_compactness_orig, sim_hist_compactness)

    mean_dist_area = sum_dist_area / repeat
    mean_dist_compactness = sum_dist_compactness / repeat
    gcore.message(_("Summary of calibrated parameters and error averaged from %s PGA runs:") % repeat)
    print "input_discount_factor=%s" % discount_factor
    print "input_compactness_mean=%s" % compactness_mean
    print "input_compactness_range=%s" % compactness_range
    print "area_distance=%s" % mean_dist_area
    print "compactness_distance=%s" % mean_dist_compactness


def run_simulation(development_start, development_end, compactness_mean, compactness_range, discount_factor, patches_file, fut_options):
    parameters = dict(patch_mean=compactness_mean, patch_range=compactness_range,
                           discount_factor=discount_factor, patch_sizes=patches_file,
                           developed=development_start)
    futures_parameters = dict(development_pressure=fut_options['development_pressure'],
                              predictors=fut_options['predictors'], n_dev_neighbourhood=fut_options['n_dev_neighbourhood'],
                              devpot_params=fut_options['devpot_params'], incentive_table=fut_options['incentive_table'],
                              num_neighbors=fut_options['num_neighbors'], seed_search=fut_options['seed_search'],
                              development_pressure_approach=fut_options['development_pressure_approach'], gamma=fut_options['gamma'],
                              scaling_factor=fut_options['scaling_factor'], num_regions=fut_options['num_regions'],
                              subregions=fut_options['subregions'], demand=fut_options['demand'],
                              output=development_end)
    parameters.update(futures_parameters)
    for not_required in ('cons_weight',):
        if options[not_required]:
            parameters.update({not_required: options[not_required]})

    gcore.run_command('r.futures.pga', flags='s', overwrite=True, **parameters)


def diff_development(development_start, development_end, subregions, development_diff):
    grast.mapcalc(exp="{res} = if({subregions} && {dev_end} && (isnull({dev_start}) ||| !{dev_start}), 1, null())".format(res=development_diff, subregions=subregions,
                  dev_end=development_end, dev_start=development_start), overwrite=True, quiet=True)


def new_development(development_end, development_diff):
    grast.mapcalc(exp="{res} = if({dev_end} > 0, 1, null())".format(res=development_diff,
                  dev_end=development_end), overwrite=True, quiet=True)


def patch_analysis(development_diff, tmp_vector_patches, threshold, output_file):
    tmp_patch_vect2 = tmp_vector_patches + '2'
    global TMP
    TMP.append(tmp_patch_vect2)
    gcore.run_command('r.to.vect', input=development_diff, output=tmp_patch_vect2, type='area', overwrite=True, quiet=True)
    gcore.run_command('v.clean', input=tmp_patch_vect2, output=tmp_vector_patches, tool='rmarea', threshold=threshold, quiet=True, overwrite=True)
    gcore.run_command('v.db.addcolumn', map=tmp_vector_patches, columns="area double precision,perimeter double precision", quiet=True)
    gcore.run_command('v.to.db', map=tmp_vector_patches, option='area', column='area', units='meters', quiet=True)
    gcore.run_command('v.to.db', map=tmp_vector_patches, option='perimeter', column='perimeter', units='meters', quiet=True)
    gcore.run_command('v.db.select', map=tmp_vector_patches, columns=['area', 'perimeter'],
                      flags='c', separator='space', file_=output_file, overwrite=True, quiet=True)


def create_histograms(input_file, hist_bins_area_orig, hist_range_area_orig, hist_bins_compactness_orig, hist_range_compactness_orig, cell_size):
    area, perimeter = np.loadtxt(fname=input_file, unpack=True)
    compact = compactness(area, perimeter)
    histogram_area, _edges = np.histogram(area / cell_size, bins=hist_bins_area_orig,
                                          range=hist_range_area_orig, density=True)
    histogram_area = histogram_area * 100
    histogram_compactness, _edges = np.histogram(compact, bins=hist_bins_compactness_orig,
                                                 range=hist_range_compactness_orig, density=True)
    histogram_compactness = histogram_compactness * 100
    return histogram_area, histogram_compactness


def write_patches_file(vector_patches, cell_size, output_file):
    gcore.run_command('v.db.select', map=vector_patches, columns='area',
                      flags='c', separator='space', file_=output_file, quiet=True)
    areas = np.loadtxt(fname=output_file)
    areas = np.round(areas / cell_size)
    np.savetxt(fname=output_file, X=areas.astype(int), fmt='%u')


def compare_histograms(hist1, hist2):
    """
    >>> hist1, edg = np.histogram(np.array([1, 1, 2, 2.5, 2.4]), bins=3, range=(0, 6))
    >>> hist2, edg = np.histogram(np.array([1, 1, 3 ]), bins=3, range=(0, 6))
    >>> compare_histograms(hist1, hist2)
    0.5
    """
    mask = np.logical_not(np.logical_or(hist1, hist2))
    hist1 = np.ma.masked_array(hist1, mask=mask)
    hist2 = np.ma.masked_array(hist2, mask=mask)
    res = 0.5 * np.sum(np.power(hist1 - hist2, 2) / (hist1.astype(float) + hist2))
    return res


def compactness(area, perimeter):
    return perimeter / (2 * np.sqrt(np.pi * area))


if __name__ == "__main__":
    options, flags = gcore.parser()
    atexit.register(cleanup)
    sys.exit(main())
