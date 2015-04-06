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


import sys
import os
import atexit
import numpy as np
import random

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

    gcore.message(_("Analyzing original patches..."))
    diff_development(dev_start, dev_end, orig_patch_diff)
    patch_analysis(orig_patch_diff, tmp_patch_vect, temp_file)
    area, compactness = np.loadtxt(fname=temp_file, unpack=True)

    # area histogram
    first_qrt, third_qrt = np.percentile(area, 25), np.percentile(area, 75)
    bin_width = 2 * (third_qrt - first_qrt) / pow(len(area), 1/3.)
    hist_bins_area_orig = int(np.ptp(area) / bin_width)
    hist_range_area_orig = (np.min(area), np.max(area))
    histogram_area_orig, _edges = np.histogram(area, bins=hist_bins_area_orig,
                                               range=hist_range_area_orig, density=True)
    # compactness histogram
    first_qrt, third_qrt = np.percentile(compactness, 25), np.percentile(compactness, 75)
    bin_width = 2 * (third_qrt - first_qrt) / pow(len(compactness), 1/3.)
    hist_bins_compactness_orig = int(np.ptp(compactness) / bin_width)
    hist_range_compactness_orig = (np.min(compactness), np.max(compactness))
    histogram_compactness_orig, _edges = np.histogram(compactness, bins=hist_bins_compactness_orig,
                                                      range=hist_range_compactness_orig, density=True)
#    import matplotlib.pyplot as plt
#    width = 0.7 * (_edges[1] - _edges[0])
#    center = (_edges[:-1] + _edges[1:]) / 2
#    plt.bar(center, histogram_compactness_orig, align='center', width=width)
#    plt.show()

    sum_dist_area = 0
    sum_dist_compactness = 0
    for i in range(repeat):
        gcore.message(_("Running FUTURES simulation {i}/{repeat}...".format(i=i, repeat=repeat)))
        run_simulation(development_start=dev_start, development_end=simulation_dev_end,
                       compactness_mean=compactness_mean, compactness_range=compactness_range, discount_factor=discount_factor)
        diff_development(dev_start, simulation_dev_end, simulation_dev_diff)
        patch_analysis(simulation_dev_diff, tmp_patch_vect, temp_file)
        sim_hist_area, sim_hist_compactness = create_histograms(temp_file, hist_bins_area_orig, hist_range_area_orig,
                                                                hist_bins_compactness_orig, hist_range_compactness_orig)
        sum_dist_area += compare_histograms(histogram_area_orig, sim_hist_area)
        sum_dist_compactness += compare_histograms(histogram_compactness_orig, sim_hist_compactness)

    gcore.message(_("Averaged error:"))
    mean_dist_area = sum_dist_area / repeat
    mean_dist_compactness = sum_dist_compactness / repeat
    print "area_distance=%s" % mean_dist_area
    print "compactness_distance=%s" % mean_dist_compactness


def run_simulation(development_start, development_end, compactness_mean, compactness_range, discount_factor):
    # here will be PGA call, r.grow is a placeholder
    gcore.run_command('r.grow', input=development_start, output=development_end, radius=4 * random.random(), new=1, old=1, overwrite=True)


def diff_development(development_start, development_end, development_diff):
    grast.mapcalc(exp="{res} = if({dev_end} && (isnull({dev_start}) ||| !{dev_start}), 1, null())".format(res=development_diff,
                  dev_end=development_end, dev_start=development_start), overwrite=True, quiet=True)


def patch_analysis(development_diff, tmp_vector_patches, output_file):
    gcore.run_command('r.to.vect', input=development_diff, output=tmp_vector_patches, type='area', flags='s', overwrite=True, quiet=True)
    gcore.run_command('v.db.addcolumn', map=tmp_vector_patches, columns="area double,compactness double", quiet=True)
    gcore.run_command('v.to.db', map=tmp_vector_patches, option='area', column='area', quiet=True)
    gcore.run_command('v.to.db', map=tmp_vector_patches, option='compact', column='compactness', quiet=True)
    gcore.run_command('v.db.select', map=tmp_vector_patches, columns=['area', 'compactness'],
                      flags='c', separator='space', file_=output_file, overwrite=True, quiet=True)


def create_histograms(input_file, hist_bins_area_orig, hist_range_area_orig, hist_bins_compactness_orig, hist_range_compactness_orig):
    area, compactness = np.loadtxt(fname=input_file, unpack=True)
    histogram_area, _edges = np.histogram(area, bins=hist_bins_area_orig,
                                          range=hist_range_area_orig, density=True)
    histogram_compactness, _edges = np.histogram(compactness, bins=hist_bins_compactness_orig,
                                                 range=hist_range_compactness_orig, density=True)
    return histogram_area, histogram_compactness


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

if __name__ == "__main__":
    options, flags = gcore.parser()
    atexit.register(cleanup)
    sys.exit(main())
