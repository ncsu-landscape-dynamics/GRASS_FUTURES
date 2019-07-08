/****************************************************************************
 *
 * MODULE:       r.futures.pga
 * AUTHOR(S):    Ross K. Meentemeyer
 *               Wenwu Tang
 *               Monica A. Dorning
 *               John B. Vogler
 *               Nik J. Cunniffe
 *               Douglas A. Shoemaker
 *               Jennifer A. Koch
 *               Vaclav Petras <wenzeslaus gmail com>
 *               Anna Petrasova
 *               (See the manual page for details and references.)
 *
 * PURPOSE:      Simulation of urban-rural landscape structure (FUTURES model)
 *
 * COPYRIGHT:    (C) 2013-2016 by Anna Petrasova and Meentemeyer et al.
 *
 *               This program is free software: you can redistribute it and/or
 *               modify it under the terms of the GNU General Public License
 *               as published by the Free Software Foundation, either version 3
 *               of the License, or (at your option) any later version.
 *
 *               This program is distributed in the hope that it will be useful,
 *               but WITHOUT ANY WARRANTY; without even the implied warranty of
 *               MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *               GNU General Public License for more details.
 *
 *               You should have received a copy of the GNU General Public
 *               License along with this program. If not, see
 *               <http://www.gnu.org/licenses/> or read the file COPYING that
 *               comes with GRASS GIS for details.
 *
 *****************************************************************************/

/**
    \file main.c
    
    The main file containing both the model code and the data handing part.
    
    The language of the code is subject to change. The goal is to use either
    C or C++, not both mixed as it is now. Update: only C is used now.
    
    Major refactoring of the code is expected.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <grass/segment.h>

#include "keyvalue.h"
#include "inputs.h"
#include "patch.h"


enum development_pressure {OCCURRENCE, GRAVITY, KERNEL};


double get_develop_probability_xy(SEGMENT *predictors, SEGMENT *devpressure,
                                  FCELL *values,
                                  struct Potential *potential_info,
                                  int region_index, int row, int col)
{
    double probability;
    int i;
    FCELL devpressure_val;
    Segment_get(devpressure, (void *)&devpressure_val, row, col);
    Segment_get(predictors, values, row, col);
    
    probability = potential_info->intercept[region_index];
    probability += potential_info->devpressure[region_index] * devpressure_val;
    for (i = 0; i < potential_info->max_predictors; i++) {
        probability += potential_info->predictors[i][region_index] * values[i];
    }
    probability = 1.0 / (1.0 + exp(-probability));
    
    return probability;
}


void compute_develop_probability(SEGMENT *probability, SEGMENT *developed,
                                 SEGMENT *predictors, SEGMENT *devpressure,
                                 SEGMENT *subregions,
                                 struct Potential *potential_info)
{
    
    int row, col;
    double prob;

    CELL developed_cell;
    CELL region_index;
    FCELL *values = G_malloc(potential_info->max_predictors * sizeof(FCELL *));
    for (row = 0; row < Rast_window_rows(); row++) {
        for (col = 0; col < Rast_window_cols(); col++) {
            Segment_get(developed, (void *)&developed_cell, row, col);
            if (developed_cell == -1) {
                Segment_get(subregions, (void *)&region_index, row, col);
                prob = get_develop_probability_xy(predictors, devpressure, values,
                                                  potential_info, region_index, row, col);
                Segment_put(probability, (void *)&prob, row, col);
            }
        }
    }
    Segment_flush(probability);
}


void initial_probability(SEGMENT *probability, SEGMENT *developed,
                         SEGMENT *predictors, SEGMENT *devpressure,
                         SEGMENT *subregions, struct SegmentMemory segmentInfo,
                         struct Potential *potential_info)
{
    if (Segment_open(probability, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segmentInfo.rows, segmentInfo.cols,
                     Rast_cell_size(FCELL_TYPE), segmentInfo.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));
    compute_develop_probability(probability, developed, predictors, devpressure,
                                subregions, potential_info);
    Segment_flush(probability);
}


void update_development_pressure(int row, int col, SEGMENT *devpressure,
                                 int neighborhood_size, double gamma, double scaling_factor,
                                 enum development_pressure devpressure_alg) {
    int i, j;
    double dist;
    double value;
    float devpressure_value;

    /* this can be precomputed */
    for (i = row - neighborhood_size; i <= row + neighborhood_size; i++) {
        for (j = col - neighborhood_size; j <= col + neighborhood_size; col++) {
            dist = get_distance(row, col, i, j);
            if (dist > neighborhood_size)
                continue;
            if (devpressure_alg == OCCURRENCE)
                value = 1;
            else if (devpressure_alg == GRAVITY)
                value = scaling_factor / pow(dist, gamma);
            else
                value = scaling_factor * exp(-2 * dist / gamma);
            Segment_get(devpressure, (void *)&devpressure_value, row, col);
            devpressure_value += value;
            Segment_put(devpressure, (void *)&devpressure_value, row, col);
            
        }
    }
}

struct Undeveloped *initialize_undeveloped(int num_subregions)
{
    struct Undeveloped *undev = (struct Undeveloped *) G_malloc(sizeof(struct Undeveloped));
    undev->max_subregions = num_subregions;
    undev->max_undeveloped = (size_t *) G_malloc(undev->max_subregions * sizeof(size_t));
    undev->num_undeveloped = (size_t *) G_calloc(undev->max_subregions, sizeof(size_t));
    undev->cell = (size_t **) G_malloc(undev->max_subregions * sizeof(size_t *));
    for (int i = 0; i < undev->max_subregions; i++){
        undev->max_undeveloped[i] = (Rast_window_rows() * Rast_window_cols()) / num_subregions;
        undev->cell[i] = (size_t *) G_malloc(undev->max_undeveloped[i] * sizeof(size_t));
    }
    return undev;
}

int main(int argc, char **argv)
{
    
    struct
    {
        struct Option
                *developed, *subregions, *predictors, *devpressure, *potentialFile,
                *demandFile, *patchFile, *output, *seed;

    } opt;

    struct
    {
        struct Flag *generateSeed;
    } flg;

    int num_predictors;
    struct KeyValueIntInt *region_map;
    //    int devDemands[MAXNUM_COUNTY][MAX_YEARS];
    struct Undeveloped *undev_cells;
    struct Demand demand_info;
    struct Potential potential_info;
    struct SegmentMemory segment_info;
    struct PatchSizes patch_info;

    G_gisinit(argv[0]);

    struct GModule *module = G_define_module();

    G_add_keyword(_("raster"));
    G_add_keyword(_("patch growing"));
    G_add_keyword(_("urban"));
    G_add_keyword(_("landscape"));
    G_add_keyword(_("modeling"));
    module->label =
            _("Simulates landuse change using FUTure Urban-Regional Environment Simulation (FUTURES).");
    module->description =
            _("Module uses Patch-Growing Algorithm (PGA) to"
              " simulate urban-rural landscape structure development.");
    
    opt.developed = G_define_standard_option(G_OPT_R_INPUT);
    opt.developed->key = "developed";
    opt.developed->required = YES;
    opt.developed->description =
            _("Raster map of developed areas (=1), undeveloped (=0) and excluded (no data)");
    opt.developed->guisection = _("Basic input");

    opt.subregions = G_define_standard_option(G_OPT_R_INPUT);
    opt.subregions->key = "subregions";
    opt.subregions->required = YES;
    opt.subregions->description = _("Raster map of subregions");
    opt.subregions->guisection = _("Basic input");

    opt.predictors = G_define_standard_option(G_OPT_R_INPUTS);
    opt.predictors->key = "predictors";
    opt.predictors->required = YES;
    opt.predictors->multiple = YES;
    opt.predictors->label = _("Names of predictor variable raster maps");
    opt.predictors->description = _("Listed in the same order as in the development potential table");
    opt.predictors->guisection = _("Potential");

    opt.devpressure = G_define_standard_option(G_OPT_R_INPUT);
    opt.devpressure->key = "development_pressure";
    opt.devpressure->required = YES;
    opt.devpressure->description =
            _("Raster map of development pressure");
    opt.devpressure->guisection = _("Development pressure");

    opt.output = G_define_standard_option(G_OPT_R_OUTPUT);
    opt.output->key = "output";
    opt.output->required = YES;
    opt.output->description =
            _("State of the development at the end of simulation");
    opt.output->guisection = _("Output");

    opt.potentialFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.potentialFile->key = "devpot_params";
    opt.potentialFile->required = YES;
    opt.potentialFile->label =
        _("Development potential parameters for each region");
    opt.potentialFile->description =
        _("Each line should contain region ID followed"
          " by parameters (intercepts, development pressure, other predictors)."
          " Values are separated by tabs."
          " First line is ignored, so it can be used for header");
    opt.potentialFile->guisection = _("Potential");
    
    opt.demandFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.demandFile->key = "demand";
    opt.demandFile->required = YES;
    opt.demandFile->description =
            _("Control file with number of cells to convert");
    opt.demandFile->guisection = _("Demand");
    
    opt.patchFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.patchFile->key = "patch_sizes";
    opt.patchFile->required = YES;
    opt.patchFile->description =
        _("File containing list of patch sizes to use");
    opt.patchFile->guisection = _("PGA");

    opt.seed = G_define_option();
    opt.seed->key = "random_seed";
    opt.seed->type = TYPE_INTEGER;
    opt.seed->required = NO;
    opt.seed->label = _("Seed for random number generator");
    opt.seed->description =
            _("The same seed can be used to obtain same results"
              " or random seed can be generated by other means.");
    opt.seed->guisection = _("Random numbers");

    flg.generateSeed = G_define_flag();
    flg.generateSeed->key = 's';
    flg.generateSeed->label =
            _("Generate random seed (result is non-deterministic)");
    flg.generateSeed->description =
            _("Automatically generates random seed for random number"
              " generator (use when you don't want to provide the seed option)");
    flg.generateSeed->guisection = _("Random numbers");

    // TODO: add mutually exclusive?
    // TODO: add flags or options to control values in series and final rasters

    // provided XOR generated
    G_option_exclusive(opt.seed, flg.generateSeed, NULL);
    G_option_required(opt.seed, flg.generateSeed, NULL);
    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    long seed_value;

    if (flg.generateSeed->answer) {
        seed_value = G_srand48_auto();
        G_message("Generated random seed (-s): %ld", seed_value);
    }
    if (opt.seed->answer) {
        seed_value = atol(opt.seed->answer);
        // this does nothing since we are not using GRASS random function
        G_srand48(seed_value);
        G_message("Read random seed from %s option: %ld",
                  opt.seed->key, seed_value);
    }
    // although GRASS random function is initialized
    // the general one must be done separately
    // TODO: replace all by GRASS random number generator?
    srand(seed_value);

    SEGMENT developed_segment;
    SEGMENT subregions_segment;
    SEGMENT devpressure_segment;
    SEGMENT predictors_segment;
    SEGMENT probability_segment;
    segment_info.rows = 64;
    segment_info.cols = 64;
    segment_info.in_memory = 4;

    num_predictors = 0;
    for (int i = 0; opt.predictors->answers[i]; i++)
        num_predictors++;
    
    //    read Subregions layer
    region_map = KeyValueIntInt_create();
    read_subregions(opt.subregions->answer, &subregions_segment, region_map);

    //    development pressure
    rast_segment_open(opt.devpressure->answer, &devpressure_segment,
                      segment_info, FCELL_TYPE);

    /* read Potential file */
    potential_info.filename = opt.potentialFile->answer;
    read_potential_file(&potential_info, region_map, num_predictors);

    /* read Demand file */
    demand_info.filename = opt.demandFile->answer;
    read_demand_file(&demand_info, region_map);
    undev_cells = initialize_undeveloped(region_map->nitems);

    /* read Patch sizes file */
    patch_info.filename = opt.patchFile->answer;
    read_patch_sizes(&patch_info);
    
    //   read developed
    G_verbose_message("Reading input rasters...");
    read_developed(opt.developed->answer, &developed_segment, &subregions_segment,
                   segment_info, undev_cells);

    /* read in predictors */
    read_predictors(opt.predictors->answers, &predictors_segment, &developed_segment,
                    segment_info, num_predictors);

    /* compute initial probability */
    initial_probability(&probability_segment, &developed_segment, &predictors_segment,
                        &devpressure_segment, &subregions_segment, segment_info, &potential_info);

    /* here do the modeling */
    enum slow_grow slow_grow_strategy = FORCE_GROW;
    
    
//    int *added_ids;
//    added_ids = (int *) G_malloc(sizeof(int) * patch_size);
//    grow_patch();
//    update_development_pressure();
//    G_free(added_ids);
    
    /* test predictors */
    int out_fd;
    double prob;
    out_fd = Rast_open_new("test_predictors", FCELL_TYPE);
    void *row_buffer = Rast_allocate_buf(FCELL_TYPE);

    FCELL *values = G_malloc(potential_info.max_predictors * sizeof(FCELL *));
    CELL reg;
    for (int row = 0; row < Rast_window_rows(); row++) {
        for (int col = 0; col < Rast_window_cols(); col++) {
            Segment_get(&subregions_segment, (void *)&reg, row, col);
            prob = get_develop_probability_xy(&predictors_segment, &devpressure_segment,
                                              values,
                                              &potential_info, reg, row, col);
            ((FCELL *) row_buffer)[col] = prob;
        }
        Rast_put_row(out_fd, row_buffer, FCELL_TYPE);
    }
    Rast_close(out_fd);
    G_free(row_buffer);
    G_free(values);
    Segment_close(&predictors_segment);

    /* write */
    int fd = Rast_open_new(opt.output->answer, CELL_TYPE);
    void *out_buffer = Rast_allocate_buf(CELL_TYPE);

    for (int row = 0; row < Rast_window_rows(); row++) {
        Segment_get_row(&developed_segment, out_buffer, row);
        Rast_put_row(fd, out_buffer, CELL_TYPE);
    }
    Rast_close(fd);
    G_free(out_buffer);
    Segment_close(&developed_segment);
    Segment_close(&subregions_segment);
    Segment_close(&devpressure_segment);
    Segment_close(&probability_segment);

    KeyValueIntInt_free(region_map);
    if (demand_info.table) {
        for (int i = 0; i < demand_info.max_subregions; i++)
            G_free(demand_info.table[i]);
        G_free(demand_info.table);
    }
    if (potential_info.predictors) {
        for (int i = 0; i < potential_info.max_predictors; i++)
            G_free(potential_info.predictors[i]);
        G_free(potential_info.predictors);
        G_free(potential_info.devpressure);
        G_free(potential_info.intercept);
    }
    if (undev_cells) {
        G_free(undev_cells->num_undeveloped);
        G_free(undev_cells->max_undeveloped);
        for (int i = 0; i < undev_cells->max_subregions; i++)
            G_free(undev_cells->cell[i]);
        G_free(undev_cells->cell);
        G_free(undev_cells);
    }

    G_free(patch_info.patch_sizes);

    return EXIT_SUCCESS;
}

