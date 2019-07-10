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
enum seed_search {RANDOM, SUITABILITY};




void get_seed(struct Undeveloped *undev_cells, int region_idx, enum seed_search method,
              int *row, int *col)
{
    int i;
    if (method == RANDOM) {
        i = (int)(G_drand48() * undev_cells->max[region_idx]);
        get_xy_from_idx(i, Rast_window_cols(), row, col);
    }
    else {
        
    }
}



double get_develop_probability_xy(struct Segments *segments,
                                  FCELL *values,
                                  struct Potential *potential_info,
                                  int region_index, int row, int col)
{
    double probability;
    int i;
    FCELL devpressure_val;
    Segment_get(&segments->devpressure_segment, (void *)&devpressure_val, row, col);
    Segment_get(&segments->predictors_segment, values, row, col);
    
    probability = potential_info->intercept[region_index];
    probability += potential_info->devpressure[region_index] * devpressure_val;
    for (i = 0; i < potential_info->max_predictors; i++) {
        probability += potential_info->predictors[i][region_index] * values[i];
    }
    probability = 1.0 / (1.0 + exp(-probability));
    return probability;
}


void compute_develop_probability(struct Segments *segments,
                                 struct Potential *potential_info)
{
    
    int row, col;
    float prob;

    CELL developed_cell;
    CELL region_index;
    FCELL *values = G_malloc(potential_info->max_predictors * sizeof(FCELL *));
    for (row = 0; row < Rast_window_rows(); row++) {
        for (col = 0; col < Rast_window_cols(); col++) {
            Segment_get(&segments->developed_segment, (void *)&developed_cell, row, col);
            if (Rast_is_null_value(&developed_cell, CELL_TYPE)) {
                Rast_set_null_value(values, 1, FCELL_TYPE);
                Segment_put(&segments->probability_segment, (void *)values, row, col);
                continue;
            }
            if (developed_cell == -1) {
                Segment_get(&segments->subregions_segment, (void *)&region_index, row, col);
                if (Rast_is_null_value(&region_index, CELL_TYPE)) {
                    Rast_set_null_value(values, 1, FCELL_TYPE);
                    Segment_put(&segments->probability_segment, (void *)values, row, col);
                    continue;
                }
                prob = get_develop_probability_xy(segments, values,
                                                  potential_info, region_index, row, col);
                Segment_put(&segments->probability_segment, (void *)&prob, row, col);
            }
        }
    }
    Segment_flush(&segments->probability_segment);
    G_free(values);
}


void initial_probability(struct Segments *segments,  struct SegmentMemory segmentInfo,
                         struct Potential *potential_info)
{
    if (Segment_open(&segments->probability_segment, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segmentInfo.rows, segmentInfo.cols,
                     Rast_cell_size(FCELL_TYPE), segmentInfo.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));
    compute_develop_probability(segments, potential_info);
    Segment_flush(&segments->probability_segment);
}


void update_development_pressure(int row, int col, struct Segments *segments,
                                 int neighborhood_size, double gamma, double scaling_factor,
                                 enum development_pressure devpressure_alg) {
    int i, j;
    double dist;
    double value;
    FCELL devpressure_value;

    /* this can be precomputed */
    for (i = row - neighborhood_size; i <= row + neighborhood_size; i++) {
        for (j = col - neighborhood_size; j <= col + neighborhood_size; j++) {
            dist = get_distance(row, col, i, j);
            if (dist > neighborhood_size)
                continue;
            if (devpressure_alg == OCCURRENCE)
                value = 1;
            else if (devpressure_alg == GRAVITY)
                value = scaling_factor / pow(dist, gamma);
            else
                value = scaling_factor * exp(-2 * dist / gamma);
            Segment_get(&segments->devpressure_segment, (void *)&devpressure_value, i, j);
            devpressure_value += value;
            Segment_put(&segments->devpressure_segment, (void *)&devpressure_value, i, j);
            
        }
    }
}

struct Undeveloped *initialize_undeveloped(int num_subregions)
{
    struct Undeveloped *undev = (struct Undeveloped *) G_malloc(sizeof(struct Undeveloped));
    undev->max_subregions = num_subregions;
    undev->max = (size_t *) G_malloc(undev->max_subregions * sizeof(size_t));
    undev->num = (size_t *) G_calloc(undev->max_subregions, sizeof(size_t));
    undev->cells = (struct UndevelopedCell **) G_malloc(undev->max_subregions * sizeof(struct UndevelopedCell *));
    for (int i = 0; i < undev->max_subregions; i++){
        undev->max[i] = (Rast_window_rows() * Rast_window_cols()) / num_subregions;
        undev->cells[i] = (struct UndevelopedCell *) G_malloc(undev->max[i] * sizeof(struct UndevelopedCell));
    }
    return undev;
}


void recompute_probabilities(struct Undeveloped *undeveloped_cells,
                             struct Segments *segments,
                             struct Potential *potential_info)
{
    int row, col, cols, rows;
    int id, i, idx, new_size;
    int region_idx;
    CELL developed;
    CELL region;
    FCELL *values;
    float probability;
    float sum;
    
    cols = Rast_window_cols();
    rows = Rast_window_rows();
    values = G_malloc(potential_info->max_predictors * sizeof(FCELL *));
    
    for (region_idx = 0; region_idx < undeveloped_cells->max_subregions; region_idx++) {
        undeveloped_cells->num[region_idx] = 0;
    }
    i = 0;
    for (row = 0; row < rows; row++) {
        for (col = 0; col < cols; col++) {
            Segment_get(&segments->developed_segment, (void *)&developed, row, col);
            if (Rast_is_null_value(&developed, CELL_TYPE))
                continue;
            if (developed != -1)
                continue;
            Segment_get(&segments->subregions_segment, (void *)&region, row, col);
            if (Rast_is_null_value(&region, CELL_TYPE))
                continue;
            
            /* realloc if needed */
            if (undeveloped_cells->num[region] >= undeveloped_cells->max[region]) {
                new_size = 2 * undeveloped_cells->max[region];
                undeveloped_cells->cells[region] = 
                        (struct UndevelopedCell *) G_realloc(undeveloped_cells->cells[region],
                                                             new_size * sizeof(struct UndevelopedCell));
                undeveloped_cells->max[region] = new_size;
            }
            
            id = get_idx_from_xy(row, col, cols);
            idx = undeveloped_cells->num[region];
            undeveloped_cells->cells[region][idx].id = id;
            /* get probability and update undevs and segment*/
            probability = get_develop_probability_xy(segments, values,
                                                     potential_info, region, row, col);
            Segment_put(&segments->probability_segment, (void *)&probability, row, col);
            undeveloped_cells->cells[region][idx].probability = probability;
            
            undeveloped_cells->num[region]++;
            
        }
    }

    for (region_idx = 0; region_idx < undeveloped_cells->max_subregions; region_idx++) {
        probability = undeveloped_cells->cells[region][0].probability;
        undeveloped_cells->cells[region][0].cumulative_probability = probability;
        for (i = 1; i < undeveloped_cells->num[region]; i++) {
            probability = undeveloped_cells->cells[region][i].probability;
            undeveloped_cells->cells[region][i].cumulative_probability += probability;
        }
        sum = undeveloped_cells->cells[region][i].cumulative_probability;
        for (i = 0; i < undeveloped_cells->num[region]; i++) {
            undeveloped_cells->cells[region][i].cumulative_probability /= sum;
        }
    }
}

int main(int argc, char **argv)
{

    struct
    {
        struct Option
                *developed, *subregions, *predictors,
                *devpressure, *nDevNeighbourhood, *devpressureApproach, *scalingFactor, *gamma,
                *potentialFile, *numNeighbors, *discountFactor,
                *demandFile, *patchFile, *output, *seed;

    } opt;

    struct
    {
        struct Flag *generateSeed;
    } flg;

    int num_predictors, num_neighbors;
    double scaling_factor, gamma;
    double discount_factor;
    enum development_pressure devpressure_alg;
    int devpressure_neighborhood;
    struct KeyValueIntInt *region_map;
    struct Undeveloped *undev_cells;
    struct Demand demand_info;
    struct Potential potential_info;
    struct SegmentMemory segment_info;
    struct PatchSizes patch_info;
    struct Segments segments;

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

    opt.nDevNeighbourhood = G_define_option();
    opt.nDevNeighbourhood->key = "n_dev_neighbourhood";
    opt.nDevNeighbourhood->type = TYPE_INTEGER;
    opt.nDevNeighbourhood->required = YES;
    opt.nDevNeighbourhood->description =
        _("Size of square used to recalculate development pressure");
    opt.nDevNeighbourhood->guisection = _("Development pressure");

    opt.devpressureApproach = G_define_option();
    opt.devpressureApproach->key = "development_pressure_approach";
    opt.devpressureApproach->type = TYPE_STRING;
    opt.devpressureApproach->required = YES;
    opt.devpressureApproach->options = "occurrence,gravity,kernel";
    opt.devpressureApproach->description =
        _("Approaches to derive development pressure");
    opt.devpressureApproach->answer = "gravity";
    opt.devpressureApproach->guisection = _("Development pressure");

    opt.gamma = G_define_option();
    opt.gamma->key = "gamma";
    opt.gamma->type = TYPE_DOUBLE;
    opt.gamma->required = YES;
    opt.gamma->description =
        _("Influence of distance between neighboring cells");
    opt.gamma->guisection = _("Development pressure");

    opt.scalingFactor = G_define_option();
    opt.scalingFactor->key = "scaling_factor";
    opt.scalingFactor->type = TYPE_DOUBLE;
    opt.scalingFactor->required = YES;
    opt.scalingFactor->description =
        _("Scaling factor of development pressure");
    opt.scalingFactor->guisection = _("Development pressure");

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

    opt.numNeighbors = G_define_option();
    opt.numNeighbors->key = "num_neighbors";
    opt.numNeighbors->type = TYPE_INTEGER;
    opt.numNeighbors->required = YES;
    opt.numNeighbors->options = "4,8";
    opt.numNeighbors->answer = "4";
    opt.numNeighbors->description =
        _("The number of neighbors to be used for patch generation (4 or 8)");
    opt.numNeighbors->guisection = _("PGA");

    opt.discountFactor = G_define_option();
    opt.discountFactor->key = "discount_factor";
    opt.discountFactor->type = TYPE_DOUBLE;
    opt.discountFactor->required = YES;
    opt.discountFactor->description = _("Discount factor of patch size");
    opt.discountFactor->guisection = _("PGA");

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
//    srand(seed_value);

    num_neighbors = atoi(opt.numNeighbors->answer);
    scaling_factor = atof(opt.scalingFactor->answer);
    gamma = atof(opt.gamma->answer);
    devpressure_neighborhood = atoi(opt.nDevNeighbourhood->answer);
    discount_factor = atof(opt.discountFactor->answer);
    if (strcmp(opt.devpressureApproach->answer, "occurrence") == 0)
        devpressure_alg = OCCURRENCE;
    else if (strcmp(opt.devpressureApproach->answer, "gravity") == 0)
        devpressure_alg = GRAVITY;
    else if (strcmp(opt.devpressureApproach->answer, "kernel") == 0)
        devpressure_alg = KERNEL;
    else
        G_fatal_error(_("Approach doesn't exist"));

    segment_info.rows = 64;
    segment_info.cols = 64;
    segment_info.in_memory = 1;

    num_predictors = 0;
    for (int i = 0; opt.predictors->answers[i]; i++)
        num_predictors++;

    //    read Subregions layer
    region_map = KeyValueIntInt_create();
    read_subregions(opt.subregions->answer, &segments, region_map);
    
    /* read Potential file */
    potential_info.filename = opt.potentialFile->answer;
    read_potential_file(&potential_info, region_map, num_predictors);

    /* read Demand file */
    demand_info.filename = opt.demandFile->answer;
    read_demand_file(&demand_info, region_map);
    undev_cells = initialize_undeveloped(region_map->nitems);

    /* read Patch sizes file */
    patch_info.filename = opt.patchFile->answer;
    read_patch_sizes(&patch_info, discount_factor);

    G_verbose_message("Reading input rasters...");
    //    development pressure
    rast_segment_open(opt.devpressure->answer, &segments.devpressure_segment,
                      segment_info, FCELL_TYPE);

    //   read developed
    read_developed(opt.developed->answer, &segments, segment_info);

    /* read in predictors */
    read_predictors(opt.predictors->answers, &segments,
                    segment_info, num_predictors);

    /* compute initial probability */
    initial_probability(&segments, segment_info, &potential_info);

    /* here do the modeling */
    enum slow_grow slow_grow_strategy = SKIP;
    
    recompute_probabilities(undev_cells, &segments, &potential_info);


    int *added_ids;
    int patch_size = 100;
    int seed_row = 352;
    int seed_col = 111;
    int step = 1;
    double alpha = 10;
    int rowx, colx;
    int found;
    added_ids = (int *) G_malloc(sizeof(int) * patch_size);
    found = grow_patch(seed_row, seed_col, added_ids, &segments,
                       num_neighbors, alpha, patch_size, step, slow_grow_strategy);
    for (int i = 0; i < found; i++) {
        get_xy_from_idx(added_ids[i], Rast_window_cols(), &rowx, &colx);
        update_development_pressure(rowx, colx, &segments, devpressure_neighborhood,
                                    gamma, scaling_factor, devpressure_alg);

    }
    G_free(added_ids);

    /* test predictors */
    int out_fd;
    double prob;
    out_fd = Rast_open_new("test_predictors", FCELL_TYPE);
    void *row_buffer = Rast_allocate_buf(FCELL_TYPE);

    FCELL *values = G_malloc(potential_info.max_predictors * sizeof(FCELL *));
    CELL reg;
    for (int row = 0; row < Rast_window_rows(); row++) {
        for (int col = 0; col < Rast_window_cols(); col++) {
            Segment_get(&segments.subregions_segment, (void *)&reg, row, col);
            if (Rast_is_null_value(&reg, CELL_TYPE))
                Rast_set_null_value(&((FCELL *) row_buffer)[col], 1, CELL_TYPE);
            else {
                prob = get_develop_probability_xy(&segments,
                                                  values,
                                                  &potential_info, reg, row, col);
                ((FCELL *) row_buffer)[col] = prob;
            }
        }
        Rast_put_row(out_fd, row_buffer, FCELL_TYPE);
    }
    Rast_close(out_fd);
    G_free(row_buffer);
    G_free(values);
    Segment_close(&segments.predictors_segment);

    /* write */
    int fd = Rast_open_new(opt.output->answer, CELL_TYPE);
    void *out_buffer = Rast_allocate_buf(CELL_TYPE);
    Segment_flush(&segments.developed_segment);
    for (int row = 0; row < Rast_window_rows(); row++) {
        Segment_get_row(&segments.developed_segment, out_buffer, row);
        Rast_put_row(fd, out_buffer, CELL_TYPE);
    }
    Rast_close(fd);
    G_free(out_buffer);
    Segment_close(&segments.developed_segment);
    Segment_close(&segments.subregions_segment);
    Segment_close(&segments.devpressure_segment);
    Segment_close(&segments.probability_segment);

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
        G_free(undev_cells->num);
        G_free(undev_cells->max);
        for (int i = 0; i < undev_cells->max_subregions; i++)
            G_free(undev_cells->cells[i]);
        G_free(undev_cells->cells);
        G_free(undev_cells);
    }

    G_free(patch_info.patch_sizes);

    return EXIT_SUCCESS;
}

