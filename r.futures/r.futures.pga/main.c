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
#include "output.h"
#include "patch.h"






enum seed_search {RANDOM, PROBABILITY};



int find_probable_seed(struct Undeveloped *undev_cells, int region)
{
    int first, last, middle;
    double p;

    p = G_drand48();
    // bisect
    first = 0;
    last = undev_cells->num[region] - 1;
    middle = (first + last) / 2;
    if (p <= undev_cells->cells[region][first].cumulative_probability)
        return 0;
    if (p >= undev_cells->cells[region][last].cumulative_probability)
        return last;
    while (first <= last) {
        if (undev_cells->cells[region][middle].cumulative_probability < p)
            first = middle + 1;
        else if (undev_cells->cells[region][middle - 1].cumulative_probability < p &&
                 undev_cells->cells[region][middle].cumulative_probability >= p) {
            return middle;
        }
        else
            last = middle - 1;
        middle = (first + last)/2;
    }
    // TODO: returning at least something but should be something more meaningful
    return 0;
}


int get_seed(struct Undeveloped *undev_cells, int region_idx, enum seed_search method,
              int *row, int *col)
{
    int i, id;
    if (method == RANDOM)
        i = (int)(G_drand48() * undev_cells->max[region_idx]);
    else
        i = find_probable_seed(undev_cells, region_idx);
    id = undev_cells->cells[region_idx][i].id;
    get_xy_from_idx(id, Rast_window_cols(), row, col);
    return i;
}



double get_develop_probability_xy(struct Segments *segments,
                                  FCELL *values,
                                  struct Potential *potential_info,
                                  int region_index, int row, int col)
{
    double probability;
    int i;
    FCELL devpressure_val;
    Segment_get(&segments->devpressure, (void *)&devpressure_val, row, col);
    Segment_get(&segments->predictors, values, row, col);
    
    probability = potential_info->intercept[region_index];
    probability += potential_info->devpressure[region_index] * devpressure_val;
    for (i = 0; i < potential_info->max_predictors; i++) {
        probability += potential_info->predictors[i][region_index] * values[i];
    }
    probability = 1.0 / (1.0 + exp(-probability));
    return probability;
}


void update_development_pressure(int row, int col, struct Segments *segments,
                                 struct DevPressure *devpressure_info) {
    int i, j;
    int cols, rows;
    double dist;
    double value;
    FCELL devpressure_value;

    cols = Rast_window_cols();
    rows = Rast_window_rows();
    /* this can be precomputed */
    for (i = row - devpressure_info->neighborhood; i <= row + devpressure_info->neighborhood; i++) {
        for (j = col - devpressure_info->neighborhood; j <= col + devpressure_info->neighborhood; j++) {
            if (i < 0 || j < 0 || i >= rows || j >= cols)
                continue;
            dist = get_distance(row, col, i, j);
            if (dist > devpressure_info->neighborhood)
                continue;
            if (devpressure_info->alg == OCCURRENCE)
                value = 1;
            else if (devpressure_info->alg == GRAVITY)
                value = devpressure_info->scaling_factor / pow(dist, devpressure_info->gamma);
            else
                value = devpressure_info->scaling_factor * exp(-2 * dist / devpressure_info->gamma);
            Segment_get(&segments->devpressure, (void *)&devpressure_value, i, j);
            if (Rast_is_null_value(&devpressure_value, FCELL_TYPE))
                continue;
            devpressure_value += value;
            Segment_put(&segments->devpressure, (void *)&devpressure_value, i, j);
            
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
    for (row = 0; row < rows; row++) {
        for (col = 0; col < cols; col++) {
            Segment_get(&segments->developed, (void *)&developed, row, col);
            if (Rast_is_null_value(&developed, CELL_TYPE))
                continue;
            if (developed != -1)
                continue;
            Segment_get(&segments->subregions, (void *)&region, row, col);
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
            undeveloped_cells->cells[region][idx].tried = 0;
            /* get probability and update undevs and segment*/
            probability = get_develop_probability_xy(segments, values,
                                                     potential_info, region, row, col);
            Segment_put(&segments->probability, (void *)&probability, row, col);
            undeveloped_cells->cells[region][idx].probability = probability;
            
            undeveloped_cells->num[region]++;
            
        }
    }

    i = 0;
    for (region_idx = 0; region_idx < undeveloped_cells->max_subregions; region_idx++) {
        probability = undeveloped_cells->cells[region_idx][0].probability;
        undeveloped_cells->cells[region_idx][0].cumulative_probability = probability;
        for (i = 1; i < undeveloped_cells->num[region_idx]; i++) {
            probability = undeveloped_cells->cells[region_idx][i].probability;
            undeveloped_cells->cells[region_idx][i].cumulative_probability = 
                    undeveloped_cells->cells[region_idx][i - 1].cumulative_probability + probability;
        }
        sum = undeveloped_cells->cells[region_idx][i - 1].cumulative_probability;
        for (i = 0; i < undeveloped_cells->num[region_idx]; i++) {
            undeveloped_cells->cells[region_idx][i].cumulative_probability /= sum;
        }
    }
}

void compute_step(struct Undeveloped *undev_cells, struct Demand *demand,
                  enum seed_search search_alg,
                  struct Segments *segments,
                  struct PatchSizes *patch_sizes, struct PatchInfo *patch_info,
                  struct DevPressure *devpressure_info, int *patch_overflow,
                  int step, int region)
{
    int i, idx;
    int n_to_convert;
    int n_done;
    int found;
    int seed_row, seed_col;
    int row, col;
    int patch_size;
    int *added_ids;
    bool force_convert_all;
    int extra;
    bool allow_already_tried_ones;
    int unsuccessful_tries;
    FCELL prob;
    CELL developed;


    added_ids = (int *) G_malloc(sizeof(int) * patch_sizes->max_patch_size);
    n_to_convert = demand->table[region][step];
    n_done = 0;
    force_convert_all = false;
    allow_already_tried_ones = false;
    unsuccessful_tries = 0;
    extra = patch_overflow[region];

    if (extra > 0) {
        if (n_to_convert - extra > 0) {
            n_to_convert -= extra;
            extra = 0;
        }
        else {
            extra -= n_to_convert;
            n_to_convert = 0;
        }
    }

    if (n_to_convert > undev_cells->num[region]) {
        G_warning("Not enough undeveloped cells (requested: %d,"
                  " available: %d). Converting all available.",
                   n_to_convert, undev_cells->num[region]);
        n_to_convert =  undev_cells->num[region];
        force_convert_all = true;
    }
    
    while (n_done < n_to_convert) {
        /* if we can't find a seed, turn off the restriction to use only untried ones */
        if (!allow_already_tried_ones && unsuccessful_tries > MAX_SEED_ITER * n_to_convert)
            allow_already_tried_ones = true;

        /* get seed's row, col and index in undev cells array */
        idx = get_seed(undev_cells, region, search_alg, &seed_row, &seed_col);
        /* skip if seed was already tried unless we switched of this check because we can't get any seed */
        if (!allow_already_tried_ones && undev_cells->cells[region][idx].tried) {
            unsuccessful_tries++;
            continue;
        }
        /* mark as tried */
        undev_cells->cells[region][idx].tried = 1;
        /* see if seed was already developed during this time step */
        Segment_get(&segments->developed, (void *)&developed, seed_row, seed_col);
        if (developed != -1) {
            unsuccessful_tries++;
            continue;
        }
        /* get probability */
        Segment_get(&segments->probability, (void *)&prob, seed_row, seed_col);
        /* challenge probability unless we need to convert all */
        if(force_convert_all || G_drand48() < prob) {
            /* ger random patch size */
            patch_size = get_patch_size(patch_sizes);
            /* grow patch and return the actual grown size which could be smaller */
            found = grow_patch(seed_row, seed_col, patch_size, step,
                               patch_info, segments, added_ids);
            /* update devpressure for every newly developed cell */
            for (i = 0; i < found; i++) {
                get_xy_from_idx(added_ids[i], Rast_window_cols(), &row, &col);
                update_development_pressure(row, col, segments, devpressure_info);
            }
            n_done += found;
        }
    }
    extra += (n_done - n_to_convert);
    patch_overflow[region] = extra;
    G_debug(2, "There are %d extra cells for next timestep", extra);
    G_free(added_ids);
}

int main(int argc, char **argv)
{

    struct
    {
        struct Option
                *developed, *subregions, *predictors,
                *devpressure, *nDevNeighbourhood, *devpressureApproach, *scalingFactor, *gamma,
                *potentialFile, *numNeighbors, *discountFactor, *seedSearch,
                *patchMean, *patchRange,
                *demandFile, *patchFile, *numSteps, *output, *outputSeries, *seed;

    } opt;

    struct
    {
        struct Flag *generateSeed;
    } flg;

    int i;
    int num_predictors;
    int num_steps;
    int region;
    int step;
    double discount_factor;
    enum seed_search search_alg;
    struct KeyValueIntInt *region_map;
    struct Undeveloped *undev_cells;
    struct Demand demand_info;
    struct Potential potential_info;
    struct SegmentMemory segment_info;
    struct PatchSizes patch_sizes;
    struct PatchInfo patch_info;
    struct DevPressure devpressure_info;
    struct Segments segments;
    int *patch_overflow;
    char *name_step;

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

    opt.outputSeries = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.outputSeries->key = "output_series";
    opt.outputSeries->required = NO;
    opt.outputSeries->label =
        _("Basename for raster maps of development generated after each step");
    opt.outputSeries->guisection = _("Output");

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

    opt.seedSearch = G_define_option();
    opt.seedSearch->key = "seed_search";
    opt.seedSearch->type = TYPE_STRING;
    opt.seedSearch->required = YES;
    opt.seedSearch->options = "random,probability";
    opt.seedSearch->answer = "probability";
    opt.seedSearch->description =
        _("The way location of a seed is determined (1: uniform distribution 2: development probability)");
    opt.seedSearch->guisection = _("PGA");
    
    opt.patchMean = G_define_option();
    opt.patchMean->key = "compactness_mean";
    opt.patchMean->type = TYPE_DOUBLE;
    opt.patchMean->required = YES;
    opt.patchMean->description =
        _("Mean value of patch compactness to control patch shapes");
    opt.patchMean->guisection = _("PGA");

    opt.patchRange = G_define_option();
    opt.patchRange->key = "compactness_range";
    opt.patchRange->type = TYPE_DOUBLE;
    opt.patchRange->required = YES;
    opt.patchRange->description =
        _("Range of patch compactness to control patch shapes");
    opt.patchRange->guisection = _("PGA");

    opt.numSteps = G_define_option();
    opt.numSteps->key = "num_steps";
    opt.numSteps->type = TYPE_INTEGER;
    opt.numSteps->required = NO;
    opt.numSteps->description =
        _("Number of steps to be simulated");
    opt.numSteps->guisection = _("Basic input");

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

    devpressure_info.scaling_factor = atof(opt.scalingFactor->answer);
    devpressure_info.gamma = atof(opt.gamma->answer);
    devpressure_info.neighborhood = atoi(opt.nDevNeighbourhood->answer);
    discount_factor = atof(opt.discountFactor->answer);
    if (strcmp(opt.devpressureApproach->answer, "occurrence") == 0)
        devpressure_info.alg = OCCURRENCE;
    else if (strcmp(opt.devpressureApproach->answer, "gravity") == 0)
        devpressure_info.alg = GRAVITY;
    else if (strcmp(opt.devpressureApproach->answer, "kernel") == 0)
        devpressure_info.alg = KERNEL;
    else
        G_fatal_error(_("Approach doesn't exist"));

    if (strcmp(opt.seedSearch->answer, "random") == 0)
        search_alg = RANDOM;
    else if (strcmp(opt.seedSearch->answer, "probability") == 0)
        search_alg = PROBABILITY;
    else
        G_fatal_error(_("Approach doesn't exist"));
    patch_info.compactness_mean = atof(opt.patchMean->answer);
    patch_info.compactness_range = atof(opt.patchRange->answer);
    patch_info.num_neighbors = atoi(opt.numNeighbors->answer);
    patch_info.strategy = SKIP;
    
    num_steps = 0;
    if (opt.numSteps->answer)
        num_steps = atoi(opt.numSteps->answer);
    
    segment_info.rows = 64;
    segment_info.cols = 64;
    segment_info.in_memory = 1;

    num_predictors = 0;
    for (i = 0; opt.predictors->answers[i]; i++)
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
 
    if (num_steps == 0)
        num_steps = demand_info.max_steps;
    
    undev_cells = initialize_undeveloped(region_map->nitems);
    patch_overflow = G_calloc(region_map->nitems, sizeof(int));

    /* read Patch sizes file */
    patch_sizes.filename = opt.patchFile->answer;
    read_patch_sizes(&patch_sizes, discount_factor);

    G_verbose_message("Reading input rasters...");
    //    development pressure
    rast_segment_open(opt.devpressure->answer, &segments.devpressure,
                      segment_info, FCELL_TYPE);

    //   read developed
    read_developed(opt.developed->answer, &segments, segment_info);

    /* read in predictors */
    read_predictors(opt.predictors->answers, &segments,
                    segment_info, num_predictors);

    /* create probability segment*/
    if (Segment_open(&segments.probability, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segment_info.rows, segment_info.cols,
                     Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));

    /* here do the modeling */

    for (step = 0; step < num_steps; step++) {
        recompute_probabilities(undev_cells, &segments, &potential_info);
        for (region = 0; region < region_map->nitems; region++) {
            compute_step(undev_cells, &demand_info, search_alg, &segments,
                         &patch_sizes, &patch_info, &devpressure_info, patch_overflow,
                         step, region);
        }
        /* export developed for that step */
        if (opt.outputSeries->answer) {
            name_step = name_for_step(opt.outputSeries->answer, step, num_steps);
            output_developed_step(&segments.developed, name_step, num_steps, true, true);
        }
    }

    /* test predictors */
    int out_fd;
    double prob;
    out_fd = Rast_open_new("test_predictors", FCELL_TYPE);
    void *row_buffer = Rast_allocate_buf(FCELL_TYPE);

    FCELL *values = G_malloc(potential_info.max_predictors * sizeof(FCELL *));
    CELL reg;
    for (int row = 0; row < Rast_window_rows(); row++) {
        for (int col = 0; col < Rast_window_cols(); col++) {
            Segment_get(&segments.subregions, (void *)&reg, row, col);
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

    /* write */
    output_developed_step(&segments.developed, opt.output->answer, num_steps, false, false);

    Segment_close(&segments.developed);
    Segment_close(&segments.subregions);
    Segment_close(&segments.devpressure);
    Segment_close(&segments.probability);
    Segment_close(&segments.predictors);

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

    G_free(patch_sizes.patch_sizes);
    G_free(patch_overflow);

    return EXIT_SUCCESS;
}

