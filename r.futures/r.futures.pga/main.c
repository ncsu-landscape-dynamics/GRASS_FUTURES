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
#include "devpressure.h"
#include "simulation.h"
#include "redistribute.h"
#include "climate.h"


struct Developables *initialize_developables(int num_subregions)
{
    struct Developables *dev = (struct Developables *) G_malloc(sizeof(struct Developables));
    dev->max_subregions = num_subregions;
    dev->max = (size_t *) G_malloc(dev->max_subregions * sizeof(size_t));
    dev->num = (size_t *) G_calloc(dev->max_subregions, sizeof(size_t));
    dev->cells = (struct DevelopableCell **) G_malloc(dev->max_subregions * sizeof(struct DevelopableCell *));
    for (int i = 0; i < dev->max_subregions; i++){
        dev->max[i] = (Rast_window_rows() * Rast_window_cols()) / num_subregions;
        dev->cells[i] = (struct DevelopableCell *) G_malloc(dev->max[i] * sizeof(struct DevelopableCell));
    }
    return dev;
}

void initialize_flood_response(struct ACDamageRelation *response_relation)
{
    response_relation->resilience_a = 1;
    response_relation->resilience_b = 0;
    response_relation->vulnerability_a = -1;
    response_relation->vulnerability_b = 0;
}

static int manage_memory(struct SegmentMemory *memory, struct Segments *segments,
                         float input_memory)
{
    int nseg, nseg_total;
    int cols, rows;
    size_t undev_size;
    size_t size;
    size_t estimate;

    memory->rows = 64;
    memory->cols = 64;
    rows = Rast_window_rows();
    cols = Rast_window_cols();

    undev_size = (sizeof(size_t) + sizeof(float) * 2 + sizeof(bool)) * rows * cols;
    estimate = undev_size;
    if (segments->use_density)
        estimate += undev_size;

    if (input_memory > 0 && undev_size > 1e9 * input_memory)
        G_warning(_("Not sufficient memory, will attempt to use more "
                    "than specified. Will need at least %d MB"), (int) (undev_size / 1.0e6));

    /* developed, subregions */
    size = sizeof(CELL) * 2;
    /* predictors, devpressure, probability */
    size += sizeof(FCELL) * 3;
    if (segments->use_weight)
        size += sizeof(FCELL);
    if (segments->use_potential_subregions)
        size += sizeof(CELL);
    /* density + capacity */
    if (segments->use_density)
        size += sizeof(FCELL) * 2;
    /* climate: HAND (F) + AC (F) + flood prob (F) + HUC (C) + adaptation (C) */
    if (segments->use_climate) {
        size += sizeof(FCELL) * 3;
        size += sizeof(CELL) * 2;
    }
    estimate += size * rows * cols;
    size *= memory->rows * memory->cols;

    nseg = (1e9 * input_memory - undev_size) / size;
    if (nseg <= 0)
        nseg = 1;
    nseg_total = (rows / memory->rows + (rows % memory->rows > 0)) *
                 (cols / memory->cols + (cols % memory->cols > 0));

    if (nseg > nseg_total || input_memory < 0)
        nseg = nseg_total;
    G_verbose_message(_("Number of segments in memory: %d of %d total"),
                      nseg, nseg_total);
    G_verbose_message(_("Estimated minimum memory footprint without using disk cache: %d MB"),
                      (int) (estimate / 1.0e6));
    return nseg;
}

int main(int argc, char **argv)
{

    struct
    {
        struct Option
                *developed, *subregions, *potentialSubregions, *predictors,
                *devpressure, *nDevNeighbourhood, *devpressureApproach, *scalingFactor, *gamma,
                *potentialFile, *numNeighbors, *discountFactor, *seedSearch,
                *patchMean, *patchRange,
                *incentivePower, *potentialWeight,
                *cellDemandFile, *populationDemandFile, *separator,
                *density, *densityCapacity, *outputDensity, *redevelopmentLag,
                *redevelopmentPotentialFile, *redistributionMatrix, *redistributionMatrixOutput,
                *HAND, *floodProbability, *depthDamageFunc, *ddf_subregions,
                *adaptiveCapacity, *HUCs, *outputAdaptation,
                *patchFile, *numSteps, *output, *outputSeries, *seed, *memory;
    } opt;

    struct
    {
        struct Flag *generateSeed;
    } flg;

    int i;
    int num_predictors;
    int num_steps;
    int nseg;
    int region;
    int HUC;
    int step;
    float memory;
    double discount_factor;
    float exponent;
    enum seed_search search_alg;
    struct RasterInputs raster_inputs;
    struct KeyValueIntInt *region_map;
    struct KeyValueIntInt *reverse_region_map;
    struct KeyValueIntInt *potential_region_map;
    struct KeyValueCharInt *predictor_map;
    struct KeyValueIntFloat *max_flood_probability_map;
    struct KeyValueIntInt *HUC_map;
    struct KeyValueIntInt * DDF_region_map;
    struct Developables *undev_cells;
    struct Developables *dev_cells;
    struct Demand demand_info;
    struct Potential potential_info;
    struct Potential redev_potential_info;
    struct SegmentMemory segment_info;
    struct PatchSizes patch_sizes;
    struct PatchInfo patch_info;
    struct DevPressure devpressure_info;
    struct Segments segments;
    struct RedistributionMatrix redistr_matrix;
    struct BBoxes bboxes;
    struct DepthDamageFunctions DDF;
    struct ACDamageRelation response_relation;
    int *patch_overflow;
    float *population_overflow;
    char *name_step;
    float leaving_population;
    bool overgrow;

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

    opt.potentialSubregions = G_define_standard_option(G_OPT_R_INPUT);
    opt.potentialSubregions->key = "subregions_potential";
    opt.potentialSubregions->required = NO;
    opt.potentialSubregions->label = _("Raster map of subregions used with potential file");
    opt.potentialSubregions->description = _("If not specified, the raster specified in subregions parameter is used");
    opt.potentialSubregions->guisection = _("Potential");

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

    opt.density = G_define_standard_option(G_OPT_R_INPUT);
    opt.density->key = "density";
    opt.density->required = NO;
    opt.density->description =
            _("Raster map of population density");
    opt.density->guisection = _("Density");

    opt.densityCapacity = G_define_standard_option(G_OPT_R_INPUT);
    opt.densityCapacity->key = "density_capacity";
    opt.densityCapacity->required = NO;
    opt.densityCapacity->description =
            _("Raster map of maximum capacity");
    opt.densityCapacity->guisection = _("Density");

    opt.redevelopmentPotentialFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.redevelopmentPotentialFile->key = "redevpot_params";
    opt.redevelopmentPotentialFile->required = NO;
    opt.redevelopmentPotentialFile->label =
        _("CSV file with redevelopment potential parameters for each region");
    opt.redevelopmentPotentialFile->description =
        _("Each line should contain region ID followed"
          " by parameters (intercepts, development pressure, other predictors).");
    opt.redevelopmentPotentialFile->guisection = _("Density");

    opt.redevelopmentLag = G_define_option();
    opt.redevelopmentLag->key = "redevelopment_lag";
    opt.redevelopmentLag->type = TYPE_INTEGER;
    opt.redevelopmentLag->required = NO;
    opt.redevelopmentLag->options = "1-";
    opt.redevelopmentLag->description =
            _("Number of steps before redevelopment can happen again in a cell developed during simulation");
    opt.redevelopmentLag->guisection = _("Density");

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

    opt.outputDensity = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.outputDensity->key = "output_density";
    opt.outputDensity->required = NO;
    opt.outputDensity->label =
        _("Basename for raster maps of density generated after each step");
    opt.outputDensity->guisection = _("Output");

    opt.potentialFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.potentialFile->key = "devpot_params";
    opt.potentialFile->required = YES;
    opt.potentialFile->label =
        _("CSV file with development potential parameters for each region");
    opt.potentialFile->description =
        _("Each line should contain region ID followed"
          " by parameters (intercepts, development pressure, other predictors).");
    opt.potentialFile->guisection = _("Potential");

    opt.cellDemandFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.cellDemandFile->key = "demand";
    opt.cellDemandFile->required = YES;
    opt.cellDemandFile->description =
            _("CSV file with number of cells to convert for each step and subregion");
    opt.cellDemandFile->guisection = _("Demand");

    opt.populationDemandFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.populationDemandFile->key = "population_demand";
    opt.populationDemandFile->required = NO;
    opt.populationDemandFile->description =
            _("CSV file with population size to accommodate");
    opt.populationDemandFile->guisection = _("Demand");

    opt.separator = G_define_standard_option(G_OPT_F_SEP);
    opt.separator->answer = "comma";
    opt.separator->description =
            _("Separator used in input CSV files");

    opt.patchFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.patchFile->key = "patch_sizes";
    opt.patchFile->required = YES;
    opt.patchFile->description =
        _("File containing list of patch sizes to use");
    opt.patchFile->guisection = _("PGA");

    opt.redistributionMatrix = G_define_standard_option(G_OPT_F_INPUT);
    opt.redistributionMatrix->key = "redistribution_matrix";
    opt.redistributionMatrix->required = NO;
    opt.redistributionMatrix->description = _("Matrix containing probabilities of moving from one subregion to another");
    opt.redistributionMatrix->guisection = _("Climate scenarios");

    opt.redistributionMatrixOutput = G_define_standard_option(G_OPT_F_OUTPUT);
    opt.redistributionMatrixOutput->key = "redistribution_output";
    opt.redistributionMatrixOutput->required = NO;
    opt.redistributionMatrixOutput->description = _("Base name for output file containing matrix of pixels moved from one subregion to another");
    opt.redistributionMatrixOutput->guisection = _("Climate scenarios");

    opt.HAND = G_define_standard_option(G_OPT_R_INPUT);
    opt.HAND->key = "hand";
    opt.HAND->required = NO;
    opt.HAND->description = _("Height Above Nearest Drainage raster");
    opt.HAND->guisection = _("Climate scenarios");

    opt.floodProbability = G_define_standard_option(G_OPT_R_INPUT);
    opt.floodProbability->key = "flood_probability";
    opt.floodProbability->required = NO;
    opt.floodProbability->description = _("Flood probability raster");
    opt.floodProbability->guisection = _("Climate scenarios");

    opt.HUCs = G_define_standard_option(G_OPT_R_INPUT);
    opt.HUCs->key = "huc";
    opt.HUCs->required = NO;
    opt.HUCs->description = _("Raster of HUCs");
    opt.HUCs->guisection = _("Climate scenarios");

    opt.adaptiveCapacity = G_define_standard_option(G_OPT_R_INPUT);
    opt.adaptiveCapacity->key = "adaptive_capacity";
    opt.adaptiveCapacity->required = NO;
    opt.adaptiveCapacity->description = _("Adaptive capacity raster");
    opt.adaptiveCapacity->guisection = _("Climate scenarios");

    opt.outputAdaptation = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.outputAdaptation->key = "output_adaptation";
    opt.outputAdaptation->required = NO;
    opt.outputAdaptation->label =
        _("Basename for raster maps of adaptation generated after each step");
    opt.outputAdaptation->guisection = _("Output");

    opt.depthDamageFunc = G_define_standard_option(G_OPT_F_INPUT);
    opt.depthDamageFunc->key = "depth_damage_functions";
    opt.depthDamageFunc->required = NO;
    opt.depthDamageFunc->description =
        _("CSV file with depth-damage function");
    opt.depthDamageFunc->guisection = _("Climate scenarios");

    opt.ddf_subregions = G_define_standard_option(G_OPT_R_INPUT);
    opt.ddf_subregions->key = "ddf_subregions";
    opt.ddf_subregions->required = NO;
    opt.ddf_subregions->description = _("Subregions raster for depth-damage functions");
    opt.ddf_subregions->guisection = _("Climate scenarios");

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

    opt.potentialWeight = G_define_standard_option(G_OPT_R_INPUT);
    opt.potentialWeight->key = "potential_weight";
    opt.potentialWeight->required = NO;
    opt.potentialWeight->label =
            _("Raster map of weights altering development potential");
    opt.potentialWeight->description =
            _("Values need to be between -1 and 1, where negative locally reduces"
              "probability and positive increases probability.");
    opt.potentialWeight->guisection = _("Scenarios");
    
    opt.incentivePower = G_define_option();
    opt.incentivePower->key = "incentive_power";
    opt.incentivePower->required = NO;
    opt.incentivePower->type = TYPE_DOUBLE;
    opt.incentivePower->answer = "1";
    opt.incentivePower->label =
        _("Exponent to transform probability values p to p^x to simulate infill vs. sprawl");
    opt.incentivePower->description =
        _("Values > 1 encourage infill, < 1 urban sprawl");
    opt.incentivePower->guisection = _("Scenarios");
    opt.incentivePower->options = "0-10";

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

    opt.memory = G_define_option();
    opt.memory->key = "memory";
    opt.memory->type = TYPE_DOUBLE;
    opt.memory->required = NO;
    opt.memory->description = _("Memory in GB");

    // TODO: add mutually exclusive?
    // TODO: add flags or options to control values in series and final rasters

    // provided XOR generated
    G_option_exclusive(opt.seed, flg.generateSeed, NULL);
    G_option_required(opt.seed, flg.generateSeed, NULL);
    G_option_requires_all(opt.density, opt.populationDemandFile, NULL);
    G_option_collective(opt.density, opt.densityCapacity, opt.outputDensity,
                        opt.redevelopmentLag, opt.redevelopmentPotentialFile, NULL);
    G_option_collective(opt.HAND, opt.redistributionMatrix, opt.populationDemandFile,
                        opt.floodProbability, opt.adaptiveCapacity, opt.HUCs,
                        opt.depthDamageFunc, NULL);
    G_option_requires(opt.outputAdaptation, opt.adaptiveCapacity, NULL);
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
    initialize_devpressure_matrix(&devpressure_info);

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
    if (opt.redevelopmentLag->answer)
        patch_info.redevelopment_lag = atoi(opt.redevelopmentLag->answer);
    
    num_steps = 0;
    if (opt.numSteps->answer)
        num_steps = atoi(opt.numSteps->answer);
    
    num_predictors = 0;
    for (i = 0; opt.predictors->answers[i]; i++)
        num_predictors++;
    
    segments.use_weight = false;
    if (opt.potentialWeight->answer) {
        segments.use_weight = true;
    }
    segments.use_potential_subregions = false;
    if (opt.potentialSubregions->answer) {
        segments.use_potential_subregions = true;
    }
    segments.use_density = false;
    if (opt.density->answer) {
        segments.use_density = true;
    }
    segments.use_climate = false;
    if (opt.HAND->answer) {
        segments.use_climate = true;
    }
    memory = -1;
    if (opt.memory->answer)
        memory = atof(opt.memory->answer);
    nseg = manage_memory(&segment_info, &segments, memory);
    segment_info.in_memory = nseg;

    potential_info.incentive_transform_size = 0;
    potential_info.incentive_transform = NULL;
    if (opt.incentivePower->answer) {
        exponent = atof(opt.incentivePower->answer);
        if (exponent !=  1)  /* 1 is no-op */
            initialize_incentive(&potential_info, exponent);
    }
    if (opt.redevelopmentPotentialFile->answer) {
        redev_potential_info.incentive_transform_size = 0;
        redev_potential_info.incentive_transform = NULL;
        if (opt.incentivePower->answer) {
            exponent = atof(opt.incentivePower->answer);
            if (exponent !=  1)  /* 1 is no-op */
                initialize_incentive(&redev_potential_info, exponent);
        }
    }

    raster_inputs.developed = opt.developed->answer;
    raster_inputs.regions = opt.subregions->answer;
    raster_inputs.devpressure = opt.devpressure->answer;
    raster_inputs.predictors = opt.predictors->answers;
    if (opt.potentialWeight->answer)
        raster_inputs.weights = opt.potentialWeight->answer;
    if (opt.potentialSubregions->answer)
        raster_inputs.potential_regions = opt.potentialSubregions->answer;
    if (opt.density->answer) {
        raster_inputs.density = opt.density->answer;
        raster_inputs.density_capacity = opt.densityCapacity->answer;
    }
    if (opt.redistributionMatrix->answer) {
        redistr_matrix.filename = opt.redistributionMatrix->answer;
        read_redistribution_matrix(&redistr_matrix);
    }
    if (opt.HAND->answer) {
        raster_inputs.HAND = opt.HAND->answer;
        raster_inputs.flood_probability = opt.floodProbability->answer;
        raster_inputs.adaptive_capacity = opt.adaptiveCapacity->answer;
        raster_inputs.HUC = opt.HUCs->answer;
        raster_inputs.DDF_regions = NULL;
        DDF.subregions_source = DDF_NONE;
        if (opt.ddf_subregions->answer) {
            /* if subregions map name used for other subregions, set NULL
               here so that no new input will be read */
            if (strcmp(opt.ddf_subregions->answer, opt.subregions->answer) == 0)
                DDF.subregions_source = DDF_DEFAULT;
            else if (opt.potentialSubregions->answer &&
                     strcmp(opt.ddf_subregions->answer, opt.potentialSubregions->answer) == 0)
                DDF.subregions_source = DDF_POTENTIAL;
            else
                DDF.subregions_source = DDF_CUSTOM;
        }
        if (DDF.subregions_source == DDF_CUSTOM)
            raster_inputs.DDF_regions = opt.ddf_subregions->answer;

        DDF.filename = opt.depthDamageFunc->answer;
        DDF.separator = G_option_to_separator(opt.separator);
        max_flood_probability_map = KeyValueIntFloat_create();
        HUC_map = KeyValueIntInt_create();
        initilize_adaptation(&segments.adaptation, &segment_info);
    }
    initialize_flood_response(&response_relation);

    //    read Subregions layer
    region_map = KeyValueIntInt_create();
    reverse_region_map = KeyValueIntInt_create();
    potential_region_map = KeyValueIntInt_create();
    predictor_map = KeyValueCharInt_create();
    DDF_region_map = KeyValueIntInt_create();
    G_verbose_message("Reading input rasters...");
    read_input_rasters(raster_inputs, &segments, segment_info, region_map,
                       reverse_region_map, potential_region_map,
                       HUC_map, max_flood_probability_map,
                       DDF_region_map);
    if (opt.HAND->answer) {
        create_bboxes(&segments.HUC, &segments.developed, &bboxes);
    }
    /* create probability segment*/
    if (Segment_open(&segments.probability, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segment_info.rows, segment_info.cols,
                     Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));

    /* fill predictor map for potential reading */
    for (i = 0; i < num_predictors; i++)
        KeyValueCharInt_set(predictor_map, raster_inputs.predictors[i], i);

    /* read Potential file */
    G_verbose_message("Reading potential file...");
    potential_info.filename = opt.potentialFile->answer;
    potential_info.separator = G_option_to_separator(opt.separator);
    read_potential_file(&potential_info,
                        opt.potentialSubregions->answer ? potential_region_map : region_map,
                        predictor_map);
    if (opt.redevelopmentPotentialFile->answer) {
        redev_potential_info.filename = opt.redevelopmentPotentialFile->answer;
        redev_potential_info.separator = G_option_to_separator(opt.separator);
        read_potential_file(&redev_potential_info,
                            opt.potentialSubregions->answer ? potential_region_map : region_map,
                            predictor_map);
    }
    /* read in predictors and aggregate to save memory */
    G_verbose_message("Reading predictors...");
    read_predictors(raster_inputs, &segments, &potential_info,
                    segment_info);

    /* read Demand file */
    G_verbose_message("Reading demand file...");
    demand_info.has_population = false;
    demand_info.cells_filename = opt.cellDemandFile->answer;
    if (opt.populationDemandFile->answer) {
        demand_info.has_population = true;
        demand_info.population_filename = opt.populationDemandFile->answer;
    }
    demand_info.separator = G_option_to_separator(opt.separator);
    read_demand_file(&demand_info, region_map);
    if (num_steps == 0)
        num_steps = demand_info.max_steps;

    /* check redistribution matrix output files */
    if (opt.redistributionMatrixOutput->answer) {
        redistr_matrix.output_basename = opt.redistributionMatrixOutput->answer;
        if (check_matrix_filenames_exist(&redistr_matrix, num_steps)) {
            if (!G_check_overwrite(argc, argv))
                G_fatal_error(_("At least one of the requested matrix output files exists. Use --o to overwrite."));
        }
    }
    if (opt.depthDamageFunc->answer) {
        /* if no DDF subregions, assume one DDF for entire area */
        if (DDF.subregions_source == DDF_NONE) {
            KeyValueIntInt_set(DDF_region_map, 1, 0);
            read_DDF_file(&DDF, DDF_region_map);
        }
        /* use subregions */
        else if (DDF.subregions_source == DDF_DEFAULT) {
            read_DDF_file(&DDF, region_map);
        }
        /* use potential subregions */
        else if (DDF.subregions_source == DDF_POTENTIAL) {
            read_DDF_file(&DDF, potential_region_map);
        }
        else {
            read_DDF_file(&DDF, DDF_region_map);
        }
    }

    /* read Patch sizes file */
    G_verbose_message("Reading patch size file...");
    patch_sizes.filename = opt.patchFile->answer;
    read_patch_sizes(&patch_sizes, region_map, discount_factor);

    undev_cells = initialize_developables(region_map->nitems);
    patch_overflow = G_calloc(region_map->nitems, sizeof(int));

    /* redevelopment */
    if (segments.use_density) {
        dev_cells = initialize_developables(region_map->nitems);
        population_overflow = G_calloc(region_map->nitems, sizeof(float));
    }
    /* here do the modeling */
    overgrow = true;
    G_verbose_message("Starting simulation...");
    leaving_population = 0;
    for (step = 0; step < num_steps; step++) {
        recompute_probabilities(undev_cells, &segments, &potential_info, false);
        if (segments.use_density)
            recompute_probabilities(dev_cells, &segments, &redev_potential_info, true);
        if (step == num_steps - 1)
            overgrow = false;
        for (region = 0; region < region_map->nitems; region++) {
            compute_step(undev_cells, dev_cells, &demand_info, search_alg, &segments,
                         &patch_sizes, &patch_info, &devpressure_info, patch_overflow,
                         population_overflow, &redistr_matrix, step, region, reverse_region_map, overgrow);
        }
        /* simulate abandonment due to climate (flooding) */
        if (segments.use_climate) {
            for (HUC = 0; HUC < HUC_map->nitems; HUC++)
                climate_step(&segments, &demand_info, &bboxes,
                             &redistr_matrix, region_map, reverse_region_map,
                             step, &leaving_population,
                             max_flood_probability_map, &DDF,
                             &response_relation, HUC);
            if (opt.redistributionMatrixOutput->answer)
                write_redistribution_matrix(&redistr_matrix, step, num_steps);
        }
        /* export developed for that step */
        if (opt.outputSeries->answer) {
            name_step = name_for_step(opt.outputSeries->answer, step, num_steps);
            output_developed_step(&segments.developed, name_step,
                                  demand_info.years[step], -1, num_steps, false, false);
        }
        /* export density for that step */
        if (opt.outputDensity->answer) {
            name_step = name_for_step(opt.outputDensity->answer, step, num_steps);
            output_step(&segments.density, &segments.developed, name_step, FCELL_TYPE);
        }
        /* export density for that step */
        if (opt.outputAdaptation->answer) {
            name_step = name_for_step(opt.outputAdaptation->answer, step, num_steps);
            output_step(&segments.adaptation, &segments.developed, name_step, CELL_TYPE);
        }
    }

    /* write */
    output_developed_step(&segments.developed, opt.output->answer,
                          demand_info.years[0], demand_info.years[step-1],
                          num_steps, false, false);

    /* close segments and free memory */
    Segment_close(&segments.developed);
    Segment_close(&segments.subregions);
    Segment_close(&segments.devpressure);
    Segment_close(&segments.probability);
    Segment_close(&segments.aggregated_predictor);
    if (opt.potentialWeight->answer) {
        Segment_close(&segments.weight);
    }
    if (opt.potentialSubregions->answer)
        Segment_close(&segments.potential_subregions);
    if (segments.use_density) {
        Segment_close(&segments.density);
        Segment_close(&segments.density_capacity);
    }
    if (segments.use_climate) {
        Segment_close(&segments.HAND);
        Segment_close(&segments.flood_probability);
        Segment_close(&segments.adaptive_capacity);
        Segment_close(&segments.HUC);
        Segment_close(&segments.adaptation);
        KeyValueIntFloat_free(max_flood_probability_map);
        KeyValueIntInt_free(HUC_map);
    }
    KeyValueIntInt_free(region_map);
    KeyValueIntInt_free(reverse_region_map);
    KeyValueCharInt_free(predictor_map);
    KeyValueIntInt_free(potential_region_map);
    KeyValueIntInt_free(DDF_region_map);
    if (demand_info.cells_table) {
        for (int i = 0; i < demand_info.max_subregions; i++)
            G_free(demand_info.cells_table[i]);
        G_free(demand_info.cells_table);
        G_free(demand_info.years);
    }
    if (demand_info.has_population) {
        for (int i = 0; i < demand_info.max_subregions; i++)
            G_free(demand_info.population_table[i]);
        G_free(demand_info.population_table);
    }
    if (potential_info.predictors) {
        for (int i = 0; i < potential_info.max_predictors; i++)
            G_free(potential_info.predictors[i]);
        G_free(potential_info.predictors);
        G_free(potential_info.devpressure);
        G_free(potential_info.intercept);
        G_free(potential_info.predictor_indices);
    }
    if (opt.redevelopmentPotentialFile->answer) {
        for (int i = 0; i < redev_potential_info.max_predictors; i++)
            G_free(redev_potential_info.predictors[i]);
        G_free(redev_potential_info.predictors);
        G_free(redev_potential_info.devpressure);
        G_free(redev_potential_info.intercept);
        G_free(redev_potential_info.predictor_indices);
    }
    for (int i = 0; i < devpressure_info.neighborhood * 2 + 1; i++)
        G_free(devpressure_info.matrix[i]);
    G_free(devpressure_info.matrix);
    if (potential_info.incentive_transform_size > 0)
        G_free(potential_info.incentive_transform);
    if (undev_cells) {
        G_free(undev_cells->num);
        G_free(undev_cells->max);
        for (int i = 0; i < undev_cells->max_subregions; i++)
            G_free(undev_cells->cells[i]);
        G_free(undev_cells->cells);
        G_free(undev_cells);
    }
    if (segments.use_density) {
        G_free(dev_cells->num);
        G_free(dev_cells->max);
        for (int i = 0; i < dev_cells->max_subregions; i++)
            G_free(dev_cells->cells[i]);
        G_free(dev_cells->cells);
        G_free(dev_cells);
        if (population_overflow)
            G_free(population_overflow);
    }
    if (opt.redistributionMatrix->answer)
        free_redistribution_matrix(&redistr_matrix);
    G_free(patch_sizes.patch_sizes);
    G_free(patch_overflow);

    return EXIT_SUCCESS;
}

