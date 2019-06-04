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


struct input
{
    const char *name;
    int fd;
    FCELL *buf;
};

struct Demand
{
    const char *filename;
    int **table;
    int max_subregions;
    int max_steps;
};

struct Potential
{
    const char *filename;
    double **predictors;
    double *intercept;
    double *devpressure;
    int max_predictors;
    int max_subregions;
};

struct SegmentMemory
{
    int rows;
    int cols;
    int in_memory;
};

struct Undeveloped
{
    int max_subregions;
    size_t *max_undeveloped;
    size_t *num_undeveloped;
    size_t **cell;
};


size_t get_idx_from_xy(int row, int col, int cols)
{
    return cols * row + col;
}

void rast_segment_open(const char *name, SEGMENT *segment, struct SegmentMemory segmentInfo, 
                       RASTER_MAP_TYPE map_type)
{
    int row;
    int rowio = Rast_open_old(name, "");
    void *raster_row = Rast_allocate_buf(map_type);

    if (Segment_open(segment, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segmentInfo.rows, segmentInfo.cols,
                     Rast_cell_size(map_type), segmentInfo.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));

    for (row = 0; row < Rast_window_rows(); row++) {
        Rast_get_row(rowio, raster_row, row, map_type);
        Segment_put_row(segment, raster_row, row);
    }
    Segment_flush(segment);
    Rast_close(rowio);          /* we won't need the raster again */
    G_free(raster_row);
}

void read_predictors(char **predictor_names, SEGMENT *segment, SEGMENT *output,
                     struct SegmentMemory segmentInfo, int ninputs)
{
    int input;
    int row;
    int col;
    int nrows = Rast_window_rows();
    int ncols = Rast_window_cols();
    size_t segment_cell_size = sizeof(FCELL) * ninputs;

    int *input_fds = G_malloc(ninputs * sizeof(int));
    /* open existing raster maps for reading */
    for (input = 0; input < ninputs; input++) {
        input_fds[input] = Rast_open_old(predictor_names[input], "");
    }
    
    if (Segment_open(segment, G_tempfile(),
                     nrows, ncols, segmentInfo.rows, segmentInfo.cols,
                     segment_cell_size, segmentInfo.in_memory) != 1)
        G_fatal_error(_("Unable to create temporary segment file"));
    
    /* allocate input buffer */
    FCELL *row_buffer = Rast_allocate_f_buf();
    FCELL *seg_buffer = G_malloc(ncols * ninputs * sizeof(FCELL));
    CELL out_mask;
    
    for (row = 0; row < nrows; row++) {
        for (input = 0; input < ninputs; input++) {
            Rast_get_f_row(input_fds[input], row_buffer, row);
            for (col = 0; col < ncols; col++) {
                    seg_buffer[col * ninputs + input] = row_buffer[col];
                    /* collect all nulls in predictors and set it in output raster */
                    if (Rast_is_null_value(&((FCELL *) row_buffer)[col], FCELL_TYPE))
                    {
                        Rast_set_c_null_value(&out_mask, 1);
                        Segment_put(output, (void *)&out_mask, row, col);
                    }
                }
        }
        if (Segment_put_row(segment, seg_buffer, row) < 1)
            G_fatal_error(_("Unable to write temporary segment file"));
    }
    Segment_flush(segment);
    Segment_flush(output);
    for (input = 0; input < ninputs; input++)
        Rast_close(input_fds[input]);
    G_free(input_fds);
    G_free(row_buffer);
    G_free(seg_buffer);
}

double get_develop_probability(SEGMENT *predictors, SEGMENT *devpressure,
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


void read_subregions(const char *subregions, SEGMENT * segment,
                    struct KeyValueIntInt *region_map)
{
    int val;

    int fd = Rast_open_old(subregions, "");
    void *buffer = Rast_allocate_c_buf();

    G_verbose_message("Reading subregions %s", subregions);
    int count_regions = 0;
    int index = 0;
    int row, col;
    int segment_rows = 64;
    int segment_cols = 64;
    int segments_in_memory = 4;
    if (Segment_open(segment, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segment_rows, segment_cols,
                     Rast_cell_size(CELL_TYPE), segments_in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));
    for (row = 0; row < Rast_window_rows(); row++) {
        Rast_get_row(fd, buffer, row, CELL_TYPE);
        void *ptr = buffer;

        for (col = 0; col < Rast_window_cols(); col++,
             ptr = G_incr_void_ptr(ptr, Rast_cell_size(CELL_TYPE))) {
            if (Rast_is_null_value(ptr, CELL_TYPE))
                ;
            //Rast_set_c_null_value(ptr, 1);
            else {
                val = *(CELL *) ptr;
                if (KeyValueIntInt_find(region_map, val, &index))
                    ; // pass
                else {
                    KeyValueIntInt_set(region_map, val, count_regions);
                    index = count_regions;
                    count_regions++;
                }
                *(CELL *) ptr = index;
            }
        }
        Segment_put_row(segment, buffer, row);
    }
    Segment_flush(segment);
    G_free(buffer);
    Rast_close(fd);
}


void read_demand_file(struct Demand *demandInfo, struct KeyValueIntInt *region_map)
{
    FILE *fp;
    if ((fp = fopen(demandInfo->filename, "r")) == NULL)
        G_fatal_error(_("Cannot open population demand file <%s>"),
                      demandInfo->filename);
    int countlines = 0;
    // Extract characters from file and store in character c 
    for (char c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n') // Increment count if this character is newline 
            countlines++; 

    rewind(fp);

    size_t buflen = 4000;
    char buf[buflen];
    if (G_getl2(buf, buflen, fp) == 0)
        G_fatal_error(_("Population demand file <%s>"
                        " contains less than one line"), demandInfo->filename);

    char **tokens;
    int ntokens;

    const char *fs = "\t";
    const char *td = "\"";

    tokens = G_tokenize2(buf, fs, td);
    ntokens = G_number_of_tokens(tokens);
    if (ntokens == 0)
        G_fatal_error("No columns in the header row");

    struct ilist *ids = G_new_ilist();
    int count;
    // skip first column which does not contain id of the region
    int i;
    for (i = 1; i < ntokens; i++) {
        G_chop(tokens[i]);
        G_ilist_add(ids, atoi(tokens[i]));
    }

    int years = 0;
    demandInfo->table = (int **) G_malloc(region_map->nitems * sizeof(int *)); 
    for (int i = 0; i < region_map->nitems; i++) {
        demandInfo->table[i] = (int *) G_malloc(countlines * sizeof(int)); 
    }
    while(G_getl2(buf, buflen, fp)) {
        if (!buf || buf[0] == '\0')
            continue;
        tokens = G_tokenize2(buf, fs, td);
        int ntokens2 = G_number_of_tokens(tokens);
        if (ntokens2 == 0)
            continue;
        if (ntokens2 != ntokens)
            G_fatal_error(_("Demand: wrong number of columns in line: %s"), buf);

        count = 0;
        int i;
        for (i = 1; i < ntokens; i++) {
            // skip first column which is the year which we ignore
            int idx;
            if (KeyValueIntInt_find(region_map, ids->value[count], &idx)) {
                G_chop(tokens[i]);
                demandInfo->table[idx][years] = atoi(tokens[i]);
            }
            count++;
        }
        // each line is a year
        years++;
    }
    demandInfo->max_subregions = region_map->nitems;
    demandInfo->max_steps = years;
    G_verbose_message("Number of steps in demand file: %d", years);
    //    if (!sParams.nSteps)
    //        sParams.nSteps = years;
    G_free_ilist(ids);
    G_free_tokens(tokens);
}


void read_potential_file(struct Potential *potentialInfo, struct KeyValueIntInt *region_map,
                         int num_predictors)
{
    FILE *fp;
    if ((fp = fopen(potentialInfo->filename, "r")) == NULL)
        G_fatal_error(_("Cannot open development potential parameters file <%s>"),
                      potentialInfo->filename);

    const char *fs = "\t";
    const char *td = "\"";

    size_t buflen = 4000;
    char buf[buflen];
    if (G_getl2(buf, buflen, fp) == 0)
        G_fatal_error(_("Development potential parameters file <%s>"
                        " contains less than one line"), potentialInfo->filename);
    potentialInfo->max_predictors = num_predictors;
    potentialInfo->intercept = (double *) G_malloc(region_map->nitems * sizeof(double)); 
    potentialInfo->devpressure = (double *) G_malloc(region_map->nitems * sizeof(double)); 
    potentialInfo->predictors = (double **) G_malloc(num_predictors * sizeof(double *)); 
    for (int i = 0; i < num_predictors; i++) {
        potentialInfo->predictors[i] = (double *) G_malloc(region_map->nitems * sizeof(double)); 
    }

    char **tokens;

    while (G_getl2(buf, buflen, fp)) {
        if (!buf || buf[0] == '\0')
            continue;
        tokens = G_tokenize2(buf, fs, td);
        int ntokens = G_number_of_tokens(tokens);
        if (ntokens == 0)
            continue;
        // id + intercept + devpressure + predictores
        if (ntokens != num_predictors + 3)
            G_fatal_error(_("Potential: wrong number of columns: %s"), buf);

        int idx;
        int id;
        double coef_intercept, coef_devpressure;
        double val;
        int j;

        G_chop(tokens[0]);
        id = atoi(tokens[0]);
        if (KeyValueIntInt_find(region_map, id, &idx)) {
            G_chop(tokens[1]);
            coef_intercept = atof(tokens[1]);
            G_chop(tokens[2]);
            coef_devpressure = atof(tokens[2]);
            potentialInfo->intercept[idx] = coef_intercept;
            potentialInfo->devpressure[idx] = coef_devpressure;
            for (j = 0; j < num_predictors; j++) {
                G_chop(tokens[j + 3]);
                val = atof(tokens[j + 3]);
                potentialInfo->predictors[j][idx] = val;
            }
        }
        // else ignoring the line with region which is not used

        G_free_tokens(tokens);
    }

    fclose(fp);
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
                *demandFile, *output, *seed;

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
    SEGMENT predictors;
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
    
    //   read developed
    G_verbose_message("Reading input rasters...");
    int rowio = Rast_open_old(opt.developed->answer, "");
    if (Segment_open(&developed_segment, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segment_info.rows, segment_info.cols,
                     Rast_cell_size(CELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));
    void *raster_row = Rast_allocate_buf(CELL_TYPE);
    int cols = Rast_window_cols();
    for (int row = 0; row < Rast_window_rows(); row++) {
        Rast_get_row(rowio, raster_row, row, CELL_TYPE);
        for (int col = 0; col < Rast_window_cols(); col++) {
            if (Rast_is_null_value(&((CELL *) raster_row)[col], CELL_TYPE))
                ;
            else {
                CELL c = ((CELL *) raster_row)[col];
                ((CELL *) raster_row)[col] = c - 1;
                /* save undeveloped cells by region */
                if(c == 0){
                    size_t idx = get_idx_from_xy(row, col, cols);
                    CELL region_index;
                    Segment_get(&subregions_segment, (void *)&region_index, row, col);
                    size_t current_num_undev = undev_cells->num_undeveloped[region_index];
                    if (current_num_undev > undev_cells->max_undeveloped[region_index])
                    {
                        size_t new_size = 2 * undev_cells->max_undeveloped[region_index];
                        undev_cells->cell[region_index] = (size_t *) G_realloc(undev_cells->cell[region_index], new_size * sizeof(size_t));
                        undev_cells->max_undeveloped[region_index] = new_size;
                    }
                    undev_cells->cell[region_index][current_num_undev] = idx;
                    undev_cells->num_undeveloped[region_index]++;
                }
            }
        }
        Segment_put_row(&developed_segment, raster_row, row);
    }
    Segment_flush(&developed_segment);
    Rast_close(rowio);
    G_free(raster_row);

    /* read in predictors */
    read_predictors(opt.predictors->answers, &predictors, &developed_segment, segment_info, num_predictors);



    

    /* here do the modeling */
    
    
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
            prob = get_develop_probability(&predictors, &devpressure_segment,
                                           values,
                                           &potential_info, reg, row, col);
            ((FCELL *) row_buffer)[col] = prob;
        }
        Rast_put_row(out_fd, row_buffer, FCELL_TYPE);
    }
    Rast_close(out_fd);
    G_free(row_buffer);
    G_free(values);
    Segment_close(&predictors);

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

    return EXIT_SUCCESS;
}

