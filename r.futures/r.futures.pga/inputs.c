/*!
   \file inputs.c

   \brief Functions to read in input files and rasters

   (C) 2016-2019 by Vaclav Petras and the GRASS Development Team

   This program is free software under the GNU General Public License
   (>=v2).  Read the file COPYING that comes with GRASS for details.

   \author Vaclav Petras
 */

#include <stdlib.h>
#include <math.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <grass/segment.h>

#include "keyvalue.h"
#include "inputs.h"


void initialize_incentive(struct Potential *potential_info, float exponent)
{
    int i;

    potential_info->incentive_transform_size = 1001;
    potential_info->incentive_transform = (float *) G_malloc(sizeof(float) *
                                                             potential_info->incentive_transform_size);
    i = 0;
    double step = 1. / (potential_info->incentive_transform_size - 1);
    while (i < potential_info->incentive_transform_size) {
        potential_info->incentive_transform[i] = pow(i * step, exponent);
        i++;
    }
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


// could be writing initial probability to speed up
void read_developed(char *filename, struct Segments *segments,
                    struct SegmentMemory segment_info)
{
    int row, col, rowio;
    void *raster_row;

    rowio = Rast_open_old(filename, "");
    if (Segment_open(&segments->developed, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segment_info.rows, segment_info.cols,
                     Rast_cell_size(CELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));
    raster_row = Rast_allocate_buf(CELL_TYPE);
    for (row = 0; row < Rast_window_rows(); row++) {
        Rast_get_row(rowio, raster_row, row, CELL_TYPE);
        for (col = 0; col < Rast_window_cols(); col++) {
            if (Rast_is_null_value(&((CELL *) raster_row)[col], CELL_TYPE))
                ;
            else {
                CELL c = ((CELL *) raster_row)[col];
                ((CELL *) raster_row)[col] = c - 1;
            }
        }
        Segment_put_row(&segments->developed, raster_row, row);
    }
    Segment_flush(&segments->developed);
    Rast_close(rowio);
    G_free(raster_row);
}





void read_predictors(char **predictor_names, struct Segments *segments,
                     struct SegmentMemory segment_info, int ninputs)
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
    
    if (Segment_open(&segments->predictors, G_tempfile(),
                     nrows, ncols, segment_info.rows, segment_info.cols,
                     segment_cell_size, segment_info.in_memory) != 1)
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
                        Segment_put(&segments->developed, (void *)&out_mask, row, col);
                    }
                }
        }
        if (Segment_put_row(&segments->predictors, seg_buffer, row) < 1)
            G_fatal_error(_("Unable to write temporary segment file"));
    }
    Segment_flush(&segments->predictors);
    Segment_flush(&segments->developed);
    for (input = 0; input < ninputs; input++)
        Rast_close(input_fds[input]);
    G_free(input_fds);
    G_free(row_buffer);
    G_free(seg_buffer);
}
void read_subregions(const char *subregions, struct Segments *segments,
                     struct SegmentMemory segment_info, struct KeyValueIntInt *region_map)
{
    int val;

    int fd = Rast_open_old(subregions, "");
    void *buffer = Rast_allocate_c_buf();

    G_verbose_message("Reading subregions %s", subregions);
    int count_regions = 0;
    int index = 0;
    int row, col;
    if (Segment_open(&segments->subregions, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segment_info.rows, segment_info.cols,
                     Rast_cell_size(CELL_TYPE), segment_info.in_memory) != 1)
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
        Segment_put_row(&segments->subregions, buffer, row);
    }
    Segment_flush(&segments->subregions);
    G_free(buffer);
    Rast_close(fd);
}

void read_weights(const char *weights, struct Segments *segments,
                  struct SegmentMemory segment_info)
{
    float val;
    int fd;
    void *buffer;
    void *ptr;
    int row, col;

    fd = Rast_open_old(weights, "");
    buffer = Rast_allocate_c_buf();

    G_verbose_message("Reading weights %s", weights);

    if (Segment_open(&segments->weight, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segment_info.rows, segment_info.cols,
                     Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));
    for (row = 0; row < Rast_window_rows(); row++) {
        Rast_get_row(fd, buffer, row, FCELL_TYPE);
        ptr = buffer;
        for (col = 0; col < Rast_window_cols(); col++,
             ptr = G_incr_void_ptr(ptr, Rast_cell_size(FCELL_TYPE))) {
            if (Rast_is_null_value(ptr, FCELL_TYPE))
                *(FCELL *) ptr = 0;
            else {
                val = *(FCELL *) ptr;
                if (val > 1)
                    val = 1;
                else if (val < -1)
                    val = -1;
                *(FCELL *) ptr = val;
            }
        }
        Segment_put_row(&segments->weight, buffer, row);
    }
    Segment_flush(&segments->weight);
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


void read_patch_sizes(struct PatchSizes *patch_info, double discount_factor)
{
    FILE *fin;
    char *size_buffer;
    int n_max_patches;
    int patch;


    patch_info->max_patches = 0;
    patch_info->max_patch_size = 0;

    G_verbose_message("Reading patch sizes...");
    fin = fopen(patch_info->filename, "rb");
    if (fin) {
        size_buffer = (char *) G_malloc(100 * sizeof(char));
        if (size_buffer) {
            /* just scan the file twice */
            n_max_patches = 0;
            while (fgets(size_buffer, 100, fin)) {
                n_max_patches++;
            }
            rewind(fin);
            if (n_max_patches) {
                patch_info->patch_sizes =
                    (int *) G_malloc(sizeof(int) * n_max_patches);
                if (patch_info->patch_sizes) {
                    while (fgets(size_buffer, 100, fin)) {
                        patch = atoi(size_buffer) * discount_factor;
                        if (patch > 0) {
                            if (patch_info->max_patch_size < patch)
                                patch_info->max_patch_size = patch;
                            patch_info->patch_sizes[patch_info->max_patches] = patch;
                            patch_info->max_patches++;
                        }
                    }
                }
            }
            free(size_buffer);
        }
        fclose(fin);
    }
}

