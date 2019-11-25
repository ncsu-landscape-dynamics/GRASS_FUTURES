/*!
   \file redistribute.c

   \brief Functions to redistribute population

   (C) 2016-2019 by Anna Petrasova and the GRASS Development Team

   This program is free software under the GNU General Public License
   (>=v2).  Read the file COPYING that comes with GRASS for details.

   \author Anna Petrasova
 */

#include <stdlib.h>
#include <stdio.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>

#include "keyvalue.h"
#include "inputs.h"
#include "redistribute.h"


/*!
 * \brief Select item based on its probability
 * 
 * https://stackoverflow.com/questions/9330394/how-to-pick-an-item-by-its-probability
 * 
 * \param probabilities array
 * \param size of array
 * \return index of array
 */
int pick_region(const float *probabilities, int size)
{
    volatile double p;
    volatile int i;
    volatile float cumulative;
    
    p = G_drand48();
    cumulative = 0;
    for (i = 0; i < size; i++) {
        cumulative += probabilities[i];
        if (p <= cumulative)
            return i;
    }
    return i;
}

/*!
 * \brief Redistribute to other county
 * \param matrix
 * \param regionID ID of region (like FIPS)
 * \return id of region where to move to, -1 if regionID not in the matrix as 'from' region
 */
int redistribute(const struct RedistributionMatrix *matrix, int regionID)
{
    int from_idx, to_idx, to_ID;
    if (!KeyValueIntInt_find(matrix->from_map, regionID, &from_idx)) {
        G_warning("Region %d is not in redistribution matrix rows", regionID);
        return -1;
    }
    to_idx = pick_region(matrix->matrix[from_idx], matrix->dim_to);
    KeyValueIntInt_find(matrix->to_map, to_idx, &to_ID);
    return to_ID;
    
}


/*!
 * \brief Reads redistribution matrix
 * 
 * Matrix is in CSV file such as:
 *      ,27642,27645, ...
 * 27642,0.22,1.01, ...   sum of line except of region code should be < 1 
 * 27645,9.12,0.28,...
 * 
 * This means probability that someone from 27642 moves to 27645 is 1.01%
 * 
 * \param matrix structure
 */
void read_redistribution_matrix(struct RedistributionMatrix *matrix)
{
    FILE *fin;
    size_t buflen = 20000;
    char *buf;
    const char *fs = ",";
    const char *td = "\"";
    int ntokens, ntokens2;
    char **tokens, **tokens2;
    int row, col;
    int new_size;

    buf = (char *) G_malloc(buflen * sizeof(char));
    
    if ((fin = fopen(matrix->filename, "r")) == NULL)
        G_fatal_error(_("Cannot open redistribution matrix file <%s>"),
                      matrix->filename);
    matrix->from_map = KeyValueIntInt_create();
    matrix->to_map = KeyValueIntInt_create();

    /* read first line */
    if (G_getl2(buf, buflen, fin) == 0)
        G_fatal_error(_("Development potential parameters file <%s>"
                        " contains less than one line"), matrix->filename);
    else {
        tokens = G_tokenize2(buf, fs, td);
        ntokens = G_number_of_tokens(tokens);
        for (col = 1; col < ntokens; col++) {
            KeyValueIntInt_set(matrix->to_map, col - 1, atoi(tokens[col]));
        }
    }

    matrix->dim_to = ntokens - 1;
    matrix->dim_from = matrix->dim_to;
    matrix->matrix = (float **) G_malloc(matrix->dim_from * sizeof(float *));
    matrix->max_dim_from = matrix->dim_from;

    row = 0;
    while (G_getl2(buf, buflen, fin)) {
        if (!buf || buf[0] == '\0')
            continue;
        tokens2 = G_tokenize2(buf, fs, td);
        ntokens2 = G_number_of_tokens(tokens2);
        if (ntokens != ntokens2)
            G_fatal_error(_("Number of fields in row in file <%s> is not consistent"), matrix->filename);
        /* allocate row */
        matrix->matrix[row] = (float *) G_malloc(matrix->dim_to * sizeof(float));
        for (col = 0; col < ntokens2; col++) {
            if (col == 0) {
                KeyValueIntInt_set(matrix->from_map, atoi(tokens2[col]), row);
//                KeyValueIntInt_set(matrix->from, row, atoi(tokens2[col]));
            }
            else {
                /* convert from % */
                matrix->matrix[row][col - 1] = atof(tokens2[col]) / 100.;
            }
        }
        row++;
        if (matrix->max_dim_from <= row) {
            new_size = 2 * matrix->max_dim_from;
            matrix->matrix = (float **) G_realloc(matrix->matrix, new_size * sizeof(float *));
            matrix->max_dim_from = new_size;
        }

        G_free_tokens(tokens2);
    }
    matrix->dim_from = row;
    G_free_tokens(tokens);
    G_free(buf);
}

void free_redistribution_matrix(struct RedistributionMatrix *matrix)
{
    int i;
    for (i = 0; i < matrix->dim_from; i++)
        G_free(matrix->matrix[i]);
    G_free(matrix->matrix);
    KeyValueIntInt_free(matrix->from_map);
    KeyValueIntInt_free(matrix->to_map);
}



void read_disturbance(const char *name, struct Segments *segments,
                      struct SegmentMemory segment_info)
{
    int row, col;
    int rows, cols;
    int fd_disturbance;
    FCELL *disturbance_row;
    FCELL *disturbance_row_in;
    bool first;

    rows = Rast_window_rows();
    cols = Rast_window_cols();

    /* open existing raster map for reading */
    fd_disturbance = Rast_open_old(name, "");

    /* Segment open */
    first = false;
    if (!segments->use_disturbance) {
        if (Segment_open(&segments->disturbance_effect, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of disturbance effect"));
        segments->use_disturbance = true;
        first = true;
    }
    disturbance_row = Rast_allocate_buf(FCELL_TYPE);
    disturbance_row_in = Rast_allocate_buf(FCELL_TYPE);
    for (row = 0; row < rows; row++) {
        /* read disturbance row */
        Rast_get_row(fd_disturbance, disturbance_row, row, FCELL_TYPE);
        if (!first)
            Segment_get_row(&segments->disturbance_effect, disturbance_row_in, row);
        for (col = 0; col < cols; col++) {
            if (Rast_is_null_value(&((FCELL *) disturbance_row)[col], FCELL_TYPE)) {
                ((FCELL *) disturbance_row)[col] = 1;
            }
            if (!first)
                ((FCELL *) disturbance_row)[col] = ((FCELL *) disturbance_row)[col] * ((FCELL *) disturbance_row_in)[col];
        }
        Segment_put_row(&segments->disturbance_effect, disturbance_row, row);
    }

    /* flush all segments */
    Segment_flush(&segments->disturbance_effect);

    /* close raster maps */
    Rast_close(fd_disturbance);

    G_free(disturbance_row);
    G_free(disturbance_row_in);
}
