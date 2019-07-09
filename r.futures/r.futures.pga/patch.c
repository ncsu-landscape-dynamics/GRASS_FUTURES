/*!
   \file patch.c

   \brief Functions to grow patches

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
#include "patch.h"


int get_patch_size(struct PatchSizes *patch_info)
{
    return patch_info->patch_sizes[(int)(G_drand48() * patch_info->max_patches)];
}

void add_neighbour(int row, int col, int seed_row, int seed_col,
                   struct CandidateNeighborsList *candidate_list,
                   SEGMENT *developed, SEGMENT *probability,
                   double alpha)
{
    int i, shouldAdd;
    double distance;
    size_t idx;
    CELL value;
    FCELL prob;
    
    Segment_get(developed, (void *)&value, row, col);
    if (Rast_is_null_value(&value, CELL_TYPE))
        return;
    if (value == -1) {
        idx = get_idx_from_xy(row, col, Rast_window_cols());
        /* need to add this cell... */
        
        /* ...either refresh its element in list if already there */
        shouldAdd = 1;
        for (i = 0; i < candidate_list->n; i++) {
            if (candidate_list->candidates[i].id == idx) {
                shouldAdd = 0;
                break;
            }
        }
        /* or add it on the end, allocating space if necessary */
        if (shouldAdd) {
            if (candidate_list->n == candidate_list->max_n) {
                candidate_list->max_n += candidate_list->block_size;
                candidate_list->candidates =  (struct CandidateNeighbor *)
                        G_realloc(candidate_list->candidates, candidate_list->max_n * sizeof(struct CandidateNeighbor));
                if (!candidate_list->candidates) {
                    G_fatal_error("Memory error in add_neighbour_if_possible()");
                }
            }
            candidate_list->candidates[candidate_list->n].id = idx;
            Segment_get(probability, (void *)&prob, row, col);
            candidate_list->candidates[candidate_list->n].potential = prob;
            distance = get_distance(seed_row, seed_col, row, col);
            candidate_list->candidates[candidate_list->n].suitability = prob / pow(distance, alpha);
            candidate_list->n++;
        }
    }
}

void add_neighbours(int row, int col, int seed_row, int seed_col,
                    struct CandidateNeighborsList *candidate_list,
                    SEGMENT *developed, SEGMENT *probability,
                    double alpha, int num_neighbors)
{
    add_neighbour(row - 1, col, seed_row, seed_col, candidate_list,
                  developed, probability, alpha);  // left
    add_neighbour(row + 1, col, seed_row, seed_col, candidate_list, 
                  developed, probability, alpha);  // right
    add_neighbour(row, col - 1, seed_row, seed_col, candidate_list, 
                  developed, probability, alpha);  // down
    add_neighbour(row, col + 1, seed_row, seed_col, candidate_list, 
                  developed, probability, alpha);  // up 
    if (num_neighbors == 8) {
        add_neighbour(row - 1, col - 1, seed_row, seed_col, candidate_list, 
                      developed, probability, alpha);
        add_neighbour(row - 1, col + 1, seed_row, seed_col, candidate_list, 
                      developed, probability, alpha);
        add_neighbour(row + 1, col - 1, seed_row, seed_col, candidate_list, 
                      developed, probability, alpha);
        add_neighbour(row + 1, col + 1, seed_row, seed_col, candidate_list, 
                      developed, probability, alpha);
    }
}

static int sort_neighbours(const void *p1, const void *p2)
{
    struct CandidateNeighbor *p1_ = (struct CandidateNeighbor *) p1;
    struct CandidateNeighbor *p2_ = (struct CandidateNeighbor *) p2;
    if (p1_->suitability > p2_->suitability) {
        return -1;
    }
    if (p2_->suitability > p1_->suitability) {
        return 1;
    }
    return 0;
}


double get_distance(int row1, int col1, int row2, int col2)
{
    return sqrt((row1 - row2) * (row1 - row2) + (col1 - col2) * (col1 - col2));
}



int grow_patch(int seed_row, int seed_col, int *added_ids,
               SEGMENT *developed, SEGMENT *probability,
               int num_neighbors, double alpha, int patch_size,
               int step, enum slow_grow strategy)
{
    int i, j, iter;
    double r, p;
    int found;
    int force, skip;
    int row, col, cols;

    struct CandidateNeighborsList candidates;
    candidates.block_size = 20;
    candidates.candidates = (struct CandidateNeighbor *) G_malloc(sizeof(struct CandidateNeighbor) * candidates.block_size);
    candidates.max_n = candidates.block_size;
    candidates.n = 0;
    
    cols = Rast_window_cols();
    found = 0;
    force = 0;
    skip = 0;

    add_neighbours(seed_row, seed_col, seed_row, seed_col,
                   &candidates, developed, probability, alpha, num_neighbors);
    iter = 0;
    while (candidates.n > 0 && found < patch_size && skip == 0) {
        i = 0;
        while (1) {
            /* challenge the candidate */
            r = G_drand48();
            p = candidates.candidates[i].potential;
            if (r < p || force) {
                /* update list of added IDs */
                added_ids[found] = candidates.candidates[i].id;
                /* update to developed */
                get_xy_from_idx(candidates.candidates[i].id, cols, &row, &col);
                Segment_put(developed, (void *)&step, row, col);
                /* remove this one from the list by copying down everything above it */
                for (j = i + 1; j < candidates.n; j++) {
                    candidates.candidates[j - 1].id = candidates.candidates[j].id;
                    candidates.candidates[j - 1].potential = candidates.candidates[j].potential;
                    candidates.candidates[j - 1].suitability = candidates.candidates[j].suitability;
                }
                /* reduce the size of the list */
                candidates.n--;
                /* find and add new candidates */
                add_neighbours(row, col, seed_row, seed_col,
                               &candidates, developed, probability, alpha, num_neighbors);
                /* sort candidates based on probability */
                qsort(candidates.candidates, candidates.n, sizeof(struct CandidateNeighbor), sort_neighbours);
                found++;
                /* restart max iterations when cell found */
                iter = 0;
                force = 0;
                break;
            }
            else {
                i++;
                if (i == candidates.n) {
                    i = 0;
                    iter++;
                    if (iter > MAX_CANDIDATE_ITER) {
                        if (strategy == FORCE_GROW) {
                            force = 1;
                        }
                        else if (strategy == SKIP) {
                            skip = 1;
                            break;
                        }
                        
                    }
                }
            }
        }
    }

    if (candidates.max_n > 0)
        G_free(candidates.candidates);

    return found;
}

