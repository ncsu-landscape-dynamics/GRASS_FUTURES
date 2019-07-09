#ifndef FUTURES_PATCH_H
#define FUTURES_PATCH_H

#include <grass/segment.h>

#include "inputs.h"


#define MAX_CANDIDATE_ITER 100

enum slow_grow { FORCE_GROW, SKIP };

struct CandidateNeighbor{
    double potential;      /* s'_i */
    double suitability;    /* s_i */
    size_t id;
    
};

struct CandidateNeighborsList
{
    int n;
    int max_n;
    int block_size;
    struct CandidateNeighbor *candidates;
    
};

int get_patch_size(struct PatchSizes *patch_info);
void add_neighbour(int row, int col, int seed_row, int seed_col,
                   struct CandidateNeighborsList *candidate_list,
                   SEGMENT *developed, SEGMENT *probability,
                   double alpha);
void add_neighbours(int row, int col, int seed_row, int seed_col,
                    struct CandidateNeighborsList *candidate_list,
                    SEGMENT *developed, SEGMENT *probability,
                    double alpha, int num_neighbors);
double get_distance(int row1, int col1, int row2, int col2);
int grow_patch(int seed_row, int seed_col, int *added_ids,
               SEGMENT *developed, SEGMENT *probability,
               int num_neighbors, double alpha, int patch_size,
               int step, enum slow_grow strategy);

#endif // FUTURES_PATCH_H
