#ifndef FUTURES_SIMULATION_H
#define FUTURES_SIMULATION_H

#include <grass/gis.h>

#include "inputs.h"
#include "patch.h"
#include "redistribute.h"


enum seed_search {RANDOM, PROBABILITY};

int find_probable_seed(struct Undeveloped *undev_cells, int region);
int get_seed(struct Undeveloped *undev_cells, int region_idx, enum seed_search method,
              int *row, int *col);
double get_develop_probability_xy(struct Segments *segments,
                                  FCELL *values,
                                  struct Potential *potential_info,
                                  int region_index, int row, int col);
void move(struct Segments *segments, const struct RedistributionMatrix *matrix,
          struct Demand *demand, struct KeyValueIntInt *region_map, int step);
void recompute_probabilities(struct Undeveloped *undeveloped_cells,
                             struct Segments *segments,
                             struct Potential *potential_info);
void compute_step(struct Undeveloped *undev_cells, struct Demand *demand,
                  enum seed_search search_alg,
                  struct Segments *segments,
                  struct PatchSizes *patch_sizes, struct PatchInfo *patch_info,
                  struct DevPressure *devpressure_info, int *patch_overflow,
                  struct RedistributionMatrix *redistr_matrix, struct KeyValueIntInt *region_map,
                  int step, int region_idx, int region_ID, bool overgrow);

#endif // FUTURES_SIMULATION_H
