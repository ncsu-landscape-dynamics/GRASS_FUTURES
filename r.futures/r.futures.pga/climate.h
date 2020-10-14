#ifndef FUTURES_CLIMATE_H
#define FUTURES_CLIMATE_H

#include <grass/segment.h>

#include "inputs.h"
#include "keyvalue.h"

float get_max_HAND(struct Segments *segments, const struct BBox *bbox, float flood_probability);
bool generate_flood(const struct KeyValueIntFloat *flood_probability_map, int region_idx, float *flood_probability);
float get_abandonment_probability(struct Segments *segments, float flood_depth, int row, int col);

#endif // FUTURES_CLIMATE_H
