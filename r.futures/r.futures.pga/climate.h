#ifndef FUTURES_CLIMATE_H
#define FUTURES_CLIMATE_H

#include <grass/segment.h>

#include "inputs.h"
#include "keyvalue.h"

float get_max_HAND(struct Segments *segments, struct BBox *bbox, float flood_probability);
bool generate_flood(struct KeyValueIntFloat *flood_probability_map, int region_idx, float *flood_probability);

#endif // FUTURES_CLIMATE_H
