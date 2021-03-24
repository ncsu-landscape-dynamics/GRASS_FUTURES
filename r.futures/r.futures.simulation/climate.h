#ifndef FUTURES_CLIMATE_H
#define FUTURES_CLIMATE_H

#include <grass/segment.h>

#include "inputs.h"
#include "map.h"

enum FloodResponse { Retreat, Adapt, Stay };

struct ACDamageRelation
{
    /* y = ax + b */
    float resilience_a;
    float resilience_b;
    float vulnerability_a;
    float vulnerability_b;
};

struct HAND_bbox_values
{
    float *array;
    uint size;
};

void initilize_adaptation(SEGMENT *adaptation,
                          const struct SegmentMemory *segment_info);
float get_max_HAND(struct Segments *segments, const struct BBox *bbox, float flood_probability,
                   struct HAND_bbox_values *HAND_bbox_vals, float percentile);
bool generate_flood(map_float_t *flood_probability_map, int region_idx, float *flood_probability);
float get_damage(struct Segments *segments, const struct DepthDamageFunctions *ddf,
                 float flood_depth, int row, int col);
enum FloodResponse flood_response(float damage, float adaptive_capacity,
                                  const struct ACDamageRelation *response);
bool is_adapted(SEGMENT *adaptation, int row, int col);
void adapt(SEGMENT *adaptation, int row, int col);

#endif // FUTURES_CLIMATE_H
