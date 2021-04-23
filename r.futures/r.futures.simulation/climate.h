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

float get_max_HAND(struct Segments *segments, const struct BBox *bbox, float flood_probability,
                   struct HAND_bbox_values *HAND_bbox_vals, float percentile);
float get_depth(SEGMENT *flood_depths, float flood_probability,
                int row, int col, FCELL *values,
                const struct FloodInputs *flood_inputs);
float get_depth_flood_level(SEGMENT *hand, float flood_level,
                            int row, int col);
bool generate_flood(map_float_t *flood_probability_map, int region_idx, float *flood_probability);
float get_damage(struct Segments *segments, const struct DepthDamageFunctions *ddf,
                 float flood_probability, float depth, int row, int col);
enum FloodResponse flood_response(float damage, float adaptive_capacity,
                                  const struct ACDamageRelation *response);
void initialize_flood_response(struct ACDamageRelation *response_relation);
bool is_adapted(SEGMENT *adaptation, float flood_probability, int row, int col);
void adapt(SEGMENT *adaptation, float flood_probability, int row, int col);
void stay(SEGMENT *adaptation, int row, int col);

#endif // FUTURES_CLIMATE_H
