/*!
   \file climate.c
   
   \brief Functions for climate scenarios (flooding)
   
   (C) 2020 by Anna Petrasova and the GRASS Development Team
   
   This program is free software under the GNU General Public License
   (>=v2).  Read the file COPYING that comes with GRASS for details.
   
   \author Anna Petrasova
 */
#include <math.h>

#include <grass/gis.h>
#include <grass/segment.h>

#include "inputs.h"
#include "keyvalue.h"

static float depth_to_damage(float depth, float height, float max_damage, float curve_shape)
{
    return pow((depth / height), curve_shape) * max_damage;
}

bool generate_flood(const struct KeyValueIntFloat *flood_probability_map, int region_idx, float *flood_probability)
{
    float max_flood_probability;
    double p;

    KeyValueIntFloat_find(flood_probability_map, region_idx, &max_flood_probability);
    if (max_flood_probability > 0) {
        p = G_drand48();
        if (p <= max_flood_probability) {
            *flood_probability = p;
            return true;
        }
    }
    return false;
}

float get_max_HAND(struct Segments *segments, const struct BBox *bbox, float flood_probability)
{
    FCELL flood_probability_value;
    FCELL HAND_value;
    FCELL max_HAND_value;
    int row, col;

    max_HAND_value = 0;
    for (row = bbox->n; row <= bbox->s; row++)
        for (col = bbox->w; col <= bbox->w; col++) {
            Segment_get(&segments->flood_probability, (void *)&flood_probability_value, row, col);
            if (flood_probability_value >= flood_probability) {
                Segment_get(&segments->HAND, (void *)&HAND_value, row, col);
                if (HAND_value > max_HAND_value) {
                    max_HAND_value = HAND_value;
                }
            }
        }
    return max_HAND_value;
}

float get_abandonment_probability(struct Segments *segments, float flood_level, int row, int col)
{
    FCELL HAND_value;
    float depth;
    float damage;

    damage = 0;
    Segment_get(&segments->HAND, (void *)&HAND_value, row, col);
    depth = flood_level - HAND_value;
    if (depth > 0) {
        damage = depth_to_damage(depth, 5, 0.6, 0.5);
        // adaptive capacity here
    }
    return damage;
}
