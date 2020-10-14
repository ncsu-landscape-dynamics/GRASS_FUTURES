/*!
   \file climate.c
   
   \brief Functions for climate scenarios (flooding)
   
   (C) 2020 by Anna Petrasova and the GRASS Development Team
   
   This program is free software under the GNU General Public License
   (>=v2).  Read the file COPYING that comes with GRASS for details.
   
   \author Anna Petrasova
 */

#include <grass/gis.h>
#include <grass/segment.h>

#include "inputs.h"
#include "keyvalue.h"


bool generate_flood(struct KeyValueIntFloat *flood_probability_map, int region_idx, float *flood_probability)
{
    float max_flood_probability;
    double p;
    
    KeyValueIntFloat_find(flood_probability_map, region_idx, &max_flood_probability);
    p = G_drand48();
    if (p <= max_flood_probability) {
        *flood_probability = p;
        return true;
    }
    return false;
}

float get_max_HAND(struct Segments *segments, struct BBox *bbox, float flood_probability)
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