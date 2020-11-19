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
#include <grass/raster.h>
#include <grass/segment.h>
#include <grass/glocale.h>

#include "inputs.h"
#include "keyvalue.h"
#include "climate.h"
#include "random.h"

void initilize_adaptation(SEGMENT *adaptation,
                          const struct SegmentMemory *segment_info)
{
    CELL *adaptation_row;
    int row;
    int rows;

    rows = Rast_window_rows();
    adaptation_row = Rast_allocate_buf(CELL_TYPE); /* uses calloc */
    if (Segment_open(adaptation, G_tempfile(), Rast_window_rows(),
                     Rast_window_cols(), segment_info->rows, segment_info->cols,
                     Rast_cell_size(CELL_TYPE), segment_info->in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map"));
    /* make sure there are zeroes */
    for (row = 0; row < rows; row++)
        Segment_put_row(adaptation, adaptation_row, row);
    Segment_flush(adaptation);
    G_free(adaptation_row);
}

void adapt(SEGMENT *adaptation, int row, int col)
{
    int adapted;

    /*  this will be level of adaptation (flood probability or depth) */
    adapted = 1;
    Segment_put(adaptation, (void *)&adapted, row, col);
}

bool is_adapted(SEGMENT *adaptation, int row, int col)
{
    int adapted;

    Segment_get(adaptation, (void *)&adapted, row, col);
    return (bool)adapted;
}

static float depth_to_damage(float depth, const struct DepthDamageFunc *func)
{
    return pow((depth / func->H), func->r) * func->M / 100;
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

float get_damage(struct Segments *segments, const struct DepthDamageFunc *func,
                 float flood_level, int row, int col)
{
    FCELL HAND_value;
    float depth;
    float damage;

    damage = 0;
    Segment_get(&segments->HAND, (void *)&HAND_value, row, col);
    depth = flood_level - HAND_value;
    if (depth > 0) {
        if (!is_adapted(&segments->adaptation, row, col))
            damage = depth_to_damage(depth, func);
    }
    return damage;
}

/*!
 * \brief Stochastically simulates flood response (retreat, adaptation, stay)
 * as a function of damage and adaptive capacity.
 * \param damage damage (0-1)
 * \param adaptive_capacity adaptive capacity (0-1)
 * \param response structure containing info about damage-AC relation
 * \return one of Retreat, Adapt, Stay
 */
enum FloodResponse flood_response(float damage, float adaptive_capacity,
                                  const struct ACDamageRelation *response)
{
    double x, y;
    double response_val;

    gauss_xy(0, 0.1, &x, &y);
    if (adaptive_capacity > 0) {
        response_val = response->resilience_a * (adaptive_capacity + y) + response->resilience_b;
        if (damage + x > response_val)
            return Retreat;
        else
            return Adapt;
    }
    else {
        response_val = response->vulnerability_a * (adaptive_capacity + y) + response->vulnerability_b;
        if (damage + x > response_val)
            return Retreat;
        else
            return Stay;
    }
}
