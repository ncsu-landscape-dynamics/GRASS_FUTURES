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

/*!
 * \brief Initialize adaptation segment
 * \param adaptation segment
 * \param segment_info
 */
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
/*!
 * \brief Adapt pixel to flooding.
 *
 * Currently, only it's 0 or 1,
 * for future there should be level of adaptation.
 *
 * \param adaptation Adaptation segment
 * \param row
 * \param col
 */
void adapt(SEGMENT *adaptation, int row, int col)
{
    int adapted;

    /*  this will be level of adaptation (flood probability or depth) */
    adapted = 1;
    Segment_put(adaptation, (void *)&adapted, row, col);
}
/*!
 * \brief Check if pixel is adapted.
 *
 * \param adaptation Adaptation segment
 * \param row
 * \param col
 * \return
 */
bool is_adapted(SEGMENT *adaptation, int row, int col)
{
    int adapted;

    Segment_get(adaptation, (void *)&adapted, row, col);
    return (bool)adapted;
}

/*!
 * \brief Converts water depth to structural damage
 * using depth-damage-function.
 *
 * \param water depth (height from bottom of structure)
 * \param func depth-damage function
 * \return structural damage (0: no damage, 1: total destruction)
 */
static float depth_to_damage(float depth, int region_idx,
                             const struct DepthDamageFunctions *ddf)
{
    int i;

    for (i = 0; i < ddf->max_levels; i++) {
        if (depth < ddf->levels[i])
            return ddf->damage[region_idx][i] / 100;
    }
    return ddf->damage[region_idx][ddf->max_levels - 1] / 100;
}

/*!
 * \brief Returns region index for which to look up
 * DDF. If DDF is the same for entire area, returns 1.
 * Doesn't check for nulls.
 *
 * \param segments Segments
 * \param ddf DepthDamageFunctions structure
 * \param row row
 * \param col col
 * \return index of a region
 */
static int get_DDF_region_index(struct Segments *segments,
                                const struct DepthDamageFunctions *ddf,
                                int row, int col) {
    CELL DDF_region_idx;

    if (ddf->use_DDF_subregions)
        Segment_get(&segments->DDF_subregions, (void *)&DDF_region_idx, row, col);
    else if (ddf->use_subregions)
        Segment_get(&segments->subregions, (void *)&DDF_region_idx, row, col);
    else if (ddf->use_potential_subregions)
        Segment_get(&segments->potential_subregions, (void *)&DDF_region_idx, row, col);
    else
        DDF_region_idx = 1;

    return DDF_region_idx;
}
/*!
 * \brief Decide if flood event occurrs and if yes,
 * what probability is the flood associated with
 * (e.g. 100yr flood = 0.01)
 *
 * \param flood_probability_map contains flood frequencies/probabilities (0-1)
 * \param region_idx region index (HUC)
 * \param flood_probability resulting flood frequency (as probability)
 * \return true if flood event occures, otherwise false
 */
bool generate_flood(const struct KeyValueIntFloat *flood_probability_map,
                    int region_idx, float *flood_probability)
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

/*!
 * \brief Get maximum Height Above Nearest Drainage within flooded area.
 *
 * This assumes flood probabilities do not correspond well to HAND,
 * therefore we take worst case inundation height.
 *
 * \param segments segments (flood probability and HAND)
 * \param bbox BBox of area where to search for it
 * \param flood_probability
 * \return max HAND value
 */
float get_max_HAND(struct Segments *segments, const struct BBox *bbox,
                   float flood_probability)
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

/*!
 * \brief Get damage caused by flooding.
 *
 * Checks if pixel was adapted (in the future adaptation
 * could be for certain level).
 * If yes, damage is 0.
 *
 * \param segments HAND and adaptation segments
 * \param ddf depth-damage functions
 * \param flood_level height of inundation
 * \param row row
 * \param col col
 * \return structural damage (0: no damage, 1: total destruction)
 */
float get_damage(struct Segments *segments, const struct DepthDamageFunctions *ddf,
                 float flood_level, int row, int col)
{
    FCELL HAND_value;
    int DDF_region_idx;
    float depth;
    float damage;

    damage = 0;
    Segment_get(&segments->HAND, (void *)&HAND_value, row, col);
    depth = flood_level - HAND_value;
    if (depth > 0) {
        if (!is_adapted(&segments->adaptation, row, col)) {
            DDF_region_idx = get_DDF_region_index(segments, ddf, row, col);
            damage = depth_to_damage(depth, DDF_region_idx, ddf);
        }
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
