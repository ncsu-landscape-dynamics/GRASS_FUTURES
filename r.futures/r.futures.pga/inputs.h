#ifndef FUTURES_INPUTS_H
#define FUTURES_INPUTS_H

#include <stdbool.h>
#include <grass/segment.h>

#include "keyvalue.h"

enum development_type {DEV_TYPE_INITIAL = 0,
                       DEV_TYPE_UNDEVELOPED = -1,
                       DEV_TYPE_ABANDONED = -2};

struct Demand
{
    const char *cells_filename;
    const char *population_filename;
    float **cells_table;
    float **population_table;
    int *years;
    int max_subregions;
    int max_steps;
    const char *separator;
    bool has_population;
};

struct Potential
{
    const char *filename;
    double **predictors;
    double *intercept;
    double *devpressure;
    int *predictor_indices;
    int max_predictors;
    int max_subregions;
    float *incentive_transform;
    int incentive_transform_size;
    const char *separator;
};

struct PatchSizes
{
    const char *filename;
    // array of patches
    int **patch_sizes;
    // array of number of patches per area
    int *patch_count;
    // maximum patch size
    int max_patch_size;
    // use single column for all regions
    bool single_column;

};

struct SegmentMemory
{
    int rows;
    int cols;
    int in_memory;
};

struct Segments
{
    SEGMENT developed;
    SEGMENT subregions;
    SEGMENT potential_subregions;
    SEGMENT devpressure;
    SEGMENT predictors;
    SEGMENT probability;
    SEGMENT weight;
    SEGMENT density;
    SEGMENT density_capacity;
    SEGMENT HAND;
    SEGMENT HUC;
    SEGMENT flood_probability;
    SEGMENT adaptive_capacity;
    bool use_weight;
    bool use_potential_subregions;
    bool use_density;
    bool use_climate;
};

struct RasterInputs
{
    const char *developed;
    const char *regions;
    const char *potential_regions;
    char **predictors;
    const char *devpressure;
    const char *weights;
    const char *density;
    const char *density_capacity;
    const char *HAND;
    const char *HUC;
    const char *flood_probability;
    const char *adaptive_capacity;
};


struct DevelopableCell
{

    size_t id;
    float probability;
    float cumulative_probability;
    bool tried;
};

struct Developables
{
    int max_subregions;
    size_t *max;
    size_t *num;
    struct DevelopableCell **cells;
};

struct BBox
{
    int n;
    int s;
    int e;
    int w;
};

struct BBoxes
{
    int n_bbox;
    int max_bbox;
    struct KeyValueIntInt *map;
    struct BBox *bbox;
};


void initialize_incentive(struct Potential *potential_info, float exponent);
void read_input_rasters(struct RasterInputs inputs, struct Segments *segments,
                        struct SegmentMemory segment_info, struct KeyValueIntInt *region_map,
                        struct KeyValueIntInt *reverse_region_map,
                        struct KeyValueIntInt *potential_region_map,
                        struct KeyValueCharInt *predictor_map, int num_predictors,
                        struct KeyValueIntInt *HUC_map,
                        struct KeyValueIntFloat *max_flood_probability_map);
void read_demand_file(struct Demand *demandInfo, struct KeyValueIntInt *region_map);
void read_potential_file(struct Potential *potentialInfo, struct KeyValueIntInt *region_map,
                         struct KeyValueCharInt *predictor_map);
void read_patch_sizes(struct PatchSizes *patch_sizes, struct KeyValueIntInt *region_map,
                      double discount_factor);
void create_bboxes(SEGMENT *raster, SEGMENT *masking, struct BBoxes *bboxes);

#endif // FUTURES_INPUTS_H
