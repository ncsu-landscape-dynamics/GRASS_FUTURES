#ifndef FUTURES_INPUTS_H
#define FUTURES_INPUTS_H

#include <grass/segment.h>

#include "keyvalue.h"

struct Demand
{
    const char *filename;
    int **table;
    int max_subregions;
    int max_steps;
};

struct Potential
{
    const char *filename;
    double **predictors;
    double *intercept;
    double *devpressure;
    int max_predictors;
    int max_subregions;
};

struct PatchSizes
{
    const char *filename;
    int max_patches;
    int *patch_sizes;
    
};

struct SegmentMemory
{
    int rows;
    int cols;
    int in_memory;
};

struct Segments
{
    SEGMENT developed_segment;
    SEGMENT subregions_segment;
    SEGMENT devpressure_segment;
    SEGMENT predictors_segment;
    SEGMENT probability_segment;
};

struct UndevelopedCell
{

    size_t id;
    float probability;
    float cumulative_probability;
};

struct Undeveloped
{
    int max_subregions;
    size_t *max;
    size_t *num;
    struct UndevelopedCell **cells;
};


void rast_segment_open(const char *name, SEGMENT *segment, struct SegmentMemory segmentInfo, 
                       RASTER_MAP_TYPE map_type);
size_t get_idx_from_xy(int row, int col, int cols);
void get_xy_from_idx(size_t idx, int cols, int *row, int *col);

void read_developed(char *filename, struct Segments *segments,
                    struct SegmentMemory segment_info);

void read_predictors(char **predictor_names, struct Segments *segments,
                     struct SegmentMemory segmentInfo, int ninputs);
void read_subregions(const char *subregions, struct Segments *segments,
                    struct KeyValueIntInt *region_map);
void read_demand_file(struct Demand *demandInfo, struct KeyValueIntInt *region_map);
void read_potential_file(struct Potential *potentialInfo, struct KeyValueIntInt *region_map,
                         int num_predictors);
void read_patch_sizes(struct PatchSizes *patch_info, double discount_factor);

#endif // FUTURES_INPUTS_H