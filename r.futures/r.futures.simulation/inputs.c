/*!
   \file inputs.c

   \brief Functions to read in input files and rasters

   (C) 2016-2019 by Vaclav Petras and the GRASS Development Team

   This program is free software under the GNU General Public License
   (>=v2).  Read the file COPYING that comes with GRASS for details.

   \author Vaclav Petras
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <grass/segment.h>

#include "inputs.h"
#include "map.h"

/*!
 * \brief Initialize arrays for transformation of probability values
 * \param potential_info
 * \param exponent
 */
void initialize_incentive(struct Potential *potential_info, float exponent)
{
    int i;

    potential_info->incentive_transform_size = 1001;
    potential_info->incentive_transform = (float *) G_malloc(sizeof(float) *
                                                             potential_info->incentive_transform_size);
    i = 0;
    double step = 1. / (potential_info->incentive_transform_size - 1);
    while (i < potential_info->incentive_transform_size) {
        potential_info->incentive_transform[i] = pow(i * step, exponent);
        i++;
    }
}

/*!
 * \brief Read input rasters into segments.
 * \param inputs
 * \param segments
 * \param segment_info
 * \param region_map
 */
void read_input_rasters(struct RasterInputs inputs, struct Segments *segments,
                        struct SegmentMemory segment_info, map_int_t *region_map,
                        map_int_t *reverse_region_map,
                        map_int_t *potential_region_map,
                        map_int_t *HUC_map,
                        map_float_t *max_flood_probability_map,
                        map_int_t *DDF_region_map)
{
    int row, col;
    int rows, cols;
    int fd_developed, fd_reg, fd_devpressure, fd_weights,
            fd_pot_reg, fd_density, fd_density_cap,
            fd_HAND, fd_flood_probability, fd_adaptive_capacity,
            fd_HUC, fd_DDF;
    int count_regions, pot_count_regions, HUC_count, DDF_count;
    int region_index, pot_region_index, HUC_index, DDF_index;
    int *region_pindex;
    int *pot_region_pindex;
    float max_flood_probability;
    CELL c;
    FCELL fc;
    bool isnull;

    CELL *developed_row;
    CELL *subregions_row;
    CELL *pot_subregions_row;
    FCELL *devpressure_row;
    FCELL *weights_row;
    FCELL *density_row;
    FCELL *density_capacity_row;
    FCELL *HAND_row;
    FCELL *flood_probability_row;
    FCELL *adaptive_capacity_row;
    CELL *HUC_row;
    CELL *DDF_row;


    rows = Rast_window_rows();
    cols = Rast_window_cols();
    count_regions = region_index = HUC_count = DDF_count = 0;
    pot_count_regions = pot_region_index = HUC_index = DDF_index = 0;
    int *pindex;
    float *pvalue;

    /* open existing raster maps for reading */
    fd_developed = Rast_open_old(inputs.developed, "");
    fd_reg = Rast_open_old(inputs.regions, "");
    if (segments->use_potential_subregions)
        fd_pot_reg = Rast_open_old(inputs.potential_regions, "");
    fd_devpressure = Rast_open_old(inputs.devpressure, "");
    if (segments->use_weight)
        fd_weights = Rast_open_old(inputs.weights, "");
    if (segments->use_density) {
        fd_density = Rast_open_old(inputs.density, "");
        fd_density_cap = Rast_open_old(inputs.density_capacity, "");
    }
    if (segments->use_climate) {
        fd_HAND = Rast_open_old(inputs.HAND, "");
        fd_flood_probability = Rast_open_old(inputs.flood_probability, "");
        fd_adaptive_capacity = Rast_open_old(inputs.adaptive_capacity, "");
        fd_HUC = Rast_open_old(inputs.HUC, "");
        if (inputs.DDF_regions)
            fd_DDF = Rast_open_old(inputs.DDF_regions, "");
    }

    /* Segment open developed */
    if (Segment_open(&segments->developed, G_tempfile(), rows,
                     cols, segment_info.rows, segment_info.cols,
                     Rast_cell_size(CELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map of development"));
    /* Segment open subregions */
    if (Segment_open(&segments->subregions, G_tempfile(), rows,
                     cols, segment_info.rows, segment_info.cols,
                     Rast_cell_size(CELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map of subregions"));
    /* Segment open development pressure */
    if (Segment_open(&segments->devpressure, G_tempfile(), rows,
                     cols, segment_info.rows, segment_info.cols,
                     Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of a raster map of development pressure"));
    /* Segment open weights */
    if (segments->use_weight)
        if (Segment_open(&segments->weight, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of weights"));
    /* Segment open potential_subregions */
    if (segments->use_potential_subregions)
        if (Segment_open(&segments->potential_subregions, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(CELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of weights"));
    /* Segment open density */
    if (segments->use_density) {
        if (Segment_open(&segments->density, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of density"));
        if (Segment_open(&segments->density_capacity, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of density capacity"));
    }
    /* Segment open HAND */
    if (segments->use_climate) {
        if (Segment_open(&segments->HAND, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of HAND"));
        if (Segment_open(&segments->flood_probability, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of flood probability"));
        if (Segment_open(&segments->adaptive_capacity, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of adaptive capacity"));
        if (Segment_open(&segments->HUC, G_tempfile(), rows,
                         cols, segment_info.rows, segment_info.cols,
                         Rast_cell_size(CELL_TYPE), segment_info.in_memory) != 1)
            G_fatal_error(_("Cannot create temporary file with segments of a raster map of HUCs"));
        if (inputs.DDF_regions)
            if (Segment_open(&segments->DDF_subregions, G_tempfile(), rows,
                             cols, segment_info.rows, segment_info.cols,
                             Rast_cell_size(CELL_TYPE), segment_info.in_memory) != 1)
                G_fatal_error(_("Cannot create temporary file with segments of a raster map of DDF subregions"));
    }
    developed_row = Rast_allocate_buf(CELL_TYPE);
    subregions_row = Rast_allocate_buf(CELL_TYPE);
    devpressure_row = Rast_allocate_buf(FCELL_TYPE);
    if (segments->use_weight)
        weights_row = Rast_allocate_buf(FCELL_TYPE);
    if (segments->use_potential_subregions)
        pot_subregions_row = Rast_allocate_buf(CELL_TYPE);
    if (segments->use_density) {
        density_row = Rast_allocate_buf(FCELL_TYPE);
        density_capacity_row = Rast_allocate_buf(FCELL_TYPE);
    }
    if (segments->use_climate) {
        HAND_row = Rast_allocate_buf(FCELL_TYPE);
        flood_probability_row = Rast_allocate_buf(FCELL_TYPE);
        adaptive_capacity_row = Rast_allocate_buf(FCELL_TYPE);
        HUC_row = Rast_allocate_buf(CELL_TYPE);
        if (inputs.DDF_regions)
            DDF_row = Rast_allocate_buf(CELL_TYPE);
    }

    for (row = 0; row < rows; row++) {
        G_percent(row, rows, 5);
        /* read developed row */
        Rast_get_row(fd_developed, developed_row, row, CELL_TYPE);
        Rast_get_row(fd_devpressure, devpressure_row, row, FCELL_TYPE);
        Rast_get_row(fd_reg, subregions_row, row, CELL_TYPE);
        if (segments->use_weight)
            Rast_get_row(fd_weights, weights_row, row, FCELL_TYPE);
        if (segments->use_potential_subregions)
            Rast_get_row(fd_pot_reg, pot_subregions_row, row, CELL_TYPE);
        if (segments->use_density) {
            Rast_get_row(fd_density, density_row, row, FCELL_TYPE);
            Rast_get_row(fd_density_cap, density_capacity_row, row, FCELL_TYPE);
        }
        if (segments->use_climate) {
            Rast_get_row(fd_HAND, HAND_row, row, FCELL_TYPE);
            Rast_get_row(fd_flood_probability, flood_probability_row, row, FCELL_TYPE);
            Rast_get_row(fd_adaptive_capacity, adaptive_capacity_row, row, FCELL_TYPE);
            Rast_get_row(fd_HUC, HUC_row, row, CELL_TYPE);
            if (inputs.DDF_regions)
                Rast_get_row(fd_DDF, DDF_row, row, CELL_TYPE);
        }
        for (col = 0; col < cols; col++) {
            isnull = false;
            /* developed */
            /* undeveloped 0 -> -1, developed 1 -> 0 */
            if (!Rast_is_null_value(&((CELL *) developed_row)[col], CELL_TYPE)) {
                c = ((CELL *) developed_row)[col];
                ((CELL *) developed_row)[col] = c - 1;
            }
            else
                isnull = true;
            /* subregions */
            if (!Rast_is_null_value(&((CELL *) subregions_row)[col], CELL_TYPE)) {
                c = ((CELL *) subregions_row)[col];
                region_pindex = map_get_int(region_map, c);
                if (!region_pindex) {
                    map_set_int(region_map, c, count_regions);
                    map_set_int(reverse_region_map, count_regions, c);
                    region_index = count_regions;
                    count_regions++;
                }
                else
                    region_index = *region_pindex;
                ((CELL *) subregions_row)[col] = region_index;
            }
            else
                isnull = true;
            if (segments->use_potential_subregions) {
                if (!Rast_is_null_value(&((CELL *) pot_subregions_row)[col], CELL_TYPE)) {
                    c = ((CELL *) pot_subregions_row)[col];
                    pot_region_pindex = map_get_int(potential_region_map, c);
                    if (!pot_region_pindex) {
                        map_set_int(potential_region_map, c, pot_count_regions);
                        pot_region_index = pot_count_regions;
                        pot_count_regions++;
                    }
                    else
                        pot_region_index = *pot_region_pindex;
                    ((CELL *) pot_subregions_row)[col] = pot_region_index;
                }
                else
                    isnull = true;
            }
            /* devpressure - just check nulls */
            if (Rast_is_null_value(&((FCELL *) devpressure_row)[col], FCELL_TYPE))
                isnull = true;
            /* density - just check nulls */
            if (segments->use_density) {
                if (Rast_is_null_value(&((FCELL *) density_row)[col], FCELL_TYPE))
                    isnull = true;
                if (Rast_is_null_value(&((FCELL *) density_capacity_row)[col], FCELL_TYPE))
                    isnull = true;
            }
            /* weights - must be in range -1, 1*/
            if (segments->use_weight) {
                if (Rast_is_null_value(&((FCELL *) weights_row)[col], FCELL_TYPE)) {
                    ((FCELL *) weights_row)[col] = 0;
                    isnull = true;
                }
                else {
                    fc = ((FCELL *) weights_row)[col];
                    if (fc > 1) {
                        G_warning(_("Probability weights are > 1, truncating..."));
                        fc = 1;
                    }
                    else if (fc < -1) {
                        fc = -1;
                        G_warning(_("Probability weights are < -1, truncating..."));
                    }
                    ((FCELL *) weights_row)[col] = fc;
                }
            }
            /* flooding */
            if (segments->use_climate) {
                if (Rast_is_null_value(&((FCELL *) HAND_row)[col], FCELL_TYPE))
                    isnull = true;
                if (Rast_is_null_value(&((FCELL *) adaptive_capacity_row)[col], FCELL_TYPE))
                    isnull = true;
                if (!Rast_is_null_value(&((CELL *) HUC_row)[col], CELL_TYPE)) {
                    c = ((CELL *) HUC_row)[col];
                    /* mapping: HUC id -> index */
                    pindex = map_get_int(HUC_map, c);
                    if (!pindex) {
                        map_set_int(HUC_map, c, HUC_count);
                        HUC_index = HUC_count;
                        HUC_count++;
                    }
                    else
                        HUC_index = *pindex;
                    ((CELL *) HUC_row)[col] = HUC_index;
                    /* save the max flood value for each HUC */
                    if (!Rast_is_null_value(&((FCELL *) flood_probability_row)[col], FCELL_TYPE)) {
                        fc = ((FCELL *) flood_probability_row)[col];
                        pvalue = map_get_int(max_flood_probability_map, HUC_index);
                        if (pvalue) {
                            max_flood_probability = *pvalue;
                            if (fc > max_flood_probability)
                                map_set_int(max_flood_probability_map, HUC_index, fc);
                        }
                        else
                            map_set_int(max_flood_probability_map, HUC_index, fc);
                    }
                }
                /* DDF subregions */
                if (inputs.DDF_regions) {
                    if (!Rast_is_null_value(&((CELL *) DDF_row)[col], CELL_TYPE)) {
                        c = ((CELL *) DDF_row)[col];
                        pindex = map_get_int(DDF_region_map, c);
                        if(!pindex) {
                            map_set_int(DDF_region_map, c, DDF_count);
                            DDF_index = DDF_count;
                            DDF_count++;
                        }
                        else
                            DDF_index = *pindex;
                        ((CELL *) DDF_row)[col] = DDF_index;
                    }
                    else
                        isnull = true;
                }
            }
            /* if in developed, subregions, devpressure or weights are any nulls
               propagate them into developed */
            if (isnull)
                Rast_set_c_null_value(&((CELL *) developed_row)[col], 1);
        }
        Segment_put_row(&segments->developed, developed_row, row);
        Segment_put_row(&segments->devpressure, devpressure_row, row);
        Segment_put_row(&segments->subregions, subregions_row, row);
        if (segments->use_weight)
            Segment_put_row(&segments->weight, weights_row, row);
        if (segments->use_potential_subregions)
            Segment_put_row(&segments->potential_subregions, pot_subregions_row, row);
        if (segments->use_density) {
            Segment_put_row(&segments->density, density_row, row);
            Segment_put_row(&segments->density_capacity, density_capacity_row, row);
        }
        if (segments->use_climate) {
            Segment_put_row(&segments->HAND, HAND_row, row);
            Segment_put_row(&segments->flood_probability, flood_probability_row, row);
            Segment_put_row(&segments->adaptive_capacity, adaptive_capacity_row, row);
            Segment_put_row(&segments->HUC, HUC_row, row);
            if (inputs.DDF_regions)
                Segment_put_row(&segments->DDF_subregions, DDF_row, row);
        }
    }
    G_percent(row, rows, 5);

    /* flush all segments */
    Segment_flush(&segments->developed);
    Segment_flush(&segments->subregions);
    Segment_flush(&segments->devpressure);
    if (segments->use_weight)
        Segment_flush(&segments->weight);
    if (segments->use_potential_subregions)
        Segment_flush(&segments->potential_subregions);
    if (segments->use_density) {
        Segment_flush(&segments->density);
        Segment_flush(&segments->density_capacity);
    }
    if (segments->use_climate) {
        Segment_flush(&segments->HAND);
        Segment_flush(&segments->flood_probability);
        Segment_flush(&segments->adaptive_capacity);
        Segment_flush(&segments->HUC);
        if (inputs.DDF_regions)
            Segment_flush(&segments->DDF_subregions);
    }
    /* close raster maps */
    Rast_close(fd_developed);
    Rast_close(fd_reg);
    Rast_close(fd_devpressure);
    if (segments->use_weight)
        Rast_close(fd_weights);
    if (segments->use_potential_subregions)
        Rast_close(fd_pot_reg);
    if (segments->use_density) {
        Rast_close(fd_density);
        Rast_close(fd_density_cap);
    }
    if (segments->use_climate) {
        Rast_close(fd_HAND);
        Rast_close(fd_flood_probability);
        Rast_close(fd_adaptive_capacity);
        Rast_close(fd_HUC);
        if (inputs.DDF_regions)
            Rast_close(fd_DDF);
    }

    G_free(developed_row);
    G_free(subregions_row);
    G_free(devpressure_row);
    if (segments->use_weight)
        G_free(weights_row);
    if (segments->use_potential_subregions)
        G_free(pot_subregions_row);
    if (segments->use_density) {
        G_free(density_row);
        G_free(density_capacity_row);
    }
    if (segments->use_climate) {
        G_free(HAND_row);
        G_free(flood_probability_row);
        G_free(adaptive_capacity_row);
        G_free(HUC_row);
        if (inputs.DDF_regions)
            G_free(DDF_row);
    }
}

/*!
 * \brief Reads predictors and aggregates them with Potential table:
 * x_1 * a + x2 * b + ...
 * Saves memory comparing to having them separately.
 *
 * \param inputs Raster inputs
 * \param segments Segments
 * \param potential Potential table
 * \param segment_info Segment memory info
 */
void read_predictors(struct RasterInputs inputs, struct Segments *segments,
                     const struct Potential *potential,
                     const struct SegmentMemory segment_info)
{
    int i;
    int row, col;
    int rows, cols;
    int *fds_predictors;
    int pred_index;
    FCELL value;
    CELL pot_index, dev_value;
    FCELL **predictor_rows;
    FCELL *aggregated_row;

    rows = Rast_window_rows();
    cols = Rast_window_cols();
    fds_predictors = G_malloc(potential->max_predictors * sizeof(int));
    for (i = 0; i < potential->max_predictors; i++) {
        fds_predictors[i] = Rast_open_old(inputs.predictors[i], "");
    }
    predictor_rows = G_malloc(potential->max_predictors * sizeof(FCELL *));
    for (i = 0; i < potential->max_predictors; i++) {
        predictor_rows[i] = Rast_allocate_buf(FCELL_TYPE);
    }
    aggregated_row = Rast_allocate_buf(FCELL_TYPE);

    /* Segment open predictors */
    if (Segment_open(&segments->aggregated_predictor, G_tempfile(), rows,
                     cols, segment_info.rows, segment_info.cols,
                     Rast_cell_size(FCELL_TYPE), segment_info.in_memory) != 1)
        G_fatal_error(_("Cannot create temporary file with segments of predictor raster maps"));

    /* read in */
    for (row = 0; row < rows; row++) {
        for (i = 0; i < potential->max_predictors; i++) {
            Rast_get_row(fds_predictors[i], predictor_rows[i], row, FCELL_TYPE);
        }
        for (col = 0; col < cols; col++) {
            ((FCELL *) aggregated_row)[col] = 0;
            Segment_get(&segments->developed, (void *)&dev_value, row, col);
            if (Rast_is_null_value(&dev_value, CELL_TYPE)) {
                continue;
            }
            for (i = 0; i < potential->max_predictors; i++) {
                /* collect all nulls in predictors and set it in output raster */
                if (Rast_is_null_value(&((FCELL *) predictor_rows[i])[col], FCELL_TYPE)) {
                    Rast_set_c_null_value(&dev_value, 1);
                    Segment_put(&segments->developed, (void *)&dev_value, row, col);
                    break;
                }
                if (segments->use_potential_subregions)
                    Segment_get(&segments->potential_subregions, (void *)&pot_index, row, col);
                else
                    Segment_get(&segments->subregions, (void *)&pot_index, row, col);
                pred_index = potential->predictor_indices[i];
                value = potential->predictors[i][pot_index] * ((FCELL *) predictor_rows[pred_index])[col];
                ((FCELL *) aggregated_row)[col] += value;
            }
        }
        Segment_put_row(&segments->aggregated_predictor, aggregated_row, row);
    }
    Segment_flush(&segments->aggregated_predictor);
    Segment_flush(&segments->developed);
    for (i = 0; i < potential->max_predictors; i++) {
        Rast_close(fds_predictors[i]);
        G_free(predictor_rows[i]);
    }
    G_free(fds_predictors);
    G_free(predictor_rows);
    G_free(aggregated_row);
}


static int _read_demand_file(FILE *fp, const char *separator,
                             float **table, int *demand_years,
                             map_int_t *region_map)
{

    size_t buflen = 4000;
    char buf[buflen];
    // read in the first row (header)
    if (G_getl2(buf, buflen, fp) == 0)
        G_fatal_error(_("Demand file"
                        " contains less than one line"));

    char **tokens;
    int ntokens;

    const char *td = "\"";

    tokens = G_tokenize2(buf, separator, td);
    ntokens = G_number_of_tokens(tokens);
    if (ntokens == 0)
        G_fatal_error("No columns in the header row");

    struct ilist *ids = G_new_ilist();
    int count;
    // skip first column which does not contain id of the region
    int i;
    for (i = 1; i < ntokens; i++) {
        G_chop(tokens[i]);
        G_ilist_add(ids, atoi(tokens[i]));
    }

    int years = 0;
    while(G_getl2(buf, buflen, fp)) {
        if (buf[0] == '\0')
            continue;

        tokens = G_tokenize2(buf, separator, td);
        int ntokens2 = G_number_of_tokens(tokens);
        if (ntokens2 != ntokens)
            G_fatal_error(_("Demand: wrong number of columns in line: %s"), buf);
        if (ntokens - 1 < map_nitems(region_map))
            G_fatal_error(_("Demand: some subregions are missing"));
        count = 0;
        int i;
        demand_years[years] = atoi(tokens[0]);
        for (i = 1; i < ntokens; i++) {
            // skip first column which is the year which we ignore
            int *idx = map_get_int(region_map, ids->value[count]);
            if (idx) {
                G_chop(tokens[i]);
                table[*idx][years] = atof(tokens[i]);
            }
            count++;
        }
        // each line is a year
        years++;
    }
    G_free_ilist(ids);
    G_free_tokens(tokens);
    return years;
}

void read_demand_file(struct Demand *demandInfo, map_int_t *region_map)
{
    FILE *fp_cell, *fp_population;
    if ((fp_cell = fopen(demandInfo->cells_filename, "r")) == NULL)
        G_fatal_error(_("Cannot open area demand file <%s>"),
                      demandInfo->cells_filename);
    int countlines = 0;
    // Extract characters from file and store in character c
    for (char c = getc(fp_cell); c != EOF; c = getc(fp_cell))
        if (c == '\n') // Increment count if this character is newline
            countlines++;
    rewind(fp_cell);

    if (demandInfo->has_population) {
        if ((fp_population = fopen(demandInfo->population_filename, "r")) == NULL)
            G_fatal_error(_("Cannot open population demand file <%s>"),
                          demandInfo->population_filename);
        int countlines2 = 0;
        for (char c = getc(fp_population); c != EOF; c = getc(fp_population))
            if (c == '\n') // Increment count if this character is newline
                countlines2++;
        if (countlines != countlines2) {
            G_fatal_error(_("Area and population demand files (<%s> and <%s>) "
                            "have different number of lines"),
                          demandInfo->cells_filename, demandInfo->population_filename);
        }
        rewind(fp_population);
    }

    demandInfo->years = (int *) G_malloc(countlines * sizeof(int));
    demandInfo->cells_table = (float **) G_malloc(map_nitems(region_map) * sizeof(float *));
    for (int i = 0; i < map_nitems(region_map); i++) {
        demandInfo->cells_table[i] = (float *) G_malloc(countlines * sizeof(float));
    }
    int num_years = _read_demand_file(fp_cell, demandInfo->separator,
                                      demandInfo->cells_table,
                                      demandInfo->years, region_map);
    demandInfo->max_subregions = map_nitems(region_map);
    demandInfo->max_steps = num_years;
    G_verbose_message("Number of steps in area demand file: %d", num_years);
    fclose(fp_cell);

    if (demandInfo->has_population) {
        int *years2 = (int *) G_malloc(countlines * sizeof(int));
        demandInfo->population_table = (float **) G_malloc(map_nitems(region_map) * sizeof(float *));
        for (int i = 0; i < map_nitems(region_map); i++) {
            demandInfo->population_table[i] = (float *) G_malloc(countlines * sizeof(float));
        }
        int num_years2 = _read_demand_file(fp_population, demandInfo->separator,
                                          demandInfo->population_table,
                                          years2, region_map);
        // check files for consistency
        if (num_years != num_years2)
            G_fatal_error(_("Area and population demand files (<%s> and <%s>) "
                            "have different number of years"),
                          demandInfo->cells_filename, demandInfo->population_filename);
        for (int i = 0; i < num_years; i++) {
            if (demandInfo->years[i] != years2[i])
                G_fatal_error(_("Area and population demand files (<%s> and <%s>) "
                                "have different years"),
                              demandInfo->cells_filename, demandInfo->population_filename);
        }
        fclose(fp_population);
        G_free(years2);
    }
}

void read_potential_file(struct Potential *potentialInfo, map_int_t *region_map,
                         map_int_t *predictor_map)
{
    FILE *fp;
    if ((fp = fopen(potentialInfo->filename, "r")) == NULL)
        G_fatal_error(_("Cannot open development potential parameters file <%s>"),
                      potentialInfo->filename);

    const char *td = "\"";
    char **tokens;
    char **header_tokens;
    int header_ntokens;
    int ntokens;
    int i;
    int *pred_idx;

    size_t buflen = 4000;
    char buf[buflen];
    if (G_getl2(buf, buflen, fp) == 0)
        G_fatal_error(_("Development potential parameters file <%s>"
                        " contains less than one line"), potentialInfo->filename);
    header_tokens = G_tokenize2(buf, potentialInfo->separator, td);
    header_ntokens = G_number_of_tokens(header_tokens);
    /* num predictors minus region id, intercept and devperssure */
    int num_predictors = header_ntokens - 3;
    if (num_predictors < 0)
        G_fatal_error(_("Incorrect header in development potential file <%s>"),
                      potentialInfo->filename);
    potentialInfo->max_predictors = num_predictors;
    potentialInfo->intercept = (double *) G_malloc(map_nitems(region_map) * sizeof(double));
    potentialInfo->devpressure = (double *) G_malloc(map_nitems(region_map) * sizeof(double));
    potentialInfo->predictors = (double **) G_malloc(num_predictors * sizeof(double *));
    potentialInfo->predictor_indices = (int *) G_malloc(num_predictors * sizeof(int));
    for (i = 0; i < num_predictors; i++) {
        potentialInfo->predictors[i] = (double *) G_malloc(map_nitems(region_map) * sizeof(double));
    }
    /* index of used predictors in columns within list of predictors */
    for (i = 0; i < num_predictors; i++) {
        pred_idx = map_get(predictor_map, header_tokens[3 + i]);
        if (pred_idx)
          potentialInfo->predictor_indices[i] = *pred_idx;
        else
            G_fatal_error(_("Specified predictor <%s> in development potential file <%s>"
                            " was not provided."), header_tokens[3 + i], potentialInfo->filename);
    }

    while (G_getl2(buf, buflen, fp)) {
        if (buf[0] == '\0')
            continue;
        tokens = G_tokenize2(buf, potentialInfo->separator, td);
        ntokens = G_number_of_tokens(tokens);
        if (ntokens == 0)
            continue;
        // id + intercept + devpressure + predictores
        if (ntokens != num_predictors + 3)
            G_fatal_error(_("Potential: wrong number of columns: %s"), buf);

        int *idx;
        int id;
        double coef_intercept, coef_devpressure;
        double val;
        int j;

        G_chop(tokens[0]);
        id = atoi(tokens[0]);
        idx = map_get_int(region_map, id);
        if (idx) {
            G_chop(tokens[1]);
            coef_intercept = atof(tokens[1]);
            G_chop(tokens[2]);
            coef_devpressure = atof(tokens[2]);
            potentialInfo->intercept[*idx] = coef_intercept;
            potentialInfo->devpressure[*idx] = coef_devpressure;
            for (j = 0; j < num_predictors; j++) {
                G_chop(tokens[j + 3]);
                val = atof(tokens[j + 3]);
                potentialInfo->predictors[j][*idx] = val;
            }
        }
        // else ignoring the line with region which is not used

        G_free_tokens(tokens);
    }

    fclose(fp);
}

void read_patch_sizes(struct PatchSizes *patch_sizes,
                      map_int_t *region_map,
                      double discount_factor)
{
    FILE *fp;
    size_t buflen = 4000;
    char buf[buflen];
    int patch;
    char** tokens;
    char** header_tokens;
    int ntokens;
    int i, j;
    int region_id;
    int *region_pid;
    const char *td = "\"";
    int num_regions;
    bool found;
    bool use_header;
    int n_max_patches;
    const char *key;
    map_iter_t iter;

    n_max_patches = 0;
    patch_sizes->max_patch_size = 0;
    fp = fopen(patch_sizes->filename, "rb");
    if (fp) {
        /* just scan the file twice */
        // scan in the header line
        if (G_getl2(buf, buflen, fp) == 0)
            G_fatal_error(_("Patch library file <%s>"
                            " contains less than one line"), patch_sizes->filename);

        header_tokens = G_tokenize2(buf, ",", td);
        num_regions = G_number_of_tokens(header_tokens);
        use_header = true;
        patch_sizes->single_column = false;
        if (num_regions == 1) {
            use_header = false;
            patch_sizes->single_column = true;
            G_verbose_message(_("Only single column detected in patch library file <%s>."
                                " It will be used for all subregions."), patch_sizes->filename);
        }
        /* Check there are enough columns for subregions in map */
        if (num_regions != 1 && num_regions < map_nitems(region_map))
            G_fatal_error(_("Patch library file <%s>"
                            " has only %d columns but there are %d subregions"), patch_sizes->filename,
                          num_regions, map_nitems(region_map));
        /* Check all subregions in map have column in the file. */
        if (use_header) {
            iter = map_iter(region_map);
            while ((key = map_next(region_map, &iter))) {
                found = false;
                for (j = 0; j < num_regions; j++)
                    if (strcmp(key, header_tokens[j]) == 0)
                        found = true;
                if (!found)
                    G_fatal_error(_("Subregion id <%s> not found in header of patch file <%s>"),
                                  key, patch_sizes->filename);
            }
        }
        // initialize patch_info->patch_count to all zero
        patch_sizes->patch_count = (int*) G_calloc(num_regions, sizeof(int));
        /* add one for the header reading above */
        if (!use_header)
            n_max_patches++;
        // take one line
        while (G_getl2(buf, buflen, fp)) {
            // process each column in row
            tokens = G_tokenize2(buf, ",", td);
            ntokens = G_number_of_tokens(tokens);
            if (ntokens != num_regions)
                G_fatal_error(_("Patch library file <%s>"
                                " has inconsistent number of columns"), patch_sizes->filename);
            n_max_patches++;
        }
        // in a 2D array
        patch_sizes->patch_sizes = (int **) G_malloc(sizeof(int * ) * num_regions);
        // malloc appropriate size for each area
        for(i = 0; i < num_regions; i++) {
            patch_sizes->patch_sizes[i] =
                    (int *) G_malloc(n_max_patches * sizeof(int));
        }
        /* read first line to skip header */
        rewind(fp);
        if (use_header)
            G_getl2(buf, buflen, fp);

        while (G_getl2(buf, buflen, fp)) {
            tokens = G_tokenize2(buf, ",", td);
            ntokens = G_number_of_tokens(tokens);
            for (i = 0; i < ntokens; i++) {
                if (strcmp(tokens[i], "") != 0 ) {
                    patch = atoi(tokens[i]) * discount_factor;
                    if (patch > 0) {
                        if (patch_sizes->max_patch_size < patch)
                            patch_sizes->max_patch_size = patch;
                        if (use_header) {
                            region_pid = map_get(region_map, header_tokens[i]);
                            if (region_pid)
                                region_id = *region_pid;
                            else
                                continue;
                        }
                        else
                            region_id = 0;
                        patch_sizes->patch_sizes[region_id][patch_sizes->patch_count[region_id]] = patch;
                        patch_sizes->patch_count[region_id]++;
                    }
                }
            }
        }
        G_free_tokens(header_tokens);
        G_free_tokens(tokens);
        fclose(fp);
    }
    /* ensure there is at least one patch in each region (of size 1) */
    for (region_id = 0; region_id < num_regions; region_id++) {
        if (patch_sizes->patch_count[region_id] == 0) {
            patch_sizes->patch_sizes[region_id][0] = 1;
            patch_sizes->patch_count[region_id]++;
        }
    }
}
/**
 * @brief Read file with Depth-Damage functions
 * Header includes inundation levels in vertical units
 * (typically meters/feet), first column is id of the region
 * and values are percentages of structural damage.
 *
 * ID,0.3,0.6,0.9
 * 101,0,15,20
 * 102,10,20,30
 * ...
 * @param ddf structure
 * @param DDF_region_map
 */
void read_DDF_file(struct DepthDamageFunctions *ddf,
                   map_int_t *DDF_region_map)
{
    FILE *fp;
    if ((fp = fopen(ddf->filename, "r")) == NULL)
        G_fatal_error(_("Cannot open file <%s> with depth-damage functions"),
                      ddf->filename);

    const char *td = "\"";
    char **tokens;
    char **header_tokens;
    int header_ntokens;
    int ntokens;
    int i;
    int num_levels;
    int *pidx;
    char *id;
    double val;
    int j;
    int nitems;
    size_t buflen = 4000;
    char buf[buflen];
    if (G_getl2(buf, buflen, fp) == 0)
        G_fatal_error(_("Depth-damage functions file <%s>"
                        " contains less than one line"), ddf->filename);
    header_tokens = G_tokenize2(buf, ddf->separator, td);
    header_ntokens = G_number_of_tokens(header_tokens);
    num_levels = header_ntokens - 1;
    if (num_levels < 0)
        G_fatal_error(_("Incorrect header in depth-damage functions file <%s>"),
                      ddf->filename);
    nitems = map_nitems(DDF_region_map);

    ddf->max_levels = num_levels;
    ddf->max_subregions = nitems;
    ddf->levels = (double *) G_malloc(num_levels * sizeof(double));
    ddf->damage = (double **) G_calloc(nitems, sizeof(double *));
    ddf->loaded = (bool *) G_malloc(nitems * sizeof(bool));
    for (i = 0; i < nitems; i++) {
        ddf->damage[i] = (double *) G_malloc(num_levels * sizeof(double));
        ddf->loaded[i] = false;
    }
    /* read inundation levels */
    for (i = 0; i < num_levels; i++) {
        ddf->levels[i] = atof(header_tokens[i + 1]);
    }

    while (G_getl2(buf, buflen, fp)) {
        if (buf[0] == '\0')
            continue;
        tokens = G_tokenize2(buf, ddf->separator, td);
        ntokens = G_number_of_tokens(tokens);
        if (ntokens == 0)
            continue;
        if (ntokens != num_levels + 1)
            G_fatal_error(_("DDF: wrong number of columns: %s"), buf);

        G_chop(tokens[0]);
        id = tokens[0];
        pidx = map_get(DDF_region_map, id);
        if (pidx) {
            for (j = 0; j < num_levels; j++) {
                G_chop(tokens[j + 1]);
                val = atof(tokens[j + 1]);
                ddf->damage[*pidx][j] = val;
            }
            ddf->loaded[*pidx] = true;
        }
        // else ignoring the line with region which is not used

        G_free_tokens(tokens);
    }
    for (i = 0; i < ddf->max_subregions; i++) {
        if (!ddf->loaded[i])
            G_fatal_error(_("DDF: not all subregions have associated DDF"));
    }

    fclose(fp);
}

/**
 * Create bounding boxes for all categories in a raster.
 * @param raster CELL map as segment
 * @param masking CELL raster map as segment containing nulls
 * @param bboxes
 */
void create_bboxes(SEGMENT *raster, SEGMENT *masking, struct BBoxes *bboxes)
{
    int rows, cols;
    int row, col;
    CELL cat;
    int *index;

    rows = Rast_window_rows();
    cols = Rast_window_cols();
    map_init(&bboxes->map);
    bboxes->max_bbox = 100;
    bboxes->n_bbox = 0;
    bboxes->bbox = (struct BBox *) G_malloc(bboxes->max_bbox * sizeof(struct BBox));
    for (row = 0; row < rows; row++)
        for (col = 0; col < cols; col++) {
            Segment_get(masking, (void *)&cat, row, col);
            if (Rast_is_null_value(&cat, CELL_TYPE))
                continue;
            Segment_get(raster, (void *)&cat, row, col);
            index = map_get_int(&bboxes->map, cat);
            if (index) {
                if (bboxes->bbox[*index].e < col)
                    bboxes->bbox[*index].e = col;
                if (bboxes->bbox[*index].w > col)
                    bboxes->bbox[*index].w = col;
                if (bboxes->bbox[*index].n > row)
                    bboxes->bbox[*index].n = row;
                if (bboxes->bbox[*index].s < row)
                    bboxes->bbox[*index].s = row;
            }
            else {
                if (bboxes->n_bbox == bboxes->max_bbox) {
                    bboxes->max_bbox *= 2;
                    bboxes->bbox =
                            (struct BBox *) G_realloc(bboxes->bbox,
                                                      bboxes->max_bbox * sizeof(struct BBox));
                }
                /* HUC idx -> bbox idx */
                map_set_int(&bboxes->map, cat, bboxes->n_bbox);
                bboxes->bbox[bboxes->n_bbox].e = col;
                bboxes->bbox[bboxes->n_bbox].w = col;
                bboxes->bbox[bboxes->n_bbox].s = row;
                bboxes->bbox[bboxes->n_bbox].n = row;
                bboxes->n_bbox++;
            }
        }
}

void update_flood_probability(const char *flood_map, struct Segments *segments,
                              map_int_t *HUC_map, map_float_t *max_flood_probability_map)
{
    int row, col;
    int fd_flood_probability;
    FCELL *flood_probability_row;
    FCELL fc;
    float *pvalue;
    float max_flood_probability;
    int HUC_index;
    map_iter_t iter;
    const char *key;

    // zero all values in HUC->max_flood map
    iter = map_iter(max_flood_probability_map);
    while ((key = map_next(max_flood_probability_map, &iter)))
        map_set(max_flood_probability_map, key, 0);

    fd_flood_probability = Rast_open_old(flood_map, "");
    flood_probability_row = Rast_allocate_buf(FCELL_TYPE);
    for (row = 0; row < Rast_window_rows(); row++) {
        Rast_get_row(fd_flood_probability, flood_probability_row, row, FCELL_TYPE);
        for (col = 0; col < Rast_window_cols(); col++) {
            if (!Rast_is_null_value(&((FCELL *) flood_probability_row)[col], FCELL_TYPE)) {
                Segment_get(&segments->HUC, (void *)&HUC_index, row, col);
                fc = ((FCELL *) flood_probability_row)[col];
                pvalue = map_get_int(max_flood_probability_map, HUC_index);
                if (pvalue) {
                    max_flood_probability = *pvalue;
                    if (fc > max_flood_probability)
                        map_set_int(max_flood_probability_map, HUC_index, fc);
                }
                else
                    map_set_int(max_flood_probability_map, HUC_index, fc);
            }
        }
        Segment_put_row(&segments->flood_probability, flood_probability_row, row);
    }
    Segment_flush(&segments->flood_probability);
    G_free(flood_probability_row);
}
