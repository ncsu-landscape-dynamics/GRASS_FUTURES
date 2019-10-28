#ifndef FUTURES_REDISTRIBUTE_H
#define FUTURES_REDISTRIBUTE_H


#include "keyvalue.h"

struct RedistributionMatrix
{
    const char *filename;
    float **matrix;
    int dim_from;
    int dim_to;
    int max_dim_from;
    struct KeyValueIntInt *from_map;  /* key is 'from' region id */
    struct KeyValueIntInt *to_map;  /* key is 'to' region index */
    //struct KeyValueIntInt *region_map;  /* key is region id */
    //struct KeyValueIntInt *index_map;  /* key is index */
};

int redistribute(const struct RedistributionMatrix *matrix, int regionID);
void read_redistribution_matrix(struct RedistributionMatrix *matrix);
void free_redistribution_matrix(struct RedistributionMatrix *matrix);
#endif // FUTURES_REDISTRIBUTE_H
