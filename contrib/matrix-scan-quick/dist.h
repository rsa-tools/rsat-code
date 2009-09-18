/***************************************************************************
 *                                                                         *
 *  dist.h
 *  
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __DIST__
#define __DIST__

#include "utils.h"

typedef struct
{
    int *data;
    int size;
    double e;
    double min;
    double max;
} values_t;

values_t *new_values(double min, double max, double e);

void free_values(values_t *values);

void values_add(values_t *values, double value);

void values_print(FILE *fout, values_t *values);

#endif
