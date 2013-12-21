/***************************************************************************
 *                                                                         *
 *  pval.h
 *
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __PVAL__
#define __PVAL__

#include "utils.h"
#include <math.h>

#ifdef __cplusplus 
extern "C" {
#endif

typedef struct
{
    double *data;
    int size;
    double w_min;
    double w_max;
} pvalues_t;

pvalues_t *new_pvalues();

void free_pvalues(pvalues_t *pvalues);

double score2pvalue(pvalues_t *pvalues, double score);

pvalues_t *read_distrib(char *filename);

#ifdef __cplusplus
}
#endif

#endif
