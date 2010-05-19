/***************************************************************************
 *                                                                         *
 *  markov.h
 *  Markov chains
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __MARKOV__
#define __MARKOV__

#include <math.h>

#include "utils.h"

typedef struct markov_s
{
    int order;
    double *S;
    double *T;
} markov_t;

markov_t *new_markov(int order);

void free_markov(markov_t *markov);

markov_t *load_markov(char *filename);

#endif
