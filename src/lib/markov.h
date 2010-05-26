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

// create a new markov model
markov_t *new_markov(int order);

// free the markov model
void free_markov(markov_t *markov);

// create a new markov model by loading oligo analysis formatted file
markov_t *load_markov(char *filename);

// return the probability of seq[pos..pos+length-1] in self
double markov_P(markov_t *self, char *seq, int pos, int length);

// print to stdout markov model
void print_markov(markov_t *self);

#endif
