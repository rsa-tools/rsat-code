/***************************************************************************
 *                                                                         *
 *  markov.h
 *  DNA Markov model
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

// create a new markov uniform model
markov_t *new_markov_uniform();

// return the probability of seq[pos..pos+length-1] in self
double markov_P(markov_t *self, char *seq, int pos, int length);

// print to stdout markov model
void print_markov(markov_t *self);

// oligo2index
int oligo2index(char *seq, int pos, int l);

// oligo2index reverse complement
int oligo2index_rc(char *seq, int pos, int l);

// index to oligo
void index2oligo(int index, int l, char *buffer);

// convert back index to oligo
void index2oligo_char(int index, int l, char *buffer);

// convert back index to oligo reverse complement
void index2oligo_rc_char(int index, int l, char *buffer);

#endif
