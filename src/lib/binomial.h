/***************************************************************************
 *                                                                         *
 *  binomial.h
 *  binomial statistics
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __BINOMIAL__
#define __BINOMIAL__

#include "utils.h"

// P[x] = P[x-1] * p (N-x+1) / (q x)
double pbinom(int n, int N, double p);


#endif
