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

#include <math.h>
#include "utils.h"

// Compute the binomial P-value P(x >= n)
// n -- number of success
// N -- number of trials
// p -- prob of success
// binomial recursive formula:
//            P[x] * p(N-x)
//   P[x+1] = -------------
//               q(x+1)
double pbinom(int n, int N, double p);

#endif
