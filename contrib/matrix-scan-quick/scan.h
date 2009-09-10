/***************************************************************************
 *                                                                         *
 *  scan.h
 *  PSSM scanning
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __SCAN__
#define __SCAN__

using namespace std;

#include <iostream> 
#include <fstream>
#include <vector>
#include <string>

#include "markov.h"
#include "matrix.h"
#include "seq.h"
#include "dist.h"

typedef struct
{
    double Wthreshold;
    int    rc;
} options_t;

// scan seq with matrix
int scan_seq(seq_t *seq, int s, Array &matrix, Markov &bg, values_t *values, double threshold, int rc);

#endif
