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


#include <iostream> 
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "utils.h"
#include "markov.h"
#include "matrix.h"
#include "seq.h"
#include "dist.h"
#include "pval.h"

typedef struct
{
    double Wthreshold;
    int    rc;
} options_t;

// scan seq with matrix
int scan_seq(FILE *fout, seq_t *seq, int s, Array &matrix, Markov &bg, values_t *values,
	     double threshold, int rc, pvalues_t *pvalues, int origin, int offset, char *matrix_name, int *scanned_pos, int first_hit, int best_hit);

#endif
