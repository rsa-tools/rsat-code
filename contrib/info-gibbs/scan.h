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

#include "sampler.h"

SITES matrix_scan(Sequences &sequences, Array &matrix, Markov &bg, Parameters &params);

#endif
