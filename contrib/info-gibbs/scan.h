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

#include "sampler.h"

SITES matrix_scan(Sequences &sequences, Array &matrix, Markov &bg, Parameters &params);

#endif
