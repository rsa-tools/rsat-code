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

SITES scan(vector<string> raw_sequences, Sequences &sequences, Array &matrix, Markov &bg, int n);

#endif
