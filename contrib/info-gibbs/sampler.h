/***************************************************************************
 *                                                                         *
 *  sampler.h
 *  motif discovery gibbs sampler
 *   
 *
 *                                                                         *
 ***************************************************************************/

#ifndef __GIBBS__
#define __GIBBS__

#include <iostream> 
#include <vector>
#include <string>
#include <cmath>
#include <ctime> 
#include <set>

using namespace std;

#include "utils.h"
#include "fasta.h"
#include "matrix.h"
#include "markov.h"

//#include "mtrand.h"

#define ALPHABET_SIZE 4
#define ALPHABET "ACGT"

extern double PSEUDO;
extern int SEED;

struct Site 
{
    int s;
    int p;

    // motif structure: m1-space-m2
    int m1;    // length of the left part of the motif (m1---m2) (!= l for dyads)
    int space; // lenght of the inner part (!= 0 for dyads)
    int m2;    // lenght of the right part (!= 0 for dyads)
    
    double score;
    
    Site(int seq=-1, int pos=-1, int leftpart=0, int spacing=0, int rightpart=0)
    {
        s = seq;
        p = pos;
        m1 = leftpart;
        space = spacing;
        m2 = rightpart;
        score = -1E300;
    }

    bool operator==(Site o)
    {
        return (s==o.s && p==o.p);
    }
};

#define SITES vector<Site>

#define LLR_SCORE     1 // log likelihood ratio
#define IC_SCORE      2 // information content
#define PQ_SCORE      3 // P/Q
#define PQ_LLR_SCORE  4 // P/Q + max LLR on each iter
#define LLR_IC_SCORE  5 // LLR + max IC on each iter

struct Parameters 
{
    int m1;      // length of the left part of the motif (m1---m2) (!= l for dyads)
    int m2;      // lenght of the right part of the motif (!= 0 for dyads)
    int minspacing;   // min space between m1 & m2
    int maxspacing;   // max space between m1 & m2
    int iter;    // number of iterations
    int n;       // a motif is composed of n sites (or words)
    float e;     // expected number of occurrences per sequence (affect n)
    int nrun;    // repeat gibbs main loop n times
    double temperature;
    bool rc;
    int dmin;
    int motifs;
    int update;
    int score_type; // LLR_SCORE or IC_SCORE
    char *title;
    bool collect;
    bool start_from_sites;
    SITES starting_sites;
    int id;
    int flanks;
    int nseq; // number of sequences (!= sequences.size() that containts rc)
    bool shift; // try to shift motif during sampling
};

// run the sampler & print the results on stdout
void run_sampler(vector<string> &raw_sequences, Sequences &sequences, Markov &markov, Parameters &params);

#endif
