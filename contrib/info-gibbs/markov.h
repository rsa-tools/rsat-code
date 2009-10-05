/***************************************************************************
 *                                                                         *
 *  markov.h
 *  Markov chain
 *    
 *  
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __MARKOV__
#define __MARKOV__

using namespace std;

#include <fstream>
#include <iostream> 
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>

#define ALPHABET_SIZE 4
#define ALPHABET "ACGT"

struct Markov
{
/*
    Stationnary
    Maximum Likelihood parameters estimation with pseudo count:
    
                       C(suffix|prefix) + pseudo
    P(suffix|prefix) = -------------------------
                         C(prefix) + N pseudo 

    N: number of suffixes
    P(suffix|prefix) (example P(C|AA))
    C(suffix|prefix): prefix-suffix count (example AAC)
    C(prefix): prefix count (example AA)

*/

    int order;
    double pseudo;     // pseudo-count (used for smoothing)
    double **T;        // transition matrix
    double *S;         // stationary vector
    int msize;         // transition matrix size (number of rows)
    double *priori;    // priori vector (pA, pC, pG, pT)
    int alphabet_size;
    
    double p;
    int prefix;
    int suffix;
    int idx;
    
    Markov(int o=-1, double p=1.0)
    {
        order = o;
        alphabet_size = ALPHABET_SIZE;
        T = NULL;
        S = NULL;
        priori = NULL;
        if (order >= 0)
            alloc(o, p);
    }
    
    void alloc(int o, double p=1.0)
    {
        pseudo = 1.0;
        order = o;
        msize = (int) pow((double) ALPHABET_SIZE, order);
        priori = new double[ALPHABET_SIZE];
        S = new double[msize];
        T = new double*[msize];
        for (int i = 0; i < msize; i++)
            T[i] = new double[ALPHABET_SIZE];

        // init values
        for (int i = 0; i < msize; i++)
        {
            S[i] = 1.0;
            for (int j = 0; j < ALPHABET_SIZE; j++)
                T[i][j] = 0.0;
        }
    }

    void dealloc()
    {
        if (order != -1)
        {
            delete []priori;
             for (int i = 0; i < msize; i++)
                delete []T[i];
            delete []T;
            delete []S;     
        }
    }

    ~Markov() 
    {
        dealloc();
    }

    Markov(const Markov& m)
    {
        order = m.order;
        pseudo = m.pseudo;
        msize = m.msize;
        alloc(order, pseudo);

        //copy values
        for (int i = 0; i < msize; i++)
        {
            S[i] = m.S[i];
            for (int j = 0; j < ALPHABET_SIZE; j++)
                T[i][j] = m.T[i][j];
        }
        for (int j = 0; j < ALPHABET_SIZE; j++)
            priori[j] = m.priori[j];
    }

    Markov& operator=(Markov &m)
    {
        dealloc();
        order = m.order;
        pseudo = m.pseudo;
        msize = m.msize;
        alloc(order, pseudo);
    
        //copy values
        for (int i = 0; i < msize; i++)
        {
            S[i] = m.S[i];
            for (int j = 0; j < ALPHABET_SIZE; j++)
                T[i][j] = m.T[i][j];
        }
        for (int j = 0; j < ALPHABET_SIZE; j++)
            priori[j] = m.priori[j];
        return m;
    }

    int word2index(int *word, int len)
    {
        idx = 0;
        int P = 1;
        for (int i = 0; i < len; i++)
        {
            idx += word[len-i-1] * P;
            P = P * ALPHABET_SIZE;
        }
        return idx;
    }

    void count(int *s, int len)
    {
        int prefix;
        int suffix;

        for (int i = 0; i < len; i++)
        {
            suffix = s[i+order]; 
            prefix = word2index(&s[i], order);

            T[prefix][suffix] += 1;
            S[prefix] += 1;
        }
    }

    double logP(int *word, int len)
    {
        if (order == 0)
        {
            p = 0.0;
            for (int i = 0; i < len; i++)
                p += log(priori[word[i]]);
        }
        else
        {
            prefix = word2index(word, order);
            p = log(S[prefix]);
            for (int i = order; i < len; i++)
            {
                suffix = word[i]; 
                prefix = word2index(&word[i-order], order);
                p += log(T[prefix][suffix]);
            }
        }
        return p;
    }

    double P(int *word, int len)
    {
        if (order == 0)
        {
            p = 1.0;
            for (int i = 0; i < len; i++)
                p *= priori[word[i]];
        }
        else
        {
            prefix = word2index(word, order);
            p = S[prefix];
            for (int i = order; i < len; i++)
            {
                suffix = word[i]; 
                prefix = word2index(&word[i-order], order);
                p *= T[prefix][suffix];
            }
        }
        return p;
    }

};

// load the model from file in MotifSampler format
// http://homes.esat.kuleuven.be/~thijs/help/help_motifsampler.html
int load_inclusive(Markov &m, string filename);

// construct a Bernoulli model from priori probs (P(A), P(C), P(G), P(T))
int bernoulli(Markov &m, double *priori);

#endif
