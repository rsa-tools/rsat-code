/***************************************************************************
 *                                                                         *
 *  matrix.h
 *
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __ARRAY__
#define __ARRAY__

#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;

#include "fasta.h"
#include "utils.h"
#include "markov.h"

struct Array 
{
    double **data;
    int I;
    int J;
    double p;
    double pseudo;
    
    Array(int dim1=0, int dim2=0, double val=0.0)
    {
        I = dim1;
        J = dim2;
        data = NULL;
        pseudo = 1.0;
        if (I > 0 && J > 0)
            alloc(I, J, val);
    }

    void alloc(int dim1, int dim2, double val=0.0)
    {
        I = dim1;
        J = dim2;
        if (I <= 0 || J <= 0)
            return;

        // alloc
        data = new double*[I];
        for (int i = 0; i < I; i++)
            data[i] = new double[J];
        
        // init
        for (int i = 0; i < I; i++){
            for (int j = 0; j < J; j++)
                data[i][j] = val;
        }
    }

    Array(const Array& a)
    {
        I = a.I;
        J = a.J;
        pseudo = a.pseudo;
        alloc(I, J);
        for (int i = 0; i < I; i++)
        {
            for (int j = 0; j < J; j++)
                data[i][j] = a.data[i][j];
        }
    }

    ~Array()
    {
        if (data != NULL)
        {
            for (int i = 0; i < I; i++)
                delete []data[i];
            delete []data;
        }
    }

    Array& operator=(Array &a)
    {
        I = a.I;
        J = a.J;
        pseudo = a.pseudo;
        alloc(I, J);
        for (int i = 0; i < I; i++)
        {
            for (int j = 0; j < J; j++)
                data[i][j] = a.data[i][j];
        }
        return a;
    }

    double *operator[](int i)
    {
        return data[i];
    }

    double **get_data()
    {
        return data;
    }
    
    string str()
    {
        stringstream buf;
        for (int i = 0; i < I; i++)
        {
            for (int j = 0; j < J; j++)
                buf << data[i][j] << " ";
            buf << endl;
        }
        return buf.str();
    }

    void transform2logfreq(Markov &markov)
    {
        for (int j = 0; j < J; j++)
        {
            double N = 0;
            for (int i = 0; i < I; i++)
                N += data[i][j];
        
            for (int i = 0; i < I; i++)
                data[i][j] = log((data[i][j] + markov.priori[i] * pseudo) / (N + pseudo));
        }
    }

    double sum(int *word)
    {
        double s = 0.0;
        for (int i = 0; i < J; i++)
            s += data[(int) word[i]][i];
        return s;
    }

    // need to call transform2logfreq before
    double logP(int *word)
    {
        p = 0.0;
        for (int i = 0; i < J; i++)
            p += data[(int) word[i]][i];
        return p;
    }

};

Array read_matrix(FILE *fp);

#endif
