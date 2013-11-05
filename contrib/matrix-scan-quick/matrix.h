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

#include "utils.h"
#include "markov.h"

struct Array 
{
    double **data;
    int I;
    int J;
    double p;
    double pseudo;
    
    Array(int dim1=0, int dim2=0, double val=0.0, double pseudo=1.0)
    {
        I = dim1;
        J = dim2;
        data = NULL;
        pseudo = pseudo;
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
        p = a.p;
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
        p = a.p;
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

    // Array reverse_complement()
    // {
    //     Array rc = Array(I, J, 0.0);
    // 
    //     int j;
    //     for (j = 0; j < J; j++)
    //     {
    //         int i;
    //         for (i = 0; i < I; i++)
    //         {
    //             int rci;
    //             if (i == 0)
    //                 rci = 3;
    //             else if (i == 1)
    //                 rci = 2;
    //             else if (i == 2)
    //                 rci = 1;
    //             else if (i == 3)
    //                 rci = 0;
    //             rc[rci][J - j - 1] = data[i][j];
    //         }
    //     }
    //     return rc;       
    // }

    double sum(int *word)
    {
        double s = 0.0;
        for (int i = 0; i < J; i++)
            s += data[(int) word[i]][i];
        return s;
    }

    // need to call transform2logfreq before
    double logP(char *word)
    {
        p = 0.0;
        for (int i = 0; i < J; i++)
            p += data[(int) word[i]][i];
        return p;
    }

};

int read_matrix(Array &matrix, char *filename, double pseudo);

#endif
