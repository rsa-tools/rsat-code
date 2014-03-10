/***************************************************************************
 *                                                                         *
 *  fasta.h
 *  simple fasta file reader
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __FASTA__
#define __FASTA__

#include <iostream> 
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#define ALPHABET_SIZE 4
#define ALPHABET "ACGT"

struct Sequence
{
    int *data;
    int len;

    Sequence(int l=0)
    {
        len = l;
        if (len > 0)
            data = new int[len];
        else
            data = NULL;
    }
    
    ~Sequence()
    {
        if (len > 0)
            delete []data;
    }

    Sequence &operator=(Sequence &s)
    {
        len = s.len;
        data = new int[len];
        for (int i = 0; i < len; i++)
        {
            data[i] = s.data[i];
        }
        return s;
    }

    int size()
    {
        return len;
    }

    int  &operator[](int i)
    {
        return data[i];
    }
};


struct Sequences
{
    int len;
    Sequence *data;

    Sequences(int l=0)
    {
        len = l;
        if (len > 0)
            data = new Sequence[len];
        else
            data = NULL;
    }

    ~Sequences()
    {
        if (len > 0)
            delete []data;
    }

    Sequences &operator=(Sequences &s)
    {
        len = s.len;
        data = new Sequence[len];
        for (int i = 0; i < len; i++){
            data[i] = s.data[i];
        }
        return s;
    }

    int size()
    {
        return len;
    }
    
    int total_size()
    {
        int total = 0;
        for (int i = 0; i < len; i++)
            total += data[i].size();
        return total;
    }

    Sequence &operator[](int i)
    {
        return data[i];
    }
};

// read the seuqneces in sequences (an empty vector)
int read_fasta(vector<string> &sequences, char *filename, bool rc);

// convert sequence to vector of int
Sequences convert_sequences(vector<string> &raw_sequences);

// compute the priori probabiliies (p(A), p(C), ...)
double *compute_priori(Sequences &sequences);


#endif
