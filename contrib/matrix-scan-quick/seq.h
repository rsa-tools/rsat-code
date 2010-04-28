/***************************************************************************
 *                                                                         *
 *  seq.h
 *  simple DNA seq manipulation
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __SEQ__
#define __SEQ__

#include "utils.h"

typedef struct
{
    char *data;
    int size;
    int msize;
    char *name; // limited to 1024 chars
} seq_t;


// create a new empty seq
seq_t *new_seq(int size);

// free the given seq
void free_seq(seq_t *seq);

// create a new reverse complement of seq
seq_t *new_seq_rc(seq_t *seq);

static inline
void seq_append_c(seq_t *seq, char c)
{
    if (seq->size >= seq->msize)
    {
        seq->msize += 8;
        seq->data = (char *) realloc(seq->data, sizeof(char) * seq->msize);
    }
    switch (c)
    {
        case 'a':
        case 'A':
            seq->data[seq->size++] = 0;
            break;

        case 'c':
        case 'C':
            seq->data[seq->size++] = 1;
            break;

        case 'g':
        case 'G':
            seq->data[seq->size++] = 2;
            break;

        case 't':
        case 'T':
            seq->data[seq->size++] = 3;
            break;

        case 'n':
        case 'N':
            seq->data[seq->size++] = -1;
            break;
    }
}

#endif
