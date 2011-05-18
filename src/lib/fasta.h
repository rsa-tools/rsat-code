/***************************************************************************
 *                                                                         *
 *  fasta.h
 *  Fasta reader
 *   
 *
 *                                                                         *
 ***************************************************************************/

#ifndef __FASTA__
#define __FASTA__

#include "utils.h"

//------------------------------------------ sequence
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

//------------------------------------------ fasta reader
typedef struct
{
    char *buffer;
    int buffze_size;
    int pos;
    FILE *fp;
} fasta_reader_t;

// create a new fasta reader
fasta_reader_t *new_fasta_reader(FILE *fp);

// destroy the given fasta reader
void free_fasta_reader(fasta_reader_t *fasta_reader);

// read and return the next sequence or NULL
seq_t *fasta_reader_next(fasta_reader_t *reader);

// restart from begining of file
void fasta_reader_rewind(fasta_reader_t *reader);

#endif
