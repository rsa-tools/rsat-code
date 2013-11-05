/***************************************************************************
 *                                                                         *
 *  cfasta.h
 *  fast fasta reader
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __CFASTA__
#define __CFASTA__

#include "seq.h"
#include "utils.h"

typedef struct
{
    char *buffer;
    int bufsize;
    int pos;
    FILE *fp;
} fasta_reader_t;

fasta_reader_t *new_fasta_reader(FILE *fp);

void free_fasta_reader(fasta_reader_t *fasta_reader);

seq_t *fasta_reader_next(fasta_reader_t *reader);


#endif
