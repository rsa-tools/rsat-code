#include "fasta.h"

#define FASTA_READER_BUFFER_SIZE 4096

//------------------------------------------ sequence
seq_t *new_seq(int size)
{
    seq_t *seq = (seq_t *) malloc(sizeof(seq_t));
    seq->size = 0;
    seq->msize = size;
    seq->data = (char *) malloc(sizeof(char) * seq->msize);
    seq->name = (char *) malloc(sizeof(char) * 1024);
    return seq;
}

void free_seq(seq_t *seq)
{
    free(seq->data);
    free(seq->name);
    free(seq);
}

seq_t *new_seq_rc(seq_t *seq)
{
    seq_t *rc = new_seq(seq->size);
    int i;
    for (i = 0; i < seq->size; i++)
        rc->data[seq->size - i - 1] = 3 - seq->data[i];
    return rc;
}

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
            seq->data[seq->size++] = 'a';
            break;

        case 'c':
        case 'C':
            seq->data[seq->size++] = 'c';
            break;

        case 'g':
        case 'G':
            seq->data[seq->size++] = 'g';
            break;

        case 't':
        case 'T':
            seq->data[seq->size++] = 't';
            break;

        case 'n':
        case 'N':
            seq->data[seq->size++] = 'n';
            break;
    }
}

//------------------------------------------ fasta
fasta_reader_t *new_fasta_reader(FILE *fp)
{
    ASSERT(fp != NULL, "invalid file");
    fasta_reader_t *fasta_reader = (fasta_reader_t *) malloc(sizeof(fasta_reader_t));
    fasta_reader->buffze_size = FASTA_READER_BUFFER_SIZE;
    fasta_reader->pos = 0;
    fasta_reader->buffer = (char *) malloc(sizeof(char) * fasta_reader->buffze_size);
    fasta_reader->fp = fp;
    return fasta_reader;
}

void free_fasta_reader(fasta_reader_t *fasta_reader)
{
    free(fasta_reader->buffer);
    free(fasta_reader);
}

static inline
char fasta_reader_getc(fasta_reader_t *reader)
{
    if (reader->pos >= reader->buffze_size) 
    {
        int read_size = fread(reader->buffer, 1, reader->buffze_size, reader->fp);
        if (read_size != reader->buffze_size)
            reader->buffer[read_size] = EOF;
        reader->pos = 0;
    }
    
    char c = reader->buffer[reader->pos];
    if (c != EOF)
        reader->pos++;
    return c;
}

seq_t *fasta_reader_next(fasta_reader_t *reader)
{
    seq_t *seq = new_seq(1000);

    // header
    char c;
    int i = 0;
    do 
    {
        c = fasta_reader_getc(reader);
        if (i < 1024 && c!= EOF && c != '\n')
            seq->name[i++] = c;
    } while (c != EOF && c != '\n');
    seq->name[i] = '\0';
    
    if (c == EOF)
        return NULL;

    // sequence
    do 
    {
        c = fasta_reader_getc(reader);
        if (c == EOF)
          break;
        if (c != '\n' && c != '>') 
            seq_append_c(seq, c);
    } while (c != '>');

    return seq;
}
