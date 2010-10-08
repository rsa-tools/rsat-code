#include "cfasta.h"

#define FASTA_READER_BUFSIZE 4096

fasta_reader_t *new_fasta_reader(FILE *fp)
{
    assert(fp != NULL);
    fasta_reader_t *fasta_reader = (fasta_reader_t *) malloc(sizeof(fasta_reader_t));
    fasta_reader->bufsize = FASTA_READER_BUFSIZE;
    fasta_reader->pos = 0;
    fasta_reader->buffer = (char *) malloc(sizeof(char) * fasta_reader->bufsize);
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
    return fgetc(reader->fp);
    // if (reader->pos >= reader->bufsize) 
    // {
    //     int read_size = fread(reader->buffer, 1, reader->bufsize, reader->fp);
    //     if (read_size != reader->bufsize)
    //         reader->buffer[read_size] = EOF;
    // 
    //     reader->pos = 0;
    // }
    // 
    // char c = reader->buffer[reader->pos];
    // if (c != EOF)
    //     reader->pos++;
    // return c;
}

seq_t *fasta_reader_next(fasta_reader_t *reader)
{
    seq_t *seq = new_seq(1000);

    // read header
    char c;
    int i = 0;
    int short_name_found = 0;
    do 
    {
        c = fasta_reader_getc(reader);
        if ((c == '\t' || c == ' ' || c == '\r') && i > 0)
            short_name_found = 1;
        if (!short_name_found && i < 1024 && \
            c!= EOF && c != '\n' && c != '\t' && c != ' ' && c != '\r' && c != '>')
            seq->name[i++] = c;
    } while (c != EOF && c != '\n');
    seq->name[i] = '\0';
    if (c == EOF)
        return NULL;

    // read sequence
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
