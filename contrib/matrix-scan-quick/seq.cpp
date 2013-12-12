#include "seq.h"

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
    {
        rc->data[seq->size - i - 1] = 3 - seq->data[i];
    }
    return rc;
}
