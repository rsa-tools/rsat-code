#include "markov.h"

markov_t *new_markov(int order)
{
    markov_t *self = malloc(sizeof(markov_t));
    self->order = order;
    int size = pow(4, order);
    self->S = malloc(sizeof(double) * size);
    self->T = malloc(sizeof(double) * (size * 4));
    int i;
    for (i = 0; i < size; i++)
        self->S[i] = 0.0;
    for (i = 0; i < size * 4; i++)
        self->T[i] = 0.0;
    return self;
}

void free_markov(markov_t *self)
{
    free(self->S);
    free(self->T);
    free(self);
}

static inline
int char2int(char c)
{
    switch (c)
    {
    case 'a':
    case 'A':
        return 0;
    case 'c':
    case 'C':
        return 1;
    case 'g':
    case 'G':
        return 2;
    break;
    case 't':
    case 'T':
        return 3;
    break;
    default:
        return -1;
    }
}

int oligo2index_char(char *seq, int pos, int l)
{
    int value = 0;
    int S = 1;
    int i;
    for (i = l - 1; i >= 0; i--) 
    {
        int v = char2int(seq[pos + i]);
        if (v == -1)
            return -1;
        value += S * v;
        S *= 4;
    }
    return value;
}

int oligo2index(int *seq, int pos, int l)
{
    int value = 0;
    int S = 1;
    int i;
    for (i = l - 1; i >= 0; i--)
    {
        if (seq[pos + i] == -1)
            return -1;
        value += S * seq[pos + i];
        S *= 4;
    }
    return value;
}

markov_t *load_markov(char *filename)
{
    FILE *fp = fopen(filename, "r");
    ENSURE(fp != NULL, "can not open file");
    // read file in oligo-analysis format
    markov_t *self = NULL;
    while (!feof(fp))
    {
        // skip comments
        int mark = getc(fp);
        ungetc(mark, fp);
        if (mark == ';' || mark == '#')
        {
            while (!feof(fp) && mark != '\n')
                mark = getc(fp);
        }
        char id[256];
        double freq;
        ENSURE(fscanf(fp, "%s\t%*s\t%d", &id, &freq) == 2, "invalid bg file");
        if (self == NULL)
        {
            int id_length = strlen(id);
            self = new_markov(id_length);
        }
        int prefix_index = oligo2index_char(id, 0, self->order);
        self->S[prefix_index] += freq;
        self->T[prefix_index + char2int(id[self->order])] = freq;
    }
    fclose(fp);
}

double markoP(markov_t *self, int *seq, int pos, int length)
{
    int prefix = oligo2index(seq, pos, length);
    double p = self->S[prefix];
    int i;
    for (i = self->order; i < length; i++)
    {
        int suffix = seq[i];
        int prefix = oligo2index(seq, i - self->order, self->order);
        p *= self->T[prefix + suffix];
    }
    return p;
}
