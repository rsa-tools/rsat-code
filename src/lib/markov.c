#include "markov.h"

markov_t *new_markov(int order)
{
    int size = pow(4, order);
    markov_t *self = malloc(sizeof(markov_t));
    self->order = order;
    self->S     = malloc(sizeof(double) * size);
    self->T     = malloc(sizeof(double) * (size * 4));
    int i;
    for (i = 0; i < size; i++)
        self->S[i] = 1e-100;
    for (i = 0; i < size * 4; i++)
        self->T[i] = 1e-100;
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
    case 't':
    case 'T':
        return 3;
    default:
        return -1;
    }
}

static inline
int int2char(int i)
{
    switch (i)
    {
    case 0:
        return 'a';
    case 1:
        return 'c';
    case 2:
        return 'g';
    case 3:
        return 't';
    default:
        return 'n';
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

int oligo2index_rc_char(char *seq, int pos, int l)
{
    int value = 0;
    int S = 1;
    int i;
    for (i = 0; i < l; i++) 
    {
        int v = char2int(seq[pos + i]);
        if (v == -1)
            return -1;
        value += S * (3 - v);
        S *= 4;
    }
    return value;
}

int oligo2index(char *seq, int pos, int l)
{
    int value = 0;
    int S = 1;
    int i;
    for (i = l - 1; i >= 0; i--) 
    {
        int v = seq[pos + i];
        if (v == -1)
            return -1;
        value += S * v;
        S *= 4;
    }
    return value;
}

int oligo2index_rc(char *seq, int pos, int l)
{
    int value = 0;
    int S = 1;
    int i;
    for (i = 0; i < l; i++) 
    {
        int v = seq[pos + i];
        if (v == -1)
            return -1;
        value += S * (3 - v);
        S *= 4;
    }
    return value;
}

void index2oligo(int index, int l, char *buffer)
{
    int S = 1;
    int i;
    buffer[l] = '\0';
    for (i = l - 1; i >= 0; i--) 
    {
        buffer[i] = (index / S) % 4;
        S *= 4;
    }
}

void index2oligo_char(int index, int l, char *buffer)
{
    int S = 1;
    int i;
    buffer[l] = '\0';
    for (i = l - 1; i >= 0; i--) 
    {
        buffer[i] = int2char((index / S) % 4);
        S *= 4;
    }
}

void index2oligo_rc_char(int index, int l, char *buffer)
{
    int S = 1;
    int i;
    buffer[l] = '\0';
    for (i = l - 1; i >= 0; i--) 
    {
        buffer[l - i - 1] = int2char(3 - (index / S) % 4);
        S *= 4;
    }
}

static
void skip_comments(FILE *fp)
{
    while (!feof(fp))
    {
        int mark = fgetc(fp);
        ungetc(mark, fp);
        if (mark == ';' || mark == '#')
        {
            while (!feof(fp) && mark != '\n')
                mark = fgetc(fp);
        }
        else
        {
            break;
        }
    }
}

static
void next_line(FILE *fp)
{
    for (;;)
    {
        int mark = fgetc(fp);
        if (feof(fp) || mark == '\n')
            break;
    }
}

markov_t *new_markov_uniform()
{
    markov_t *self = new_markov(0);
    self->S[0] = 0.25;
    self->S[1] = 0.25;
    self->S[2] = 0.25;
    self->S[3] = 0.25;
    return self;
}

markov_t *load_markov(char *filename)
{
    // open file
    FILE *fp = fopen(filename, "r");
    ENSURE(fp != NULL, "can not open file");
    markov_t *self = NULL;

    // read file
    for (;;)
    {
        skip_comments(fp);
        if (feof(fp))
            break;
        char id[256];
        float freq;
        int n = fscanf(fp, "%s\t%*s\t%f", id, &freq);
        ENSURE(n == 2, "invalid bg file");
        next_line(fp);
        //INFO("%s %.3f", id, freq);
        if (self == NULL)
        {
            int id_length = strlen(id) - 1;
            self = new_markov(id_length);
        }
        int prefix_index = oligo2index_char(id, 0, self->order);
        self->S[prefix_index] += freq;
        self->T[4 * prefix_index + char2int(id[self->order])] = freq;
    }
    fclose(fp);

    // compute S & T
    int size = pow(4, self->order);
    int i;
    // S
    double sum = 0;
    for (i = 0; i < size; i++)
        sum += self->S[i];
    for (i = 0; i < size; i++)
        self->S[i] = self->S[i] / sum;
    // T
    for (i = 0; i < size * 4; i += 4)
    {
        int j;
        double sum = 0.0;
        for (j = i; j < i + 4; j++)
            sum += self->T[j];
        for (j = i; j < i + 4; j++)
            self->T[j] = self->T[j] / sum;
    }

    //
    return self;
}

void print_markov(markov_t *self)
{
    int size = pow(4, self->order);

    // S
    printf("S\n");
    int i;
    for (i = 0; i < size; i++)
        printf("%.3f\n", self->S[i]);

    // T
    printf("T\n");
    for (i = 0; i < size * 4; i++)
        printf("%.3f\n", self->T[i]);
}

double markov_P(markov_t *self, char *seq, int pos, int length)
{
    // bernoulli
    if (self->order == 0)
    {
        double p = 1.0;
        int i;
        for (i = 0; i < length; i++)
        {
            int prefix = oligo2index(seq, i, 1);
            if (prefix == -1)
                return 0.0;
            p *= self->S[prefix];
        }
        return p;
    }

    // markov order >= 1
    int prefix = oligo2index(seq, pos, length - 1);
    if (prefix == -1)
        return 0.0;

    double p = self->S[prefix];
    // INFO("p=%.3f", p);
    int i;
    for (i = self->order; i < length; i++)
    {
        int suffix = seq[i];
        int prefix = oligo2index(seq, i - self->order, self->order);
        if (suffix == -1 || prefix == -1)
            return 0.0;
        p *= self->T[4 * prefix + suffix];
    }
    // INFO("p=%.3f", p);
    return p;
}
