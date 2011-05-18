#include "count.h"

// spaced motif
//      l
// -----------
// ***-----***
//  m  sp   m
//
// E-value = stats->ntests * stats->pv[i]
// expected occurrences = stats->N * stats->p[i]

static inline
int count_array_size(int l)
{
    int i;
    int size = 1;
    for (i = 0; i < l; i++)
        size *= 4;
    return size;
}

long *new_count_array(int l)
{
    int size = count_array_size(l);
    long *array = malloc(sizeof(long) * size);
    int i;
    for (i = 0; i < size; i++)
        array[i] = 0;
    return array;
}

void init_last_position_array(long *array, int l)
{
    int i;
    int size = count_array_size(l);
    for (i = 0; i < size; i++)
        array[i] = -l;
}

static
int oligo2index_full(seq_t *seq, int pos, int l, int sp)
{
    if (sp > 0)
    {
        int m = (l - sp) / 2;
        int S = count_array_size(m);
        return S * oligo2index(seq->data, pos, m) + \
                   oligo2index(seq->data, pos + m + sp, m);
    }
    else
    {
        return oligo2index(seq->data, pos, l);
    }
}

static
int oligo2index_rc_full(seq_t *seq, int pos, int l, int sp)
{
    if (sp > 0)
    {
        int m = (l - sp) / 2;
        int S = count_array_size(m);
        return oligo2index_rc(seq->data, pos, m) + \
               S * oligo2index_rc(seq->data, pos + m + sp, m);
    }
    else
    {
        return oligo2index_rc(seq->data, pos, l);
    }
}
void count_occ(count_t *count, seq_t *seq)
{
    int index   = -1;
    int l       = count->oligo_length;
    int sp      = count->spacing;
    int rc      = count->rc;
    int noov    = count->noov;

    // INFO("olio length=%d", l);
    // INFO("olio spacing=%d", sp);

    if (noov)
        init_last_position_array(count->last_position, l);

    int i;
    for (i = 0; i < seq->size - l + 1; i++) 
    {
        index = oligo2index_full(seq, i, l, sp);
        // INFO("index=%d", index);
        if (rc)
        {
            long index_rc = oligo2index_rc_full(seq, i, l, sp);
            count->palindromic[index] = index == index_rc;
            index = MIN(index, index_rc);
        }

        // invalid position
        if (index == -1)
            continue;

        // increment position counter
        count->position_count++;
            
        // overlapping occurrences
        if (noov)
        {
            if (count->last_position[index] + l - 1 >= i)
            {
                //overlapping_occ[index]++;
                continue;
            }
            count->last_position[index] = i;
        }

        // count
        count->count_table[index]++;
        count->occ_count++;
    }
}

count_t *new_count(int l, int sp, int rc, int noov)
{
    int indexed_length = l;
    if (sp != -1)
        indexed_length = l - sp;
    count_t *count = (count_t *) malloc(sizeof(count_t));
    count->size             = count_array_size(indexed_length);
    count->count_table      = new_count_array(indexed_length);
    count->last_position    = new_count_array(indexed_length);
    count->palindromic      = new_count_array(indexed_length);
    count->position_count   = 0;
    count->occ_count        = 0;
    count->oligo_length     = l;
    count->monomer_length   = (l - sp) / 2;
    count->spacing          = sp;
    count->rc               = rc;
    count->noov             = noov;
    return count;
}

void free_count(count_t *count)
{
    free(count->count_table);
    free(count->last_position);
    free(count->palindromic);
    free(count);
}

stats_t *new_stats(count_t *count, markov_t *bg)
{
    stats_t *stats = (stats_t *) malloc(sizeof(stats_t));
    stats->size     = count->size;
    stats->p_table  = malloc(sizeof(double) * stats->size);
    stats->pv_table = malloc(sizeof(double) * stats->size);
    stats->ntests   = 0;
    stats->N        = count->position_count;

    // computes stats
    char oligo[16];
    int i;
    for (i = 0; i < count->size; i++)
    {
        long n = count->count_table[i];
        if (n == 0)
            continue;

        index2oligo(i, count->oligo_length, oligo);
        double p = markov_P(bg, oligo, 0, count->oligo_length);
        //INFO("oligo=%s p=%G", name, p);
        // IF BG CONSTRUCTED ON 1STR
        // if (count->rc && !count->palindromic[i])
        //     p *= 2;
        //n_exp = N * p;
        stats->p_table[i]  = p;
        stats->pv_table[i] = pbinom(n, stats->N, p);
        stats->ntests += 1;
        // ev = pv * count->test_count;
        // sig = -log10(ev);
    }
    return stats;
}

void free_stats(stats_t *stats)
{
    free(stats->p_table);
    free(stats->pv_table);
    free(stats);
}
