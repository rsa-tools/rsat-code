#include "count.h"

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

// static
// int oligo2index(seq_t *seq, int pos, int l)
// {
//     int value = 0;
//     int S = 1;
//     int i;
//     for (i = l - 1; i >= 0; i--) 
//     {
//         if (seq->data[pos + i] == -1)
//             return -1;
//         value += S * seq->data[pos + i];
//         S *= 4;
//     }
//     return value;
// }
// 
// static
// int oligo2index_rc(seq_t *seq, int pos, int l)
// {
//     int value = 0;
//     int S = 1;
//     int i;
//     for (i = 0; i < l; i++) 
//     {
//         if (seq->data[pos + i] == -1)
//             return -1;
//         value += S * (3 - seq->data[pos + i]);
//         S *= 4;
//     }
//     return value;
// }

static
int oligo2index(seq_t *seq, int pos, int l, int sp)
{
    if (sp > 0)
    {
        int m = (l - sp) / 2;
        int S = count_array_size(m);
        return S * oligo2index_char(seq->data, pos, m) + \
                   oligo2index_char(seq->data, pos + m + sp, m);
    }
    else
    {
        return oligo2index_char(seq->data, pos, l);
    }
}

static
int oligo2index_rc(seq_t *seq, int pos, int l, int sp)
{
    if (sp > 0)
    {
        int m = (l - sp) / 2;
        int S = count_array_size(m);
        return oligo2index_rc_char(seq->data, pos, m) + \
               S * oligo2index_rc_char(seq->data, pos + m + sp, m);
    }
    else
    {
        return oligo2index_rc_char(seq->data, pos, l);
    }
}
void count_occ(long *count_table, long *last_position, int l, int sp, seq_t *seq, int rc, int noov)
{
    if (noov)
        init_last_position_array(last_position, l);

    int index = -1;
    int position_count = 0;
    int total_count = 0;

    int i;
    for (i = 0; i < seq->size - l + 1; i++) 
    {
        index = oligo2index(seq, i, l, sp);
        if (rc)
            index = MIN(index, oligo2index_rc(seq, i, l, sp));

        // invalid position
        if (index == -1)
            continue;

        // increment position counter
        position_count++;
            
        // overlapping occurrences
        if (noov)
        {
            if (last_position[index] + l - 1 >= i)
            {
                //overlapping_occ[index]++;
                continue;
            }
            last_position[index] = i;
        }

        // count
        count_table[index]++;
        total_count++;
    }
}
