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
void count_occ(count_t *count, int l, int sp, seq_t *seq, int rc, int noov)
{
    if (noov)
        init_last_position_array(count->last_position, l);

    int index = -1;

    int i;
    for (i = 0; i < seq->size - l + 1; i++) 
    {
        index = oligo2index_full(seq, i, l, sp);
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
        count->total_count++;
    }
}

// void print_header(FILE *fp, int l, int sp, int noov, long *oov_occ, int argc, char *argv[])
// {
//     // if (VERBOSITY >= 1)
//     // {
//     // 
//     //     // print command line        
//     //     fprintf(output_fp, "; ");
//     //     int i;
//     //     for (i = 0; i < argc; i++) 
//     //     {
//     //         fprintf(output_fp, "%s ", argv[i]);            
//     //     }
//     //     fprintf(output_fp, "\n");
//     //     fprintf(output_fp,
//     //         "; oligomer length               %d\n", oligo_length);
//     // 
//     //     fprintf(output_fp, 
//     //         "; column headers\n"
//     //         ";    1    seq    oligomer sequence\n"
//     //         ";    2    id     oligomer identifier\n"
//     //         ";    3    freq   relative frequencies (occurrences per position)\n" 
//     //         ";    4    occ    occurrences\n"
//     //     );
//     // 
//     //     if (noov) 
//     //     {
//     //         fprintf(output_fp, 
//     //             ";    5    ovl_occ    overlapping occurrences\n"
//     //         );
//     //     }
//     // }
// 
//     // header
//     // if (overlapping_occ)
//     //     fprintf(output_fp, "#seq\tidentifier\tobserved_freq\tocc\tovl_occ\n");
//     // else
//     fprintf(fp, "#seq\tidentifier\tobserved_freq\tocc\n");
// }

// void print_body(FILE *fp, count_t *count)
// {
//     fprintf(fp, "#seq\tidentifier\tobserved_freq\tocc\n");
// }

count_t *new_count(int l)
{
    count_t *count = (count_t *) malloc(sizeof(count_t));
    count->size             = count_array_size(l);
    count->count_table      = new_count_array(l);
    count->last_position    = new_count_array(l);
    count->palindromic      = new_count_array(l);
    count->position_count   = 0;
    count->total_count      = 0;
    return count;
}

void free_count(count_t *count)
{
    free(count->count_table);
    free(count->last_position);
    free(count->palindromic);
    free(count);
}
