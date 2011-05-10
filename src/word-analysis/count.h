/***************************************************************************
 *                                                                         *
 *  count.h
 *  
 *   
 *
 *                                                                         *
 ***************************************************************************/
#ifndef __COUNT__
#define __COUNT__

#include "utils.h"
#include "fasta.h"
#include "markov.h"

typedef struct count_s
{
    long *count_table;
    long *last_position;
    int position_count;
    int total_count;
} count_t;


// new count data
count_t *new_count(int l);

// free count data
void free_count(count_t *count);

// count array size
int count_size(int l);

// count oligo occurrences in given sequence
// last_position: used by noov to store last occurrence 
// l: total oligo length
// sp: spacing (> 0 for dyads)
// rc: use the 2 strands of seq
// noov: discard overlapping occurrences
void count_occ(count_t *count, int l, int sp, seq_t *seq, int rc, int noov);

#endif
