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
#include "binomial.h"

// count table
typedef struct count_s
{
    long *count_table;   // occurence count table
    long *last_position; // last occurrence position table (used internally)
    long *palindromic;   // palyndromic status table
    int size;            // table size
    int position_count;  // total number of scanned positions
    int occ_count;       // total number of occurrences
    int oligo_length;
    int monomer_length;
    int spacing;
    int rc;
    int noov;
} count_t;

// new count data
count_t *new_count(int l, int sp, int rc, int noov);

// free count data
void free_count(count_t *count);

// count oligo occurrences in given sequence
// last_position: used by noov to store last occurrence 
// l: total oligo length
// sp: spacing (> 0 for dyads)
// rc: use the 2 strands of seq
// noov: discard overlapping occurrences
void count_occ(count_t *count, seq_t *seq);

// stats table
typedef struct stats_s
{
    double *p_table;    // expected frequency table
    double *pv_table;   // Binomial P-value table
    long size;          // table size
    long ntests;        // number of oligo tested
    long N;             // total number of scanned positions
} stats_t;

// new stats data
stats_t *new_stats(count_t *count, markov_t *bg);

// free stats data
void free_stats(stats_t *stats);

#endif
