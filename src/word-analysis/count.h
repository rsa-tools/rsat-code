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


// count oligo occurrences in given sequence
// last_position: used by noov to store last occurrence 
// l: total oligo length
// sp: spacing (> 0 for dyads)
// rc: use the 2 strands of seq
// noov: discard overlapping occurrences
void count_occ(long *count_table, long *last_position, int l, int sp, seq_t *seq, int rc, int noov);

#endif
