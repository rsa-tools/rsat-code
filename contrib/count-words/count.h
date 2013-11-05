#ifndef __COUNT__
#define __COUNT__

#include "utils.h"

void count_in_file(FILE *input_fp, FILE *output_fp, int oligo_length, \
        int spacing, int add_rc, int noov, int grouprc, int argc, char *argv[], int header);

#endif
