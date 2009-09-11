/***************************************************************************
 *                                                                         *
 *  quick-scan
 *  
 *   
 *
 *                                                                         *
 ***************************************************************************/

using namespace std;

#include <iostream> 
#include <vector>
#include <string>
#include <cstring>

#include <string.h>

#include "utils.h"
#include "markov.h"
#include "matrix.h"
#include "scan.h"
#include "cfasta.h"
#include "seq.h"
#include "dist.h"

int VERSION = 200909;
char *COMMAND_LINE;

/*
 *
 * USAGE & HELP
 *
 */
void usage()
{
    printf("usage: matrix-scan-quick -m matrix -i seq.fa\n\n");
}

void help()
{
    printf(
"\n"
"NAME\n"
"        matrix-scan-quick\n"
"\n"
"VERSION\n"
"        200909\n"
"\n"
"AUTHOR\n"
"        Matthieu Defrance <defrance@bigre.ulb.ac.be>\n"
"\n"
"DESCRIPTION\n"
"        Faster and limited version of matrix-scan.\n"
"        This program takes as input a matrix and a sequence file, \n"
"        and returns either the matching positions (default), or the\n"
"        full distribution of scores observed in the whole sequence\n"
"        (option -distrib).\n"
"\n"
"CATEGORY\n"
"        sequences\n"
"        pattern matching\n"
"        PSSM\n"
"        \n"
"USAGE        \n"
"        matrix-scan-quick -m matrix [-i sequences] [-bgfile background] [-h | --help]\n"
"\n"
"INPUT FORMATS\n"
"  Sequence file\n"
"    Only sequences in FASTA format are supported.\n"
"\n"
"  Matrix file\n"
"    Only the tab format is supported.\n"
"    see convert-matrix for details.\n"
"\n"
"  Background file\n"
"    Only the INCLUSive format is supported.\n"
"    see convert-background-model for details.\n"
"\n"
"OUTPUT FORMAT\n"
"  Matches (default):\n"
"    The output is a tab-delimited file, with one row per match.\n"
"\n"
"  Distribution (option -distrib):\n"
"    The output is a tab-delimited file, with one row per weight\n"
"    score.\n"
"\n"
"SCORING SCHEME\n"
"    See matrix-scan -h for details.\n"
"\n"
"ARGUMENTS\n"
"    -h, --help            show this help message and exit.\n"
"\n"
"    -i #                  read sequence from # (must be in FASTA format).\n"
"                          if not specified, the standard input is used.\n"
"\n"
"    -m #                  read the matrix # (must be in tab format).\n"
" \n"
"    -bgfile #             use # as background model (must be in INCLUSive format).\n"
"                          by default an equiprobable model is used.\n"
"\n"
"    -2str                 scan both DNA strands.\n"
"\n"
"    -1str                 scan only one DNA strand.\n"
"\n"
"    -t #                  only capture site with a score >= #.\n"
"\n"
"    -distrib              output the weight score distribution.\n"
"\n"
"    -e #                  precision parameter for the -distrib option.\n"
"                          use -e 0.1 to compute the distribution with 1 decimal.\n"
"\n"
   );
}

/*
 *
 * MAIN
 *
 */
int main(int argc, char *argv[])
{
    VERBOSITY = 0;

    char *seqfile = NULL;
    char *bgfile  = NULL;
    char *matfile = NULL;
    int distrib = 0;
    int rc = TRUE;
    double precision = 0.1;
    double theshold = -1000.0;

    // // construct command line string
    // string cmdline = "";
    // for (int i = 0; i < argc; i++)
    // {
    //     cmdline += argv[i];
    //     cmdline += " ";
    //     COMMAND_LINE = (char *) strdup(cmdline.c_str());
    // }

    int i;
    for (i = 1; i < argc; i++) 
    {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) 
        {
            help();
            exit(0);
        } else if (strcmp(argv[i], "--version") == 0) 
        {
            printf("%d\n", VERSION);
            exit(0);
        } else if (strcmp(argv[i], "-v") == 0) 
        {
            ASSERT(argc > i + 1, "-v requires a nummber (0, 1 or 2)");
            VERBOSITY = atoi(argv[++i]);
            ASSERT(VERBOSITY >= 0 && VERBOSITY <= 2, "invalid verbosity level (should be 0, 1 or 2)");
        } 
        else if (strcmp(argv[i], "-1str") == 0) 
        {
            rc = FALSE;
        } 
        else if (strcmp(argv[i], "-2str") == 0) 
        {
            rc = TRUE;
        } 
        else if (strcmp(argv[i], "-i") == 0) 
        {
            ASSERT(argc > i + 1, "-i requires a filename");
            seqfile = argv[++i];
        } 
        else if (strcmp(argv[i], "-m") == 0) 
        {
            ASSERT(argc > i + 1, "-m requires a filename");
            matfile = argv[++i];
        } 
        else if (strcmp(argv[i], "-bgfile") == 0) 
        {
            ASSERT(argc > i + 1, "-bgfile requires a filename");
            bgfile= argv[++i];
        } 
        else if (strcmp(argv[i], "-distrib") == 0) 
        {
            distrib = TRUE;
        } 
        else if (strcmp(argv[i], "-e") == 0) 
        {
            ASSERT(argc > i + 1, "-e requires a number");
            precision = atof(argv[++i]);
            ASSERT(precision >= 0.0001 && precision <= 10, "invalid precision");
        } 
        else if (strcmp(argv[i], "-t") == 0) 
        {
            ASSERT(argc > i + 1, "-t requires a number");
            theshold = atof(argv[++i]);
        }
        else
        {
            WARNING("invalid option %s", argv[i]);
        }
    }

    if (argc <= 1)
    {
        usage();
        return 0;
    }

    // set bg model
    Markov markov;
    if (bgfile != NULL)
    {   
        if (!load_inclusive(markov, bgfile))
            ERROR("can not load bg model");
    }
    else
    {
        double priori[4] = {0.25, 0.25, 0.25, 0.25};
        bernoulli(markov, priori);
    }

    if (matfile == NULL)
    {
        ERROR("You should specify at least a matrix file and a DNA sequence file");
    }
    Array matrix;
    read_matrix(matrix, matfile);
    matrix.transform2logfreq(markov);

    // values (distrib)
    values_t *values = NULL;
    if (distrib)
        values = new_values(-1000, 1000, precision);

    // sequences
    FILE *fp;
    if (seqfile == NULL)
        fp = stdin;
    else
        fp = fopen(seqfile, "r");

    if (fp == NULL)
        ERROR("unable to open '%s'", seqfile);

    fasta_reader_t *reader = new_fasta_reader(fp);

    if (!distrib)
        fprintf(stdout, "#seq_id\tft_type\tft_name\tstrand\tstart\tend\tsequence\tweight\n");
    
    // scan all sequences
    int s = 1;
    while (1)
    {
        seq_t *seq = fasta_reader_next(reader);
        if (seq == NULL)
            break;

        scan_seq(seq, s++, matrix, markov, values, theshold, rc);
        free_seq(seq);
    }
    
    if (distrib)
        values_print(values);

    //scan(raw_sequences, sequences, matrix, markov, rc);

    return 0;
}
