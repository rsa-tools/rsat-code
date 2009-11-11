/***************************************************************************
 *                                                                         *
 *  info-gibbs
 *  Gibbs sampler for DNA motif discovery based on information content
 *   
 *
 *                                                                         *
 ***************************************************************************/

//  Defrance, M. and van Helden, J. info-gibbs: a motif discovery algorithm
//  that directly optimizes information content during sampling,
//  Bioinformatics. 2009;25:2715-2722.

using namespace std;

#include <iostream> 
#include <vector>
#include <string>
#include <getopt.h>
#include <cstring>

#include <string.h>

#include "utils.h"
#include "fasta.h"
#include "markov.h"
#include "sampler.h"
#include "scan.h"

int VERSION = 200911;
char *COMMAND_LINE;

/*
 *
 * USAGE & HELP
 *
 */
void usage()
{
    printf("usage: info-gibbs -w motif_width -e expected_occ_per_seq -i seq.fa\n\n");
}

void help()
{
    printf(
"\n"
"NAME\n"
"        info-gibbs\n"
"\n"
"VERSION\n"
"        ");
    printf("%d\n", VERSION);
    printf(
"\n"
"AUTHOR\n"
"        Matthieu Defrance <defrance@bigre.ulb.ac.be>\n"
"\n"
"DESCRIPTION\n"
"        Gibbs sampling algorithm for motifs discovery.\n"
"        Searches for highly conserved motifs in a set of DNA sequences.\n"
"        Convervation is based on the motif information content (Relative Entropy).\n"
"\n"
"CATEGORY\n"
"        sequences\n"
"        motif discovery\n"
"        \n"
"USAGE        \n"
"        info-gibbs -l motiflength [-i inputfile] [-h | --help]\n"
"\n"
"ARGUMENTS\n"
"  GENERAL OPTIONS\n"
"    -h, --help            show this help message and exit\n"
"\n"
"    -v #, --verbosity=#   set verbosity to level #\n"
"                              0 no verbosity\n"
"                              1 low verbosity\n"
"                              2 high verbosity\n"
"                              3 maximal verbosity + debug + trace\n"
"\n"
"    -i #, --input=#       read sequence from # (must be in FASTA format)\n"
"                          if not specified, the standard input is used\n"
"\n"
"    -w #, --width=#       set the motif width to #\n"
"                          when the option dyad is used # represents the length monad1 + monad2\n"
"                          EXAMPLE: --length=7\n"
"\n"
"    --maxspacing=#        set maximal spacing between motif monad to # (only for dyadic motif).\n"
"                          in this case the parameter -l corresponds to length of monad1 + monad2 (without spacing)\n"
"\n"
"    --minspacing=#        set minimal spacing between motif monad to # (only for dyadic motif).\n"
"                          in this case the parameter -l corresponds to length of monad1 + monad2 (without spacing)\n"
"\n"
"    -s #, --strand=#      search in foward strand + or in both strands +-\n"
"                          EXAMPLE: --strand=+\n"
"\n"
"    -n #, --iter=#        maximum number of Gibbs sampling iterations\n"
"\n"
"    --sites=#             number of motif occurrences that are expected to\n"
"                          be found (incompatible with -e)\n"
"\n"
"    -e #, --mean_sps=#    mean number of motif occurrences (sites) expected per sequence\n"
"                          that are expected to be found (incompatible with -w)\n"
"                          DEFAULT: 1\n"
"\n"
"    -m #, --motifs=#      number of motifs to extract (one by default)\n"
"\n"
"    -b #, --bgfile=#      use # predefined INCLUSive background model\n"
"                          [http://homes.esat.kuleuven.be/~thijs/help/help_motifsampler.html#background]    \n"
"                          EXAMPLE --bgfile=mybgfile\n"
"\n"
"    -d #, --dmin=#        set minimal distance between 2 motif occurrences to #\n"
"\n"
"    -t #                  set the temperature (should be in range [0.6 1.4])\n"
"                          DEFAULT: 1.0\n"
"\n"
"    -r #  --nrun=#        try to run the Gibbs sampling seach # times\n"
"\n"
"    --finalcycle          try to collect the N best sites using their weight scores \n"
"\n"
"    --sigmatrix=#         start sampling form sites collected by scanning the sequences with sig matrix # \n"
"\n"
"    --rseed=#             set random seed to #\n"
"\n"
"    --title=#             add title # to output\n"
"\n"
"    -V, --version         print version\n"
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
    Parameters params;

    params.iter     = 2000;    // iterations
    params.n        = 0;       // a motif is composed of w sites (or words)
    params.e        = 1.0;     // expected number of motif occurrences per sequence
    params.nrun     = 4;       // run gibbs main loop n times
    params.rc       = true;    // search also on reverse strand
    params.dmin     = 0;       // minimal distance between 2 sites
    params.motifs   = 1;       // number of motifs t find
    params.m1       = 6;       // left part motifs length
    params.m2       = 0;       // right part motifs length
    params.update   = 1;
    params.start_from_sites = false;
    params.minspacing   = 0;
    params.maxspacing   = 0;
    params.temperature  = 1.0;
    params.score_type = LLR_IC_SCORE;
    params.title = (char *) "";
    params.finalcycle = false;

    int optchar;
    int l;
    char *strand;    // "+-" "+"
    char *seqfile = NULL;
    char *bgfile  = NULL;
    char *matfile = NULL;

    if (argc <= 1)
    {
        usage();
        return 0;
    }

    // construct command line string
    string cmdline = "";
    for (int i = 0; i < argc; i++)
    {
        cmdline += argv[i];
        cmdline += " ";
        COMMAND_LINE = (char *) strdup(cmdline.c_str());
    }

    // options
    static const char *optString =  "l:g:G:w:i:n:t:r:hs:v:b:u:d:m:e:VLIQp:FM:";
    static const struct option longOpts[] = {
        { "input",       required_argument, NULL, 'i' },
        { "strand",      required_argument, NULL, 's' },
        { "width",        required_argument, NULL, 'w' },
        { "maxspacing",  required_argument, NULL, 'G' },
        { "minspacing",  required_argument, NULL, 'g' },
        { "sites",       required_argument, NULL, 'S' },
        { "mean_sps",    required_argument, NULL, 'e' },
        { "nrun",        required_argument, NULL, 'r' },
        { "iter",        required_argument, NULL, 'n' },
        { "temperature", required_argument, NULL, 't' },
        { "bgfile",      required_argument, NULL, 'b' },
        { "verbose",     required_argument, NULL, 'v' },
        { "update",      required_argument, NULL, 'u' },
        { "dmin",        required_argument, NULL, 'd' },
        { "motifs",      required_argument, NULL, 'm' },
        { "rseed",       required_argument, NULL, 'z' },
        { "pseudo",      required_argument, NULL, 'p' },
        { "title",       required_argument, NULL, 'T' },
        { "sigmatrix",      required_argument, NULL, 'M' },

        { "llr",         no_argument,       NULL, 'L' },
        { "pqllr",       no_argument,       NULL, 'X' },
        { "llric",       no_argument,       NULL, 'C' },
        { "ic",          no_argument,       NULL, 'I' },
        { "pq",          no_argument,       NULL, 'Q' },
        { "help",        no_argument,       NULL, 'h' },
        { "finalcycle",  no_argument,       NULL, 'F' },
        { "version",     no_argument,       NULL, 'V' },
        { NULL,          no_argument,       NULL, 0   }
    };
    int longIndex;
    
    while ((optchar = getopt_long(argc, argv, optString, longOpts, &longIndex )) != -1)
    {   
        switch (optchar) 
        { 
            case 'h':
            help();
            return 0;
            break;

            case 'V':
            printf("%d\n", VERSION);
            return 0;
            break;
            
            case 'w': 
            l = atoi(optarg);
            CHECK_VALUE(l, 1, 100, "invalid value for length (should be between 0 and 100)")
            if (params.maxspacing != 0)
            {
                params.m1 = l / 2;
                params.m2 = l - params.m1;
            }
            else
            {
                params.m1 = l;
                params.m2 = 0;
            }
            break; 

            case 'G': 
            params.maxspacing = atoi(optarg);
            l = params.m1;
            params.m1 = l / 2;
            params.m2 = l - params.m1;
            CHECK_VALUE(params.maxspacing, params.minspacing, 30, "invalid number spacing (should be between minspacing and 30)")
            break; 

            case 'g': 
            params.minspacing = atoi(optarg);
            if (params.maxspacing != 0)
                CHECK_VALUE(params.minspacing, 0, params.maxspacing, "invalid number spacing (should be between 0 and maxspacing)")
            else
                CHECK_VALUE(params.minspacing, 0, 30, "invalid number spacing (should be between 0 and 30)")
            break; 

            case 'S': 
            params.n = atoi(optarg);
            break; 

            case 'e': 
            params.e = atof(optarg);
            CHECK_VALUE(params.e, 0.001, 1000, "invalid value for e (should be between 0.001 and 1000)")
            break; 

            case 'r': 
            params.nrun = atoi(optarg);
            CHECK_VALUE(params.nrun, 1, 10000, "invalid number of run (should be between 1 and 10000)")
            break; 

            case 'u': 
            UPDATE = atoi(optarg);
            params.update = atoi(optarg);
            break; 

            case 't': 
            CHECK_VALUE(params.temperature, 0.1, 2.0, "invalid value for temperature (should be between 0.1 and 2.0)")
            params.temperature = atof(optarg);
            break; 

            case 'n':
            params.iter = atoi(optarg);
            CHECK_VALUE(params.iter, 0, 1000000, "invalid number of iterations (should be between 0 and 1.000.000)")
            break; 
            
            case 'i':
            seqfile = (char *) strdup(optarg);
            break;

            case 'M':
            matfile = (char *) strdup(optarg);
            break;

            case 'T':
            params.title = (char *) strdup(optarg);
            break;

            case 'b':
            bgfile = (char *) strdup(optarg);
            break;

            case 'm':
            params.motifs = atoi(optarg);
            CHECK_VALUE(params.motifs, 1, 10, "invalid number of motifs (should be between 1 and 10)")
            break;

            case 'v':
            VERBOSITY = atoi(optarg);
            CHECK_VALUE(VERBOSITY, 0, 3, "invalid verbosity level (should be 0, 1, 2 or 3)")
            break;

            case 'z':
            SEED = atoi(optarg);
            break;

            case 'd':
            params.dmin = atoi(optarg);
            CHECK_VALUE(params.dmin, 0, 10000000, "invalid distance")
            break;

            case 'F':
            params.finalcycle = true;
            break;

            case 'L':
            params.score_type = LLR_SCORE;
            break;

            case 'C':
            params.score_type = LLR_IC_SCORE;
            break;

            case 'X':
            params.score_type = PQ_LLR_SCORE;
            break;
            
            case 'I':
            params.score_type = IC_SCORE;
            break;

            case 'Q':
            params.score_type = PQ_SCORE;
            break;

            case 'p':
            PSEUDO = atof(optarg);
            break;

            case 's':
            strand = (char *) strdup(optarg);
                if (strcmp(strand, "+-") == 0)
                    params.rc = true;
                else
                    params.rc = false;
            break;

            default:
            usage();
            return 0;
        }
    }

    // read the sequences
    if (seqfile == NULL)
    {
        ERROR("must provide at least a fasta file");
    }
    
    vector<string> raw_sequences;
    if (read_fasta(raw_sequences, seqfile, params.rc) == 0)
    {
            ERROR("can not load sequences");        
    }
    
    Sequences sequences = convert_sequences(raw_sequences);

    // set expected number of sites
    if (params.n == 0)
    {
        if (params.rc == true)
            params.n = (int) (0.5 * sequences.len * params.e);
        else
            params.n = (int) (sequences.len * params.e);
    }

    // set bg model
    Markov markov;
    if (bgfile != NULL)
    {
        if (load_inclusive(markov, bgfile) == 0)
        {
            ERROR("can not load bg model");
        }
    }
    else
    {
        //double priori[4] = {0.25, 0.25, 0.25, 0.25};
        bernoulli(markov, compute_priori(sequences));
    }

    // read optional input sig matrix
    if (matfile != NULL)
    {
        Array matrix;
        read_matrix(matrix, matfile);
        params.start_from_sites = true;
        params.m1 = matrix.J;
        params.m2 = 0;
        params.minspacing   = 0;
        params.maxspacing   = 0;        
        SITES sites = scan(raw_sequences, sequences, matrix, markov, params.n);        
        params.starting_sites = sites;
    }

    // run gibbs sampler
    run_sampler(raw_sequences, sequences, markov, params);

    return 0;
}
