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

#include <iostream> 
#include <vector>
#include <string>
#include <cstring>

using namespace std;

#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include "utils.h"
#include "fasta.h"
#include "markov.h"
#include "sampler.h"
#include "scan.h"

int VERSION = 20140213;
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
"        info-gibbs -w motif_width [-i inputfile] [-h | --help]\n"
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
"                          that are expected to be found (incompatible with --sites)\n"
"                          DEFAULT: 1\n"
"\n"
"    --zoops               try to find 0 or 1 site per sequence\n"
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
"    --collect             try to collect the N best sites using their weight scores \n"
"                          (during the collection --dmin and --zoops are not taken into account)\n"
"\n"
"    --seedmatrix=#        start sampling form sites collected by scanning the sequences with matrix #\n"
"\n"
"    --seedmatrix_sites=#  when using seed matrix specify the number of sites for each matrix (n1,n2,n3) \n"
"\n"
"    --flanks=#            when using --seedmatrix add extra # positions around the matrix\n"
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
 * parse seedmatrix_sites option (seed)
 *
 */
int parse_seedmatrix_sites(char *arg, int *values)
{
    int count = 0;
    char *tk = NULL;
    do
    {
        ASSERT(count < 16, "too much matrices");
        tk = strtok (arg, ",");
        if (tk != NULL)
        {
            //DEBUG("tk: %s", tk);
            values[count++] = atoi(tk);
        }
        arg = NULL;
    } while (tk != NULL);
    return count;
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
    params.iter             = 2000;    // iterations
    params.n                = 0;       // a motif is composed of w sites (or words)
    params.e                = 1.0;     // expected number of motif occurrences per sequence
    params.nrun             = 4;       // run gibbs main loop n times
    params.rc               = true;    // search also on reverse strand
    params.dmin             = 0;       // minimal distance between 2 sites
    params.motifs           = 1;       // number of motifs t find
    params.m1               = 6;       // left part motifs length
    params.m2               = 0;       // right part motifs length
    params.update           = 1;
    params.start_from_sites = false;
    params.minspacing       = 0;
    params.maxspacing       = 0;
    params.temperature      = 1.0;
    params.score_type       = LLR_IC_SCORE;
    params.title            = (char *) "";
    params.collect          = false;
    params.id               = 1;
    params.flanks           = 0;
    params.nseq             = 0;
    params.shift            = true;

    int optchar                = 0;
    int  l                     = 8;
    int  w_is_set              = FALSE;
    char *strand               = NULL;    // "+-" "+"
    char *seqfile              = NULL;
    char *bgfile               = NULL;
    char *matfile              = NULL;
    int  seedmatrix_sites[16];
    int  seedmatrix_sites_count = 0;

    for (int i = 0; i < 16; i++)
    {
        seedmatrix_sites[i] = -1;
    }

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
    static const char *optString =  "l:g:G:w:i:n:t:r:hs:v:b:u:d:m:e:VLIQp:FM:Y:";
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
        { "dmin",        required_argument, NULL, 'd' },
        { "motifs",      required_argument, NULL, 'm' },
        { "rseed",       required_argument, NULL, 'z' },
        { "pseudo",      required_argument, NULL, 'p' },
        { "title",       required_argument, NULL, 'T' },
        { "seedmatrix",   required_argument, NULL, 'M' },
        { "seedmatrix_sites",  required_argument, NULL, 'Y' },
        { "flanks",      required_argument, NULL, 'K' },
        { "llr",         no_argument,       NULL, 'L' },
        { "pqllr",       no_argument,       NULL, 'X' },
        { "llric",       no_argument,       NULL, 'C' },
        { "ic",          no_argument,       NULL, 'I' },
        { "pq",          no_argument,       NULL, 'Q' },
        { "help",        no_argument,       NULL, 'h' },
        { "collect",     no_argument,       NULL, 'F' },
        { "zoops",       no_argument,       NULL, 'Z' },

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
            if (matfile != NULL)
                ERROR("-w option is incompatible with the --seedmatrix option");
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
            w_is_set = TRUE;
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

            case 'Y': 
                seedmatrix_sites_count = parse_seedmatrix_sites(optarg, seedmatrix_sites);
            break; 

            case 'e': 
            params.e = atof(optarg);
            CHECK_VALUE(params.e, 0.001, 1000, "invalid value for e (should be between 0.001 and 1000)")
            break; 

            case 'Z': 
            if (params.e >= 1.0)
                params.e = 1.0;
            params.dmin = 10000000;
            CHECK_VALUE(params.e, 0.001, 1000, "invalid value for e (should be between 0.001 and 1000)")
            break; 

            case 'r': 
            params.nrun = atoi(optarg);
            CHECK_VALUE(params.nrun, 1, 10000, "invalid number of run (should be between 1 and 10000)")
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
            if (w_is_set)
                ERROR("-w option is incompatible with the --seedmatrix option");

            matfile = (char *) strdup(optarg);
            break;

            case 'K': 
            params.flanks = atoi(optarg);
            CHECK_VALUE(params.flanks, 0, 100, "invalid number of run (should be between 0 and 100)")
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
            params.collect = true;
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

    // monitor start & end time
    time_t rawtime;
    time(&rawtime);
    struct tm * start_time;
    start_time = localtime(&rawtime);

    VERBOSE1("running info-gibbs version %d\n", VERSION);

    // read the sequences
    // if (seqfile == NULL)
    // {
    //     ERROR("must provide at least a fasta file");
    // }
    
    // read sequences
    vector<string> raw_sequences;
    if (read_fasta(raw_sequences, seqfile, params.rc) == 0)
    {
        ERROR("can not load sequences");        
    }
    
    Sequences sequences = convert_sequences(raw_sequences);
    if (params.rc == true)
        params.nseq = sequences.size() / 2;
    else
        params.nseq = sequences.size();

    // set expected number of sites
    if (params.n == 0 && params.rc == true)
        params.n = (int) (sequences.len * params.e / 2.0);
    else if (params.n == 0)
        params.n = (int) (sequences.len * params.e);

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
        //bernoulli(markov, priori);
        bernoulli(markov, compute_priori(sequences));
    }

    // read optional input sig matrices
    if (matfile != NULL)
    {
        params.shift = false;
        if (params.dmin != 0)
            ERROR("option --seedmatrix is incompatible with --zoops or --dmin=#")
        // open file
        FILE *fp = fopen(matfile, "r");
        if (fp == NULL)
            ERROR("unable to open '%s'", matfile);
        
        while (1)
        {
            Array matrix = read_matrix(fp);
            if (matrix.J == 0)
                break;
            matrix.transform2logfreq(markov);
            params.start_from_sites = true;
            if (seedmatrix_sites_count > 0)
            {
                if (params.id - 1 < 16 && seedmatrix_sites[params.id - 1] != -1)
                    params.n = seedmatrix_sites[params.id - 1];
                else
                    ERROR("invalid seedmatrix_sites provided value");
            }
            params.m1 = matrix.J + params.flanks * 2;
            params.m2 = 0;
            params.minspacing   = 0;
            params.maxspacing   = 0;
            SITES sites = matrix_scan(sequences, matrix, markov, params);
            params.starting_sites = sites;
            run_sampler(raw_sequences, sequences, markov, params);
            params.id           += 1;
        }
    }
    else
    {
        // run gibbs sampler
        run_sampler(raw_sequences, sequences, markov, params);
    }

    //DEBUG("END info-gibbs");
    //VERBOSE1("end info-gibbs version %d\n", VERSION);
    
    // time info
    char time_buffer[256];
    time(&rawtime);
    struct tm * end_time;
    end_time = localtime(&rawtime);
    if (VERBOSITY >= 1)
    {
        strftime (time_buffer, 256, "%Y_%m_%d.%H%M%S", start_time);
        printf("; Job started %s\n", time_buffer);
        strftime (time_buffer, 256, "%Y_%m_%d.%H%M%S", end_time);
        printf("; Job done    %s\n", time_buffer);
    }
    return 0;
}
