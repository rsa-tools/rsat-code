//
//
//  word-analysis
//  
//
//  
// 

#include "utils.h"
#include "count.h"
#include "markov.h"
#include "binomial.h"


int VERSION = 20110510;

// ===========================================================================
// =                            usage & help
// ===========================================================================
void usage(char *progname)
{
    printf("usage: %s -l length [-i inputfile] [-h]\n", progname);
}

void help(char *progname)
{
    printf(
"NAME\n"
"        word-analysis\n"   
"\n"
"AUTHOR\n"
"        Matthieu Defrance\n"
"\n"
"DESCRIPTION\n"
"        Calculates oligomer frequencies in a set of sequences,\n"
"        and detects overrepresented oligomers.\n"
"\n"
"CATEGORY\n"
"        sequences\n"
"        pattern discovery\n"
"\n"
"USAGE\n"
"        word-analysis -l length [-i inputfile]\n"
"\n"
"ARGUMENTS\n"
"    INPUT OPTIONS\n"
"        --version        print version\n"
"        -v #             change verbosity level (0, 1, 2)\n"
"        -l #             set oligomer length to # (monad size when using dyads)\n"
"        -bg #            load the background model from # (oligo-analysis format)\n"
"        -2str            add reverse complement\n"
"        -1str            do not add reverse complement\n"
"        -noov            do not allow overlapping occurrences\n"
// "        -grouprc         group reverse complement with the direct sequence\n"
// "        -nogrouprc       do not group reverse complement with the direct sequence\n"
"\n"
"\n"

    );
}

// ===========================================================================
// =                            main
// ===========================================================================
// int parse_spacing(char *arg, int *values)
// {
//     int count = 0;
//     char *tk = NULL;
//     do
//     {
//         tk = strtok (arg, "-");
//         if (tk != NULL)
//         {
//             //DEBUG("tk: ", tk);
//             values[count++] = atoi(tk);
//         }
//         arg = NULL;
//     } while (count != 2 && tk != NULL);
//     return count;
// }

int main(int argc, char *argv[])
{
    // default options
    char *input_filename    = NULL;
    char *output_filename   = NULL;
    char *bg_filename       = NULL;
    int rc                  = TRUE;
    int noov                = FALSE;
    int oligo_length        = 1;
    // int spacing = -1;
    // int spacing_range[2]    = {-1, -1};
    // int grouprc = TRUE;
    // 
    // // options
    // if (argc == 1) 
    // {
    //     usage(argv[0]);
    //     exit(0);
    // }
    // 

    // parse options
    int i;
    for (i = 1; i < argc; i++) 
    {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) 
        {
            help(argv[0]);
            exit(0);
        } 
        else if (strcmp(argv[i], "--version") == 0) 
        {
            printf("%d\n", VERSION);
            exit(0);
        } 
        else if (strcmp(argv[i], "-bg") == 0) 
        {
            ENSURE(argc > i + 1, "-v requires a number (0, 1 or 2)");
            bg_filename = argv[++i];
        } 
    //     else if (strcmp(argv[i], "-v") == 0) 
    //     {
    //         ASSERT(argc > i + 1, "-v requires a nummber (0, 1 or 2)");
    //         VERBOSITY = atoi(argv[++i]);
    //         ASSERT(VERBOSITY >= 0 && VERBOSITY <= 2, "invalid verbosity level (should be 0, 1 or 2)");
    //     } 
        else if (strcmp(argv[i], "-1str") == 0) 
        {
            rc = FALSE;
        } 
        else if (strcmp(argv[i], "-2str") == 0) 
        {
            rc = TRUE;
        } 
        else if (strcmp(argv[i], "-noov") == 0) 
        {
            noov = TRUE;
        } 
        else if (strcmp(argv[i], "-l") == 0)
        {
            ASSERT(argc > i + 1, "-l requires a nummber");
            oligo_length = atoi(argv[++i]);
            ASSERT(oligo_length >= 1 && oligo_length <= 14, "invalid oligo length");
        } 
        // else if (strcmp(argv[i], "-sp") == 0)
        // {
        //     ASSERT(argc > i + 1, "-sp requires an argument");
        //     parse_spacing(argv[++i], spacing_tab);
        //     spacing = spacing_tab[0];
        //     //printf("spacing =%d %d", spacing_tab[0], spacing_tab[1]);
        // } 
        else if (strcmp(argv[i], "-i") == 0) 
        {
            ASSERT(argc > i + 1, "-i requires a string");
            input_filename = argv[++i];
        } 
        else if (strcmp(argv[i], "-o") == 0) 
        {
            ASSERT(argc > i + 1, "-o requires a string");
            output_filename = argv[++i];
        } 
        else 
        {
            FATAL_ERROR("invalid option %s", argv[i]);
        }
    }

    //markov_t *m = load_markov("../2nt_upstream-noorf_Escherichia_coli_K12-noov-2str.freq");
    markov_t *bg = NULL;
    if (bg_filename)
        bg = load_markov(bg_filename);
    else
        bg = new_markov_uniform();
    //print_markov(bg);
    // double x = markov_P(bg, "acg", 0, 3);
    // INFO("x=%G", x);
    FILE *input_fp = stdin;
    if (input_filename) 
    {
        input_fp = fopen(input_filename, "r");
        //ENSURE(input_fp != NULL, "can not read from file '%s'", input_filename);
        ENSURE(input_fp != NULL, "can not read from file");
    }
    FILE *output_fp = stdout;
    if (output_filename) 
    {
       output_fp = fopen(output_filename, "w");
       //ENSURE(output_fp != NULL, "can not write to file '%s'", output_filename);
       ENSURE(output_fp != NULL, "can not write to file");
    }

    // read fasta file & compute count table
    count_t *count = new_count(oligo_length);
    fasta_reader_t *reader = new_fasta_reader(input_fp);
    while  (TRUE)
    {
        seq_t *seq = fasta_reader_next(reader);
        if (seq == NULL)
            break;
        count_occ(count, oligo_length, 0, seq, rc, noov);
        free_seq(seq);
    }
    // 
    // fprintf(output_fp, 
    //     "; column headers\n"
    //     ";    1    seq    oligomer sequence\n"
    //     ";    2    id     oligomer identifier\n"
    //     ";    3    freq   relative frequencies (occurrences per position)\n" 
    //     ";    4    occ    occurrences\n"
    // );

    fprintf(output_fp, "#seq\tid\tocc\tocc_P\tocc_E\tocc_sig\n");
    
    // binomial stats on count table
    char name[16];
    char name_rc[16];
    char oligo[16];
    char id[128];
    for (i = 0; i < count->size; i++)
    {
        long n = count->count_table[i];
        if (n == 0)
            continue;
        index2oligo(i, oligo_length, oligo);
        index2oligo_char(i, oligo_length, name);
        long N = count->position_count;
        double p = markov_P(bg, oligo, 0, oligo_length);
        if (rc && !count->palindromic[i])
            p *= 2;
        //printf("%ld %ld %f\n", n, N, p);

        double pv = pbinom(n, N, p);
        double ev = pv * count->test_count;
        double sig = -log10(ev);
        if (rc)
        {
            index2oligo_rc_char(i, oligo_length, name_rc);
            sprintf(id, "%s|%s", name, name_rc);
        }
        else
        {
            sprintf(id, "%s", name);
        }
        fprintf(output_fp, "%s\t%s\t%ld\t%G\t%G\t%G\n", name, id, n, pv, ev, sig);
    }

    // free data
    free_fasta_reader(reader);
    free_count(count);

    // fflush(output_fp);
    if (input_filename)
        fclose(input_fp);
    if (output_filename)
        fclose(output_fp);
    return 0;
}


// count_occ(&count, oligo_length, int sp, seq_t *seq, rc, noov);
// 
// count_in_file(input_fp, output_fp, oligo_length, spacing, add_rc, noov, grouprc, argc, argv, 1);
// // if (spacing != -1)
// // {
// //     for (spacing = spacing_tab[0] + 1; spacing <= spacing_tab[1]; spacing++)
// //     {
// //         fseek(input_fp, SEEK_SET, 0);
// //         count_in_file(input_fp, output_fp, oligo_length, spacing, add_rc, noov, grouprc, argc, argv, 0);
// //     }
// // }
// 

