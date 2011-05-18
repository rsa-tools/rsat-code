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
#include <math.h>

int VERSION = 20110518;

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
"AUTHORS\n"
"        Jacques van Helden\n"
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
"        --version        print version information.\n"
"        -v #             change verbosity level (0, 1, 2).\n"
"        -l #             set oligomer length to # (monad size when using dyads).\n"
"        -expfreq #       load the background model from # (oligo-analysis format).\n"
"        -2str            oligonucleotide occurrences found on both stands are summed.\n"
"        -1str            inactivates the summation of occurrences on both strands.\n"
"        -noov            do not allow overlapping occurrences.\n"
"        -count           only repport oligo count.\n"

// "        -grouprc         group reverse complement with the direct sequence\n"
// "        -nogrouprc       do not group reverse complement with the direct sequence\n"
"\n"
"\n"

    );
}

// ===========================================================================
// =                            main
// ===========================================================================
int parse_spacing(char *arg, int *values)
{
    int count = 0;
    char *tk = NULL;
    do
    {
        tk = strtok(arg, "-");
        if (tk != NULL)
        {
            //DEBUG("tk: ", tk);
            values[count++] = atoi(tk);
        }
        arg = NULL;
    } while (count != 2 && tk != NULL);
    return count;
}

count_t *count_in_file(FILE *fp, int oligo_length, int spacing, int rc, int noov)
{
    // read fasta file & compute count table
    count_t *count = new_count(oligo_length, spacing, rc, noov);
    fasta_reader_t *reader = new_fasta_reader(fp);
    while  (TRUE)
    {
        // INFO("%c", fgetc(input_fp));
        seq_t *seq = fasta_reader_next(reader);
        if (seq == NULL)
            break;
        count_occ(count, seq);
        free_seq(seq);
    }
    free_fasta_reader(reader);
    return count;
}

void construct_name(int rc, int l, int sp, int index, char *buffer)
{
    if (sp == -1)
    {
        if (rc)
            index2oligo_rc_char(index, l, buffer);
        else
            index2oligo_char(index, l, buffer);

    }
    else
    {
        char buffer1[128];
        char buffer2[128];
        int m = (l - sp) / 2;
        if (rc)
        {
            index2oligo_rc_char(index, l - sp, buffer1);
            index2oligo_rc_char(index, l - sp, buffer2);
        }
        else
        {
            index2oligo_char(index, l - sp, buffer1);
            index2oligo_char(index, l - sp, buffer2);
        }
        buffer1[m] = '\0';
        buffer2[2 * m] = '\0';
        sprintf(buffer, "%sn{%d}%s", &buffer1[0], sp, &buffer2[m]);
    }
}

void construct_id(int rc, char *name, char *name_rc, char *buffer)
{
    if (rc)
        sprintf(buffer, "%s|%s", name, name_rc);
    else
        sprintf(buffer, "%s", name);
}

void write_count(FILE *output_fp, count_t *count)
{
    // write header
    fprintf(output_fp, "#seq\tid\tobserved_freq\tocc\n");
    
    // binomial stats on count table
    char name[32];
    char name_rc[32];
    char id[128];
    long n;
    int i;
    for (i = 0; i < count->size; i++)
    {
        n = count->count_table[i];
        if (n == 0)
            continue;
        double f = n / (double) count->position_count;
        construct_name(0, count->oligo_length, count->spacing, i, name);
        construct_name(1, count->oligo_length, count->spacing, i, name_rc);
        construct_id(count->rc, name, name_rc, id);
        fprintf(output_fp, "%s\t%s\t%.13f\t%ld\n", name, id, f, n);
    }
}

// void write_count(FILE *output_fp, count_t *count, int count_only)
// {
//     // write header
//     if (count_only)
//         fprintf(output_fp, "#seq\tid\tobserved_freq\tocc\n");
//     else
//         fprintf(output_fp, "#seq\tid\texp_freq\tocc\texp_occ\tocc_P\tocc_E\tocc_sig\n");
//     
//     // binomial stats on count table
//     char name[16];
//     char name_rc[16];
//     char oligo[16];
//     char id[128];
//     long N      = 0;
//     long n      = 0;
//     double p    = 1.0;
//     double pv   = 1.0;
//     double ev   = 1.0;
//     double sig  = 0.0;
//     double freq = 0.0;
//     long n_exp  = 1.0;
// 
//     for (i = 0; i < count->size; i++)
//     {
//         n = count->count_table[i];
//         if (n == 0)
//             continue;
//         N = count->position_count;
//         freq = n / (double) N;
//         index2oligo(i, oligo_length, oligo);
//         index2oligo_char(i, oligo_length, name);
//         if (rc)
//         {
//             index2oligo_rc_char(i, oligo_length, name_rc);
//             sprintf(id, "%s|%s", name, name_rc);
//         }
//         else
//         {
//             sprintf(id, "%s", name);
//         }
//         if (!count_only)
//         {
//             p = markov_P(bg, oligo, 0, oligo_length);
//             //INFO("oligo=%s p=%G", name, p);
//             if (rc && !count->palindromic[i])
//                 p *= 2;
//             n_exp = N * p;
//             // if (FALSE && n <= n_exp)
//             // {
//             //     pv = 1.0;
//             //     ev = 1000.0;
//             //     sig = -1000.0;
//             // }
//             // else
//             // {
//                 pv = pbinom(n, N, p);
//                 ev = pv * count->test_count;
//                 sig = -log10(ev);
//             // }
//         }
//         if (count_only)
//             fprintf(output_fp, "%s\t%s\t%.13f\t%ld\n", name, id, freq, n);
//         else
//             fprintf(output_fp, "%s\t%s\t%.13f\t%ld\t%ld\t%G\t%G\t%G\n", 
//                     name, id, p, n, n_exp, pv, ev, sig);
//     }
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
    int count_only          = FALSE;
    int spacing             = -1;
    //int spacing_range[2]    = {-1, -1};
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
        else if (strcmp(argv[i], "-expfreq") == 0) 
        {
            ENSURE(argc > i + 1, "-v requires a number (0, 1 or 2)");
            bg_filename = argv[++i];
        } 
        else if (strcmp(argv[i], "-count") == 0) 
        {
            count_only = TRUE;
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
        else if (strcmp(argv[i], "-sp") == 0)
        {
            ASSERT(argc > i + 1, "-sp requires an argument");
            //parse_spacing(argv[++i], spacing_range);
            //spacing = spacing_range[0];
            spacing = atoi(argv[++i]);
            //printf("spacing =%d %d", spacing_range[0], spacing_range[1]);
        } 
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

    FILE *output_fp = stdout;
    if (output_filename) 
    {
       output_fp = fopen(output_filename, "w");
       //ENSURE(output_fp != NULL, "can not write to file '%s'", output_filename);
       ENSURE(output_fp != NULL, "can not write to file");
    }

    FILE *input_fp = stdin;
    if (input_filename) 
    {
        input_fp = fopen(input_filename, "r");
        //ENSURE(input_fp != NULL, "can not read from file '%s'", input_filename);
        ENSURE(input_fp != NULL, "can not read from file");
    }

    if (spacing != -1)
    {
        oligo_length = oligo_length * 2 + spacing;
    }

    // read fasta file & compute count table
    count_t *count = count_in_file(input_fp, oligo_length, spacing, rc, noov);
    if (count_only)
    {
        write_count(output_fp, count);
        free_count(count);
    }

    // fflush(output_fp);
    if (input_filename)
        fclose(input_fp);
    if (output_filename)
        fclose(output_fp);
    return 0;
}
