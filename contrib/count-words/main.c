// 
//  main.c
//  count-words
//  
//
//  
// 

#include <time.h>
#include "utils.h"
#include "count.h"

int VERSION = 20100518;

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
"        count-words\n"   
"\n"
"AUTHOR\n"
"        Matthieu Defrance\n"
"\n"
"DESCRIPTION\n"
"        calculates oligomer frequencies from a set of sequences\n"
"\n"
"CATEGORY\n"
"        sequences\n"
"        pattern discovery\n"
"\n"
"USAGE\n"
"        count-words -l length -i inputfile\n"
"\n"
"ARGUMENTS\n"
"    INPUT OPTIONS\n"
"        --version        print version\n"
"        -v #             change verbosity level (0, 1, 2)\n"
"        -l #             set oligomer length to # (monad size when using dyads)\n"
"        -i #             input filename\n"
"        -2str            add reverse complement\n"
"        -1str            do not add reverse complement\n"
"        -sp #-#          spacing between the two parts of the dyads from # to # \n"
"        -noov            do not allow overlapping occurrences\n"
"        -grouprc         group reverse complement with the direct sequence\n"
"        -nogrouprc       do not group reverse complement with the direct sequence\n"
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
        tk = strtok (arg, "-");
        if (tk != NULL)
        {
            //DEBUG("tk: ", tk);
            values[count++] = atoi(tk);
        }
        arg = NULL;
    } while (count != 2 && tk != NULL);
    return count;
}

int main(int argc, char *argv[])
{
    char *input_filename = NULL;
    char *output_filename = NULL;
    int add_rc = TRUE;
    int noov = FALSE;
    int oligo_length = 1;
    int spacing = -1;
    int spacing_tab[2] = {0, 0};
    int grouprc = TRUE;
    
    // options
    if (argc == 1) 
    {
        usage(argv[0]);
        exit(0);
    }

    int i;
    for (i = 1; i < argc; i++) 
    {
        if (strcmp(argv[i], "--help") == 0) 
        {
            help(argv[0]);
            exit(0);
        } 
        else if (strcmp(argv[i], "--version") == 0) 
        {
            printf("%d\n", VERSION);
            exit(0);
        } 
        else if (strcmp(argv[i], "-v") == 0) 
        {
            ASSERT(argc > i + 1, "-v requires a nummber (0, 1 or 2)");
            VERBOSITY = atoi(argv[++i]);
            ASSERT(VERBOSITY >= 0 && VERBOSITY <= 2, "invalid verbosity level (should be 0, 1 or 2)");
        } 
        else if (strcmp(argv[i], "-h") == 0) 
        {
            help(argv[0]);
            exit(0);
        } 
        else if (strcmp(argv[i], "-1str") == 0) 
        {
            add_rc = FALSE;
        } 
        else if (strcmp(argv[i], "-2str") == 0) 
        {
            add_rc = TRUE;
        } 
        else if (strcmp(argv[i], "-noov") == 0) 
        {
            noov = TRUE;
        } 
        else if (strcmp(argv[i], "-grouprc") == 0) 
        {
            grouprc = TRUE;
        } 
        else if (strcmp(argv[i], "-nogrouprc") == 0) 
        {
            grouprc = FALSE;
        } 
        else if (strcmp(argv[i], "-l") == 0)
        {
            ASSERT(argc > i + 1, "-l requires a nummber");
            oligo_length = atoi(argv[++i]);
        } 
        else if (strcmp(argv[i], "-sp") == 0)
        {
            ASSERT(argc > i + 1, "-sp requires an argument");
            parse_spacing(argv[++i], spacing_tab);
            spacing = spacing_tab[0];
            //printf("spacing =%d %d", spacing_tab[0], spacing_tab[1]);
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
            ERROR("invalid option %s", argv[i]);
        }
    }

    // monitor start & end time
    time_t rawtime;
    time(&rawtime);
    struct tm * start_time;
    start_time = localtime(&rawtime);

    FILE *input_fp = NULL; //stdin;
    FILE *output_fp = stdout;
    if (input_filename) 
    {
        input_fp = fopen(input_filename, "r");
    }
    if (!input_fp) 
    {
       ERROR("can not read from file '%s'", input_filename);
    }
    if (output_filename) 
    {
       output_fp = fopen(output_filename, "w");
       if (!output_fp) 
       {
           ERROR("can not write to file '%s'", output_filename);
       }
    }
    
    count_in_file(input_fp, output_fp, oligo_length, spacing, add_rc, noov, grouprc, argc, argv, 1);

    if (spacing != -1)
    {
        for (spacing = spacing_tab[0] + 1; spacing <= spacing_tab[1]; spacing++)
        {
            fseek(input_fp, SEEK_SET, 0);
            count_in_file(input_fp, output_fp, oligo_length, spacing, add_rc, noov, grouprc, argc, argv, 0);
        }
    }

    fflush(output_fp);
    if (input_filename)
        fclose(input_fp);
    if (output_filename)
        fclose(output_fp);

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
