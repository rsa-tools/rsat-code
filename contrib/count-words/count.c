// 
//  count.c
//  count-words
//  
//  Created by Matthieu on 2008-11-04.
// 

#include "count.h"

#define ALPHABET_SIZE 4
static long position_count = 0;
static long total_count = 0;

// ===========================================================================
// =                            Fast fasta reader
// ===========================================================================
#define MAIN_BUFFER_SIZE 100000
static char main_buffer[MAIN_BUFFER_SIZE];
static int main_buffer_pos = 0;

inline char get_char(FILE *fp)
{
    if (main_buffer_pos >= MAIN_BUFFER_SIZE) 
    {
        int read_size = fread(main_buffer, 1, MAIN_BUFFER_SIZE, fp);
        if (read_size != MAIN_BUFFER_SIZE)
            main_buffer[read_size] = EOF;
        main_buffer_pos = 0;
    }
    return main_buffer[main_buffer_pos++];
}

int fasta_next(string_buffer_t *buffer, FILE *fp)
{
    // read header
    char c;
    do 
    {
        c = get_char(fp);
    } while (c != EOF && c != '\n');

    // read sequence
    int i = 0;
    do 
    {
        c = get_char(fp);
        if (c == EOF)
          break;
        if (c != '\n' && c != '>') 
        {
          STRING_BUFFER_REALLOC(buffer, i + 1);
          buffer->data[i++] = c;
          buffer->current_size = i;
          // string_buffer_set_char(buffer, c, i++);
        }

    } while (c != '>');

    STRING_BUFFER_REALLOC(buffer, i + 1);
    buffer->data[i] = '\0';

    if (c == '>')
        return TRUE;
    else
        return FALSE;
}
// ===========================================================================
// =                            Count oligos
// ===========================================================================
long *new_count_array(int oligo_length)
{
    int i;
    int size = 1;
    for (i = 0; i < oligo_length; i++)
        size *= ALPHABET_SIZE;
    
    long *array = malloc(sizeof(long) * size);
    for (i = 0; i < size; i++)
        array[i] = 0;

    return array;
}

void init_last_position_array(long *array, int oligo_length)
{
    int i;
    int size = 1;
    for (i = 0; i < oligo_length; i++)
        size *= ALPHABET_SIZE;
    for (i = 0; i < size; i++)
        array[i] = -oligo_length;
}

inline int oligo2int(char *string, int pos, int length)
{
    int value = 0;
    int S = 1;
    int i;
    for (i = length - 1; i >= 0; i--) 
    {
        switch (string[pos + i]) 
        {
        case 'A':
        case 'a':
            value += S * 0;
            break;
        case 'C':
        case 'c':
            value += S * 1;
            break;
        case 'G':
        case 'g':
            value += S * 2;
            break;
        case 'T':
        case 't':
            value += S * 3;
            break;
        default:
            return -1;
        }
        S *= ALPHABET_SIZE;
    }
    return value;
}

inline int oligo2int_rc(char *string, int pos, int length)
{
    int value = 0;
    int S = 1;
    int i;
    for (i = 0; i < length; i++) 
    {
        switch (string[pos + i]) 
        {
        case 'A':
        case 'a':
            value += S * 3;
            break;
        case 'C':
        case 'c':
            value += S * 2;
            break;
        case 'G':
        case 'g':
            value += S * 1;
            break;
        case 'T':
        case 't':
            value += S * 0;
            break;
        default:
            return -1;
        }
        S *= ALPHABET_SIZE;
    }
    return value;
}

void count_occ(long *count_table, long *last_position, long *overlapping_occ, char *string, \
                int oligo_length, int add_rc, int noov, int grouprc)
{
    if (noov)
        init_last_position_array(last_position, oligo_length);

    int i;
    int index = -1;
    int index_f = -1;  // forward strand
    int index_r = -1;  // reverse strand
    
    int string_size = (int) strlen(string);
    for (i = 0; i < string_size - oligo_length + 1; i++) 
    {

            // compute index
            index = index_f = oligo2int(string, i, oligo_length);
            if (add_rc) 
            {
                index_r = oligo2int_rc(string, i, oligo_length);
                index = MIN(index, index_r);
            }

            if (index == -1) // bad position
                continue;

            // increment position counter
            position_count++;
            
            // overlapping occurrences
            if (noov && last_position[index] + oligo_length - 1 >= i) 
            {
                overlapping_occ[index]++;
                if (add_rc && !grouprc)
                    overlapping_occ[MAX(index_r, index_f)]++;
                continue;
            }
            if (noov)
                last_position[index] = i;

            // count
            count_table[index]++;
            total_count++;

             // count on other strand when occurrences are not grouped
            if (add_rc && !grouprc)
                count_table[MAX(index_r, index_f)]++;
            
    }
}

void print_header(FILE *output_fp, int oligo_length, int noov, int argc, char *argv[])
{
    if (VERBOSITY < 1)
        return;
        
    // print command line        
    fprintf(output_fp, "; ");
    int i;
    for (i = 0; i < argc; i++) 
    {
        fprintf(output_fp, "%s ", argv[i]);            
    }
    fprintf(output_fp, "\n");
    fprintf(output_fp,
        "; oligomer length              	%d\n", oligo_length);

    fprintf(output_fp, 
        "; column headers\n"
        ";    1    seq    oligomer sequence\n"
        ";    2    id     oligomer identifier\n"
        ";    3    freq   relative frequencies (occurrences per position)\n" 
        ";    4    occ    occurrences\n"
    );

    if (noov) 
    {
        fprintf(output_fp, 
            ";    5    ovl_occ    overlapping occurrences\n"
        );
    }
}

void print_count_array(FILE *output_fp, long *count_array, long *overlapping_occ, int oligo_length, int add_rc)
{
    char letter[ALPHABET_SIZE] = "acgt";
    int i, k;
    int size = 1;
    for (i = 0; i < oligo_length; i++)
        size *= ALPHABET_SIZE;    
    
    char oligo_buffer[oligo_length + 1];
    oligo_buffer[oligo_length] = '\0';
    for (i = 0; i < oligo_length; i++)
        oligo_buffer[i] = letter[0];
    // header
    if (overlapping_occ)
        fprintf(output_fp, "#seq\tidentifier\tobserved_freq\tocc\tovl_occ\n");
    else
        fprintf(output_fp, "#seq\tidentifier\tobserved_freq\tocc\n");

    k = 0;
    for (i = 0; i < size; i++) 
    {
        if (count_array[i] != 0) 
        {
            // special case for overlapping_occ
            if (overlapping_occ) 
            {
                fprintf(output_fp, "%s\t%s\t%.13f\t%d\t%d\n", \
                oligo_buffer, oligo_buffer, count_array[i] / (double) position_count, (int) count_array[i], (int) overlapping_occ[i]);
            } 
            else 
            {
                fprintf(output_fp, "%s\t%s\t%.13f\t%d\n", \
                oligo_buffer, oligo_buffer, count_array[i] / (double) position_count, (int) count_array[i]);
            }
        }
        for (k = oligo_length - 1; k >= 0; k--) 
        {
            if (oligo_buffer[k] == 't') 
            {
                oligo_buffer[k] = 'a';
            } 
            else if (oligo_buffer[k] == 'a') 
            {
                oligo_buffer[k] = 'c';
                break;
            } 
            else if (oligo_buffer[k] == 'c') 
            {
                oligo_buffer[k] = 'g';
                break;
            } 
            else if (oligo_buffer[k] == 'g') 
            {
                oligo_buffer[k] = 't';
                break;
            }
        }
    }
}

void count_in_file(FILE *input_fp, FILE *output_fp, int oligo_length, int add_rc, \
    int noov, int grouprc, int argc, char *argv[])
{
    long *count = new_count_array(oligo_length);
    long *last_position = NULL;
    long *overlapping_occ = NULL;
    position_count = 0;
    total_count = 0;
 
    if (noov) 
    {
        last_position = new_count_array(oligo_length);
        overlapping_occ = new_count_array(oligo_length);
    }
        
    string_buffer_t *buffer = new_string_buffer();
    int end = FALSE;
    do 
    {
        end = fasta_next(buffer, input_fp);
        count_occ(count, last_position, overlapping_occ, buffer->data, oligo_length, add_rc, noov, grouprc);
    } while (end != FALSE);

    print_header(output_fp, oligo_length, noov, argc, argv);
    print_count_array(output_fp, count, overlapping_occ, oligo_length, add_rc);
    free(count);
    free(last_position);
    free_string_buffer(buffer);
}
