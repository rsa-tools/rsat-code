#include "matrix.h"
#include "utils.h"

/*
  restrictive matrix reader
  only well formatted matrices in tab format can be read
  base order must be A C G T
  SHOULD BE REWRITTEN IN A LESS RESTRICTIVE WAY

  example of valid format:
; motif.length                  6       
; motif.conservation            0.91    
; 
a |     0.030   0.030   0.030   0.910   0.910   0.030
c |     0.030   0.910   0.910   0.030   0.030   0.030
g |     0.910   0.030   0.030   0.030   0.030   0.910
t |     0.030   0.030   0.030   0.030   0.030   0.030
*/
int read_matrix(Array &matrix, char *filename, double pseudo)
{
    // open file
    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
        ERROR("unable to open '%s'", filename);

    // read data
    float p[4][256];
    char base[4];
    int l = 0;
    int c = -1;
    while (l < 4 && !feof(fp))
    {
        if (fscanf(fp, "%c", &base[l]) < 0)
            return 0;

        if (base[l] == ';') // skip comments
        {
            char buffer[1024];
            if (fgets(buffer, 1024, fp) == NULL)
                return 0;
            continue;
        }

        if (fscanf(fp, " | ") < 0)
            return 0;
        c = -1;
        while (++c < 256)
        {
            float v;
            if (fscanf(fp, "%f", &v) != 1)
                break;
            p[l][c] = v;
        }
        l++;
    }

    // check validity
    ASSERT(base[0] == 'a' || base[0] == 'A', "invalid matrix format");
    ASSERT(base[1] == 'c' || base[1] == 'C', "invalid matrix format");
    ASSERT(base[2] == 'g' || base[2] == 'G', "invalid matrix format");
    ASSERT(base[3] == 't' || base[3] == 'T', "invalid matrix format");

    // convert to matrix
    matrix.alloc(l, c);
    matrix.pseudo = pseudo;

    int i;
    for (i = 0; i < l; i++)
    {
        int j;
        for (j = 0; j < c; j++)
            matrix[i][j] = p[i][j];
    }

    // security check
    while (!feof(fp))
    {
        if (fscanf(fp, "%c", &base[0]) < 0)
            return 0;

        if (base[0] == ';' || base[0] == ' ') // skip comments
        {
            char buffer[1024];
            if (fgets(buffer, 1024, fp) == NULL)
                return 0;
            continue;
        }
        else if (base[0] == '\n')
        {
            
        }
        else
        {
            WARNING("only first matrix in %s is used", filename);
            break;;
        }
    }

    // close stream
    fclose(fp);
    
    return 1;
}
