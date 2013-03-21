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
Array read_matrix(FILE *fp)
{
    // end of stream ?
    if (feof(fp) != 0)
        return Array(0,0);
        
    // read data
    float p[4][256];
    char base[4];
    int l = 0;
    int c = -1;

    while (!feof(fp))
    {
        int junk = getc(fp);
        ungetc(junk, fp);
        if (junk == '/')
        {
            while (!feof(fp) && junk != '\n')
            {
                junk = getc(fp);
            }
            
            break;
        }
        
        if (fscanf(fp, "%c", &base[l]) == EOF)
            break;

        if (base[l] == ';') // skip comments
        {
            char buffer[1024];
            if (fgets(buffer, 1024, fp) == NULL)
                break;
            if (feof(fp))
                return Array(0, 0);
            else
                continue;
        }

        if (fscanf(fp, " | ") == EOF)
            break;
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
    ASSERT(c >= 1, "invalid matrix format");
    Array matrix = Array(l, c);
    int i;
    for (i = 0; i < l; i++)
    {
        int j;
        for (j = 0; j < c; j++)
            matrix[i][j] = p[i][j];
    }
    
    return matrix;
}
