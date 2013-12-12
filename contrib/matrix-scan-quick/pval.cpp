#include "utils.h"
#include "pval.h"

pvalues_t *new_pvalues()
{
    pvalues_t *pvalues = (pvalues_t *) malloc(sizeof(pvalues_t));
    pvalues->w_min = 1E300;
    pvalues->w_max = -1E300;
    pvalues->size = 0;
    pvalues->data = (double *) malloc(sizeof(double) * 1);
    return pvalues;
}

void free_pvalues(pvalues_t *pvalues)
{
    if (pvalues == NULL)
        return;
    free(pvalues->data);
    free(pvalues);
}

double score2pvalue(pvalues_t *pvalues, double score)
{
    if (pvalues == NULL)
        return 0.0;
    double index  = ((score - (double) pvalues->w_min) / (double) (pvalues->w_max - pvalues->w_min)) * ((double) pvalues->size);
    index = MIN(MAX(round(index) - 1, 0), pvalues->size - 1);
    return pvalues->data[(int) index];
}

/*
  score distrib file reader (generated with matrix-distrib)

#weight proba   cum     Pval    ln_Pval log_P   sig
-5.9    NA      0.0e+00 1.0e+00 0.000   0       0.000
-5.8    5.6e-02 5.6e-02 1.0e+00 0.000   0       0.000

*/
pvalues_t *read_distrib(char *filename)
{
    // open file
    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
        ERROR("unable to open '%s'", filename);

    if (fscanf(fp, "#%*s\t%*s%*s\t%*s%*s\t%*s\t%*s\n") < 0)
        ERROR("unable read '%s'", filename);
        

    float weight, Pval;
    pvalues_t *table = new_pvalues();
    
    while (!feof(fp))
    {
        int r = fscanf(fp, "%G\t%*s\t%*G\t%G\t%*G\t%*G\t%*G\n", &weight, &Pval);
        if (r != 2)
            break;
        table->w_min = MIN(table->w_min, weight);
        table->w_max = MAX(table->w_max, weight);
        table->data = (double *) realloc(table->data, sizeof(double) * (table->size + 1));
        table->data[table->size] = Pval;
        table->size += 1;
        //printf("w=%G pval=%G\n", weight, Pval);
    }

    if (table->size == 0)
        ERROR("invalid distrib file '%s", filename);

    // close stream
    fclose(fp);

    return table;
}
