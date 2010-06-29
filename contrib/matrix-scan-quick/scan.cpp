#include "scan.h"

static inline
int word_is_valid(seq_t *seq, int start, int l)
{
    int i;
    for (i = start; i < start + l; i++)
    {
        if (seq->data[i] < 0)
            return 0;
    }
    return 1;
}

static
void set_buffer(char *buffer, char *seq, int i, int l)
{
    static char BASES[4] = {'A', 'C', 'G', 'T'};
    buffer[l] = '\0';
    int k;
    for (k = 0; k < l; k++)
    {
        buffer[k] = BASES[(int) seq[i + k]];
    }
}

int scan_seq(FILE *fout, seq_t *seq,  int s, Array &matrix, Markov &bg, values_t *values, 
            double threshold, int rc, pvalues_t *pvalues, int origin, char *matrix_name, int *scanned_pos)
{
    char buffer[256];
    int l = matrix.J;
    ASSERT(l < 256, "invalid matrix size");
    int a = 0;
    int b = 0;
    seq_t *seqrc = NULL;
    if (rc)
        seqrc = new_seq_rc(seq);

    int maxpos = seq->size - l;
    int i;
    for (i = 0; i <= maxpos; i++)
    {
        if (!word_is_valid(seq, i, l))
            continue;
        double W;
        if (bg.order == 0)
            W = matrix.logP(&seq->data[i]) - bg.logPBernoulli(&seq->data[i], l);
        else
            W = matrix.logP(&seq->data[i]) - bg.logP(&seq->data[i], l);

        // position
        if (origin == -1) // start
            a = i + 1;
        else if (origin == 0) // center
            a = i - seq->size / 2;
        else // end
            a = i - seq->size;
        b = a + l - 1;

        (*scanned_pos) += 1;

        double Pval = score2pvalue(pvalues, W);
        if ((pvalues == NULL && W >= threshold) || (pvalues != NULL && Pval <= threshold))
        {
            if (values != NULL)
            {
                values_add(values, W);
            }
            else
            {
                //const char *seqstr = "?";
                set_buffer(buffer, seq->data, i, l);
                fprintf(fout, "%s\t%s\t%s\t%c\t%d\t%d\t%s\t%G", seq->name, "site", matrix_name, 'D', a, b, buffer, W);
                if (pvalues != NULL)
                    fprintf(fout, "\t%G", Pval);
                fprintf(fout, "\n");
            }
        }

        if (!rc)
            continue;

        double Wrc;
        if (bg.order == 0)
            Wrc = matrix.logP(&seqrc->data[maxpos - i]) - bg.logPBernoulli(&seqrc->data[maxpos - i], l);
        else
            Wrc = matrix.logP(&seqrc->data[maxpos - i]) - bg.logP(&seqrc->data[maxpos - i], l);

        double Pval_rc = score2pvalue(pvalues, Wrc);
        if ((pvalues == NULL && Wrc >= threshold) || (pvalues != NULL && Pval_rc <= threshold))
        {
            if (values != NULL)
            {
                values_add(values, Wrc);
            }
            else
            {
                //const char *seqrcstr = "?";
                set_buffer(buffer, seqrc->data, seq->size - i - l, l);
                fprintf(fout, "%s\t%s\t%s\t%c\t%d\t%d\t%s\t%G", seq->name, "site", matrix_name, 'R', a, b, buffer, Wrc);
                if (pvalues != NULL)
                    fprintf(fout, "\t%G", Pval_rc);
                fprintf(fout, "\n");
            }
        }
    }
    
    if (rc)
        free_seq(seqrc);
    
    return 1;
}
