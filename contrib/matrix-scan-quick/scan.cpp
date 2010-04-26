#include "scan.h"

// int scan(vector<string> raw_sequences, Sequences &sequences, Array &matrix, Markov &bg, bool rc)
// {
//     double threshold = 0.0;
//     int s = 0;
//     int l = matrix.J;
// 
//     printf("#seq_id\tft_type\tft_name\tstrand\tstart\tend\tsequence\tweight\n");
// 
//     for (s = 0; s < sequences.size(); s++)
//     {
//         int maxpos = sequences[s].size() - l;
//         int i;
//         for (i = 0; i <= maxpos; i++)
//         {
//             double W = matrix.logP(&sequences[s][i]) - bg.logP(&sequences[s][i], l);
//             if (W > threshold)
//             {
//                 const char *seq = raw_sequences[s].substr(i, l).c_str();
//                 printf("%d\t%s\t%s\t%c\t%d\t%d\t%s\t%.3f\n", s + 1, "site", "matrix", 'D', i + 1, i + l, seq, W);
//             }
//         }
//     }
//     return 1;
// }

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

// int scan_seq_distrib(seq_t *seq, Array &matrix, Markov &bg, values_t *values)
// {
//     double threshold = -1000.0;
//     int l = matrix.J;
// 
//     int maxpos = seq->size - l;
//     int i;
//     for (i = 0; i <= maxpos; i++)
//     {
//         if (!word_is_valid(seq, i, l))
//             continue;
//         double W = matrix.logP(&seq->data[i]) - bg.logP(&seq->data[i], l);
//         // if (W >= threshold)
//         //     values_add(values, W);
//         // if (W > 0)
//         fprintf(fout, "W=%.3f\n", W);
//     }
//     return 1;
// }



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

int scan_seq(FILE *fout, seq_t *seq, int s, Array &matrix, Markov &bg, values_t *values, double threshold, int rc, pvalues_t *pvalues, int origin)
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

        if (W < threshold)
            continue;

        double Pval = score2pvalue(pvalues, W);

        // position
        if (origin == -1) // start
            a = i + 1;
        else if (origin == 0) // center
            a = i - seq->size / 2;
        else // end
            a = i - seq->size;
        b = a + l - 1;

        if (values != NULL)
        {
            values_add(values, W);
        }
        else
        {
            //const char *seqstr = "?";
            set_buffer(buffer, seq->data, i, l);
            fprintf(fout, "%d\t%s\t%s\t%c\t%d\t%d\t%s\t%G\t%G\t%G\t%G\n", s, "site", "matrix", 'D', a, b, buffer, W, 0.0, 0.0, Pval);
        }

        if (!rc)
            continue;

        double Wrc;
        if (bg.order == 0)
            Wrc = matrix.logP(&seqrc->data[maxpos - i]) - bg.logPBernoulli(&seqrc->data[maxpos - i], l);
        else
            Wrc = matrix.logP(&seqrc->data[maxpos - i]) - bg.logP(&seqrc->data[maxpos - i], l);

        if (Wrc < threshold)
            continue;

        double Pval_rc = score2pvalue(pvalues, Wrc);

        if (values != NULL)
        {
            values_add(values, Wrc);
        }
        else
        {
            //const char *seqrcstr = "?";
            set_buffer(buffer, seqrc->data, seq->size - i - l, l);
            fprintf(fout, "%d\t%s\t%s\t%c\t%d\t%d\t%s\t%G\t%G\t%G\t%G\n", s, "site", "matrix", 'R', a, b, buffer, Wrc, 0.0, 0.0, Pval_rc);
        }
    }
    
    if (rc)
        free_seq(seqrc);
    
    return 1;
}






