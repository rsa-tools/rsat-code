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
//         fprintf(stdout, "W=%.3f\n", W);
//     }
//     return 1;
// }

int scan_seq(seq_t *seq, int s, Array &matrix, Markov &bg, values_t *values, double threshold, int rc)
{
    int l = matrix.J;
    seq_t *seqrc = NULL;
    if (rc)
        seqrc = new_seq_rc(seq);
    int maxpos = seq->size - l;
    int i;
    for (i = 0; i <= maxpos; i++)
    {
        if (!word_is_valid(seq, i, l))
            continue;
        double W = matrix.logP(&seq->data[i]) - bg.logP(&seq->data[i], l);
        if (W < threshold)
            continue;

        if (values != NULL)
        {
            values_add(values, W);
        }
        else
        {
            const char *seqstr = "?";
            fprintf(stdout, "%d\t%s\t%s\t%c\t%d\t%d\t%s\t%.3f\n", s, "site", "matrix", 'D', i + 1, i + l, seqstr, W);
        }

        if (!rc)
            continue;

        double Wrc = matrix.logP(&seqrc->data[maxpos - i]) - bg.logP(&seqrc->data[maxpos - i], l);
        if (Wrc < threshold)
            continue;
        if (values != NULL)
        {
            values_add(values, Wrc);
        }
        else
        {
            const char *seqrcstr = "?";
            fprintf(stdout, "%d\t%s\t%s\t%c\t%d\t%d\t%s\t%.3f\n", s, "site", "matrix", 'R', i + 1, i + l, seqrcstr, Wrc);
        }
    }
    
    if (rc)
        free_seq(seqrc);
    
    return 1;
}






