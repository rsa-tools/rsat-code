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

int scan_seq(FILE *fout,         // output file
	     seq_t *seq,         // sequence to scan
	     int s,              
	     Array &matrix,      // position-specific scoring matrix
	     Markov &bg,         // background model
	     values_t *values,   
	     double threshold,   // shreshold on either weight score or P-value
	     int rc,             // if 1, also scan reverse complementary sequence
	     pvalues_t *pvalues, 
	     int origin,         // reference for the position 0
	     int offset,         // offset added to each position
	     char *matrix_name,  // name of the position-specific scoring matrix, printed in 3rd column
	     int *scanned_pos,    
	     int first_hit,      // if 1, only print the first hit per sequence
	     int best_hit        // if 1, only print the best hit per sequence
	     )
{
    char buffer[256];
    int l = matrix.J;
    ASSERT(l < 256, "invalid matrix size");
    int a = 0;
    int b = 0;
    seq_t *seqrc = NULL;

    // Best values for option -best_hit_per_seq
    double bestW = -9999999999; // Best weight score for current sequence
    double bestPval = 1;
    int bestA = 0;
    int bestB = 0;
    char bestS;
    //    char *bestWline;            // line to print for the best score of the current sequence
      
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

        a = a - offset;
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

		// Store best hit if option -best_hit_per_seq has been activated
		if (best_hit)
		{
		  if (W > bestW) {
		    bestW = W;
		    bestA = a;
		    bestB = b;
		    bestS = 'D';
		    if (pvalues != NULL)
		      bestPval = Pval;

		    // snprintf(bestWline, "%s\t%s\t%s\t%c\t%d\t%d\t%s\t%G", seq->name, "site", matrix_name, 'D', a, b, buffer, W);
		    // if (pvalues != NULL)
		    //   snprintf(bestWline, "\t%G", Pval);
		    // snprintf(bestWline, "\n");
		  }

		// print current hit
		} else
		{
                  fprintf(fout, "%s\t%s\t%s\t%c\t%d\t%d\t%s\t%G", seq->name, "site", matrix_name, 'D', a, b, buffer, W);
                  if (pvalues != NULL)
                    fprintf(fout, "\t%G", Pval);
                  fprintf(fout, "\n");
		}
            }

            if (first_hit)
                break;
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

		if (best_hit)
		{
		  if (Wrc > bestW) {
		    bestW = Wrc;
		    bestA = a;
		    bestB = b;
		    bestS = 'R';
		    if (pvalues != NULL)
		      bestPval = Pval_rc;
		    
		    // snprintf(bestWline, "%s\t%s\t%s\t%c\t%d\t%d\t%s\t%G", seq->name, "site", matrix_name, 'R', a, b, buffer, Wrc);
		    // if (pvalues != NULL)
		    //   snprintf(bestWline, "\t%G", Pval_rc);
		    // snprintf(bestWline, "\n");
		  }
		    
		} else
		{

                  fprintf(fout, "%s\t%s\t%s\t%c\t%d\t%d\t%s\t%G", seq->name, "site", matrix_name, 'R', a, b, buffer, Wrc);
                  if (pvalues != NULL)
                      fprintf(fout, "\t%G", Pval_rc);
                  fprintf(fout, "\n");
		}
            }
            
            if (first_hit)
                break;
        }
    }
    
    if (rc)
        free_seq(seqrc);

    if (best_hit)
    {
      fprintf(fout, "%s\t%s\t%s\t%c\t%d\t%d\t%s\t%G", seq->name, "site", matrix_name, bestS, bestA, bestB, buffer, bestW);
      if (pvalues != NULL)
	fprintf(fout, "\t%G", bestPval);
      fprintf(fout, "\n");
    }
    
    return 1;
}
