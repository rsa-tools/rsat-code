#include "scan.h"
#include "sampler.h"

SITES empty_sites(int n, int l)
{
    SITES motif;
    for (int i = 0; i < n; i++)
    {
        Site s = Site(-1,-1,l);
        s.score = -1E300;
        motif.push_back(s);
    }
    return motif;
}

static inline
void try_add_sites(SITES &sites, int s, int p, double score)
{
    // find site with min score
    double min = 1E300;
    int minpos = -1;
    int i;
    
    // find an open slot
    for (i = 0; i < (int) sites.size(); i++)
    {
        if (sites[i].score == -1E300)
        {
            min = sites[i].score;
            minpos = i;
            break;
        }
    }

    // find a slot with minimal score
    if (minpos == -1)
    {
        for (i = 0; i < (int) sites.size(); i++)
        {
            if (sites[i].score < min)
            {
                min = sites[i].score;
                minpos = i;
            }
        }
    }
    
    // update site if needed
    if (score > sites[minpos].score)
    {
        //DEBUG("ADD s=%d p=%d i=%d score=%G %G", s, p, minpos, score, sites[minpos].score);
        sites[minpos].s = s;
        sites[minpos].p = p;
        sites[minpos].score = score;
        
    }
}

void check_sites(SITES &sites)
{
    int i;
    for (i = 0; i < (int) sites.size(); i++)
    {
        if (sites[i].p == -1)
        {
            //DEBUG("i = %d", i);
            ERROR("invalid motif constructed from input matrix");
        }
    }
}

static inline
int word_is_valid(int *seq, int start, int l)
{
    int i;
    for (i = start; i < start + l; i++)
    {
        if (seq[i] < 0)
            return 0;
    }
    return 1;
}

SITES matrix_scan(Sequences &sequences, Array &matrix, Markov &bg, Parameters &params)
{
    int l = matrix.J;
    SITES sites = empty_sites(params.n, l + params.flanks * 2);
    for (int s = 0; s < sequences.size(); s++)
    {
        int maxpos = sequences[s].size() - (l + params.flanks);        
        for (int i = params.flanks; i <= maxpos; i++)
        {
            if (!word_is_valid(sequences[s].data, i - params.flanks, l + params.flanks * 2))
                continue;            
            //double W = matrix.sum(&sequences[s].data[i]);
            double W = matrix.logP(&sequences[s].data[i]) - bg.logP(&sequences[s].data[i], l);
            try_add_sites(sites, s, i - params.flanks, W);
            //DEBUG("i=%d, w=%f", i, W);
        }
    }
    
    check_sites(sites);
    return sites;
}
