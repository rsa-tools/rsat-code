#include "scan.h"
#include "sampler.h"

SITES empty_sites(int n, int l)
{
    SITES motif;
    for (int i = 0; i < n; i++)
        motif.push_back(Site(-1,-1,l));
    return motif;
}

static inline
void try_add_sites(SITES &sites, int s, int p, double score)
{
    // find site with min score
    double min = 1E300;
    int minpos = -1;
    int i;
    for (i = 0; i < (int) sites.size(); i++)
    {
        if (sites[i].score < min)
        {
            min = sites[i].score;
            minpos = i;
        }
    }
    
    // update site if needed
    if (score > sites[minpos].score)
    {
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
            ERROR("invalid motif constructed from input matrix");
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

SITES scan(vector<string> raw_sequences, Sequences &sequences, Array &matrix, Markov &bg, int n)
{
    int l = matrix.J;
    SITES sites = empty_sites(n, l);
    
    for (int s = 0; s < sequences.size(); s++)
    {
        int maxpos = sequences[s].size() - l;        
        for (int i = 0; i <= maxpos; i++)
        {
            if (!word_is_valid(sequences[s].data, i, l))
                continue;            
            double W = matrix.sum(&sequences[s].data[i]);
            try_add_sites(sites, s, i, W);
        }
    }
    
    check_sites(sites);
    return sites;
}
