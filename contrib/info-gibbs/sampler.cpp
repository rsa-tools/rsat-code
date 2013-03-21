#include "sampler.h"
#include "scan.h"

int SEED = time(NULL);
double PSEUDO = 1.0;

//MTRand_closed mt(SEED);
//#define RAND mt()

/***************************************************************************
 *                                                                         *
 *  LLR (log likelihood ratio) table
 *                                                                         *
 ***************************************************************************/
/*
 *   llr =   Sum_k log P(w^k | M) - log P(w^k | B)
 *   llr = Sum_i Sum_k log P(w_i^k | M) - Sum_k log P(w^k | B)
 *         ----------------------------   ----------------
 *               alpha                        beta
 *   
 *   llr =~ 1/N Sum_i Sum_j C_ij /n log C_ij /n - Sum_k log P(w^k | B)
 */
struct LLR
{
    int l;
    int n;
    double *p;
    double **T; 
    double beta;
    double alpha;
    double P;
    int i,j,k;
    Array Tab;
    //Markov markov;
    Sequence *sequences; // avoid sequences copy
    double **logp;

    LLR(SITES &motif, int r, Array &matrix, int length, 
        Markov &markov_model, Sequences &seqs, double **logprob, double pseudo = PSEUDO)
    {
        l         = length;
        n         = motif.size();
        p         = markov_model.priori;
        sequences = seqs.data; // sequences store only data !
        logp      = logprob;
        
        // compute the matrix
        for (i = 0; i < l; i++)
        {
            for (j = 0; j < matrix.J; j++)
                matrix[i][j] = 0;
            for (k = 0; k < n; k++)
            {
                if (k == r)
                    continue;
                if (motif[k].space > 0 && i >= motif[k].m1)
                    j = sequences[motif[k].s][motif[k].p+motif[k].space+i];
                else
                    j = sequences[motif[k].s][motif[k].p+i];
                matrix[i][j] += 1;
            }
        }

        // compute alpha
        int J;
        Tab.alloc(l, matrix.J);

        for (i = 0; i < l; i++)
        {
            for (j = 0; j < matrix.J; j++)
            {
                alpha = 0.0;
                for (k = 0; k < n; k++)
                {
                    if (k == r)
                        continue;
                    if (motif[k].space > 0 && i >= motif[k].m1)
                        J = sequences[motif[k].s][motif[k].p+motif[k].space+i];
                    else
                        J = sequences[motif[k].s][motif[k].p+i];
                    if (j == J)
                        alpha += log ( (matrix[i][J] + 1 + p[J] * pseudo) / (n + pseudo) );
                    else
                        alpha += log ( (matrix[i][J] + p[J] * pseudo) / (n + pseudo) );                 
                }
                // the new letter contribution
                alpha += log ( (matrix[i][j] + 1 + p[j] * pseudo) / (n + pseudo) );
                Tab.data[i][j] = alpha;
            }           
        }
        T = Tab.get_data();

        // compute beta
        beta = 0.0;
        for (k = 0; k < n; k++)
        {
            if (k == r)
                continue;
            beta += logp[motif[k].s][motif[k].p];
        }
    }

    double llr(Site &site)
    {
        alpha = 0.0;
        for (i = 0; i < l; i++)
        {
            if (site.space > 0 && i >= site.m1)
                j = sequences[site.s][site.p+site.space+i];
            else
                j = sequences[site.s][site.p+i];

            alpha += T[i][j];
        }
        return alpha - beta - logp[site.s][site.p];
    }
};

void count_matrix(Array &matrix, SITES &motif, Sequences &sequences)
{
    int i, j, k;
    int n = (int) motif.size();
    for (i = 0; i < matrix.I; i++)
    {
        for (j = 0; j < matrix.J; j++)
        {
            matrix[i][j] = 0;
        }
        
        for (k = 0; k < n; k++)
        {
            j = (int) sequences[motif[k].s][motif[k].p+i];
            matrix[i][j] += 1;
        }
    }
}

/**
 *  Convert a list of words to a frequency matrix using priori probabilities p
 *  m_{b,i} = ( f_{b,i} + p_i ) / (N + 1)
 */
inline void freq_matrix(Array &matrix, SITES &motif, Sequences &sequences, \
                        Markov &markov, double pseudo=PSEUDO, int r=-1)
{
    int i,j,k;
    int n = (int) motif.size();
    for (i = 0; i < matrix.I; i++)
    {
        for (j = 0; j < matrix.J; j++)
        {
            matrix[i][j] = 0;
        }

        for (k = 0; k < n; k++)
        {
            if (r == k) // remove 
                continue;
            if (motif[k].space > 0 && i >= motif[k].m1) // right dyad part m1---m2
                j = (int) sequences[motif[k].s][motif[k].p+motif[k].space+i];
            else
                j = (int) sequences[motif[k].s][motif[k].p+i];

            matrix[i][j] += 1;
        }
        for (j = 0; j < matrix.J; j++)
            matrix[i][j] = (matrix[i][j] + markov.priori[j] * pseudo) / (n+pseudo);
    }
}

/*
    Relative Entropy (information content)[Hertz 1999]
    Bernoulli
    IC = \sum_{i=1}^w \sum_{b=1}^4 f_{b,i} ln (f_{b,i} / p_b)

    matrix -- frequency matrix
    p      -- priori probability [0.25, 0.25, 0.25, 0.25]
*/
double IC_Bernoulli(SITES &motif, Sequences &sequences, Array &matrix, \
                    Markov &markov, double pseudo=PSEUDO)
{
    // matrix
    freq_matrix(matrix, motif, sequences, markov, pseudo);

    // iseq
    double iseq = 0.0;
    for (int i = 0; i < matrix.I; i++)
    {
        for (int j = 0; j < matrix.J; j++)
        {
            if (matrix[i][j] != 0.0)
                iseq += matrix[i][j] * log(matrix[i][j] / markov.priori[j]);
        }
    }
    return iseq;
}

// do not allow spaces
double PQ_Bernoulli(SITES &motif, Sequences &sequences, Array &matrix, \
                    Markov &markov, double pseudo=PSEUDO, int r = 0)
{
    // matrix
    freq_matrix(matrix, motif, sequences, markov, pseudo, r);
    int j;
    double p = 1.0;
    double q = 1.0;
    for (int i = 0; i < matrix.I; i++)
    {
        j = (int) sequences[motif[r].s][motif[r].p+i];
        if (matrix[i][j] != 0.0)
            p *= matrix[i][j];

        q *= markov.priori[j];
    }
    return p/q;
}

// PM[start..start+l-1] w
double P_M(Array &matrix, int w, int l, int start) 
{
    int N = (int) pow((double)ALPHABET_SIZE, (double)l-1);
    double p = 1.0;
    for (int i = 0; i < l; i++) 
    {
        int j = (int) (w / N);
        w = w % N;
        N = (int) (N / ALPHABET_SIZE);
        p = p * matrix[i+start][j];
    }
    return p;
}

/*
    Information content
    Markov
    See Paper
    IC = A - (X + Y)
*/
double IC_Markov(SITES &motif, Sequences &sequences, Array &matrix, \
                 Markov &markov, double pseudo = PSEUDO)
{
    // special Bernoulli case
    if (markov.order == 0)
        return IC_Bernoulli(motif, sequences, matrix, markov, pseudo);

    // matrix
    freq_matrix(matrix, motif, sequences, markov, pseudo);

    int I = matrix.I; // length
    int J = matrix.J; // ALPHABET_SIZE
    int order = markov.order;
    double A = 0.0;
    //A = sum[i in (0..l-1)] sum[j in (ALPHABET.indices)] m[i][j] * ln m[i][j]
    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < J; j++) 
        {
            A += matrix[i][j] * log(matrix[i][j]);
        }
    }

    // X = sum[wp in (0..S.lastIndex)] P_mi(wp, order, 0) * ln S[wp]
    double X = 0.0;
    int nwp = (int) pow((double) J, (double) order); // 4^order (size of S, stationnay vector)
    for (int wp = 0; wp < nwp; wp++) 
    {
        X += P_M(matrix, wp, order, 0) * log(markov.S[wp]);
    }

    //Y = sum[i in (0..l-1-order)] sum[ws in ALPHABET]
    // sum[wp in (0..S.lastIndex)] P_mi(wp*ALPHABET.size+ws, order+1, i) * ln T[wp][ws]    
    double Y = 0.0;
    for (int i = 0; i < I-order; i++) 
    {
        for (int wp = 0; wp < nwp; wp++) 
        {
            for (int ws = 0; ws < ALPHABET_SIZE; ws++) 
            {
                Y += P_M(matrix, wp * ALPHABET_SIZE + ws, order + 1, i) * log(markov.T[wp][ws]);
            }
        }
    }
    return A - (X + Y);
}

/*
  Log likelihood ratio
             P(w|M)
  \sum_w log -------
             P(w|Bg)
*/
double llr_Bernoulli(SITES &motif, Sequences &sequences, Array &matrix, \
                     Markov &markov, double pseudo = PSEUDO)
{
    // matrix
    freq_matrix(matrix, motif, sequences, markov, pseudo);
    // llr
    int n = (int) motif.size();
    int j;
    double s = 0.0;
    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < matrix.I; i++)
        {
            j = (int) sequences[motif[k].s][motif[k].p+i];
            s += log(matrix[i][j] / markov.priori[j]);
        }
    }
    return s;
}

double llr(SITES &motif, Sequences &sequences, Array &matrix, Markov &markov, double pseudo=PSEUDO)
{
    // matrix
    freq_matrix(matrix, motif, sequences, markov, pseudo);
    // llr
    int n = (int) motif.size();
    int j;
    double s = 0.0;
    for (int k = 0; k < n; k++)
    {
        double logP_m = 0.0;
        for (int i = 0; i < matrix.I; i++)
        {
            j = (int) sequences[motif[k].s][motif[k].p+i];
            logP_m += log(matrix[i][j]);
        }
        s += logP_m - markov.logP(&sequences[motif[k].s][motif[k].p], matrix.I);
    }
    return s;
}

/***************************************************************************
 *                                                                         *
 *  SITES = all available motif positions
 *                                                                         *
 ***************************************************************************/
void print_sites(SITES &sites)
{
    for (unsigned int i = 0; i < sites.size(); i++)
    {
        printf("%d %d\n", sites[i].s + 1, sites[i].p + 1);
    }
}

SITES mask_motif(SITES &sites, SITES &motif)
{
    SITES masked_sites;
    for (int i = 0; i < (int) sites.size(); i++)
    {
        bool remove_site = false;
        for (int j = 0; j < (int) motif.size(); j++)
        {
            if (sites[i].s == motif[j].s && sites[i].p == motif[j].p)
            {
                remove_site = true;
            }
        }
        if (!remove_site)
            masked_sites.push_back(sites[i]);
    }
    return masked_sites;
}

SITES all_sites(Sequences &sequences, int l, int m1, int m2)
{
    SITES sites;
    for (int s = 0; s < sequences.size(); s++)
    {
        int len = sequences[s].size();
        for (int p = 0; p < len-l+1; p++)
        {
            bool is_valid_site = true;
            
            // remove AAAAAAAA, .....
            int max_base_occ = 8;
            for (int base = 0; base < ALPHABET_SIZE; base++)
            {
                int c = 0;
                for (int y = 0; y < l; y++){
                    if (sequences[s][p+y] == base)
                        c++;
                    else
                        c = 0;
                    if (c >= max_base_occ)
                    {
                        is_valid_site = false;
                        break;
                    }
                }
            }            
            
            for (int x = 0; x < l; x++)
            {
                if (sequences[s][p+x] == -1)
                {
                    is_valid_site = false;
                    break;
                }
            }
            if (is_valid_site)
            {
                sites.push_back(Site(s,p,m1,l-m1-m2,m2));
            }
        }
    }
    return sites;
}

SITES remove_neighbours(SITES &allsites, SITES &motif, int nseq, int dmin=0, int r=-1)
{
    if (dmin == 0)
        return allsites;
    
    SITES sites = allsites;
    for (int k = 0; k < (int) motif.size(); k++)
    {
        if (k == r)
            continue;
        for (int i = 0; i < (int) sites.size(); i++)
        {
            if (ABS(sites[i].s - motif[k].s) != 0 && ABS(sites[i].s - motif[k].s) != nseq)
                continue;
            if (ABS(sites[i].p - motif[k].p) <= dmin)
                sites[i].p = -1; //invalidate site
        }
    }
    SITES new_sites;
    for (int i = 0; i < (int) sites.size(); i++)
    {
        if (sites[i].p != -1) // invalid site ?
        {
            new_sites.push_back(sites[i]);
        }
    }
    return new_sites;
}

/***************************************************************************
 *                                                                         *
 *  MOTIF
 *                                                                         *
 ***************************************************************************/
void print_motif(SITES &motif, vector<string> &raw_sequences, Sequences &sequences, int l, double ic, double llr, bool rc=false)
{
    int n = motif.size();

    printf("; log likelihood ratio          %.3f\n", llr);
    printf("; information content           %.3f\n", ic);
    printf("; motif width                   %d\n", l);
    printf("; sites                         %d\n", (int) motif.size());
    printf("; (seq and pos start at 1)\n");
    printf("; seq\tstrand\tpos\tsite\n");

    //sites
    int seqsize = raw_sequences.size();
    int seq;
    char strand_label;

    for (int i = 0; i < n; i++)
    {
        int s = motif[i].s;
        int p = motif[i].p;

        if (rc == true && s >= seqsize / 2)
        {
            strand_label = '-';
            seq = s - seqsize / 2 + 1;
        }
        else
        {
            strand_label = '+';
            seq = s + 1;
        }
        const char *word = raw_sequences[s].substr(p,l).c_str();
        printf("; %d\t%c\t%d\t%s\n", seq, strand_label, p + 1, word);
    }

    //matrix
    Array matrix = Array(l, ALPHABET_SIZE);
    count_matrix(matrix, motif, sequences);

    for (int j = 0; j < ALPHABET_SIZE; j++)
    {
        printf("%c | ", ALPHABET[j]);
        for (int i = 0; i < l; i++)
        {
            printf("\t%d", (int) matrix[i][j]);
        }
        printf("\n");
    }
    printf("//\n");
}

bool inline is_in_sites(Site &site, SITES &sites)
{
    for (unsigned int i = 0; i < sites.size(); i++)
    {
        if (site.p == sites[i].p && site.s == sites[i].s)
        {
            return true;
        }
    }
    return false;
}

SITES random_motif(SITES &allsites, int n, int nseq, int dmin = 0)
{
    SITES motif;
    int j = 0;
    SITES sites = allsites;
    VERBOSE2("generating random motif\n");
    while ((int) motif.size() < n && j++ < n*2 && (int) sites.size() > 0)
    {
        int i = (int) (RAND * sites.size());
        if (!is_in_sites(sites[i], motif))
        {
            motif.push_back(sites[i]);
            sites = remove_neighbours(sites, motif, nseq, dmin);
        }
    }
    if ((int) motif.size() != n)
        WARNING("number of sites is smaller than required");
    return motif;
}

/***************************************************************************
 *                                                                         *
 *  SPEEDUP STRUCTURES
 *                                                                         *
 ***************************************************************************/
struct Is_a_site 
{
    bool **cache;
    Sequences *seqs;

    Is_a_site(Sequences &sequences, SITES &sites)
    {
        int S = sequences.size();
        seqs  = &sequences;
        cache = new bool*[S];
        // alloc & set to false
        for (int s = 0; s < S; s++)
        {
            int len = sequences[s].size();
            cache[s] = new bool[len];
            for (int p=0; p<len; p++)
                cache[s][p] = false;
        }
        // init with sites
        for (int k = 0; k < (int) sites.size(); k++)
        {
            cache[sites[k].s][sites[k].p] = true;
        }
    }
    
    bool is_valid(Site &site)
    {
        if (site.s < 0 or site.s >= seqs->size())
        {
            return false;
        }
        else if (site.p < 0 or site.s >= seqs->data[site.s].size())
        {
            return false;
        }
        else
        {
            return cache[site.s][site.p];
        }
    }
};
 
struct SamplingData 
{
    int l;
    int nsites;          // max number of sites
    int S;               // number of sequences
    double **logp;       // log P(word) cache for each p
    double *cdf;         // cumulative dist function 
    Site *sampled_sites; // sites (corresponds to cdf)

    SamplingData(Sequences &sequences, SITES &sites, Markov &markov, int l)
    {
        S = (int) sequences.size();
        nsites = (int) sites.size();
        cdf = new double[nsites];
        sampled_sites = new Site[nsites];
        logp = new double*[S];

        // alloc logP
        for (int s = 0; s < S; s++)
        {
            int len = sequences[s].size();
            logp[s] = new double[len];
        }

        for (int i = 0; i < nsites; i++)
        {
            int s = sites[i].s;
            int p = sites[i].p; 
            logp[s][p] = markov.logP(&sequences[s][p], l);
        }
    }
};

/***************************************************************************
 *                                                                         *
 *  SHIFTING
 *                                                                         *
 ***************************************************************************/
SITES shifted(SITES &motif, SITES &sites, int delta, Is_a_site &cache)
{
    if (delta == 0)
        return motif;

    SITES shifted_motif = motif;
    for (unsigned int i = 0; i < motif.size(); i++)
    {
        Site site = Site(motif[i].s, motif[i].p+delta, motif[i].m1, motif[i].space, motif[i].m2);
        if (!cache.is_valid(site))
        {
            return motif;
        }
        shifted_motif[i] = site;
    }
    
    return shifted_motif;
}

SITES shift(SITES &motif, Sequences &sequences, SITES &sites, Array &matrix, Markov &markov, Is_a_site &cache)
{
    SITES best_motif;
    SITES current_motif;
    double best_ic = 0.0;
    double current_ic;

    for (int delta = -1; delta <= 1; delta++)
    {
        current_motif = shifted(motif, sites, delta, cache);
        current_ic = IC_Markov(current_motif, sequences, matrix, markov);
        if (current_ic > best_ic)
        {
            best_motif = current_motif;
            best_ic = current_ic;
        }
    }
    return best_motif;
}

/***************************************************************************
 *                                                                         *
 *  SAMPLING
 *                                                                         *
 ***************************************************************************/
int UPDATE = 0;
void sample_update(SITES &allsites, Sequences &sequences, SITES &motif, Array &matrix, Markov &markov,\
                    int l, double temperature, SamplingData &data, int nseq, int dmin=0, int score_type=LLR_SCORE)
{
    int n = motif.size();
    int r = 0;

    double beta = 1.0 / temperature;
    double *cdf = data.cdf;
    double S = 0.0;
    double val = 0.0;
    int i;

    if (score_type == IC_SCORE)
        beta = beta * n * ((n + exp(PSEUDO)) / n);

    // choose words
    r = (int) (RAND * n);

    Site *sampled_sites = data.sampled_sites;
    LLR llr_table = LLR(motif, r, matrix, l, markov, sequences, data.logp, PSEUDO);

    // // set update
    // if (UPDATE == -1)
    //     UPDATE = n;
    
    SITES sites = remove_neighbours(allsites, motif, nseq, dmin, r);
    //SITES &sites = allsites;
    // DEBUG("allsites=%d sites=%d motif=%d\n", (int) allsites.size(), (int) sites.size(), (int) motif.size());
    // DEBUG("motif r=%d", r);
    // print_sites(motif);
    //DEBUG("sites");
    //print_sites(sites);

    int nsites = sites.size();
    double s = 0.0;

    SITES tmp_motif;
    if (score_type == IC_SCORE || score_type == PQ_SCORE || score_type == PQ_LLR_SCORE || score_type == LLR_IC_SCORE)
        tmp_motif = motif;
    
    for (i = 0; i < nsites; i++)
    {
        if (score_type == LLR_SCORE || score_type == LLR_IC_SCORE) 
        {
            //with alpha table
            s = exp(beta * llr_table.llr(sites[i]));
        } 
        else if (score_type == PQ_SCORE || score_type == PQ_LLR_SCORE) 
        {
            tmp_motif[r] = sites[i];
            s = PQ_Bernoulli(tmp_motif, sequences, matrix, markov, PSEUDO, r);
            //DEBUG("s=%.3f", s);
        } 
        else if (score_type == IC_SCORE)
        {
            tmp_motif[r] = sites[i];
            //s = exp(beta * llr(tmp_motif, sequences, matrix, markov, PSEUDO));
            s = exp(beta * IC_Markov(tmp_motif, sequences, matrix, markov, PSEUDO));
        }

        if (s >= 1e300) // avoid math overflow
            s = 1e300;
        S += s;
        sampled_sites[i] = sites[i];
        cdf[i] = S;
        
    }
    if (S >= 1.7e308)
    {
        DEBUG("WARNING math overflow error\n");
        return;
    }

    // choose new site
    // random choose
    val = RAND * S;
    for (i = 0; i < nsites; i++)
    {
        if (cdf[i] > val)
            break;
    }

    // update motif
    Site new_site = sampled_sites[i];
    if (!is_in_sites(new_site, motif))
    {
        motif[r] = new_site;
    }
}

/***************************************************************************
 *                                                                         *
 *  FIND ONE MOTIF
 *                                                                         *
 ***************************************************************************/
struct Result 
{
    SITES motif;
    double ic;
    double llr;
    double avg_llr;
    double avg_ic;
    int l;
};

Result find_one_motif(vector<string> &raw_sequences, Sequences & sequences, SITES &sites, Markov &markov, Parameters &params)
{
    int l = params.m1 + params.m2;
    int n = params.n;
    int max_iter = params.iter;
    double temperature = params.temperature;
    int n_run = params.nrun;
    //bool rc = params.rc;
    //int dmin = params.dmin;
    
    double current_llr     = 0.0;
    double current_ic      = 0.0;
    double best_ic         = 0.0;
    double best_llr        = 0.0;
    double best_ic_in_run  = 0.0;
    double best_llr_in_run = 0.0;
    double avg_llr         = 0.0;
    double avg_ic          = 0.0;

    Array matrix = Array(l, markov.alphabet_size);
    SITES motif;
    SITES best_motif;
    SITES best_motif_in_run;

    // init efficiency only dedicated data
    SamplingData data = SamplingData(sequences, sites, markov, l);

    Is_a_site sites_cache = Is_a_site(sequences, sites);
    FILE *trace = NULL;

    if ((int) sites.size() < n)
    {
        ERROR("parameters do not allow to find a motif");
    }

    for (int run = 0; run < n_run; run++)
    {
        VERBOSE1("processing run %d/%d\n", run+1, n_run);
        NEW_TRACE(trace, run+1);
        //DEBUG("sites=%d n=%d", (int) allsites.size(), n)
        if (params.start_from_sites)
            motif = params.starting_sites;
        else
            motif = random_motif(sites, n, params.nseq, params.dmin);

        current_ic = IC_Markov(motif, sequences, matrix, markov);
        if (max_iter == 0)
        {
            if (current_ic > best_ic)
            {
                best_ic = current_ic;
                if (params.score_type == IC_SCORE || params.score_type == LLR_IC_SCORE)
                    best_motif = motif;
            }
            if (current_llr > best_llr)
            {
                best_llr = current_llr;
                if (params.score_type == LLR_SCORE || params.score_type == PQ_LLR_SCORE)
                    best_motif = motif;
            }
        }
        int iter=-1;
        //VERBOSE2("sampling\n");
        TRACE(trace, "iter\tic\tic_max\tllr\tbest_llr\n");
        while (++iter < max_iter)
        {
            sample_update(sites, sequences, motif, matrix, markov, l, temperature, data, params.nseq, params.dmin, params.score_type);
            if (params.shift)
                motif = shift(motif, sequences, sites, matrix, markov, sites_cache);
            current_llr = llr(motif, sequences, matrix, markov);
            current_ic = IC_Markov(motif, sequences, matrix, markov);
            best_ic_in_run = max(current_ic, best_ic_in_run);
            best_llr_in_run = max(current_llr, best_llr_in_run);

            if (current_ic > best_ic)
            {
                best_ic = current_ic;
                if (params.score_type == IC_SCORE || params.score_type == LLR_IC_SCORE)
                    best_motif = motif;
            }
            if (current_llr > best_llr && params.score_type != PQ_SCORE)
            {
                best_llr = current_llr;
                if (params.score_type == LLR_SCORE || params.score_type == PQ_LLR_SCORE)
                    best_motif = motif;
            }

            VERBOSE3("[%d] ic=%.2f bestic=%.2f llr=%.2f bestllr=%.2f\n", iter, current_ic, best_ic_in_run, current_llr, best_llr_in_run);
            TRACE(trace, "%i\t%.3f\t%.3f\t%.3f\t%.3f\n", iter, current_ic, best_ic_in_run, current_llr, best_llr_in_run);
            
        }
        CLOSE_TRACE(trace);
        avg_llr += best_llr_in_run;
        avg_ic += best_ic_in_run;
        
        if (params.score_type == PQ_SCORE && current_llr >= best_llr)
        {
            best_llr = current_llr;
            best_motif = motif;
        }
        
    }
    Result result;
    result.ic = best_ic;
    result.llr = best_llr;
    result.motif = best_motif;
    result.l = l;
    result.avg_llr = avg_llr / n_run;
    result.avg_ic = avg_ic / n_run;
    return result;
}

/***************************************************************************
 *                                                                         *
 *  FINAL CYCLE
 *                                                                         *
 ***************************************************************************/
// log P(site|M) - log P(site|Bg)
double Score(Site &site, SITES &motif, Sequences &sequences, Array &matrix, Markov &markov, double pseudo=PSEUDO)
{
    // matrix
    freq_matrix(matrix, motif, sequences, markov, pseudo);
    // score
    int j;
    double s = 0.0;
    double logP_m = 0.0;
    for (int i = 0; i < matrix.I; i++)
    {
        j = (int) sequences[site.s][site.p+i];
        logP_m += log(matrix[i][j]);
    }
    s += logP_m - markov.logP(&sequences[site.s][site.p], matrix.I);
    return s;
}

Result collect_sites(Result result, SITES &sites, Sequences &sequences, Markov &markov, Parameters &params)
{
    Array matrix = Array(params.m1 + params.m2, markov.alphabet_size);
    count_matrix(matrix, result.motif, sequences);
    matrix.transform2logfreq(markov);
    params.minspacing   = 0;
    params.maxspacing   = 0;
    params.flanks = 0;

    Result new_result;
    new_result.l = params.m1 + params.m2;
    new_result.motif = matrix_scan(sequences, matrix, markov, params);
    new_result.llr = llr(result.motif, sequences, matrix, markov);
    new_result.ic  = IC_Markov(result.motif, sequences, matrix, markov);
    //print_sites(result.motif);
    return result;
}

/***************************************************************************
 *                                                                         *
 *  MAIN GIBBS
 *                                                                         *
 ***************************************************************************/
void run_sampler(vector<string> &raw_sequences, Sequences &sequences, Markov &markov, Parameters &params)
{
    
    //clock_t start_clock = clock();
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    
    vector<Result> all_results;
    // init random number generator
    srand(SEED);
    // all motifs
    for (int i = 0; i < params.motifs; i++)
    {
        double best_ic = 0.0;
        double best_llr = 0.0;
        Result best_result;
        Result result;
        int l = 0;
        VERBOSE1("starting to sample motif %d.(%d/%d)\n", params.id, i+1, params.motifs);
        // all spacing
        for (int spacing = params.minspacing; spacing <= params.maxspacing; spacing++)
        {
            //VERBOSE1("spacing is set to %d\n", spacing);
            l = params.m1 + spacing + params.m2;
            SITES sites = all_sites(sequences, l, params.m1, params.m2);            
            for (int m = 0; m < (int) all_results.size(); m++)
                sites = mask_motif(sites, all_results[m].motif);

            VERBOSE2("starting [find_one_motif]");
            result = find_one_motif(raw_sequences, sequences, sites, markov, params);
            result.l = l;
            VERBOSE2("end [find_one_motif]");

            // print_motif(result.motif, raw_sequences, sequences, result.l, result.ic, result.llr, params.rc);
            if (params.collect)
            {
                // final cycle
                VERBOSE1("ic before final cycle=%.3f\n", result.ic);
                result = collect_sites(result, sites, sequences, markov, params);
                VERBOSE1("ic after final cycle=%.3f\n", result.ic);
            }

            if (result.ic > best_ic)
            {
                best_ic = result.ic;
                if (params.score_type == IC_SCORE || params.score_type == LLR_IC_SCORE)
                    best_result = result;
            }
            
            if (result.llr > best_llr)
            {
                best_llr = result.llr;
                if (params.score_type == LLR_SCORE || params.score_type == PQ_SCORE || params.score_type == PQ_LLR_SCORE)
                    best_result = result;
            }
        }
        
        all_results.push_back(best_result);
    }


    if (params.id == 1)
    {
        printf("; info-gibbs %d\n", VERSION);
        printf("; %s\n", COMMAND_LINE);
        printf("; title                         %s\n", params.title);
        printf("; started at                    %s", asctime (timeinfo));
        printf("; random seed                   %d\n", SEED);
        printf("; number of runs                %d\n", params.nrun);
        printf("; number of iterations          %d\n", params.iter);
        printf("; sequences                     %d\n", (int) params.nseq);
        printf("; total size in bp              %d\n", (int) sequences.total_size());
        printf("; expected motif occurrences    %d\n", params.n);
        printf("; prior                         a:%.3f|c:%.3f|g:%.3f|t:%.3f\n", markov.priori[0], markov.priori[1], 
                                                                        markov.priori[2], markov.priori[3]);
        printf("; motifs fo find                %d\n", params.motifs);
        printf(";\n");
    }
    
    //clock_t end_clock = clock();
    //printf("; elapsed time (in seconds)     %.3f\n", (end_clock - start_clock) / (double) CLOCKS_PER_SEC);

    for (int m = 0 ; m < (int) all_results.size(); m++)
    {
        printf("; motif                         %d.%d\n", params.id, m + 1);
        printf("; avg.llr                       %.3f\n", all_results[m].avg_llr);
        printf("; avg.ic                        %.3f\n", all_results[m].avg_ic);
        print_motif(all_results[m].motif, raw_sequences, sequences, all_results[m].l, all_results[m].ic, all_results[m].llr, params.rc);
    }
}
