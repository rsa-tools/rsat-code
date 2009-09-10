#include "markov.h"
#include "utils.h"

#define SNF               1
#define OLIGO_FREQUENCY   2
#define TRANSITION_MATRIX 3

int bernoulli(Markov &m, double *priori)
{
    m.dealloc();
    m.alloc(0);
    m.S[0] = 1.0;
    for (int i = 0; i < ALPHABET_SIZE; i++)
    {
        m.priori[i] = priori[i];
        m.logpriori[i] = log(priori[i]);
        m.T[0][i] = priori[i];
    }
    return 1;
}

// --
// #INCLUSive Background Model v1.0
// #
// #Order = 1
// #Organism = athaliana
// #Sequences = /users/sista/thijs/scratch/DNA/intergenic.tfa
// #Path = 
// #
// 
// #snf
// 0.3449  0.1581  0.1556  0.3414  
// 
// #oligo frequency 
// 0.3449
// 0.1581
// 0.1556
// 0.3414
// 
// #transition matrix
// 0.3911  0.1516  0.1482  0.3091  
// 0.3760  0.1703  0.1268  0.3269  
// 0.3630  0.1389  0.1671  0.3311  
// 0.2756  0.1678  0.1711  0.3855  
// --

int load_inclusive(Markov &m, string filename)
{
    m.dealloc();

    ifstream f;
    f.open(filename.c_str());
    if (!f) 
    {
        ERROR("unable to open '%s'", filename.c_str());
        return 0;
    }
    string line;
    char const *buffer;
    int order = -1;
    int step = 0;
    int s = 0;
    int t = 0;
    float val[4];
    
    getline(f, line);
    if (strncmp(line.c_str(), "#INCLUSive", 10) != 0)
    {
        f.close();
        return 0;
    }
    
    while (getline(f, line))
    {
        buffer = line.c_str();
        if (line[0] == '#')
        {
            if (sscanf(buffer, "%*s %*s %d", &order) == 1)
                m.alloc(order);
            else if (line == "#snf")
                step = SNF;
            else if (line == "#oligo frequency")
                step = OLIGO_FREQUENCY;
            else if (line == "#transition matrix")
                step = TRANSITION_MATRIX;
        }
        else
        {
            if (order == -1)
                break;
            if (step == SNF)
            {
                if (sscanf(buffer, "%f %f %f %f", &val[0], &val[1], &val[2], &val[3]) == 4)
                {
                    m.priori[0] = (double) val[0];
                    m.priori[1] = (double) val[1];
                    m.priori[2] = (double) val[2];
                    m.priori[3] = (double) val[3];

                    m.logpriori[0] = log(m.priori[0]);
                    m.logpriori[1] = log(m.priori[1]);
                    m.logpriori[2] = log(m.priori[2]);
                    m.logpriori[3] = log(m.priori[3]);
                }
            }
            else if (step == OLIGO_FREQUENCY)
            {
                 if (s < m.msize)
                 {
                    sscanf(buffer, "%f", &val[0]);
                    m.S[s] = (double) val[0];
                    s++;
                 }
            }
            else if (step == TRANSITION_MATRIX)
            {
                if (t < m.msize)
                {
                    if (sscanf(buffer, "%f %f %f %f", &val[0], &val[1], &val[2], &val[3]) == 4)
                    {
                        m.T[t][0] = (double) val[0];
                        m.T[t][1] = (double) val[1];
                        m.T[t][2] = (double) val[2];
                        m.T[t][3] = (double) val[3];                                        
                        t++;
                    }
                }
            }
        }
    }
    f.close();

    // special case when order=0
    if (m.order == 0)
        m.S[0] = 1.0;

    if (m.order == -1)
        return 0;
    else
        return 1;
}
