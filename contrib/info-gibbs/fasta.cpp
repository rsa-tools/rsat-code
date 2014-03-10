#include "fasta.h"

/***************************************************************************
 *                                                                         *
 *  CONVERT & PARSE SEQUENCES
 *                                                                         *
 ***************************************************************************/

Sequences convert_sequences(vector<string> &raw_sequences)
{
    Sequences sequences = Sequences((int) raw_sequences.size());
    for (int s = 0; s < sequences.size(); s++)
    {
        int len = (int) raw_sequences[s].size();
        sequences[s].data = new int[len];
        sequences[s].len = len;
        for (int i = 0; i < sequences[s].size(); i++)
        {
            switch (raw_sequences[s][i])
            {
                case 'A':
                sequences[s][i] = 0;
                break;
                
                case 'C':
                sequences[s][i] = 1;                
                break;

                case 'G':
                sequences[s][i] = 2;
                break;

                case 'T':
                sequences[s][i] = 3;
                break;
                
                default:
                sequences[s][i] = -1;
            }
        }
    }
    return sequences;
}

// compute priori probabilities (pA, pC, pG, pT)
double *compute_priori(Sequences &sequences)
{
    double *priori = new double[ALPHABET_SIZE]; 
    int total = 0;
    
    for (int i = 0; i < ALPHABET_SIZE; i++)
    {
        priori[i] = 0;
    }
    
    for (int s = 0; s < sequences.size(); s++)
    {
        for (int p = 0; p < sequences[s].size(); p++)
        {
            if (sequences[s][p] != -1)
            {
                priori[sequences[s][p]]++;
                total++;
            }
        }
    }
        
    for (int i = 0; i < ALPHABET_SIZE; i++)
    {
        priori[i] /= total;
    }
    return priori;
}

/***************************************************************************
 *                                                                         *
 *  CONVERT & PARSE SEQUENCES
 *                                                                         *
 ***************************************************************************/
string reverse_complement(string s)
{
    string rc = string(s);
    int l = s.size();

    for (int i = 0; i < l; i++)
    {
        if (s[i] == 'A')
            rc[l-i-1] = 'T';
        else if (s[i] == 'T')
            rc[l-i-1] = 'A';
        else if (s[i] == 'C')
            rc[l-i-1] = 'G';
        else if (s[i] == 'G')
            rc[l-i-1] = 'C';
        else
            rc[l-i-1] = s[i];
    }
    return rc;
}

int read_fasta(vector<string> &sequences, char *filename, bool rc=false)
{
    vector<string> titles;
    int i = 0;

    // read sequences
    istream *f;
    ifstream ifn;

    if (filename == NULL)
    {
        f = &cin;
    }
    else
    {
        ifn.open(filename);
        f = &ifn;
    }

    if (!f)
        return 0;

    string line;
    while (getline(*f, line))
    {
        if (line[0] == '>')
        {
            sequences.push_back(string());
            titles.push_back(string(line.substr(1)));
        }
        else
        {
            if (sequences.size() > 0)
                sequences.back() += line;
        }
    }
    
    if (filename != NULL)
        ifn.close();

    // convert sequences
    int n = sequences.size();
    for (int k = 0; k < n; k++)
    {
        // uppercase
        int length = sequences[k].length();
        for (i = 0; i < length; i++)
            sequences[k][i] = toupper(sequences[k][i]);
            
        // reverse complement
        if (rc) 
        {
            sequences.push_back(reverse_complement(sequences[k]));
            titles.push_back(titles[k] + " rc");        
        }
    }
    
    if (n == 0)
        return 0;
    else
        return 1;
}
