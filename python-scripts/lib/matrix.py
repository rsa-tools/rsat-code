"""
:NAME:
matrix

:DESC:
PSSM (Position Specific Scoring Matrix) utilities.
A matrix is defined as a bi-dimensional array of flaot.
Dimension 0 defines position.
Dimension 1 defines letter in alphabet.

"""

import os
import tempfile
import bisect
import random
from math import log
from copy import copy

PSEUDO = 1.0

LETTER2J = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3, 'N' : 4}
ALPHABET = ['A', 'C', 'G', 'T']

SEQLOGO = '/Users/matthieu/Recherche/Software/weblogo/seqlogo'


########################################
#                                      #
# SEQLOGO
#                                      #
########################################
def words2weblogo(words, outdir='.', bin=SEQLOGO, type='PNG', title='', filename='logo'):
    junk, seqfilename = tempfile.mkstemp()
    f = open(seqfilename, 'w')
    f.write('\n'.join(words))
    f.close()
    print '%s -F %s -Y -c  -n -t "%s" -f %s > %s/%s.%s' % (bin, type, title, seqfilename, outdir, filename, type)
    os.system('%s -F %s -Y -c  -n -t "%s" -f %s > %s/%s.%s' % (bin, type, title, seqfilename, outdir, filename, type.lower()))
    os.remove(seqfilename)

########################################
#                                      #
# MATRIX
#                                      #
########################################
def words2countmatrix(words, priori):
    """
    Convert a list of words to a simple count matrix
    """
    w = len(words[0])
    m = [ [0.0] * 4 for i in range(w) ]
    N = len(words)
    for i in range(w):
        for s in range(N):
            try:
                letter = words[s][i]
            except IndexError:
                letter = 'N'
            if letter == 'A':
                m[i][0] = m[i][0] + 1.0
            elif letter == 'C':
                m[i][1] = m[i][1] + 1.0
            elif letter == 'G':
                m[i][2] = m[i][2] + 1.0  
            elif letter == 'T':
                m[i][3] = m[i][3] + 1.0
            elif letter == 'N':
                m[i][0] = m[i][0] + priori[0]
                m[i][1] = m[i][1] + priori[1]
                m[i][2] = m[i][2] + priori[2]
                m[i][3] = m[i][3] + priori[3]
    return m

def words2matrix(words, priori, pseudo=None):
    """
    Convert a list of words to a frequency matrix using priori probabilities p
    m_{b,i} = ( f_{b,i} + p_i ) / ( N + 1) [Hertz 1999]
    """
    pseudo = pseudo or PSEUDO
    c = words2countmatrix(words, priori)
    w, B = len(c), len(c[0])
    f = [ [0.0] * len(ALPHABET) for i in range(w) ]
    
    for i in range(w):
        N = float(sum(c[i]))
        for b in range(B):
            f[i][b] = (c[i][b] + priori[b] * pseudo) / (N + pseudo)
    return f

########################################
#                                      #
# MATRIX INPUT/OUTPUT
#                                      #
########################################
def matrix2txt(matrix, title=''):
    w = len(matrix)

    str = [ '; %s' % (title) ]
    for b in range(4):
        if b == 0:
            line = [ 'A|' ]
        if b == 1:
            line = [ 'C|' ]
        if b == 2:
            line = [ 'G|' ]
        if b == 3:
            line = [ 'T|' ]

        for i in range(w):
            line += ['%.2f' % matrix[i][b] ]
        str += [ '   '.join(line) ]

    return '\n'.join(str)

def matrix2tab(matrix, title='', count=False):
    w = len(matrix)

    if title != '':
        str = [ '; %s' % (title) ]
    else:
        str = []

    for b in range(len(ALPHABET)):
        if b == 0:
            line = [ 'a |' ]
        if b == 1:
            line = [ 'c |' ]
        if b == 2:
            line = [ 'g |' ]
        if b == 3:
            line = [ 't |' ]

        for i in range(w):
            if count:
                line += ['%d' % matrix[i][b] ]
            else:
                line += ['%.3f' % matrix[i][b] ]
        str += [ '\t'.join(line) ]

    return '\n'.join(str)

def tab2matrix(f):
    m = None
    for line in f:
        line = line.strip()
        if line.startswith(';') or line == '':
            continue
        if line.startswith('//'):
            break

        elements = line.split('\t')
        letter = elements[0].strip()[0].upper()
        
        if (elements[1] == '|'):
            elements = elements[2:]
        else:
            elements = elements[1:]

        if m is None:
            l = len(elements)
            m = [ [0.0] * len(ALPHABET) for i in range(l) ]
        J = LETTER2J[letter]
        for i in range(len(elements)):
            m[i][J] = float(elements[i])

    return m

def tab2matrices(f):
    l = []
    while True:
        m = tab2matrix(f)
        if m == None:
            return l
        else:
            l.append(m)

########################################
#                                      #
# RANDOM SITES
#                                      #
########################################
def wchoice(l, frequencies):
    """
    l -- list
    frequences -- associated unnormalized frequencies
    return a weighted choice function

    """
    assert len(l) == len(frequencies)
    S = 0.0
    cdf = []
    for f in frequencies:
        S += f
        cdf += [ S ]
    return lambda : l[bisect.bisect(cdf, random.random() * S)]

def random_site_generator(matrix):
    choose = []
    length = len(matrix)
    for i in range(length):
        choose += [ wchoice(ALPHABET, matrix[i]) ]
        
    while 1:
        site = ['N'] * length
        for i in range(length):
            site[i] = choose[i]()
        yield ''.join(site)


########################################
#                                      #
# CONSENSUS
#                                      #
########################################
def consensus_05(matrix, fmin=0.5):
    w = len(matrix)

    str = []
    for i in range(w):
        if matrix[i][0] > fmin:
            str += [ 'A' ]
        elif matrix[i][1] > fmin:
            str += [ 'C' ]
        elif matrix[i][2] > fmin:
            str += [ 'G' ]
        elif matrix[i][3] > fmin:
            str += [ 'T' ]
        else:
            str += [ 'n' ]

    return ''.join(str)

'''
A                       (Adenine) 
C                       (Cytosine)
G                       (Guanine)
T                       (Thymine)
R       = A or G        (puRines)
Y       = C or T        (pYrimidines)
W       = A or T        (Weak hydrogen bonding)
S       = G or C        (Strong hydrogen bonding)
M       = A or C        (aMino group at common position)
K       = G or T        (Keto group at common position)
H       = A, C or T     (not G)
B       = G, C or T     (not A)
V       = G, A, C       (not T)
D       = G, A or T     (not C)
N       = G, A, C or T  (aNy)
'''
CONSENSUS = { 'A' : 'A',
              'C' : 'C',
              'G' : 'G',
              'T' : 'T',
              'AG' : 'R',
              'CT' : 'Y',
              'AT' : 'W',
              'CG' : 'S',
              'AC' : 'M',
              'GT' : 'K',
              'ACT' : 'H',
              'CGT' : 'B',
              'ACG' : 'V',
              'AGT' : 'D',
              'ACGT' : 'N'
              }
EXPLAIN_CONSENSUS = dict([(v,k) for k,v in CONSENSUS.items()])

def explain_consensus(w):
    str = []
    for letter in w:
        letter = EXPLAIN_CONSENSUS[letter]
        if len(letter) == 4:
            letter = 'N' 
        elif len(letter) != 1:
            letter = '[' + letter + ']'
        str.append(letter)
    return ''.join(str)

def consensus(matrix, priori, mask=False):
    w = len(matrix)

    str = []
    for i in range(w):
        letter = ''
        for b in range(4):
            if matrix[i][b] != 0.0 and log(matrix[i][b]  / priori[b]) >= 0:
                letter = letter + 'ACGT'[b]
        if mask and (len(letter) == 4 or len(letter) == 3 or len(letter) == 2 or len(letter) == 0):
            letter = 'n'
        else:
            letter = CONSENSUS[letter]

        #if len(letter) != 1:
        #    letter = '[' + letter + ']'

        str.append(letter)
    return ''.join(str)

########################################
#                                      #
# INFORMATION CONTENT / LLR
#                                      #
########################################
def logo_IC(matrix):
    '''
    seqlogo information content
    
    '''
    w = len(matrix)
    S = 0.0
    for i in range(w):
        s = 0
        for b in range(len(ALPHABET)):
            if matrix[i][b] != 0.0:
                s += -matrix[i][b] * (log(matrix[i][b]) / log(2))
        S += (2.0 - s)
    return  S

def Q(word, matrix):
    w = len(matrix)
    P = 1.0
    for i in range(w):
        P *= matrix[i][LETTER2J[word[i]]]
    return P

def Iseq(matrix, priori):
    """
    Relative Entropy [Hertz 1999]
    Iseq = \sum_{i=1}^w \sum_{b=1}^4 f_{b,i} ln f_{b,i} / p_b

    matrix -- frequency matrix
    p      -- priori probability [0.25, 0.25, 0.25, 0.25]
    """
    I = 0.0
    w, B = len(matrix), len(matrix[0])
    for i in range(w):
        for b in range(B):
            if matrix[i][b] != 0.0:
                I += matrix[i][b] * log(matrix[i][b] / priori[b])

    return I

def llr(words, matrix, priori):
    S = 0.0
    for w in words:
        for i in range(len(w)):
            j = LETTER2J[w[i]]
            S += log ( matrix[i][j] / priori[j] )
    return S

def llr_markov(words, matrix, mm):
    S = 0.0
    for w in words:
        for i in range(len(w)):
            if w[i] == 'N':
                continue
            S += log ( matrix[i][LETTER2J[w[i]]] )

        S += -log(mm.P(w))
    return S

def Iseq_markov(words, matrix, mm):
    return llr_markov(words, matrix, mm) / len(words) / len(words[0])


class LLR(object):
    """
    llr =   Sum_k log P(w^k | M) - log P(w^k | B)
    llr = Sum_i Sum_k log P(w_i^k | M) - Sum_k P(w^k | B)
          ----------------------------   ----------------
                alpha                        beta
    
    llr =~ 1/N Sum_i Sum_j C_ij /n log C_ij /n - Sum_k log P(w^k | B)
    """

    def __init__(self, words, mm, pseudo=None):
        pseudo = pseudo or PSEUDO
        l = len(words[0])
        p = mm.priori
        self.l = l
        self.N = len(words)
        C = words2countmatrix(words, p)

        #
        # Compute alpha
        #
        T = [ [0.0]*5 for i in range(l) ]
        for i in range(l):
            N = sum(C[i])
            for j in range(4):
                alpha = 0.0
                for w in words:
                    if w[i] == 'N':
                        continue
                    J = LETTER2J[w[i]]
                    if j == J :
                        alpha += log ( (C[i][J] + 1 + p[J] * pseudo) / (N + 1 + pseudo) )
                    else:
                        alpha += log ( (C[i][J] + p[J] * pseudo) / (N + 1 + pseudo) )

                # add prob for the new letter
                alpha += log ( (C[i][j] + 1 + p[j] * pseudo) / (N + 1 + pseudo) )
                T[i][j] = alpha
            # SPECIAL 'N' CASE
            alpha = 0.0
            for w in words:
                if w[i] == 'N':
                    continue
                J  = LETTER2J[w[i]]
                alpha += log ( (C[i][J] + p[J] + p[J] * pseudo) / (N + 1 + pseudo) )
            T[i][4] = alpha

        #
        # Compute beta
        #
        beta = 0.0
        for w in words:
            #beta += log(mm.P(w))
            beta += mm.logP(w)

        self.T = T
        self.beta = beta
        self.mm = mm
        
        #print T
        #print beta

    def llr(self, word):
        alpha = 0
        i = 0
        while i < self.l: 
            alpha += self.T[i][LETTER2J[word[i]]]
            i += 1
        #print alpha, self.beta
        return alpha - self.beta - self.mm.logP(word)

    def IC(self, word):
        return self.llr(word) / float(self.N + 1) / self.l


#def llr(words, mm, alpha=0.1):
#    matrix = words2matrix(words, mm.priori, alpha)
#    #print matrix2tab(matrix)
#    return llr_markov(words, matrix, mm)


########################################
#                                      #
# TESTS
#                                      #
########################################
'''
def test0():
    mm = markov.MM(0)
    mm.learn( [ 'ATGCTGGGCTAGGATGCGTGAGGCGTGGATACCGATTCG' * 10])
    priori = mm.priori
    words = ['AATTGAT', 'AANGTTA', 'GNAAAAA']
    print '->', llr(words , mm, 0.0001)
    print '--'
    
    T = LLR(['AATTGAT', 'AANGTTA'], mm, 0.0001)
    print T.llr('GNAAAAA')
'''

def test_matrix():
    markovmodel = markov.MM(0)
    markovmodel.learn(['ATGCGCTCGGGCAGGCGATGGCTTGG'])

    mymatrix = [ [0.2, 0.2, 0.3, 0.3], [0.2, 0.2, 0.3, 0.3] ]
    mycountmatrix = [ [2, 2, 4, 4], [4, 2, 4, 2] ]
    priori = [0.25, 0.25, 0.25, 0.25]
    words = ['AAA', 'TTT']

    print words2countmatrix(['AAA', 'TTT'], priori)
    print words2matrix(['AAA', 'TTT'], priori, 0.2)
    print matrix2tab(mymatrix, title='test matrix')
    print tab2matrix('simple.mat')
    print explain_consensus('ACCWNATC')
    print consensus(mycountmatrix, priori, False)

    print Q('ATT', mymatrix)
    print Iseq(mymatrix, priori)

    matrix = words2matrix(words, priori, 0.2)
    print llr_markov(words, matrix, markovmodel)
    print Iseq_markov(words, matrix, markovmodel)
    
    #myllr = LLR(words, markovmodel, 0.1)
    #print myllr.llr('AAA')
    #print myllr.IC('AAA')    


if __name__ == '__main__':
    test_matrix()


