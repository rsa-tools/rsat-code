import os
import tempfile
from math import log
from copy import copy

from numpy import array, zeros, shape, log as Nlog, sum as Nsum
#from Numeric import array, zeros, shape, transpose, log as Nlog, sum

PSEUDO = 1 
SEQLOGO = '/Users/matthieu/Recherche/Software/weblogo/seqlogo'

def words2weblogo(words, outdir='.', bin=SEQLOGO, type='PNG', title='', filename='logo'):
    junk, seqfilename = tempfile.mkstemp()
    f = open(seqfilename, 'w')
    f.write('\n'.join(words))
    f.close()
    print '%s -F %s -Y -c  -n -t "%s" -f %s > %s/%s.%s' % (bin, type, title, seqfilename, outdir, filename, type)
    os.system('%s -F %s -Y -c  -n -t "%s" -f %s > %s/%s.%s' % (bin, type, title, seqfilename, outdir, filename, type.lower()))
    os.remove(seqfilename)
    

def words2countmatrix(words, priori):
    """
    Convert a list of words to a simple count matrix
    """
    w = len(words[0])
    m = zeros( (w, 4), 'f')
    N = len(words)
    for i in range(w):
        for s in range(N):
            letter = words[s][i]
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


def words2countmatrix_bis(words):
    """
    Convert a list of words to a simple count matrix
    """
    w = len(words[0])
    m = zeros( (w, 4), 'f')
    N = len(words)
    for i in range(w):
        for s in range(N):
            letter = words[s][i]
            if letter == 'A':
                m[i][0] = m[i][0] + 1.0
            elif letter == 'C':
                m[i][1] = m[i][1] + 1.0
            elif letter == 'G':
                m[i][2] = m[i][2] + 1.0  
            elif letter == 'T':
                m[i][3] = m[i][3] + 1.0

    return m


def words2matrix(words, priori, pseudo=PSEUDO):
    """
    Convert a list of words to a frequency matrix using priori probabilities p
    m_{b,i} = ( f_{b,i} + p_i ) / ( N + 1) [Hertz 1999]
    """

    c = words2countmatrix(words, priori)
    w, B = shape(c)
    f = zeros( (w, B), 'f')
    
    for i in range(w):
        N = float(sum(c[i]))
        for b in range(B):
            f[i][b] = (c[i][b] + priori[b] * pseudo) / (N + pseudo)
    return f


def matrix2tab_slides(matrix, title=''):
    w = shape(matrix)[0]

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


def matrix2tab(matrix, title=''):
    w = shape(matrix)[0]

    str = [ '; %s' % (title) ]
    for b in range(4):
        if b == 0:
            line = [ 'a |' ]
        if b == 1:
            line = [ 'c |' ]
        if b == 2:
            line = [ 'g |' ]
        if b == 3:
            line = [ 't |' ]

        for i in range(w):
            line += ['%.4f' % matrix[i][b] ]
        str += [ '\t'.join(line) ]

    return '\n'.join(str)


def consensus_05(matrix, fmin=0.5):
    w = shape(matrix)[0]

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

def consensus(matrix, priori, mask=True):
    w = shape(matrix)[0]

    str = []
    for i in range(w):
        letter = ''
        for b in range(4):
            if matrix[i][b] != 0.0 and log(matrix[i][b]  / priori[b]) >= 0:
                letter = letter + 'ACGT'[b]
        if mask and len(letter) == 4 or len(letter) == 3 or len(letter) == 2 or len(letter) == 0:
            letter = 'n'
        else:
            letter = CONSENSUS[letter]

        #if len(letter) != 1:
        #    letter = '[' + letter + ']'

        str.append(letter)
    return ''.join(str)


def info_content(matrix):
    w,B = matrix.shape
    S = 0.0
    for i in range(w):
        s = 0
        for b in range(B):
            if matrix[b][i] != 0.0:
                s += -matrix[b][i] * (log(matrix[b][i]) / log(2))
        S += (2.0 - s)
    return  S


def Q(word, matrix):
    w = matrix.shape[0]
    P = 1.0
    for i in range(w):
        if word[i] == 'A':
            P *= matrix[i][0]
        elif word[i] == 'C':
            P *= matrix[i][1]
        elif word[i] == 'G':
            P *= matrix[i][2]
        elif word[i] == 'T':
            P *= matrix[i][3]
    return P


def Iseq(matrix, priori):
    """
    Relative Entropy [Hertz 1999]
    (~ log-likelihood ratio / -N)
    Iseq = \sum_{i=1}^w \sum_{b=1}^4 f_{b,i} ln f_{b,i} / p_b

    matrix -- frequency matrix
    p      -- priori probability [0.25, 0.25, 0.25, 0.25]
    """

    return sum(sum(matrix * Nlog( matrix / priori)))

def Iseq_ref(matrix, priori):
    """
    Relative Entropy [Hertz 1999]
    (~ log-likelihood ratio / -N)
    Iseq = \sum_{i=1}^w \sum_{b=1}^4 f_{b,i} ln f_{b,i} / p_b

    matrix -- frequency matrix
    p      -- priori probability [0.25, 0.25, 0.25, 0.25]
    """
    I = 0.0
    w, B = matrix.shape
    for i in range(w):
        for b in range(B):
            if matrix[i][b] != 0.0:
                I += matrix[i][b] * log(matrix[i][b] / priori[b])

    return I


def llr_old(words, matrix, priori):
    """
    llr [MEME]
    llr = log (Pr(sites | motif) / Pr(sites | back)) 
    Approximated by Iseq * N
    """
    return Iseq(matrix, priori) * (len(words))


def llr_ref(words, matrix, priori):
    S = 0.0
    for w in words:
        for i in range(len(w)):
            letter = w[i]
            if letter == 'A':
                S += log ( matrix[i][0] / priori[0] )
            elif letter == 'C':
                S += log ( matrix[i][1] / priori[1] )
            elif letter == 'G':
                S += log ( matrix[i][2] / priori[2] )            
            elif letter == 'T':
                S += log ( matrix[i][3] / priori[3] )
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

LETTER2J = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3, 'N' : 4}

class LLR_old(object):
    """
    llr =   Sum_k log P(w^k | M) - log P(w^k | B)
    llr = Sum_i Sum_k log P(w_i^k | M) - Sum_k P(w^k | B)
          ----------------------------   ----------------
                alpha                        beta
    
    llr =~ 1/N Sum_i Sum_j C_ij /n log C_ij /n - Sum_k log P(w^k | B)
    """
    def __init__(self, words, mm, pseudo=PSEUDO):
        l = len(words[0])
        p = mm.priori
        self.l = l
        self.N = len(words)
        C = words2countmatrix(words, priori)

        #
        # Compute alpha
        #
        T = [ [0.0]*4 for i in range(l) ]
        

        for i in range(l):
            N = sum(C[i])
            for j in range(4):
                alpha = 0.0
                for w in words:
                    if w[i] == 'N':
                        continue
                    J = LETTER2J[w[i]]
                    if j == J:
                        alpha += log ( (C[i][J] + 1 + p[J] * pseudo) / (N + 1 + pseudo) )
                    else:
                        alpha += log ( (C[i][J] + p[J] * pseudo) / (N + 1 + pseudo) )                            
                alpha += log ( (C[i][j] + 1 + p[j] * pseudo) / (N + 1 + pseudo) )

                T[i][j] = alpha

        #
        # Compute beta
        #
        beta = 0.0
        for w in words:
            beta += log(mm.P(w))


        self.T = T
        self.beta = beta
        self.mm = mm
        
    def llr(self, word):
        alpha = 0
        for i in range(self.l):
            if word[i] != 'N':
                alpha += self.T[i][LETTER2J[word[i]]]

        return alpha - self.beta - log(self.mm.P(word))
        

    def IC(self, word):
        return self.llr(word) / float(self.N + 1) / self.l



class LLR(object):
    """
    llr =   Sum_k log P(w^k | M) - log P(w^k | B)
    llr = Sum_i Sum_k log P(w_i^k | M) - Sum_k P(w^k | B)
          ----------------------------   ----------------
                alpha                        beta
    
    llr =~ 1/N Sum_i Sum_j C_ij /n log C_ij /n - Sum_k log P(w^k | B)
    """
    def __init__(self, words, mm, pseudo=PSEUDO):
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
            beta += log(mm.P(w))


        self.T = T
        self.beta = beta
        self.mm = mm
        
    def llr(self, word):
        alpha = 0
        i = 0
        while i < self.l: 
            alpha += self.T[i][LETTER2J[word[i]]]
            i += 1

        return alpha - self.beta - log(self.mm.P_cached(word))


    def IC(self, word):
        return self.llr(word) / float(self.N + 1) / self.l


def llr(words, mm, alpha=0.1):
    matrix = words2matrix(words, mm.priori, alpha)
    #print matrix2tab(matrix)
    return llr_markov(words, matrix, mm)


if __name__ == '__main__':

    
    import markov
    mm = markov.MM(0)
    mm.learn( [ ('ATGCTGGGCTAGGATGCGTGAGGCGTGGATACCGATTCG') * 10])
    priori = mm.priori
    #print priori
    
    #print llr(['AAAA', 'AAAT', 'NNNN'], mm, 0.0001)

    #sys.exit()
    
    
    words = ['AATTGAT', 'AANGTTA', 'GNAAAAA']

    print '->', llr(words , mm, 0.0001)
    print '--'
    
    T = LLR(['AATTGAT', 'AANGTTA'], mm, 0.0001)
    print T.llr('GNAAAAA')





