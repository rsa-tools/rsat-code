"""
:NAME:
markov

:DESC:
Simple markov chain utilities

:HELP:
m = MM(1) # order=1
m.learn(['ATGGGATT', 'GTTGGTGGGTTGCT])
print m.P('ATGT')


"""

import gzip
from math import log
from dna import reverse_complement

ALPHABET = 'ACGT'

ALPHABET_SIZE = len(ALPHABET)
ALPHABET2index = dict(zip(list(ALPHABET), range(ALPHABET_SIZE)))


def sequences2bernoulli(s, pseudo=1.0):
    c = 0.0
    A = 0
    C = 0
    G = 0
    T = 0

    for site in s:
        c += len(site)
        A += site.count('A')
        C += site.count('C')
        G += site.count('G')
        T += site.count('T')
    return [(A+pseudo)/(c+pseudo*4), (C+pseudo)/(c+pseudo*4), (G+pseudo)/(c+pseudo*4), (T+pseudo)/(c+pseudo*4)]


class MM(object):
    """
    -Stationnary
    -Maximum Likelyhood parameters estimation with pseudo count
                       C(suffix|prefix) + pseudo
    P(suffix|prefix) = -------------------------
                         C(prefix) + N pseudo 

    N: number of suffixes
    P(suffix|prefix) (example P(C|AA))
    C(suffix|prefix): prefix-suffix count (example AAC)
    C(prefix): prefix count (example AA)
    """

    def __init__(self, order=0, pseudo=1.0):
        self.order = order
        self.T = {}
        self.S = {}
        self.Tpseudo = {} # prob when count = 0
        self.Spseudo = 0  # prob when count = 0
        self.pseudo = float(pseudo) #pseudo count
        self.priori = [1.0/ALPHABET_SIZE] * ALPHABET_SIZE

    def _count(self, sequence):
        order = self.order
        T = self.T
        S = self.S
        for i in range(len(sequence) - order):
            prefix = sequence[i:i+order]
            
            # skip N
            if prefix.find('N') != -1 or sequence[i+order] == 'N':
                continue
            
            t = T.get(prefix, {})
            t[sequence[i+order]] = t.get(sequence[i+order], 0) + 1
            T[prefix] = t
            S[prefix] = S.get(prefix, 0) + 1

        #end of sequence
        #if order >= 1:
        #    prefix = sequence[-order:]
        #    S[prefix] = S.get(prefix, 0) + 1


    def learn(self, sequences):
        """
        sequences -- list of strings
        """
        assert(type(sequences) is list)
        self.T = T = {}
        self.S = S = {}
        
        for s in sequences:
            self._count(s)
        self.freq()

        #bernoulli
        self.priori = sequences2bernoulli(sequences)
    
    def freq(self):
        T = self.T
        S = self.S
        pseudo = self.pseudo
        order = self.order
        Tpseudo = self.Tpseudo
        #compute frequencies
        count = float(sum(S.values()))
        self.Spseudo = pseudo / (count + pseudo * ALPHABET_SIZE**order)
        for prefix in S:
            S[prefix] = (S[prefix] + pseudo) / (count + pseudo * ALPHABET_SIZE**order)
        for prefix in T:
            c = float(sum(T[prefix].values()))
            Tpseudo[prefix] = pseudo / (c + pseudo * ALPHABET_SIZE) 
            for suffix in T[prefix]:
                T[prefix][suffix] = (T[prefix][suffix] + pseudo) / (c + pseudo * ALPHABET_SIZE) 

    def logP(self, word):
        return log(self.P(word))

    def P_bernoulli(self, word):
        T = self.T
        p = 1.0
        for letter in word:
            if letter != 'N':
                p *= T[''][letter]
        return p
                
    def P(self, word):
        if self.order == 0:
            return self.P_bernoulli(word)
        order = self.order
        T = self.T
        assert(len(word) >= order +1)

        p = self.S.get(word[:order], self.Spseudo)
        for i in range(order, len(word)):
            try:
                p *= T[word[i-order:i]][word[i]]
            except KeyError:
                try:
                    p *= self.Tpseudo[word[i-order:i]]
                except KeyError:
                    p *= self.priori[ ALPHABET2index[word[i]] ]

        return p



def oligo2MM(filename):
    """Load data from oligo-analysis formated file (can be gzipped)
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)

    rc = 0
    priori = {'A': 0, 'C' : 0, 'G' : 0, 'T' : 0} 
    i = 0
    for line in f:
        if line.startswith('#') or line.startswith(';'):
            if line.find('grouped by pairs of reverse complements') > 0:
                rc = 1
            continue
        elements = line.strip().split()
        w, freq, count = elements[0], float(elements[2]), int(elements[3])
        w = w.upper()

        # choose markov order
        if i == 0:
            mm = MM(len(w) - 1, pseudo=0.0)
            #mm.order = len(w) - 1
        i += 1

        if rc:
            wrc = reverse_complement(w)
            prefix = wrc[:-1]
            if w != wrc:
                freq = freq / 2.0
                mm.S[prefix] = mm.S.get(prefix, 0) + freq
                mm.T[prefix] = mm.T.get(prefix, {})
                mm.T[prefix][wrc[-1]] = freq
                for letter in wrc:
                    priori[letter] += freq
                #priori[prefix] += freq

        prefix = w[:-1]
        mm.S[prefix] = mm.S.get(prefix, 0) + freq
        mm.T[prefix] = mm.T.get(prefix, {})
        mm.T[prefix][w[-1]] = freq
        #priori
        for letter in w:
            priori[letter] += freq
        #priori[prefix] += freq

    S = float(sum(priori.values()))
    mm.priori = [priori[b] / S for b in ALPHABET]
    mm.freq()
    #print mm.priori
    #print mm.order
    #print mm.S
    #print mm.T
    return mm        


class MMError(MM):
    """
    Markov Model for motif like AANNNTT

    """
    def __P(self, word):
        """return P(word) in model
        """
        # MM0
        if self.order == 0:
            p = 1.0
            for letter in word:
                if letter != 'N':
                    p *= self.T[''][letter]
            return p
        # MM >= 1
        words = Bio.sequence.IUPAC2list(word)
        return sum([super(MMError, self).P(w) for w in words]) 

    def set_NExtension(self, NExtension):
        self.NExtension = NExtension

    def P(self, word, NExtension=(1,1)):
        """
        NExtension can be set with set_NExtension
        """
        if word.find('N') == -1:
            return super(MMError, self).P(word)
        else:
            NExtension = NExtension or self.NExtension
            return sum([self.__P(iw) for iw in Bio.sequence.extendn(word, NExtension)])


def test():
    sequences = [ ('A' * 20 + 'T' *10) * 10]
    #print sequences
    mm1 = MM(1)
    mm1._count('ATGGGTTT')
    #mm1.freq()

    mm1 = MM(1)
    mm1.learn(sequences)
    print mm1.T
    print mm1.S
    
    word = 'AATT'
    print word , mm1.P(word)
    print word , mm1.logP(word)
   
    word = 'AAAA'
    print word , mm1.P(word) 

    mm0 = MM(0)
    mm0.learn(sequences)
    mm0.P_bernoulli('ATT')
    
    mm0.load_oligo_file('oligo.freq')
    
    mm = oligo2MM('oligo.freq')
    print mm.P('ATGTGGTGTTC')
    
if __name__ == '__main__':
    test()


    
    
    
    
    
            