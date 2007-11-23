
import cache

def sequences2bernoulli(s):
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
    return [A/c, C/c, G/c, T/c]



class MM(object):

    def __init__(self, order=1, pseudo=0.01):
        self.order = order
        self.T = {}
        self.S = {}
        
        self.pseudo = pseudo
        self.priori = [0.25] * 4

        self.P_cached = cache.MemoryCache(self.P)


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
        if order >= 1:
            prefix = sequence[-order:]
            S[prefix] = S.get(prefix, 0) + 1


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
        #compute frequencies
        count = float(sum(S.values()))
        for prefix in S:
            S[prefix] /= count
        for prefix in T:
            c = float(sum(T[prefix].values()))
            for suffix in T[prefix]:
                T[prefix][suffix] /= c


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
        #pseudo = self.pseudo
        assert(len(word) >= order +1)

        p = self.S.get(word[:order], self.pseudo * order)

        for i in range(order, len(word)):
            try:
                p *= T[word[i-order:i]][word[i]]
            except KeyError:
                p *= self.pseudo

        return p
            

def oligo2MM(filename):
    f = open(filename)
    f.readline() #skip first line

    priori = {'A': 0, 'C' : 0, 'G' : 0, 'T' : 0} 
    i = 0
    for line in f:
        w,junk,freq = line.strip().split()
        w = w.upper()
        freq = float(freq)
        # choose markov order
        if i == 0:
            mm = MM(len(w)-1)
        i += 1

        prefix = w[:-1]
        mm.S[prefix] = mm.S.get(prefix, 0) + freq
        #mm.S[w[1:]] = mm.S.get(w[1:], 0) + freq        
        
        mm.T[prefix] = mm.T.get(prefix, {})
        mm.T[prefix][w[-1]] = freq

        #priori
        
        for letter in w:
            priori[letter] += freq

    S = float(sum(priori.values()))

    mm.priori = [priori[b] / S for b in 'ACGT']
    mm.freq()
    return mm
    
    
def test():
    sequences = [ ('A' * 20 + 'T' *10) * 10]
    print sequences
    mm1 = MM(1)
    
    
    mm1.build(sequences)
    print mm1.T
    print mm1.S
    
    word = 'AATT'
    print word , mm1.P(word)
   
    word = 'AAAA'
    print word , mm1.P(word) 
    
if __name__ == '__main__':
    #test()
    mm = oligo2MM('bg.freq')
    print mm.T
    print mm.S
    print mm.priori

    #mm2 = MM(1)


    
    
    
    
    
            