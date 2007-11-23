import random
from copy import copy
from matrix import *


class Mask(set):
    pass

class Motif(object):

    def __init__(self, sequences, bg, spl=None):
        self.sequences = sequences
        self.bg = bg # Background
        self.spl = spl or set() #sequence, position, length
        self._fitness = None
        self._matrix = None
        self._words = None
        self._consensus = None
        self._N = 0 #number of words
        self._l = None
        
    def words(self):
        if self._words == None:
            words = []
            for s,p,l in self.spl:
                words += [ self.sequences[s][p:p+l] ]
            self._words = words
        return self._words

    def word(self, spl):
        s,p,l = spl
        return self.sequences[s][p:p+l]

    def consensus(self):
        if self._consensus == None:
            self._consensus = consensus(self.matrix(), self.bg.priori)
            #self._consensus = consensus_05(self.matrix(), 0.8)
        return self._consensus

    def l(self):
        if self._l == None:
            self._l = len(self._words[0])
        return self._l
        
    def fitness(self):
        if self._fitness == None:
            words = self.words()
            #self._fitness = Iseq_ref(self.matrix(), self.bg.priori) / float(len(words[0]))
            self._fitness = Iseq(self.matrix(), self.bg.priori) / float(len(words[0]))
        return self._fitness

    def matrix(self):
        if self._matrix == None:
            self._matrix = words2matrix(self.words(), self.bg.priori)        
        return self._matrix

    def N(self):
        return len(self.spl)

    def __contains__(self, spl):
        return spl in self.spl

    def pop(self):
        self._fitness = None
        self._words = None
        self._matrix = None
        self._consensus = None
        self._N -= 1
        return self.spl.pop()
        
    def add(self, spl):
        if self._l == 0:
            self._l = spl[-1]
        self.spl.add(spl)
        self._fitness = None
        self._words = None
        self._matrix = None
        self._N += 1
        self._consensus = None

    def remove(self, spl):
        self.spl.remove(spl)    
        self._fitness = None
        self._words = None
        self._matrix = None
        self._N -= 1
        self._consensus = None        

    def discard(self, spl):
        self.spl.discard(spl)    
        self._fitness = None
        self._words = None
        self._matrix = None
        self._consensus = None
        self._N -= 1


    #def update(self):
    #    self.spl.sort()
    #    self._fitness = None
    #    self.fitness()

    def __eq__(self, other):
        return self.sequences == other.sequences and self.spl == other.spl

    def __cmp__(self, other):
        #print self.fitness(), other.fitness(), cmp(self.fitness(), other.fitness())
        if self.spl == other.spl:
            return 0
        else:
            return cmp(self.fitness(), other.fitness())
        
    def copy(self):
        return Motif(self.sequences, self.bg, copy(self.spl))        

    def str(self):
        title = [ 'IC=%.3f' % (self.fitness()) ]
        title += ['N=%d' % len(self.spl) ]
        title += ['Consensus=%s' % self.consensus() ]    
        return matrix2tab(self.matrix(), '\t'.join(title))

    def simple_str(self):
        return 'IC=%.3f %s N=%d' % (self.fitness(), self.consensus(), self.N())

    def full_str(self):
        words = self.words()
        matrix = self.matrix()
        str = [';']
        str += ['; prior                 \t' + 'a:%f|c:%f|g:%f|t:%f' % tuple(self.bg.priori) ]
        str += ['; total.information     \t%f ' % (Iseq_ref(self.matrix(), self.bg.priori)) ]
        str += ['; information.per.column\t%f' % (Iseq_ref(self.matrix(), self.bg.priori) / self.l()) ]
        str += ['; sites                 \t%d' % len(self.spl) ]
        str += ['; consensus.IUPAC       \t%s' % self.consensus() ]
        str += [';' ]

        str += ['; seq\tpos\tword']
        for s,p,l in self.spl:
            str += [ '; %i\t%i' % (s,p) + '\t%s' % self.sequences[s][p:p+l] ]
        str += [ ';' ]
        #str += [ matrix2tab(self.matrix(), '') ]
        str += [ matrix2tab(words2countmatrix(self.words(), self.bg.priori)) ]
        return '\n'.join(str)
        

def shifted(i, x):
    newI = Motif(i.sequences, i.bg)
    for spl in i.spl:
        s,p,l = spl
        if p+x <= len(i.sequences[s])-l and p+x >= 0:
            newI.add( (s,p+x,l) )
        else:
            newI.add( (s,p,l) )

    return newI
        
        
def shift(motif):
    best = motif
    for x in [-2,-1,1,2]:
        m = shifted(motif, x)
        if motif.fitness() > best.fitness():
            best = m
    return best


def choose_word(i, l):
    s = random.randint(0, len(i.sequences)-1)
    
    if len(i.sequences[s]) >= l:
        p = random.randint(0, len(i.sequences[s])-l)
    else:
        return None

    return s,p,l

def init(sequences, bg, mask, parameters):
    l = parameters['l']
    words = parameters['words']

    m = Motif(sequences, bg)
    while(m.N()) < words:
        w = choose_word(m, l)
        if w != None and w not in mask:
            m.add(w)
    return m


def words2motif(words, sequences, bg, parameters):
    """
    construct motif with given set of words
    """
    N = parameters['words']
    m = Motif(sequences, bg)
    for word in set(words):
        l = len(word)
        found = []
        for s in range(len(sequences)):
            sequence = sequences[s]
            start = 0
            while True: 
                i = sequence.find(word, start)
                if i == -1:
                    break
                start = i+1
                found += [(s,i,l)]
            
        count = words.count(word)
        for spl in random.sample(found, min(len(found), count)):
            m.add(spl)
    
    while(m.N()) < N:
        spl = choose_word(m, l)
        m.add(spl)
    return m
            



        