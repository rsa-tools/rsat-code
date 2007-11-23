import random
from math import *
import Core.UI.cli as cli
from motif import *
import rand
import matrix
import cache

#from pylab import *

def final(motif, options):
    l = options['l']

    # sort positions given PM/PB
    id = []
    prob = []
    for s in range(len(motif.sequences)):
        sequence = motif.sequences[s]
        for p in range(0, len(sequence) - l + 1):
            id.append( (s,p,l) )
            word = motif.word( (s,p,l) )
            PM = matrix.Q(word, motif.matrix())
            PB = motif.bg.P(word)
            prob.append( PM / PB )

    L = zip(prob, id)
    L.sort(reverse=True)

    # compute IC
    m = Motif(motif.sequences, motif.bg)
    F = []
    max = -1
    for junk, spl in L[:100]:
        m.add(spl)
        f = m.fitness()
        if f >= max:
            max = f
            best = m.copy()
        F += [f]

    #plot(range(1,len(F)+1), F)
    #show()
    return best
