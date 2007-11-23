import random
from math import *

from Core.Types.list import SortedList
import Core.UI.cli as cli

from motif import *
import matrix
import cache
#from pylab import *



def sample(motif, spl, e, N, parameters):
    l = parameters['l']
    n = float(motif.N())

    newmotif = Motif(motif.sequences, motif.bg)

    W = 0.8
    Z = W / (1-W)
    m = len(spl)
    pm = (m + e * Z) / (N + N * Z)
    pbg = (N-m + (N-e) * Z) / (N + N*Z)

    p = 1.0*e / N 
    L = range(len(motif.sequences))
    random.shuffle(L)
    for s in L:
        sequence = motif.sequences[s]
        pos = range(0, len(sequence) - l + 1)
        for i in pos:
            word = motif.word( (s,i,l) )
    
            Q = matrix.Q(word, motif.matrix())
            P = motif.bg.P(word)
            
            #p = pm
            pmotif = Q * p
            pbackground = P *(1-p)

            #print pmotif, pbackground
            if random.random() < pmotif / (pmotif + pbackground) and newmotif.N() < parameters['maxwords']:
                newmotif.add((s,i,l))

    if len(newmotif.spl) < parameters['minwords']:
        return motif

    return newmotif

def gibbs(sequences, motif, bg, parameters):
    
    SHIFT_COUNT = parameters['shift']
    tab = parameters['tab']
    info = cli.Info(parameters['iter'], 1)
    p = parameters['p']
    #motif = init(sequences, bg, parameters)
    ibest = 0
    #bestList = SortedList(1, False)
    best = motif.copy()
    count = 0
    k = 0
    
    while count < parameters['iter']: #  and count - ibest < 100:
        count += 1
        
        if SHIFT_COUNT != 0 and count % SHIFT_COUNT == 0:
            motif = shift(motif)

        motif = sample(motif, motif.spl, parameters['e'], parameters['N'], parameters)

        if motif.fitness() > best.fitness():
            #if motif.N() == parameters['words'] or parameters['alpha'] == 0.0:
            best = motif.copy()
            ibest = count
            
        info('%d best=%.2f %s N=%d | IC=%.2f %s N=%d' % (count, best.fitness(), best.consensus(), best.N(), motif.fitness(), motif.consensus(), motif.N()) )
        tab.append( dict(bestIC=best.fitness(), IC=motif.fitness(), word=motif.consensus(), 
                                bestWord=best.consensus(), i=count, time=parameters['clock'].get_value(), N=motif.N(), bestN=best.N() ))

    return best



if __name__ == '__main__':
    pass
