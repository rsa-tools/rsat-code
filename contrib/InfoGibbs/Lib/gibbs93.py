import random
from math import *
#from pylab import *

from Core.Types.list import SortedList
import Core.UI.cli as cli
from motif import *
import rand
import matrix
import cache

def sample(motif, spl, mask, options, EM=False):
    l = options['l']
    alpha = options['alpha']
    T = LLR(motif.words(), motif.bg)
    #myCache = cache.MemoryCache(T.llr)

    prob93 = []
    prob = []
    id = []
    for s in range(len(motif.sequences)):
        sequence = motif.sequences[s]
        for p in range(0, len(sequence) - l + 1):
            id.append( (s,p) )
            word = motif.word( (s,p,l) )
            if (s,p,l) in motif or (s,p,l) in mask or (not EM and (s,p,l) == spl):
                prob.append(-100000000)
            else:
                prob.append( T.llr(word) )

            if options['sampling'] == '93':
                Q = matrix.Q(word, motif.matrix())
                P = motif.bg.P(word)
                prob93.append( Q / P ) 
            
    prob = [ exp(x*alpha) for x in prob]

    if options['sampling'] == 'IC':
        prob = probIC
    elif options['sampling'] == 'llr':
        prob = prob
    elif options['sampling'] == '93':
        prob = prob93
 
    if EM:
        idx = prob.index(max(prob))    
        s, newp = id[idx]
    else:
        s,newp = rand.wchoice(id, prob)()
        
    if not (s,newp,l) in motif:
        motif.add( (s,newp,l) )
    else:
        motif.add( spl )

def update(motif, junk):
    return motif.pop()

def EM(motif, n, mask, parameters):
    #final EM
    count = 0
    param = {}
    param.update(parameters)

    for i in range(n):
        count += 1
        s,p,l = update(motif, k)
        sample(motif, (s,p,l), mask, param, True)
        motif = shift(motif)
    return motif

def gibbs(sequences, motif, bg, mask, parameters):
    SHIFT_COUNT = parameters['shift']
    EM_CYCLE = parameters['EM']
    tab = parameters['tab']
    info = parameters['info']
    nrun = parameters['nrun']
    #
    # INIT
    #
    #motif = init(sequences, bg, parameters)

    ibest = 0
    best = motif.copy()
    if parameters['best'] is None:
        parameters['best'] = best

    count = 0
    k = 0
    while count < parameters['iter'] :#and count - ibest < 100:
        count += 1
        if SHIFT_COUNT != 0 and count % SHIFT_COUNT == 0:
            motif = shift(motif)

        s,p,l = update(motif, k)
        k = (k + 1) % len(motif.sequences)
    
        sample(motif, (s,p,l), mask, parameters)

        if motif.fitness() > best.fitness():
            best = motif.copy()
            ibest = count

            if best.fitness() > parameters['best'].fitness():
                 parameters['best'] = best

        info('[%d] %d best=%.2f %s N=%d | IC=%.2f %s N=%d | IC=%.2f %s N=%d' % (nrun, count, parameters['best'].fitness(), parameters['best'].consensus(), parameters['best'].N(), 
            best.fitness(), best.consensus(), best.N(), motif.fitness(), motif.consensus(), motif.N()) )
        tab.append( dict(bestIC=best.fitness(), IC=motif.fitness(), word=motif.consensus(), 
                            bestWord=best.consensus(), i=count, time=parameters['clock'].get_value(), N=motif.N(), bestN=best.N() ))


    #EM
    best = EM(best, EM_CYCLE, mask, parameters)
    info('[%d] %d best=%.2f %s N=%d | IC=%.2f %s N=%d | IC=%.2f %s N=%d' % (nrun, count, parameters['best'].fitness(), parameters['best'].consensus(), parameters['best'].N(), 
            best.fitness(), best.consensus(), best.N(), motif.fitness(), motif.consensus(), motif.N()) )    
    tab.append( dict(bestIC=best.fitness(), IC=motif.fitness(), word=motif.consensus(), 
                            bestWord=best.consensus(), i=count, time=parameters['clock'].get_value(), N=motif.N(), bestN=best.N() ))

    return best



if __name__ == '__main__':
    pass
