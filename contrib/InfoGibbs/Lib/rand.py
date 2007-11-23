import random
import bisect
#from pylab import *
#from gamdlib import *

def word_prob(word, priori):
    P = 1.0
    for letter in word:
        if letter == 'A':
            P *= priori[0]
        elif letter == 'C':
            P *= priori[1]
        elif letter == 'G':
            P *= priori[2]
        elif letter == 'T':
            P *= priori[3]
    return P
        

#def normalize(t):
#    tmax = float(max(t))
#    tmin = min(t)
#    return [ (i - tmin) / (tmax - tmin) for i in t ]


def wchoice(l, frequencies):
    """
    l -- list
    frequences -- associated unnormalized frequencies
    return a weighted choice function

    """
    assert(len(l) == len(frequencies))
    S = 0.0
    cdf = []
    for f in frequencies:
        S += f
        cdf += [ S ]
    #plot(range(len(cdf)), [i/float(S) for i in cdf])
    #print cdf, S
    return lambda : l[bisect.bisect(cdf, random.random() * S)]
