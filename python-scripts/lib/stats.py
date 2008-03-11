"""
:NAME:
stats

:DESC:
Simple statistical utilities
    -mean
    -var (variance)
    -pmf (probability mass function)
    -cdf (cumulative distribution function)
    -table (count table)
    -...
    
"""
# IMPORTS
#-------------------------------------------------------------------------------
from math import *
#-------------------------------------------------------------------------------
# IMPORTS

# CODE
#-------------------------------------------------------------------------------
def pmf_sum(d):
    """
    should return 1.0
    """
    return sum(d.values())


def pmf(values):
    """
    values -- list of values
    Return probability mass function (pmf)
    dict key=value, value=prob
    
    """
    count = {}
    for v in values:
        count[v] = count.get(v, 0) + 1

    N = float(len(values))
    for v in count:
        count[v] /= N

    return count


def pmf_tail(pmf, pvalue, lowerValue=0):
    """
    keep only distribution corresponding to tail >= pvalue
    group other values on 0
    
    """
    v = list(pfm.keys())
    v.sort()
    p = map(pmf.get, v)

    i = len(p)
    pv = 0.0
    while i >  0 and pv < pvalue:
        i -= 1    
        pv += p[i]
    
    if v[i] <= 0:
        raise Error

    dr = {lowerValue: 1.0 - pv + p[i]}
    for j in range(i+1, len(v)):
        dr[v[j]] = p[j]

    return dr
 
    
def cdf(pmf):
    cdf = {}
    cp = 0.0
    values = list(pmf.keys())
    values.sort()
    for v in values:
        p = pmf[v]
        cp += p
        cdf[v] = cp
    return cdf


def cdf_upper(pmf):
    """
    p(x) : P[X>=x]
    """
    cdf = {}
    cp = 1.0
    values = list(pmf.keys())
    values.sort()
    for v in values:
        p = pmf[v]
        cdf[v] = cp
        cp -= p
    return cdf


def probs(pmf, a=None, b=None, lowerProb=None):
    """
    return list of prob in range (a,b+1) 
    default : a,v = min(d.keys()), max(d.keys())
    
    """
    a = a or min(pmf.keys())
    b = b or max(pmf.keys())
    x = []
    for i in range(a, b+1):
        x.append(d.get(i, 0.0))
    return x


########################################
#                                      #
#         Generic functions
#                                      #
########################################

def multi(x, func, *args):
    """
    do y = func(x, *args) for each x
    return list of y if x is list or y 
    """
    if type(x) == list:
        return [func(i, *args) for i in x]    
    else:
        return func(x, *args)


def mean(values):
    return sum(values) / float(len(values))


def var(values):
    m = mean(values)
    return sum( [(x-m)**2 for x in values] ) / float(len(values)-1)


def mean_var(values):
    """
    return tuple mean and var
    """
    return mean(values), var(values)


def table(v, binFunc=None):
    """Count table
    v       -- values
    binFuc  -- binarize function default=round
    return tuple (values, counts)
    """
    binFunc = binFunc or (lambda x: round(x))
    v = map(binFunc, v)
    #v = map(float, v)
    t = {}
    for i in v:
        t[i] = t.get(i, 0) + 1
    l = t.items()
    l.sort()
    return [i[0] for i in l], [i[1] for i in l]

########################################
#                                      #
#         
#                                      #
########################################

def thdist(x_range, step, func, *params):
    """x,y for theorical distribution
    step -- number of point per unit
    example thdist( (-1.0,1.0), 2, norm, a, b) 
    return x,y
    """
    x_range = [i/float(step) for i in range(int(x_range[0]*step), int(x_range[1]*step+1))]
    return x_range, multi(x_range, func, *params)



#-------------------------------------------------------------------------------
# CODE
