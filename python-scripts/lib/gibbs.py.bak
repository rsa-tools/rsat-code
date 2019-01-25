"""
:NAME:
gibbs

:DESC:
Gibbs sampling algorithm (site sampling)

"""
import random
import bisect
from copy import copy
from math import *

import cli
import matrix

########################################
#
# UTILS
#
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


########################################
#
# TRACE
#
########################################
#cli.TRACE = 0
DIR = ''
traceic = cli.Trace('ic')
#trace2 = Trace('score')
#traceT = Trace('temp')
#tracetaboo = Trace('taboo')

########################################
#
# MOTIF __str__
#
########################################
def motif_sites2str(motif, strand='+'):
    str = []
    nsequences = len(motif.sequences)
    sites_list = list(motif)
    sites_list.sort()
    for (s,p) in sites_list:
        if strand == '+-' and s >= nsequences / 2:
            strandlabel = '-'
            seq = s - nsequences / 2 + 1
        else:
            strandlabel = '+'
            seq = s+1
        if hasattr(motif, 'm1'): # dyad
            word = motif.word((s,p))
            word = word[:motif.m1] + '-' * motif.spacing + word[-motif.m2:]
        else:
            word = motif.word((s,p))
        str += [ '; %d\t%s\t%d\t%s' % (seq,strandlabel,p, word) ]
    return str


def motif2fullstr(motif, sequences, bg, strand='+'):
    nsites = len(motif)
    nsequences = len(sequences)
    words = Motif.words(motif)
    mymatrix = motif.matrix()
    consensus = motif.consensus()
    l = len(mymatrix)
    IC = matrix.Iseq(mymatrix, bg.priori)

    str = []
    str += ['; prior                        \t' + 'a:%.4f|c:%.4f|g:%.4f|t:%.4f' % tuple(bg.priori) ]
    #str += ['; log.likelihood.ratio         \t%.3f' % (matrix.llr(words, mymatrix, bg.priori)) ]
    str += ['; total.information            \t%.3f' % (IC) ]
    str += ['; information.per.column       \t%.3f' % (IC / l) ]
    str += ['; sites                        \t%d' % len(motif) ]
    str += ['; consensus.IUPAC              \t%s' % consensus ]
    str += [';' ]
    #str += ['; seq   \tsequence (starting from 1)']
    #str += ['; strand\tforward/reverse (+/-)']
    #str += ['; pos   \tposition in sequence (starting from 0)' ]
    #str += [';#seq\tstrand\tpos\tword' ]
    str += motif_sites2str(motif, strand)
    str += [ matrix.matrix2tab(matrix.words2countmatrix(words, bg.priori), count=True) ]
    return '\n'.join(str)


########################################
#
# MOTIF / SITE
#
########################################
class Motif(set):
    def __init__(self, l, sequences, bg, *args):
        self.l = l
        self.sequences = sequences
        self.bg = bg
        self.dmin=0 #mininal distance between 2 motifs
        super(Motif, self).__init__(*args)

    def __copy__(self):
        m = Motif(self.l, self.sequences, self.bg, iter(self))
        m.dmin = self.dmin
        return m

    def copy(self):
        return self.__copy__()

    def words(self):
        return [self.sequences[s][p:p+self.l] for (s,p) in self]
        
    def word(self, site):
        (s,p) = site
        return self.sequences[s][p:p+self.l]

    def matrix(self):
        return matrix.words2matrix(self.words(), self.bg.priori)

    def consensus(self):
        return matrix.consensus(self.matrix(), self.bg.priori)

    def IC(self):
        return matrix.Iseq(self.matrix(), self.bg.priori) # / float(len(words[0]))

    def __str__(self):
        mymatrix = self.matrix()
        consensus = matrix.consensus(mymatrix, self.bg.priori)
        l = len(mymatrix)
        IC = matrix.Iseq(mymatrix, self.bg.priori)
        return 'IC=%.3f\t%s' % (IC, consensus)


class Dyad(Motif):
    def __init__(self, m1, m2, spacing, *args, **kwargs):
        self.l = m1 + spacing + m2
        self.m1 = m1
        self.m2 = m2
        self.spacing = spacing
        super(Dyad, self).__init__(self.l, *args, **kwargs)
    
    def __copy__(self):
        m = Dyad(self.m1, self.m2, self.spacing, self.sequences, self.bg, iter(self))
        m.dmin = self.dmin
        return m

    def words(self):
        return [ self.sequences[s][p:p+self.m1] + self.sequences[s][p+self.l-self.m2:p+self.l] for (s,p) in self ]        
    
    def word(self, site):
        (s,p) = site
        return self.sequences[s][p:p+self.m1] + self.sequences[s][p+self.l-self.m2:p+self.l]

    def consensus(self):
        dyadconsensus = super(Dyad, self).consensus()
        return dyadconsensus[0:self.m1] + ('n{%s}' % self.spacing) + dyadconsensus[self.m1:]


'''
class Site(object):
    __slots__ = ['s', 'p', 'l']

    def __init__(self, s, p, l):
        self.s = s
        self.p = p
        self.l = l
    
    def __str__(self):
        return '%d %d %d' % (self.s, self.p, self.l)
'''


########################################
#
# MOTIF UTILS
#
########################################
def shifted(motif, sites, delta):
    new_motif = motif.copy()
    new_motif.clear()

    for (s,p) in motif:
        if (s,p+delta) in sites:
            new_motif.add( (s,p+delta) )
        else:
            new_motif.add( (s,p) )
    return new_motif


def shift(motif, sites):
    best = motif
    for delta in [-1, 1]:
        new_motif = shifted(motif, sites, delta)
        if new_motif.IC() > best.IC():
            best = new_motif
    return best
        

def find_occurrences(motif, max=100):
    '''
    return a sorted list of motif occurrences (sites)
    sort positions given PM/PB
    '''
    ids = []
    probs = []
    sequences = motif.sequences
    l = motif.l
    mymatrix = motif.matrix()
    for s in range(len(sequences)):
        for p in range(0, len(sequences[s]) - l + 1):
            ids += [ (s,p) ]
            word = motif.word( (s,p) )
            probs += [ matrix.Q(word, mymatrix) / motif.bg.P(word) ]

    L = zip(probs, ids)
    L.sort(reverse=True)
    return [i[1] for i in L[:max]]


########################################
#
# SCORING
#
########################################
def score_IC(motif, sites):
    words = motif.words()
    probs = [0.0] * len(sites)
    N = len(words) + 1
    beta = 1 / T

    for i in range(len(sites)):
        word = motif.word(sites[i])
        mymatrix = matrix.words2matrix(words + [word], motif.bg.priori)
        probs[i] = exp(beta * matrix.Iseq(mymatrix, motif.bg.priori) * N)
    return probs


def score_93(motif, sites):
    mymatrix = motif.matrix()
    probs = [0.0] * len(sites)
    for i in range(len(sites)):
        #word = sequences[site.s][site.p:site.p+site.l]
        word = motif.word(sites[i])
        probs[i] = matrix.Q(word, mymatrix) / motif.bg.P(word)
    return probs

I = 0

#T = 0.90

def score_llr(motif, sequences, bg, sites):
    global I, T
    words = motif.words()
    llr = matrix.LLR(words, bg).llr
    #sequences[site.s][site.p:site.p+m1] + sequences[site.s][site.p+site.l-m2:site.p+site.l]
    #return [exp(llr(sequences[site.s][site.p:site.p+site.l])) for site in sites]
    #return [exp(llr(motif.word(site))) for site in sites]
    I += 1
    #T = 1.6 - motif.IC() / 6.0
    #T = 1.0 - 0.2 * (motif.IC() / 6.0)
    #T = 0.8

    beta = 1 / T
    #beta = 1
    l =  [exp(beta*llr(motif.word(site))) for site in sites]
    #print max(l) / min(l)
    #print l
    #print len(l)
    #traceT.writeln(T, True)
    #trace2.write_values(l, True)
    return l



########################################
#
# ALGORITHM
#
########################################
def final_cycle(motif):
    m = motif.copy()
    m.clear()
    F = []

    min_occ = len(motif)
    max_occ = len(motif) * 4

    occurrences = find_occurrences(motif, max_occ)
    bestIC = 0.0
    
    for site in occurrences:
        m.add(site)
        ic = m.IC()
        F += [ic]

        if len(m) >= min_occ and ic >= bestIC:
            bestIC = ic
            best = m.copy()

    #print F
    #import pylab
    #pylab.plot(range(1,len(F)+1), F)
    #pylab.show()
    return best

maxchoosen = 0

def sample(motif, sites, score, EM):
    global allsites, maxchoosen
    sites_list = list(sites)
    #print list(motif)
    #print sorted(sites_list)

    if score == 'llr':
        probs = score_llr(motif, motif.sequences, motif.bg, sites_list) 
    elif score == 'IC':
        probs = score_IC(motif, sites_list) 
    elif score == '93':
        probs = score_93(motif, sites_list) 

    if EM:
        imax = probs.index(max(probs))
        newsite = sites_list[imax]
    else:
        #Gibbs
        #imax = probs.index(max(probs))
        #maxite = sites_list[imax]

        #newsite = wchoice(sites_list, probs)
        #print probs
        #if sorted(list(maxite)) == sorted(list(newsite)):
        #    maxchoosen += 1
        return wchoice(sites_list, probs)
    return newsite


def random_site(sites):
    return random.choice(list(sites))

def run_EM(motif, sites, max_iteration, score, EM):
    bestIC = motif.IC()
    for i in range(max_iteration):
        site = motif.pop()
        #site = random_site(motif)
        #motif.remove(site)
        #sites.remove(site)

        noov_sites = no_overlapping(sites, motif)
        newite = sample(motif, noov_sites, score, EM)

        if newite not in motif:
            motif.add(newite)
        else:
            motif.add(site)


        if motif.IC() <= bestIC:
            break
        bestIC = motif.IC()
            
        #print 
        #print motif.IC()
    return motif


maxIC = 2.0

#T = 0.95
def run_gibbs(motif, sites, max_iteration, score, EM):
    global taboo, tabooCount, T, allsites, ICs, maxchoosen, maxIC

    N = len(motif)
    iteration = 0
    best_ic = 0
    bestmotif = motif

    info = cli.Info(max_iteration, verbosity=4)
    allsites = sites
    while iteration < max_iteration:
        iteration += 1

        #sites.difference_update(motif)
        lastIC = motif.IC()
        site = random_site(motif)
        motif.remove(site)        

        noov_sites = no_overlapping(sites, motif)
        #newite = sample(motif, noov_sites, score, EM)
        f = sample(motif, noov_sites, score, EM)

        newite = f()
        #for x in range(len(motif)):
        #    motif.pop()


        #sites.add(site)
        #T = T - 0.001
        if newite == site or newite in motif:
            motif.add(site)
            #EM = 0
        else:
            #EM = 1
            #T = 0.9
            #tracetaboo.writeln('%d\t%f' % (iteration, tabooCount/float(iteration)))
            if newite not in motif:
                motif.add(newite)
            else:
                motif.add(site)

        #multi update 
        for x in range(N-1):
            oldsite = motif.pop()
            newsite2 = f()
            if newsite2 not in motif:
                motif.add(newsite2)
            else:
                motif.add(oldsite)

        #if  random.random() < 0.01:
        motif = shift(motif, sites)

        #current_llr = matrix.llr(motif.words(), motif.matrix(), motif.bg.priori)
        current_llr = 0.0
        ic = motif.IC()
        if ic >= best_ic:
            bestmotif = motif.copy()
            best_ic = ic

        
        #def PT(IC, T):
        #    return exp(IC/T) / sum([exp(IC/T) for IC in ICs])


        IC = ic #* len(motif)           

        '''

        T1 = 0.6
        T2 = 0.95

        PT1T2 = 1 / (1+PT(IC, T1) / PT(IC, T2)) 
        #print PT1T2
        if T == T1 and random.random() < PT1T2:
            T = T2
        elif T == T2 and random.random() < (1-PT1T2):
            T = T1
        '''


        #traceic.writeln('%.2f\t%.2f\t%.3f' % (ic, best_ic, T), 1)
        info('ic=%.2f bestIC=%.2f t=%.4f llr=%.3f' % (ic, best_ic, T, current_llr/N) )

    return bestmotif


########################################
#
# INIT
#
########################################
def all_sites(sequences, bg, l):
    sites = set()
    for s in range(len(sequences)):
        sequence = sequences[s]
        for p in range(0, len(sequence) - l + 1):
            site = (s,p)
            sites.add(site)
    return sites


def no_overlapping(sites, motif):
    '''
    return no overlapping sites (distance(motif1, motif2) >= motif.dim)
    '''
    if motif.dmin == 0:
        return sites

    mysites = copy(sites)
    for (s,p) in motif:
        pmin = max(0, p-motif.dmin+1)
        pmax = min(len(motif.sequences[s])-1, p+motif.dmin-1)
        for pos in range(pmin, pmax+1):
            mysites.discard( (s,pos) )
    return mysites


#def is_overlapping(motif, site, l):
#    print list(motif)
#    print site
#    (s,p) = site
#    for d in range(1, l+1):
#        if (s,p-d) in motif:
#            return True
#    return False

    
def random_sites(motif, sites, N):
    '''
    Add random sites to motif
    '''
    i = 0
    while len(motif) < N:
        i += 1
        if i > N*10: #avoid infinite loop
            break
        sites = no_overlapping(sites, motif)
        if len(sites) == 0:
            break
        mysite = random_site(sites)
        motif.add(mysite)
    return motif


def matrix2motif(m, l, N, sequences, bg):
    motif = Motif(l, sequences, bg)
    for site in find_occurrences(motif)[:N]:
        motif.add( site )
    return motif


########################################
#
# RUN
#
########################################
def run_trials(startmotif, sequences, bg, sites, N, gibbsiterations, EMiterations, trials, score):
    global traceic, ICs

    TRACE_PREFIX = 'ic'

    bestIC = 0.0
    bestmotif = startmotif
    info = cli.Info(trials, 1)

    for t in range(trials):
        ICs = set()
        info('run=%d' % (t+1))
        
        traceic = cli.Trace('%s.%d.ic' % (TRACE_PREFIX, t+1) )
        traceic.writeln('#iter\tic\tbestic\tT')
        m = startmotif.copy()
        
        #random motif
        if len(startmotif) == 0:
            random_sites(m, sites, N)

        gibbsmotif = run_gibbs(m, sites, gibbsiterations, score, False)
        EMmotif = run_gibbs(gibbsmotif, sites, EMiterations, score, True)
        ic = EMmotif.IC()
        if ic >= bestIC:
            bestIC = ic
            bestmotif = EMmotif.copy()
        
    return bestmotif


def run_dyad(sequences, startmotif, bg, length, N, gibbsiterations, EMiterations, nmotifs=1, ntrials=1, score='llr', strand='+', spacing=None, finalcycle=False):
    motifs = []

    for i in range(nmotifs):
        bestmotif = None
        bestIC = 0.0
        for space in range(spacing[0], spacing[1]+1):
            sites = all_sites(sequences, bg, length*2+space)
            [ sites.difference_update(m) for m in motifs ] #mask already found motif
            startmotif = Dyad(length, length, space, sequences, bg)
            motif = run_trials(startmotif, sequences, bg, sites, N, gibbsiterations, EMiterations, ntrials, score)

            if motif.IC() > bestIC:
                bestIC = motif.IC()
                bestmotif = motif

        if finalcycle:
            bestmotif = final_cycle(bestmotif)
        motifs += [ bestmotif ]

    return [motif2fullstr(motif, sequences, bg, strand) for motif in motifs]


def run(sequences, startmotif, bg, length, N, gibbsiterations, EMiterations, nmotifs=1, ntrials=1, score='llr', strand='+', dmin=0, finalcycle=False):
    #bg = markov.MM(0)
    #bg.learn(sequences) 

    sites = all_sites(sequences, bg, length)
    startmotif = startmotif or Motif(length, sequences, bg)
    startmotif.dmin = dmin

    motifs = []
    for i in range(nmotifs):
        bestmotif = run_trials(startmotif, sequences, bg, sites, N, gibbsiterations, EMiterations, ntrials, score)
        if finalcycle:
            bestmotif = final_cycle(bestmotif)
        motifs += [ bestmotif ]
        sites.difference_update(bestmotif) #mask found motif

    return [motif2fullstr(motif, sequences, bg, strand) for motif in motifs]




########################################
#
# TESTS
#
########################################

def test_gibbs():
    import markov
    l = 3
    sequences = ['CATTGG', 'TTAGTG']
    bg = markov.MM(0)
    bg.learn(['ATGCGCTCGGGCAGGCGATGGCTTGG']) 
    #sp1 = set([(0,1), (1,2)])
    #sp2 = set([(0,2), (1,1)])

    m = set()

    m.add((1,2))
    m.add((0,2))    
    m.add((0,1))    
    #print m
    
    sites = all_sites(sequences, bg, l)
    bestmotif = run_gibbs(m, sequences, bg, sites, 10)

    print motif2str(bestmotif, sequences, bg, '+-')


if __name__ == '__main__':
    test_gibbs()
    
    
