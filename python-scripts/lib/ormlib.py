"""
:NAME:
ormlib

:DESC:
ORM library

:BUG:
count_tree Ntab end of seq (PARTIAL FIX WHEN spacing=(1:1) )

:HELP:
H -> hash table key=oligonucleotide value=list of occurrence position
N -> scanned bases (count) per position for each word size

"""

import sys
from math import log10
import bisect
import gzip
#import sets

import cli
import dna
import markov
import ST
from dna import reverse_complement, find_location
from dist import pbinom, ppois, pbinom_right_left, pbinom_left
from tab import load_oligo_file

def list_unique(seq):
    keys = {}
    for e in seq:
        keys[e] = 1
    return keys.keys()

##
#
# cache
#
##
#pbinom_right_left_cached = Core.cache.MemoryCache(pbinom_right_left)
#pbinom_cached = Core.cache.MemoryCache(pbinom)
#ppois_cached = Core.cache.MemoryCache(ppois)

##############################################################################
#                                                                            #
#                                 COUNT WORDS
#                                                                            #
##############################################################################
def count_words_hash(sequences, l, searchLocation, strand = '+-', overlap=False):
    """Count each word of length l in sequences
    l               -- oligonucleotide length
    searchLocation  -- location tuple example (-200,-1)
    strand          -- + or +- 
    overlap         -- allow auto-overlapping

    return N, H
    """
    location = find_location(sequences)
    H = {} #hash table key=oligonucleotide value=list of occurrence position
    N = {} #scanned bases per position for each word size
    N[l] = [0] * (searchLocation[1] - searchLocation[0] + 1) 
    scannedPositions = 0
    scannedWords = 0
    info = cli.Info(len(sequences), 1, 1)

    for s in sequences:
        info('Counting words in [%+05d:%+05d]' % (searchLocation[0], searchLocation[1]))
        dna = s.sequence
        HS = {} # for current sequence

        a,b = max(searchLocation[0], s.location[0]), min(searchLocation[1], s.location[1] - l+1)
        
        for I in range(a, b+1):
            scannedPositions += 1
            i = I - s.location[0]

            if dna[i:i+l].find('N') >= 0:
                continue

            if strand == '+':
                w = dna[i:i+l]
            elif strand == '+-':
                wf = dna[i:i+l]
                wr =  reverse_complement(dna[i:i+l])
                w = min(wf, wr)

            N[l][I-searchLocation[0]] += 1

            if not overlap and HS.get(w, [a-l])[-1] + l > I:
                continue

            H.setdefault(w, []).append( I )
            HS.setdefault(w, []).append( I )
            scannedWords += 1
    return dict(N=N, H=H, scannedPositions=scannedPositions, scannedWords=scannedWords)


def count_dyads_hash(sequences, l, spacing, searchLocation, strand = '+-', overlap=False):
    """Count each dyad of length l in sequences
    l               -- oligonucleotide length
    spacing         -- spacing
    searchLocation  -- location tuple example (-200,-1)
    strand          -- + or +- 
    overlap         -- allow auto-overlapping

    """
    lmonad = l
    l = 2*lmonad + spacing

    location = find_location(sequences)
    H = {}
    N = {}
    N[l] = [0] * (searchLocation[1] - searchLocation[0] + 1) 
    scannedPositions = 0
    scannedWords = 0
    info = cli.Info(len(sequences), 1, 1)

    for s in sequences:
        info('Counting words in [%+05d:%+05d]' % (searchLocation[0], searchLocation[1]))
        dna = s.sequence
        HS = {}
        a,b = max(searchLocation[0], s.location[0]), min(searchLocation[1], s.location[1] - l+1)

        for I in range(a, b+1):
            scannedPositions += 1
            i = I - s.location[0]
            if dna[i:i+l].find('N') >= 0 :
                continue

            if strand == '+':
                w = dna[i:i+lmonad] + 'N' * spacing + dna[i+lmonad+spacing:i+2*lmonad+spacing]
            elif strand == '+-':
                wf = dna[i:i+lmonad] + 'N' * spacing + dna[i+lmonad+spacing:i+2*lmonad+spacing]
                wr =  reverse_complement(wf)
                w = min(wf, wr)

            N[l][I-searchLocation[0]] += 1
            if not overlap and HS.get(w, [a-l])[-1] + l > I:
                continue

            H.setdefault(w, []).append( I )
            HS.setdefault(w, []).append( I )
            scannedWords += 1
    return dict(N=N, H=H, scannedPositions=scannedPositions, scannedWords=scannedWords)


def count_words_tree(sequences, l, searchLocation, strand = '+-', overlap=False, error=0, spacing=(1,1)):
    """Count each word of length l in sequences
    l               -- oligonucleotide length
    searchLocation  -- location tuple example (-200,-1)
    strand          -- + or +- 
    overlap         -- allow auto-overlapping

    return N, H
    """
    location = find_location(sequences)
    #H = {} #hash table key=oligonucleotide value=list of occurrence position
    N = {} #scanned base count per position fo each word size
    N[l] = [0] * (searchLocation[1] - searchLocation[0] + 1) 
    scannedPositions = 0
    scannedWords = 0
    info = cli.Info(len(sequences), 1, 1)

    #
    # construct SuffixTree
    #
    st = ST.SuffixTree(maxDepth=l, overlapping=overlap, maxIUPAC=error, NExtension=spacing, storePosition=True)

    for s in sequences:
        info('Counting words in [%+05d:%+05d]' % (searchLocation[0], searchLocation[1]))
        a,b = max(searchLocation[0], s.location[0]), min(searchLocation[1], s.location[1] - l+1)       
        dna = s.get_dna((a,b+l+1))
        st.add_dna(dna, shift=a)
        for I in range(a,b+1):
            i = I - s.location[0]
            w = dna[i:i+l]
            if w.find('N') >= 0:
                continue
            N[l][I-searchLocation[0]] += 1
 
        #@DEBUG
        #ST.display(st.root, maxDepth=6, full=1)
          
    #
    # Count
    #

    #@DEBUG
    #keys = st.extract(minLength=l, maxLength=l).keys()
    #keys.sort()
    #print '\n'.join(keys)

    C = st.extract(minLength=l, maxLength=l)

    if strand == '+-':
        H = ST.get_positions_two_strands(C, overlap)
    else:
        H = ST.get_positions(C)

    return dict(N=N, H=H, scannedPositions=scannedPositions, scannedWords=scannedWords)


def count_words(sequences, l, searchLocation, strand = '+-', overlap=False, params=None):
    if params['count'] == 'hash' and params['dyad']:
        all = dict(N={}, H={}, scannedWords=0, scannedPositions=0)
        for spacing in range(params['spacing'][0], params['spacing'][1]+1):
            r = count_dyads_hash(sequences, l, spacing, searchLocation, strand=strand, overlap=overlap)
            all['N'].update(r['N'])
            all['H'].update(r['H'])
            all['scannedWords'] += r['scannedWords']
            all['scannedPositions'] += r['scannedPositions']
        return all
        
    elif params['count'] == 'hash' and not params['dyad']:
        return count_words_hash(sequences, l, searchLocation, strand=strand, overlap=overlap)
    elif params['count'] == 'tree':
        return count_words_tree(sequences, l, searchLocation, strand=strand, overlap=overlap, error=params['error'], spacing=params['spacing'])



##############################################################################
#                                                                            #
#                                 EXTRACT WINDOWS
#                                                                            #
##############################################################################
def get_window(location, a, b, size):
    """
    return window centered on a,b of width size (or at least)
    """
    center = a + (b-a) / 2.0
    start = max(location[0], center - int(size / 2.0))
    start = min(start, location[1] - size + 1)
    return int(start), int(start + size - 1)


def extract_windows(sequences, bg, location=None, wl=None, params=None):
    location = location or bg.location
    r = count_words(sequences, bg.l, location, strand=bg.strand, overlap=bg.overlap, params=params)

    R = []    
    H = r['H']
    N = r['N']
    #print H
    wl = wl or H.keys()    
    info = cli.Info(len(wl), 1, 1)

    for w in wl:
        info('Extracting windows')

        # if params['heuristic'] == 'score':
        #     x = extract_windows_score(N, H, bg, location=location, w=w, params=params)
        # else:
        extractor = Extractor(N, H, bg, location=location, w=w, params=params)
        x = extractor.run()
        R.extend(x)

    #print pbinom_cached.stats()
    return {'R': R, 'scannedWords' : r['scannedWords']}


class Extractor:
    def __init__(self, N, H, bg, location, w, params):
        R = []
        l = H[w]
        l.sort()

        self.MIN_WIDTH = params['width'][0]
        self.MAX_WIDTH = params['width'][1]
        
        self.MIN_OCC = params['occ'][0]
        self.MAX_OCC = params['occ'][1]

        self.l = l
        self.R = R
        self.N = N
        self.H = H
        self.bg = bg
        self.location = location
        self.w = w
        self.params = params

    def __next(self, a, b, obsOcc):
        width = b-a+1    
        if width < self.MIN_WIDTH or width > self.MAX_WIDTH:
            return
        if obsOcc < self.MIN_OCC or obsOcc > self.MAX_WIDTH:
            return
    
        n = sum(self.N[len(self.w)][a-self.location[0]:b-self.location[0]+1])
        try:
            obsFreq = obsOcc / float(n)
        except:
            cli.warning('n error')

        expFreq = self.bg.freq(self.w, (a,b))
        expOcc = expFreq * n

        #pv = ppois(obsOcc, expOcc)
        #pv = ppois_cached(obsOcc, expOcc)
        #pv = pbinom_right_left_cached(obsOcc, n, expFreq)
        if self.params['under']:
            pv = pbinom_left(obsOcc, n, expFreq)
        else:
            pv = pbinom(obsOcc, n, expFreq)

        ev = 1.0
        label = '%s|%s' % (self.w, reverse_complement(self.w))
        w = self.w
        
        spaces = self.w.count('N')
        if spaces >= 1: 
            label = label.replace('N'*spaces, 'n{%d}' % spaces)
            w = self.w.replace('N'*spaces, 'n{%d}' % spaces)
        self.R.append( [w, label, obsFreq, expFreq, obsOcc, expOcc, pv, ev, -log10(ev), a, b, b-a+1, 0, n, 0, 0] ) 

    def run(self):
        l = self.l
        params = self.params

        if len(l) < params['occ'][0] or len(l) > params['occ'][1] :
            return []    

        #ALL
        if params['window'] == 'none':
            a,b = self.location
            obsOcc = len(l)
            self.__next(a,b,obsOcc)
    
        # FIXED SIZE
        elif params['window_group'] and params['center'] != None:
            window_l = params['window_group']
            center = params['center']
            while True:
                a = center - window_l / 2
                b = center + window_l / 2
                obsOcc = bisect.bisect_right(l, b) - bisect.bisect_left(l, a)
                self.__next(a,b,obsOcc)
                window_l *= 2
                a = center - window_l / 2
                b = center + window_l / 2
                if window_l > self.location[1] - self.location[0]:
                    break
        
        # elif params['window_group']:
        #     window_l = params['window_group']
        #     for a in range(self.location[0], self.location[1]-window_l+1+1, window_l):
        #         for b in range(a+window_l-1, self.location[1]+1, window_l):
        #             print (a,b)
        #             obsOcc = bisect.bisect_right(l, b) - bisect.bisect_left(l, a)
        #             self.__next(a,b,obsOcc)

        elif params['window'] != 'variable' and params['window'] != 'none':
            window_with = int(params['window'])
            for a in range(self.location[0], self.location[1]-window_with+1+1, window_with):
                b = a + window_with - 1
                obsOcc = bisect.bisect_right(l, b) - bisect.bisect_left(l, a)
                self.__next(a,b,obsOcc)


        # SLICES
        elif params['window'] == 'variable':
            if len(l) == 1:
                return []

            nwin = min(params['slices'], len(l)-1)
            step = (len(l)-1) / float(nwin)

            for I in range(nwin):
                for J in range(I+1, nwin+1):
                    i,j = int(round(I*step)), int(round(J*step))
                    a,b = l[i], l[j]
                    obsOcc = j - i + 1
                    self.__next(a,b,obsOcc)
                    
        # # ALL
        # elif params['heuristic'] == 'all':
        #     if len(l) == 1:
        #         return []
        #         
        #     #L = list(sets.Set(l))
        #     L = list_unique(l)
        #     L.sort()
        # 
        #     for i in range(len(L)-1):
        #         for j in range(i+1, len(L)):
        #             a,b = L[i], L[j]
        #             if b-a+1 > params['width'][1]:
        #                 break
        #             obsOcc = l.index(b) - l.index(a) + l.count(b)
        #             self.__next(a,b,obsOcc)

        R = self.R
        H = self.H
        for r in R:
            r[7] = r[6] * len(H) * len(R) #ev
            try:
                r[8] = -log10(r[7]) #sig
            except:
                r[8] = float("inf")

            r[12] = len(R)

        return R

    def full_window(self):
        pass
        """
        #test size
        if b-a+1 < params['width'][0] :
            a,b = get_window(location, a, b, int(params['width'][0]))
            
            #only allow keep full window 
            breakloop = 0
            continueloop = 0
            for x in range(len(l)):
                if l[x] >= a and l[x] <= l[i] and x < i:
                    breakloop = 1
                    break

                if l[x] <= b and l[x] >= l[j] and x > j:
                    continueloop = 1
                    break

            if breakloop:
                break
            if continueloop:
                continue
        """              
                    

def extract_windows_score(N, H, bg, location, w, params):
    """
    """
    R = []
    ratio = params['ratio']
    # convert position to new ref
    l = [p-location[0] for p in H[w]]

    if len(l) < params['occ'][0] or len(l) > params['occ'][1]:
        return []

    #extract windows
    mu = bg.mu(w, location)[:-len(w)+1]

    #print w
    mu = [mu[i]*N[len(w)][i] for i in range(len(mu))]
    alpha = [x * ratio for x in mu]
    
    for a,b,obsOcc in score(l, alpha, mu):
        a,b = a+location[0], b+location[0]
        if  b-a+1 < params['width'][0] or b-a+1 > params['width'][1]:
            continue

        #test occ
        if obsOcc < params['occ'][0] or obsOcc > params['occ'][1]:
            continue
                        
        n = sum(N[len(w)][a-location[0]:b-location[0]+1])

        try:
            obsFreq = obsOcc / float(n)
        except:
            cli.warning('n error')
            pass

        expFreq = bg.freq(w, (a,b))
        expOcc = expFreq * n

        pv = Stats.dist.ppois(obsOcc, expOcc)
        #pv = Stats.dist.pbinom(obsOcc, n, expFreq)
        ev = 1.0

        label = '%s|%s' % (w, reverse_complement(w))
        R.append( [w, label, obsFreq, expFreq, obsOcc, expOcc, pv, ev, -log10(ev), a, b, b-a+1, 0, n, 0, 0] ) 
        #cli.info(R[-1])           
    

    for r in R:
        r[7] = r[6] * len(H) * len(R)
        r[8] = -log10(r[7])
        r[12] = len(R)

    return R


def score(l, alpha, mu):
    """
    l     -- list of positions (p)

    alpha -- Prob in target model
    mu    -- Prob in bacground model

    return score tab

    """

    k = [0] * len(mu)
    for p in l:
        k[p] += 1

    minSize = 2
    maxSize = 1000
    minCount = 1
    maximum = 0.0
    a = -1
    b = -1
    count = 0
    regions = []

    #
    S = 0
    for i in range(len(mu)):
        s = k[i] * log( alpha[i] / mu[i] ) + ( mu[i] - alpha[i] )
        S = max(0.0, S + s)

        # extract regions
        if S > 0.0 and a == -1:
            a = i
            count += 1
            
        if S > maximum:
            maximum = S
            b = i

        if S == 0 or i == len(l)-1:
            if a != -1 and (b-a+1) >= minSize and (b-a+1) <= maxSize and count >= minCount:
                regions += [ (a,b+1, count) ]
                maximum = 0.0
                count = 0
            a = -1            

    return regions

                                     
##############################################################################
#                                                                            #
#                                 SELECT
#                                                                            #
##############################################################################
def convert_thresholds(defaults, MIN, MAX, columnHeader, columnType):
    """
    defaults -- dict key tuple (min,max)
    return dict key=param, value a tuple (min,max)
    """
    t = {}
    
    for k in defaults:
        t[k] = list(defaults[k])

    for colname, th in MIN:
        if t.has_key(colname):
            type = columnType[columnHeader.index(colname)]
            t[colname][0] = type(th)

    for colname, th in MAX:
        if t.has_key(colname):
            type = columnType[columnHeader.index(colname)]
            t[colname][1] = type(th)

    return t


def select(R, thresholds, columnHeader):
    for colname in thresholds:
        c = columnHeader.index(colname)
        f = lambda x: 1
        if thresholds[colname][0] != None and thresholds[colname][1] != None:
            f = lambda x: x[c] >= thresholds[colname][0] and x[c] <= thresholds[colname][1]
        elif thresholds[colname][0] != None:
            f = lambda x: x[c] >= thresholds[colname][0]             
        elif thresholds[colname][1] != None:
            f = lambda x: x[c] <= thresholds[colname][1]
           
        R = [i for i in R if f(i)]
    return R
    

def sort(R, criteria, columnHeader):
    l = []
    for c in criteria:
        if c[0] == '+':
            growing = True
        else:
            growing = False
        
        i = columnHeader.index(c[1:])
        l.append( (i, growing) )

    def str_cmp(a, b):
        if a < b:
            return -1
        if a > b:
            return 1
        return 0        

    def mycmp(a, b):
        for i,G in l:
            if type(a[i]) == str:
                acmp = str_cmp
            else:
                acmp = cmp

            if G and acmp(a[i], b[i]) != 0:
                return acmp(a[i], b[i])
            if not G and acmp(b[i], a[i]) != 0:
                return acmp(b[i], a[i])
        return 0
    R.sort( mycmp )
    return R


def window_rank(R, columnHeader):
    c = columnHeader.index('w_rank')
    c2 = columnHeader.index('seq')
    C = {}
    for i in range(len(R)):
        C[R[i][0]] = C.get(R[i][c2], 0) + 1
        R[i][c] = C[R[i][c2]]
    return R


def rank(R, columnHeader):
    c = columnHeader.index('rank')
    for i in range(len(R)):
        R[i][c] = i+1
    return R


def format_header(header, **kwargs):
    return header % kwargs


def format_output(l, columnHeader, columnHeaderChar, formatRow):
    s = []
    s += [columnHeaderChar + '\t'.join(columnHeader)]
    for i in l:
        s += [formatRow % tuple(i) ]        
    return '\n'.join(s)


##############################################################################
#                                                                            #
#                                 BACKGROUND
#                                                                            #
##############################################################################

def smooth(tab, windowWidth=500):
    """
    Smooth given table by averaging value in a sliding window of width windowSize
    @param windowSize: int

    """
    smoothTab = [0.0] * len(tab)
    for j in range(len(tab)):
        a, b = window(j, windowWidth, 0, len(tab))
        smoothTab[j] = sum(tab[a:b+1]) / float(len(tab[a:b+1]))
    return smoothTab


def window(center, width, a, b):
    """
    @param center: int window center
    @param width: int window width
    
    @Return: [start,end] window of width width  in [a,b]
    """
    start = max(a, center - int(width / 2.0))
    start = min(start, b - width)
    return [start, start + width - 1]


class MM(markov.MMError):
    
    def __getitem__(self, key):
        if self.strand == '+-':
            wrc = reverse_complement(key)
            if key != wrc:
                if self.dyad:
                    return self.P(key[:self.monad]) * self.P(key[-self.monad:]) + self.P(wrc[:self.monad]) * self.P(wrc[-self.monad:])
                else:    
                    return self.P(key) + self.P(wrc)
            else:
                if self.dyad:
                    return self.P(key[:self.monad]) * self.P(key[-self.monad:])
                else:    
                    return self.P(key)
        else:
            if self.dyad:
                return self.P(key[:self.monad]) * self.P(key[-self.monad:])
            else:    
                return self.P(key)

    def get(self, key, default):
        return self.__getitem__(key)


class Bg(dict):
    """
    key: window location tuple
    value: H (word count)

    """
    def __init__(self, location=None, W=None, step=None, l=2, strand='', overlap=False, params=None):
        """
        W -- window size
        
        """
        self.W = W
        self.step = step or W
        self.l = l
        self.intervals = []
        self.location = location
        self.strand = strand
        self.overlap = overlap
        self.params = params

    def __repr__(self):
        return 'BG location=[%+05d:%+05d] l=%d strand=%s overlap=%d W=%d' \
        % (self.location[0], self.location[1], self.l, self.strand, self.overlap, self.W)

    def build(self, sequences):
        info = cli.Info((self.location[1] - self.l + 1 - self.location[0]) / self.step + 1) 
        for i in range(self.location[0], self.location[1] - self.l + 1, self.step):
            loc = (i, i+self.W-1)
            info('building BG pos=%i' % i)
            r = count_words(sequences, self.l, loc, strand=self.strand, overlap=self.overlap, params=self.params)
            H = r['H']
            N = r['N']

            #compute frequency
            for k in H:
                H[k] = len(H[k]) / float(sum(N[len(k)]))
            self[loc] = H
        sorted_keys = self.keys()
        sorted_keys.sort()
        self.intervals = sorted_keys

    def build_from_oligo_file(self, filename, location=(-100000, 100000)):
        self.W = location[1] - location[0] +1
        self.location = location
        F = load_oligo_file(filename) 
        self[location] = F
        self.intervals.append(location)

    def build_markov_from_oligo_file(self, filename, location=(-100000, 100000)):
        self.W = location[1] - location[0] +1
        self.location = location
        mm = MM()
        mm.load_oligo_file(filename) 
        mm.overlap = self.overlap
        mm.strand = self.strand
        mm.set_NExtension(self.params['spacing'])
        self[location] = mm
        self.intervals.append(location)

    def build_markov(self, sequences, order=1):
        info = cli.Info((self.location[1] - self.l + 1 - self.location[0]) / self.step + 1) 
        for i in range(self.location[0], self.location[1] - self.l + 1, self.step):
            info('building Markov BG pos=%i' % i)
            loc = (i, i+self.W-1)
            mm = MM(order)
            mm.dyad = self.params['dyad']
            mm.monad = self.l
            mm.overlap = self.overlap
            mm.strand = self.strand
            mm.set_NExtension(self.params['spacing'])
            mm.learn([s.get_dna(loc) for s in sequences])
            self[loc] = mm

        sorted_keys = self.keys()
        sorted_keys.sort()
        self.intervals = sorted_keys

    def freq(self, word, location):
        #find intervals to use
        a,b = location
        inLocation = False
        f = 0
        n = 0

        for interval in self.intervals:
            if a >= interval[0] and a <= interval[1]:
                inLocation = True
            if inLocation:
                F = self[interval]
                x = min(b, interval[1]) + 1 - max(a, interval[0])
                #change this ?
                f += F.get(word, 0.0) * x
                if F.get(word, None) == None:
                    cli.warning('can not find key="%s" in BG' % word)
                n += x
            if b <= interval[1]:
                break

        return f / float(n)

    def mu(self, word, location):
        """
        score scanning
        """
        mu_tab = [0] * (location[1] -location[0] + 1)

        #find intervals to use
        a,b = location
        inLocation = False

        f = 0
        n = 0
        for interval in self.intervals:
            if a >= interval[0] and a <= interval[1]:
                inLocation = True
            
            if inLocation:
                F = self[interval]

                i,j = max(a, interval[0]), min(b, interval[1])
                x = j-i+1
                #change this ?
                f = F.get(word, 0.0)
                if not F.get(word, 0):
                    cli.warning('can not find key="%s" in BG' % word)
                mu_tab[i-a:i-a+x] = [f] * x
            if b <= interval[1]:
                break

        return smooth(mu_tab, 200)

