"""
:NAME:
ST

:DESC:
suffix tree for dna sequence indexation using suffix tree
nodes: dict key: IUPAC + '$' for data {sequence, [positions]}
transition are labelled by IUPAC code

code	description
A	Adenine
C	Cytosine
G	Guanine
T	Thymine
U	Uracil
R	Purine (A or G)
Y	Pyrimidine (C, T, or U)
M	C or A
K	T, U, or G
W	T, U, or A
S	C or G
B	C, T, U, or G (not A)
D	A, T, U, or G (not C)
H	A, T, U, or C (not G)
V	A, C, or G (not T, not U)
N	Any base (A, C, G, T, or U)

suffix tree is reffered as st

motif occurrence location [a,b] length w=(b-a+1)

foword strand  position=i

reverse strand position=L-1-b = L-w-a

"""
import sys

import cli
import dna

MIN_INT = -sys.maxint +1

A = 0
C = 1
G = 2
T = 3
N = 4
R = 5
Y = 6
M = 7
K = 8
W = 9
S = 10
B = 11
D = 12
H = 13
V = 14

IUPAC_LENGTH = 5
#IUPAC_LENGTH = 15

def base2IUPAC(base, full=0):
    #IUPAC = [base]
    if base == 'A':
        IUPAC = [A]
    elif base == 'C':
        IUPAC = [C]
    elif base == 'G':
        IUPAC = [G]
    elif base == 'T':
        IUPAC = [T]        
    else:
        IUPAC = []

    IUPAC.append(N) #N

    if IUPAC_LENGTH == 5:
        return IUPAC

    if base == 'A' or base == 'G':
        IUPAC.append(R) #R
    if base == 'C' or base == 'T':
        IUPAC.append(Y) #Y
    if base == 'C' or base == 'A':
        IUPAC.append(M) #M
    if base == 'G' or base == 'T':
        IUPAC.append(K) #K
    if base == 'A' or base == 'T':
        IUPAC.append(W) #W
    if base == 'C' or base == 'G':
        IUPAC.append(S) #S
    if base == 'C' or base == 'G' or base == 'T':
        IUPAC.append(10) #B
    if base == 'A' or base == 'G' or base == 'T':
        IUPAC.append(11) #D
    if base == 'A' or base == 'C' or base == 'T':
        IUPAC.append(12) #H
    if base == 'A' or base == 'C' or base == 'G':
        IUPAC.append(13) #V

    return IUPAC

def IUPACcode2str(x):
    if x == A:
        return 'A'
    elif x == C:
        return 'C'
    elif x == G:
        return 'G'
    elif x == T:
        return 'T'
    elif x == N:
        return 'N'
    elif x == R:
        return 'R'
    elif x == Y:
        return 'Y'
    elif x == M:
        return 'M'
    elif x == K:
        return 'K'
    elif x == W:
        return 'W'
    elif x == S:
        return 'S'
    elif x == B:
        return 'B'
    elif x == D:
        return 'D'
    elif x == H:
        return 'H'
    elif x == V:
        return 'V'
    else:
        return '?'


def str2IUPACcode(x):
    if x == 'A':
        return A
    elif x == 'C':
        return C
    elif x == 'G':
        return G
    elif x == 'T':
        return T
    elif x == 'N':
        return N
    elif x == 'R':
        return R
    elif x == 'Y':
        return Y
    elif x == 'M':
        return M
    elif x == 'K':
        return K
    elif x == 'W':
        return W
    elif x == 'S':
        return S
    elif x == 'B':
        return B
    elif x == 'D':
        return D
    elif x == 'H':
        return H
    elif x == 'V':
        return V
    else:
        return N


##############################################################################
#                                                                            #
#                                 PRINT
#                                                                            #
##############################################################################
def display(node, maxDepth=3, full=False):
    #print '\n'
    for code in range(IUPAC_LENGTH):
        subnode = node.childs[code]
        if subnode:
            PREFIX = '|  ' * node.depth + '+--'
            if full:
                print PREFIX + IUPACcode2str(code) + ' [C%d,I%d,N%d,D%d]' % (subnode.count, subnode.IUPACcount, subnode.NCount, subnode.depth)
            else:
                print PREFIX + IUPACcode2str(code) + '[%d]' % (subnode.count)
            display(subnode, maxDepth, full)
    #print '|  ' * node.depth


def print_count(c):
    l = [ (k,v) for k,v in c.items() ]
    l.sort()
    for w,wcount in l:
        wrc = dna.reverse_complement(w)
        print '%s|%s %4d' % (w, wrc, wcount)


def count(node, path, c, minLength, maxLength):
    """
    c -- dict key=work value=count
    """
    if node is None:
        return

    if len(path) >= minLength and len(path) <= maxLength:
        if path[-1] != '-':
            c[path] = node.count

    for code in range(IUPAC_LENGTH):
        X = IUPACcode2str(code)
        count(node.childs[code], path+X, c, minLength, maxLength)


##############################################################################
#                                                                            #
#                                 EXTRACT
#                                                                            #
##############################################################################
def extract(node, path, c, minLength, maxLength):
    """
    c -- dict key=work value={seq1 : [positions seq1], seq2 : [positions seq2]} 
    """
    if node is None:
        return

    if len(path) >= minLength and len(path) <= maxLength:
        if path[-1] != '-':
            c[path] = node.positions

    for code in range(IUPAC_LENGTH):
        X = IUPACcode2str(code)
        extract(node.childs[code], path+X, c, minLength, maxLength)


def get_positions_two_strands(c, overlap=False):
    positions = {}
    for wf in c:
        wrc = dna.reverse_complement(wf)
        w = min(wf, wrc)
        if positions.has_key(w):
            continue

        l = []
        f = c[wf]
        for i in f:
            l += [j[0] for j in f[i]]

        if wf != wrc:
            r = c.get(wrc, {})
            for i in r:
                l += [j[0] for j in r[i]]

        l.sort()
        positions[w] = l
        
    return positions


def get_positions(c):

    positions = {}
    for w in c:
        l = []
        f = c[w]
        for i in f:
            l += [j[0] for j in f[i]]

        l.sort()
        positions[w] = l

    return positions


def group_rc(c):
    """
    TODO : grouprc and overlapping on both strands
    """
    groupedc = {}
    for w in c:
        wrc = dna.reverse_complement(w)
        x = min(w, wrc)
        if groupedc.has_key(x):
            continue
        groupedc[x] = c[w]
        if w != wrc and c.has_key(wrc):
            groupedc[x] += c[wrc]
        else:
            groupedc[x] *= 2            
    return groupedc


##############################################################################
#                                                                            #
#                                 NODE
#                                                                            #
##############################################################################
class Node(object):
    """
    code -- label for the link  parent -> node
    """
    __slots__ = ['code', 'depth', 'childs', 'count', 'lastPosition', 'IUPACcount', 'positions', 'NCount']

    def __init__(self, code=None, maxChild=IUPAC_LENGTH, depth=0, IUPACcount=0):
        self.code = code
        self.depth = depth
        self.childs = [None] * maxChild

        self.count = 0
        self.lastPosition = (0, MIN_INT,MIN_INT) #(seq, start, end)
        self.IUPACcount = IUPACcount
        self.positions = {}
        self.NCount = 0


def get_node(node, code):
    if node.childs[code] == None: 
        xnode = Node(code, IUPAC_LENGTH, node.depth+1, node.IUPACcount)
        node.childs[code] = xnode
        if code >= N:
            xnode.IUPACcount += 1
    else:
        xnode = node.childs[code]
    return xnode


##############################################################################
#                                                                            #
#                           SUFFIX TREE
#                                                                            #
##############################################################################
class SuffixTree(object):
    """
    """
    def __init__(self, maxDepth=4, overlapping=False, maxIUPAC=1, NExtension=(1,1), storePosition=False):
        assert(NExtension[0] >= 1)
        assert(NExtension[1] >= NExtension[0])
        #parameters
        self.maxDepth = maxDepth
        self.overlapping = overlapping
        self.maxIUPAC = maxIUPAC
        self.NExtension = NExtension 

        #data
        self.root = Node(None, IUPAC_LENGTH, 0)
        self.nodeCount = 1
        #self.n = 0 #number of scanned motifs
        self.N = 0 #number of scanned positions

        self.s = 0 #number of sequences
        self.storePosition = storePosition


    def add_word(self, dna, i, s, shift=0):
        """
        """
        overlapping = self.overlapping
        nodes = [self.root]

        jmax = self.maxDepth + self.maxIUPAC * self.NExtension[1]

        for j in range(min(len(dna)-i, jmax)):
            #self.n += 1
            newnodes = []
            IUPAC = base2IUPAC(dna[i+j], 1)

            for node in nodes:
                for code in IUPAC:

                    # new N
                    if code == N and node.code != N: 

                        #can not end with N
                        if node.depth +1 >= self.maxDepth:
                            continue
                        if node.depth == 0:
                            continue

                    # extend N
                    if code == N and node.code == N: 
                        #N size constraint
                        if node.NCount >= self.NExtension[1]:
                            continue
                        node.NCount += 1
                        newnodes.append(node)
                        continue              

                    #branch after N
                    if code != N and node.code == N:
                        if node.NCount < self.NExtension[0]:
                           continue


                    #IUPAC count
                    if code >= N:
                        if node.IUPACcount >= self.maxIUPAC: # B1 > 
                            continue

                    ##
                    # get or create node
                    #
                    xnode = get_node(node, code)
                    #
                    ##    

                    ##
                    # NCount IUPACcount
                    #
                    if code == N:
                        xnode.NCount = 1
                    else:
                        xnode.NCount = node.NCount
                    #
                    ##


                    ##
                    # count
                    #
                    if overlapping or s !=  xnode.lastPosition[0] or i+shift > xnode.lastPosition[2]:
                        xnode.lastPosition = (s, i+shift,i+j+shift)
                        xnode.count += 1
                        if self.storePosition and xnode.depth == self.maxDepth:
                            xnode.positions.setdefault(s, []).append( (i+shift, i+j+shift) )

                    #
                    ##

                    if xnode.depth < self.maxDepth:
                        newnodes.append(xnode)

            nodes = newnodes
            if len(nodes) == 0:
                return

    def add_dna(self, dna, shift=0):
        #for i in range(1):
        for i in range(len(dna)):
            self.add_word(dna, i, self.s, shift)
            self.N += 1
        self.s += 1

    def add_sequences(self, sequences):
        info = cli.Info(len(sequences))
        for i in range(len(sequences)):
            info('Proccessing sequence %s' % sequences[i].id)
            self.add_dna(sequences[i].sequence)

    def __repr__(self):
        return 'NExtension=%d:%d maxIUPAC=%d' % (self.NExtension[0], self.NExtension[1], self.maxIUPAC)
            

    def count(self, minLength=1, maxLength=None):
        maxLength = maxLength or self.maxDepth
        c = {}
        count(self.root, '', c, minLength, maxLength)
        return c

    def extract(self, minLength=1, maxLength=None):
        maxLength = maxLength or self.maxDepth
        c = {}
        extract(self.root, '', c, minLength, maxLength)
        return c

    def get(self, w):
        node = self.root
        for letter in w:
            code = str2IUPACcode(letter)
            node = node.childs[code]
        return node

##############################################################################
#                                                                            #
#                                 TESTS
#                                                                            #
##############################################################################
def realtest(N=2):
    #filename = 'E2F.fa'
    filename = 'Test/MM0.fa'    
    sequences = dna.fasta2sequences(filename)
    print sequences
    st = SuffixTree(maxDepth=N, overlapping=True, maxIUPAC=N, NExtension=(1,1), storePosition=0)
    info = cli.Info(len(sequences))
    for i in range(len(sequences)):
        info('Processing sequence %s' % sequences[i].id)
        st.add_dna(sequences[i].sequence)

    #display(st.root)
    #count = group_rc(st.count(minLength=N, maxLength=N))
    #print_count(count)
    return st


def test(N=2):
    n = 5
    st = SuffixTree(maxDepth=n, overlapping=True, maxIUPAC=1, NExtension=(6,6), storePosition=1)
    dna = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    print dna
    st.add_dna(dna)
    display(st.root, full=1)
    print st.extract(minLength=n, maxLength=n)
    keys = st.extract(minLength=n, maxLength=n).keys()
    print keys
    
    #print st.nodeCount
    #count = group_rc(st.count(minLength=N, maxLength=N))
    #print_count(count)


if __name__ == '__main__':
    from Core.Devel.utils import memory_usage
    m = memory_usage()
    #str = test(4)
    x = realtest(int(sys.argv[1]))
    print 'memory usage : %.3fM' % (memory_usage() - m)

