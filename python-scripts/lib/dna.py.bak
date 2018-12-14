"""
:NAME:
bio

:DESC:
Utilities for dna sequences manipulation

:HELP:
A sequence is stored as a simple string 'ACTGGGT'
A set of sequences is stored as a list of sequences ['ATGGAT', 'TGGTTGGT']

reverse_complement(sequence)
read_fasta(input, strand='+-')
write_fasta(output, sequences, labels)

"""

import sys
import time
import re
import gzip

LETTER2J = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3, 'N' : 4}
BASES = ['A', 'C', 'G', 'T']
COMPLEMENT = { 
       'A' : 'T',
       'T' : 'A',
       'C' : 'G',
       'G' : 'C',
       'n' : 'n',
       'N' : 'N',
       'R' : 'Y',
       'Y' : 'R',
       'M' : 'K',
       'K' : 'M',
       'W' : 'W',
       'S' : 'S',
       'B' : 'V',
       'V' : 'B',
       'D' : 'H',
       'H' : 'D',
}

'''
def write_data(f, data):
    opened = False
    if type(f) is str:
        f = open(f, 'w')
        opened = True
    f.write(data)
    f.flush()
    if opened:
        f.close()
'''

def reverse_complement(sequence): 
    ns = [COMPLEMENT.get(b, b) for b in sequence] 
    ns.reverse()
    return ''.join(ns)

##############################################################################
#                                                                            #
#                       FASTA READER/WRITER
#                                                                            #
##############################################################################
def read_fasta_full(input, strand='+'):
    opened = False
    if type(input) is str:
        input = open(input, mode='rU')
        opened = True
    f = FastaReader(input)
    sequences = []
    names = []
    for name, seq in f:
        seq = seq.upper()
        if strand.find('+') > -1:
            sequences += [ seq ]
            names += [ name ]
        if strand.find('-') > -1:
            sequences += [ reverse_complement(seq) ]
            names += [ name + ' rc' ]
    if opened:
        input.close()
    return (names, sequences)

def read_fasta(input, strand='+'):
    return read_fasta_full(input, strand)[1]

def write_fasta(f, sequences, labels=None):
    opened = False
    if type(f) is str:
        f = open(f, 'w')
        opened = True
    writer = FastaWriter(f)
    for i in range(len(sequences)):
        if labels:
            title = labels[i]
        else:
            title = str(i+1)
        writer.append(title, sequences[i])
    if opened:
        f.close()

class Writer(object):
    """
    Generic writer class
    """
    LINE_WIDTH = 80 #line width for output

    def append(self, title, dna=None):
        """Append a sequence to file
        @arg title:    title string
        @arg dna:      dna string
        
        implement this method
        """
        pass

class FastaWriter(Writer):
    """
    Fasta class to write in FASTA format
    """
    LINE_WIDTH = 80 #line width for output

    def __init__(self, f):
        """
        @arg f: file
        """
        Writer.__init__(self)
        self.out = f

    def append(self, title, dna=None):
        """Append a sequence to FASTA file
        @arg title:    title string
        @arg dna:      dna string

        """
        self.out.write('>%s\n' % title)
        for i in range(0, len(dna), self.__class__.LINE_WIDTH):
            self.out.write('%s\n' % dna[i:i+self.__class__.LINE_WIDTH])
        self.out.flush()

class Reader(object):
    """Generic class to read FASTA with builtin iterator
    """

    def read(self, n=1):
        """
        Read n (or less if eof) sequences from fasta file default=1
        @arg n: number of sequences to read [default=1]

        @return: list of tuple (title, sequence)
        """
        l = []
        count = 0
        for id, sequence in self:
            if count >= n:
                break
            count += 1
            l.append( (id, sequence) )
        return l

    def __iter__(self):
        return self

    def next(self, decode=False):
        """
        Read next sequence in file

        @param decode:  if True decode title
        @return (title, sequence) 
        """
        raise NotImplementedError

class FastaReader(Reader):
    """Fasta Reader
    Use this one
    """
    def __init__(self, f):
        Reader.__init__(self)
        self.f = f
        self.last = ''
        self.sequence = []

    def next(self):
        self.item = False
        self.eof  = False

        while 1:
            if self.last:
                self.l = self.last
            else:
                line = self.f.readline()
                if line == '':
                    self.eof = True
                    break
                else:
                    self.l = line.strip()
            self.last = ''
            if self.l == '':
                continue
            elif self.l[0] == '#':
                continue
            elif self.l[0] == '>':
                if not self.item:
                    self.item  = True
                    self.sequence = []
                    self.title = self.l[1:]
                else:
                    self.last = self.l
                    break
            else:
                self.sequence.append(self.l)
        if self.item:
            return self.title, re.sub(r'[^ACGTUNacgtun]', '', ''.join(self.sequence))
        if self.eof:
            raise StopIteration

########################################
#                                      #
#           DEPRECATED CODE
#                                      #
########################################
class Sequences(list):
    def __init__(self, l=[], name=''):
        super(Sequences, self).__init__(l)
        self.name = name
        self.info  = {}
        self.idx   = None
        self.location = find_location(self)
        
    def get_location(self):
        self.location = find_location(self)
        return self.location
        
    def __repr__(self):
        self.get_location()
        return '%s [%d sequences] [%+05d:%+05d]' % (self.name, len(self), self.location[0], self.location[1])

    def create_idx(self):
        self.idx =  types.list.index(self, lambda s: s.id)
        return self.idx

    def change_locations(self, new_origins):
        """update each sequence location
        """
        assert(len(self)==len(new_origins))
        for i in range(len(self)):
            self[i].location  = (self[i].location[0] - new_origins[i], self[i].location[1] - new_origins[i])

class Sequence:
    """Simple sequence object
    """
    def __init__(self, id=None, location=None, sequence=''):
        """
        id        -- string ex NM_000043
        location  -- (start, end)  start, end: int position (relative to TSS) default=(0, len(sequence)-1)
                     start and end are included
        sequence  -- DNA string 'ATGC...' default=''
        NOTE: sequence can be automaticaly reshaped
        """

        self.id       = id
        self.location = location or (0,len(sequence)-1) #relative location
        self.sequence = sequence.upper()

    def reverse_complement(self):
        self.sequence = reverse_complement(self.sequence)

    def __len__(self):
        return self.location[1] - self.location[0] + 1

    def get_dna(self, location=None):
        location = location or self.location
        if location[0] > self.location[1] or location[1] < self.location[0]:
            return ''
        if self.location[0] >= location[0]:
            a = 0
        else:
            a = location[0]-self.location[0]
        if self.location[1] <= location[1]:
            return self.sequence[a:]
        else:
            b = self.location[1] - location[1]
            return self.sequence[a:-b]

def is_valid_location(start, end, min, max):
    if start >= end:
        return False
    if start < min or start > max:
        return False
    if end < min or end > max:
        return False
    return True

def find_location(sequences):
    if len(sequences) > 0:
        return (min([i.location[0] for i in sequences]), max([i.location[1] for i in sequences]))
    else:
        return 0,0

def fasta2sequences(uf, location=None, rightPosition=None, leftPosition=None, centerPosition=None):
    """
    rightPosition    -- right bound
    leftPosition     -- left bound    
    centertPosition  -- center  

    """
    try:
        if type(uf) is str:
            if uf.endswith('.gz'):
                f = gzip.open(uf)
            else:
                f = open(uf)
        else:
            f = uf

        l = []
        for title, dna in FastaReader(f):
            dna = dna.upper()
            loc = location
            if not loc and rightPosition != None:
                loc = (rightPosition - len(dna) + 1, rightPosition)
            if not loc and leftPosition != None:
                loc = (leftPosition, leftPosition + len(dna) - 1)
            if not loc and centerPosition != None:
                loc = (1 + centerPosition - int(round(len(dna) / 2.0, 0)), centerPosition + len(dna) / 2) 

            l.append( Sequence(id, loc, dna) )
        r = Sequences(l)
        return r
    except:
        sys.stderr.write('Error while reading input sequence')
        raise
