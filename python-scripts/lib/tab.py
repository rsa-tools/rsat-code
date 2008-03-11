"""
:NAME:
tab

:DESC:
Tab delimited text file reader/writer

:TODO:
COLUMN_TYPE
COLUMN_FORMAT
COLUMN_HEADER
COLUMN_HEADER_HELP

"""

import os
import string
import time
import re
import gzip

COMMENT_CHAR = ';'
TAB_HEADER_CHAR = '#'

class Timer:
    def __init__(self):
        self.start    = time.time()
        self.duration = 0.0
        self.stopped  = False
        
    def stop(self):
        if not self.stopped:
            self.duration = time.time() - self.start
            self.stopped = True

    def get_value(self):
        return time.time() - self.start    

    def __repr__(self):
        self.duration = time.time() - self.start
        gm = time.gmtime(self.duration)
        gms = time.strftime('%Hh%Mm%Ss', gm)
        return '%.3fs (%s)' % (self.duration, gms)


def magic_type(string):
    try:
        int(string)
        return int
    except:
        pass
    try:
        float(string)
        return float
    except:
        pass
    return str
        

class Tab(object):
    """
    .header : list
    .data : list of lists
    """

    def __init__(self, header):
        self.header = header
        self.data = []
        self.index_header()
        self.timer = Timer()
        self.comments = None

    def index_header(self):
        #build index
        self.index = {}
        for i in range(len(self.header)):
            self.index[self.header[i]] = i        

    def convert_data(self):
        if len(self.data) == 0:
            return
        firstline = self.data[0]
        types = []
        for item in firstline:
            types += [ magic_type(item) ]
        for i in range(len(self.data)):
            self.data[i] = [types[j](self.data[i][j]) for j in range(len(self.data[i]))]

        self.types = types

    def get_value(self, line, column):
        return self.data[line][self.index[column]]

    def append(self, data):
        line = [None] * len(self.header)
        for k,v in data.items():
            line[self.index[k]] = v
        self.data.append(line)

    def filter(self, thresholds):
        """Filter (in place)
        thresholds -- {'col' : (min, max), 'col2' : (None, max), ... }
        convert    -- function to apply to column before filtering

        """
        for colname in thresholds:
            c = self.index[colname]

            f = lambda x: 1
            if thresholds[colname][0] != None and thresholds[colname][1] != None:
                f = lambda x: x[c] >= thresholds[colname][0] and x[c] <= thresholds[colname][1]
            elif thresholds[colname][0] != None:
                f = lambda x: x[c] >= thresholds[colname][0]
            elif thresholds[colname][1] != None:
                f = lambda x: x[c] <= thresholds[colname][1]

            return [ i for i in self.data if f(i) ]

    def select(self, columns):
        """Select a list of column (inplace)
        columns -- list of column names

        """
        for i in range(len(self.data)):
            self.data[i] = [self.data[i][self.index[c]] for c in columns]

        self.header = columns
        self.index_header()

    def __len__(self):
        return len(self.data)

    def insert_column(self, name, position=0, values=None):
        for i in range(len(self.data)):
            self.data[i].insert(position, values[i])

        self.header.insert(position, name)
        self.index_header()

    def sort(self, criteria):
        """Sort Tab
        criteria -- list of key Example ['+col1', '-col3']
        
        """
        if len(criteria) == 0 or criteria == ['']:
            return

        l = []
        for c in criteria:
            if c[0] == '+':
                growing = True
            else:
                growing = False

            i = self.index[c[1:]]
            l.append( (i, growing) )

        def mycmp(a, b):
            for i,G in l:
                if type(a[i]) == str:
                    acmp = Core.types.string.cmp 
                else:
                    acmp = cmp

                if G and acmp(a[i], b[i]) != 0:
                    return acmp(a[i], b[i])
                if not G and acmp(b[i], a[i]) != 0:
                    return acmp(b[i], a[i])
            return 0

        self.data.sort( mycmp )

    def rank(self):
        c = self.index['rank']
        for i in range(len(self.data)):
            self.data[i][c] = i+1

    def get_column(self, id):
        c = self.index[id]
        type = magic_type(str(self.data[0][c]))
        return [type(l[c]) for l in self.data]            

    ##
    # OUTPUT
    ##
    def to_txt(self, comments=None, help=None, addTimeInfo=False):
        comments = comments or self.comments
        content = []
        #comments
        if addTimeInfo:
            content += [ COMMENT_CHAR + ' date\t%(date)s' % {'date' : time.ctime()} ]
            content += [ COMMENT_CHAR + ' runningTime\t%s' % str(self.timer) ]
        if comments:
            for k,v in comments.items():
                content += [ COMMENT_CHAR + ' %s\t%s'% (str(k), str(v)) ]

        #header help
        if help:
            content += [ help ]
            #for i in range(len(self.header)):
            #    print >> f, '; column headers'
            #    print >> f, COMMENT_CHAR, '\t%d\t\t\t%s' % (i+1, headerHelp[self.header[i]])


        #header
        content += [ TAB_HEADER_CHAR + '\t'.join(self.header) ]

        #data
        for line in self.data:
            content += [ '\t'.join([str(i) for i in line]) ]

        return '\n'.join(content)

    def to_html(self):
        header = self.header
        content = []
        content += ['<table id="tx" class="sortable">']
        content += ['<tr class="header"><td>' + '</td><td>'.join([str(label) for label in header]) + '</td></tr>' ]
        for data in self.data:
            content += [ '<tr><td>' + '</td><td>'.join([str(data[i]) for i in range(len(header))]) + '</td></tr>' ]
        content += ['</table>' ]

        return '\n'.join(content)

##
#
# INPUT
#
##
def read(filename):
    f = open(filename)

    comments = []
    header = []
    data = []
    last = None
    for line in f:
        if line.startswith(COMMENT_CHAR):
            comments += [line[1:].strip()]
            last = 'comment'
        elif line.startswith(TAB_HEADER_CHAR):
            header = map(string.strip, line[1:].split('\t'))
            last = 'header'
        else:
            if last == 'comment':
                header = map(string.strip, comments.pop().split('\t'))
            if header:
                data += [ map(string.strip, line.split('\t')) ]
                last = 'data'
        
    t = Tab(header)
    t.data = data
    t.comments = dict( [ (line.split('\t')[0], '\t'.join(line.split('\t')[1:])) for line in comments])

    t.convert_data()

    return t

    #return dict(comments=comments, header=header, data=data, filename=filename)


def create_filename(id, params, ext=None):
    p = []
    for k,v in params.items():
        v = str(v).replace(',', '+')
        p += [ '%s_%s' % (k, v)]
    if ext != None:
        name = '%s,%s%s' % (id, ','.join(p), ext)
    else:
        name = '%s,%s' % (id, ','.join(p))
    name = name.replace(' ', '')
    name = name.replace('(', '')
    name = name.replace(')', '')
    return name


def filename2id(filename):
    id = os.path.basename(filename).split('.')[0]
    return id


def filename2easyreading(filename):
    str = '.'.join(os.path.basename(filename).split('.')[:-1])
    return str.replace(',', ' ').replace('_', '=')


def decode_filename(filename):
    str = '.'.join(os.path.basename(filename).split('.')[:-1])
    l = str.split(',')
    id = l[0]

    d = {}
    for opt in l[1:]:
        k, v = opt.split('_')
        d[k] = v
    return id, d


def convert_thresholds(MIN, MAX, columnHeader, columnType, defaults=None):
    """
    COLUMN_HEADER = ['seq', 'identifier', 'obs_freq', 'exp_freq']
    COLUMN_TYPE   = [str  ,  str        , float     , float     ]
    thresholds = convert_thresholds(options.min, options.max, COLUMN_HEADER, COLUMN_TYPE)
    """
    defaults = defaults or {}
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


def string2filter(s):
    """
    string format : Se:(0.4,None)&pv:(0.001,0.1)'
    return dict that can be used for filtering
    """
    #filter
    d = {}
    for filter in s.split('&'):
        id, values = filter.split(':')
        min, max = values.strip().strip('()').split(',')
        try:
            min = float(min)
        except:
            min = None
        try:
            max = float(max)
        except:
            max = None

        d[id]=(min,max)
    return d


# oligo-analysis 
DYAD = re.compile('([ACGT]+)N\{(\d+)\}([ACGT]+)')
def load_oligo_file(filename):
    """Load data from oligo-analysis formated file (can be gziped)

    """
    F = {}

    if filename.endswith('.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)

    i = 0
    for line in f:
        if line.startswith('#') or line.startswith(';'):
            continue
        elements = line.strip().split()
        w, freq = elements[0], elements[2]  
        w = w.upper()
        dyad = DYAD.findall(w)
        if dyad:
            w = dyad[0][0] + 'N' * int(dyad[0][1]) + dyad[0][2]
        freq = float(freq)
        F[w] = freq

    return F


def load_assemblies(filename):
    if filename.endswith('.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)

    i = 0
    assemblies = None
    for line in f:
        line = line.strip()
        if line.startswith(';assembly') or line.startswith('; Isolated patterns'):
            if assemblies == None:
                assemblies = []
            else:
                assemblies += [assembly]
            assembly = []


        if line.startswith('#') or line.startswith(';'):
            continue
        if line == '':
            continue
        elements = line.strip().split()
        w = elements[0].upper().replace('.', 'N')
        assembly += [w]
    return assemblies


if __name__ == '__main__':
    t = read('test.tab')
    data = t.data
    print '\n'.join(t.comments)
    #print t['header']
    #print data[0:3]
