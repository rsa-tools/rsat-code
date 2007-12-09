#!/usr/bin/env python
import os, sys, glob

[sys.path.insert(1,os.path.join(sys.path[0], path)) for path in glob.glob1(os.path.abspath(sys.path[0]),'*.zlib')]

import optparse
import time

import Core.Sys.timer
import Core.Persistence.op as op
import Core.UI.cli as cli
import Core.types
import Core.utils
import Bio.sequence

import ORMLib.ormlib as ormlib

DEBUG = 1
cli.UPDATE = 1
########################################
#
# OUTPUT HEADER
#
########################################

HEADER = ''';
; ORM
; Detection of locally overrepresented motifs
; %(command)s
; %(date)s
; running time                     %(runningTime)s
; oligomer length                  %(l)d
; strand                           %(strand)s
; window width                     %(windowWidth)s
; bg window width                  %(bgwindowWidth)s            
; Input file                       %(inputfilename)s
; nb of sequences                  %(numberOfSequences)d
; scanned region                   [%(start)+05d:%(end)+05d]
; scanned words                    %(scannedWords)d
; lower thresholds                 %(lowerThreshold)s
; upper thresholds                 %(upperThreshold)s
; sorting criteria                 %(sort)s
;
; column headers
;    1    seq                oligomer sequence
;    2    identifier         oligomer identifier
;    3    obs_freq           observed relative frequency
;    4    exp_freq           expected relative frequency
;    5    occ                observed occurrences
;    6    exp_occ            expected occurrences
;    7    occ_P              occurrence probability
;    8    occ_E              E-value for occurrences
;    9    occ_sig            occurrence significance
;    10   start              occurrence start location
;    11   end                occurrence end location
;    12   width              window width          
;    13   n_win              number of analysed windows (for a given word)
;    14   n_pos              number of analysed positions (for a given window) 
;    15   w_rank             window rank
;    16   rank               occurrence rank (based on given criteria)'''

COLUMN_HEADER_CHAR = ';'
COLUMN_HEADER = ['seq', 'identifier', 'obs_freq', 'exp_freq', 'occ', 'exp_occ', 'occ_P', 'occ_E', 'occ_sig', 'start', 'end', 'width', 'n_win', 'n_pos', 'w_rank', 'rank']
COLUMN_TYPE   = [str  ,  str        , float     , float     , int  , int      , float  , float  , float    , int    , int  , int   , int     , int,     int,      int]
FORMAT_ROW = '%s\t%s\t%.2e\t%.2e\t%d\t%.2f\t%.2e\t%.2e\t%.2f\t%+05d\t%+05d\t%d\t%d\t%d\t%d\t%d'


########################################
#
# COMMAND LINE
#
########################################

description = '''ORM is a command line software that searches for locally overrepresented motifs in a set of DNA sequences.
[ Matthieu Defrance defrance@scmbb.ulb.ac.be ]
'''
version = '2007.12'
usage = """
        %prog -b bgmodel [-i inputfile]
                [-o outputfile] 
                [--max parameter #][--min parameter #]
                [-h | --help]

        Example : %prog  -i seq.fasta --markov=2 -v2 --max rank 10
"""

parser = optparse.OptionParser(usage=usage, version=version, description=description)
parser.add_option("-s", "--strand", dest="strand", help="strand + or +- (default=+-)", choices=['+', '+-'], default='+-')
parser.add_option("-p", "--overlap", dest="overlap", help="allow overlapping word occurences", action="store_true", default=False)
parser.add_option("-l", "--length", dest="l", help="oligomer length", action="store", type="int", default=6)

#BG
parser.add_option("-m", "--markov", dest="markov", help="use markov model of order n", action="store", type="int", default=-1)

parser.add_option("--bgfile", action="store", dest="bgfile", type="string", default=None)
parser.add_option("--bgoligo", action="store", dest="bgoligo", type="string", default=None)
parser.add_option("--bgoligomarkov", action="store", dest="bgoligomarkov", type="string", default=None)

parser.add_option("--max", action="append", dest="max", type="string", nargs=2, help="upper threshold example: --max w_win 1", default=[])
parser.add_option("--min", action="append", dest="min", type="string", nargs=2, help="lower threshold example: --min occ_sig 0", default=[])
parser.add_option("--sort", action="append", dest="sort", type="string", default=[], help="sort ouput with this criteria name. Prefix the name with a + (or -) for growing (or not) order")
parser.add_option("-v", "--verbosity", action="store", dest="verbosity", type="int", default=0)

parser.add_option("--window", action="store", dest="window", help="use fixed window width", type='int', default=None)
parser.add_option("--bgwindow", action="store", dest="bgwindow", help="widow size in bg model", type='int', default=None)

parser.add_option("--location", dest="location", help="region to scan (example -2000:-1)", action="store", type="string", default=None)
parser.add_option("-o", "--output", action="store", dest="output", default=sys.stdout)
parser.add_option("-i", "--input", action="store", dest="input", default=sys.stdin)
parser.add_option("-r", "--ratio", action="store", dest="ratio", type="float", default=2.0)
parser.add_option("--heuristic", action="store", dest="heuristic", type="string", default='slices', help='heutistic for extracting windows (form slower to faster): all, slices, score (default=slices)')
parser.add_option("--count", action="store", dest="count", type="string", default='hash', help='method for counting oligos (hash or tree)')
parser.add_option("--spacing", dest="spacing", help="spacing range (example 0:10)", action="store", type="string", default='1:1')
parser.add_option("-e", "--error", action="store", dest="error", type="int", default=0)
parser.add_option("--slices", action="store", dest="slices", type="int", help="number of slices to use (only when --heuristic=slice)", default=10)
parser.add_option("--right", dest="right", action="store", type="int", help="use sequence right bound as reference (ok for upstream)", default=None)
parser.add_option( "--left", dest="left", action="store", type="int", help="use sequence left bound as reference (ok for downstream)", default=None)

(options, args) = parser.parse_args()
cli.VERBOSITY = options.verbosity


########################################
#
# RUN FUNCTION
#
########################################

def run(args, options):
    timer = Core.Sys.timer.Timer()
    spacing = ( int(options.spacing.split(':')[0]), int(options.spacing.split(':')[1]) )


    ##
    #
    # input sequence
    #
    ##
    if options.left == None and options.right == None:
        options.right = -1

    s = Bio.sequence.fasta2sequences(options.input, rightPosition=options.right, leftPosition=options.left)

    # location
    if options.location:
        location = ( int(options.location.split(':')[0]), int(options.location.split(':')[1]) ) 
    else:
        location = s.location
        #location = bg.location

    if not Core.types.interval.is_included(location, s.location):
        cli.warning('Sequences are shorter then the given location')
        location = s.location

    ##
    #
    # background
    #
    ##
    options.bgwindow = options.bgwindow or 1000000
    if options.bgwindow > s.location[1] - s.location[0] + 1:
        options.bgwindow = s.location[1] - s.location[0] + 1

    if options.bgfile:
        bg = op.oload(options.bgfile)
    else:
        bgparams = dict(spacing=spacing, count=options.count, error=options.error)
        bg = ormlib.Bg(location=location, W=options.bgwindow, l=options.l, strand=options.strand, overlap=options.overlap, params=bgparams)
        if options.markov >= 0:
            bg.build_markov(s, order=options.markov)
        elif options.bgoligo:
            bg.build_from_oligo_file(options.bgoligo, location)
        elif options.bgoligomarkov:
            bg.build_markov_from_oligo_file(options.bgoligomarkov, location)            
        else:
            bg.build(s)

    if options.window > s.location[1] - s.location[0] + 1:
        options.window = s.location[1] - s.location[0] + 1

    if not Core.types.interval.is_included(location, bg.location):
        cli.error('Given location do not match Background model location')
        return


    ##
    # Thresholds
    defaults = {'width' : (2, 1000000),
                'occ'   : (1, 100000),
                'occ_sig'   : (None, None),
                'occ_E' : (None, None),
                'occ_P' : (None, None),
                'start' : (None, None),
                'end'   : (None, None)
                }
    thresholds = ormlib.convert_thresholds(defaults, options.min, options.max, COLUMN_HEADER, COLUMN_TYPE)

    params = {
              'ratio' : options.ratio,
              'spacing' : spacing,
              'binomial' : 1,
              'error' : options.error,
              'spacing' : spacing,
              'heuristic' : options.heuristic,
              'count'     : options.count,
              'slices'    : options.slices,
              'window' : options.window
              }
    params.update(thresholds)


    ##
    #
    # Run ORM
    #
    ##
    R = ormlib.extract_windows(s, bg, location=location, params=params)
    r = R['R']

    r = ormlib.select(r, thresholds, COLUMN_HEADER)
    options.sort = options.sort or ['-occ_sig']
    r = ormlib.sort(r, options.sort, COLUMN_HEADER)
    r = ormlib.window_rank(r, COLUMN_HEADER)
    thresholds = ormlib.convert_thresholds({'w_rank' : (None, None)}, options.min, options.max, COLUMN_HEADER, COLUMN_TYPE)
    r = ormlib.select(r, thresholds, COLUMN_HEADER)
    r = ormlib.rank(r, COLUMN_HEADER)
    thresholds = ormlib.convert_thresholds({'rank' : (None, None)}, options.min, options.max, COLUMN_HEADER, COLUMN_TYPE)
    r = ormlib.select(r, thresholds, COLUMN_HEADER)


    if options.window != None:
        windowWidth = '%d' % options.window
    else:
        windowWidth = 'variable'

    ##
    #
    # output
    #
    ##
    output = Core.utils.uopen(options.output, mode='w')


    print >> output, ormlib.format_header(HEADER, l=options.l, 
                                          inputfilename=str(options.input), 
                                          numberOfSequences=len(s), 
                                          start=location[0], end=location[1],
                                          command=' '.join(sys.argv),
                                          windowWidth = windowWidth,
                                          bgwindowWidth = bg.W,
                                          strand=options.strand,
                                          lowerThreshold = str(options.min),
                                          upperThreshold = str(options.max),
                                          sort = str(options.sort),
                                          date = time.ctime(),
                                          runningTime = str(timer),
                                          scannedWords = R['scannedWords'],
                                          )

    print >> output, ormlib.format_output(r, COLUMN_HEADER, COLUMN_HEADER_CHAR, FORMAT_ROW)


########################################
#
# 
#
########################################

if __name__ == '__main__':
    try:
        #if len(sys.argv) == 1:
        #    parser.print_help()
        #else:
        run(args, options)
        #import Core.Devel.devel as devel
        #devel.prof('run(args, options)')
    except:
        #import traceback
        #traceback.print_exc(file=sys.stdout)
        #cli.error('Error while running orm')
        #parser.error('while running')
        if DEBUG:
            raise
