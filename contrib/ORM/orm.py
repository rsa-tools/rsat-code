#!/usr/bin/env python
'''NAME
        ORM

VERSION
        %(version)s

AUTHOR
        Matthieu Defrance (defrance@scmbb.ulb.ac.be)

DESCRIPTION
        compute oligomer frequencies in a set of sequences,
        and detects locally overrepresented oligomers.

CATEGORY
        sequences
        pattern discovery'''
VERSION = '2007.12.11'
USAGE = '''orm.py -l OLIGOMER_LENGTH [-i inputfile]
                [-o outputfile] 
                [--max parameter #][--min parameter #]
                [-h | --help]

        Example : %prog -i seq.fasta -l 6 --markov=2 -v2 --max rank 10
'''
import os, sys, glob
[sys.path.insert(1,os.path.join(sys.path[0], path)) for path in glob.glob1(os.path.abspath(sys.path[0]),'*.zlib')]

import optparse
import time
import Core.timer
import Core.Persistence.op as op
import Core.UI.cli as cli
import Core.types
import Core.IO.utils
import Bio.sequence
from Bio.scripthelper import *
import ORMLib.ormlib as ormlib

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
; version                          %(version)s
; date                             %(date)s
; running time                     %(runningTime)s
; oligomer length                  %(l)d
; strand                           %(strand)s
; window width                     %(windowWidth)s
; bg window width                  %(bgwindowWidth)s            
; input file                       %(inputfilename)s
; bg file                          %(bgfile)s   
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
parser = optparse.OptionParser(usage=USAGE, description=globals()['__doc__'] % {'version': VERSION}, version=VERSION, add_help_option=1, formatter=Formatter())

parser.add_option("-l", "--length", action="store", dest="l", type="int", default=None, help="set oligomer length to L. REQUIRED ARGUMENTS!\nEXAMPLE: --length=7",)
parser.add_option("-i", "--input", action="store", dest="input", default=sys.stdin, metavar="FILE", help="read sequence from FILE (must be in FASTA format)\nif not specified, the standard input is used.",)
parser.add_option("-o", "--output", action="store", dest="output", default=sys.stdout, metavar="FILE", help="output results to FILE\nif not specified, the standard output is used.")
parser.add_option("-s", "--strand", dest="strand", choices=['+', '+-'], default='+-', help="search in foward strand (+) or in both strands (+-)\nDEFAULT:%default\nEXAMPLE: --strand=+",)
parser.add_option("-p", "--overlap", dest="overlap", action="store_true", default=False, help="allow overlapping word occurences\nEXAMPLE: --overlap")

#BG
parser.add_option("-m", "--markov", dest="markov", help="use a Markov model of order ORDER as background model. Markov model parameters are computed from the input sequences\nEXAMPLE: --markov=2 (Markov chain of order 2 generated with input sequences)", action="store", type="int", metavar="ORDER", default=-1)

parser.add_option("--bgfile", action="store", dest="bgfile", type="string", help="background file to use\nEXAMPLE --bgfile=mybgfile", default=None)
parser.add_option("--bgoligo", action="store", dest="bgoligo", type="string", help="use oligo-analysis background model\nEXAMPLE: --bgolio=myfile.gz", default=None)
parser.add_option("--bgoligomarkov", action="store", dest="bgoligomarkov", type="string", help="use a Markovian background model loaded from an oligo-analysis file\nEXAMPLE: --bgoligomarkov=myfile.gz", default=None)

parser.add_option("--max", action="append", dest="max", type="string", nargs=2, help="limit output to items that have PARAM <= VALUE\nEXAMPLE: --max rank 10\nSupported parameters: "+ str(','.join(COLUMN_HEADER)), metavar="PARAM VALUE", default=[])
parser.add_option("--min", action="append", dest="min", type="string", nargs=2, help="limit output to items that have PARAM >= VALUE\nEXAMPLE: --min occ_sig 0\nSupported parameters:" + str(','.join(COLUMN_HEADER)), metavar="PARAM VALUE", default=[])
parser.add_option("--sort", action="append", dest="sort", type="string", default=[], help="sort ouput according to parameter PARAM in growing order (+) or inverse (-)\nEXAMPLE: --sort=+label\nSupported parameters:" + str(','.join(COLUMN_HEADER)), metavar="[+][-]PARAM")

parser.add_option("--window", action="store", dest="window", help="use fixed window width of length LENGTH EXAMPLE: --window=20 (use a window of length 20)", type='int', metavar="LENGTH", default=None)
parser.add_option("--windowgroup", action="store", dest="window_group", help="use fixed window width of length LENGTH, 2*LENGTH, ... ,N*LENGTH\nEXAMPLE: --window=20 (use a window of length 20)", type='int', metavar="LENGTH", default=None)

parser.add_option("--bgwindow", action="store", dest="bgwindow", help="use a widow size of length LENGTH in background model EXAMPLE: --bgwindow=200 (use a background window of length 200)", type='int', metavar="LENGTH", default=None)
parser.add_option("--slices", action="store", dest="slices", type="int", help="number of slices to use\nThis option is relevant only when heuristic is set to slices\nDEFAULT: %default", default=10)
parser.add_option("--right", dest="right", action="store", type="int", help="use sequence right bound position as reference position POSITION. This should be used for upstream sequences.\nEXAMPLE: --right=-1 (use right bound of input sequence as position -1)", metavar="POSITION", default=None)
parser.add_option( "--left", dest="left", action="store", type="int", help=" use sequence left bound position as reference position POSITION. This should be used for downstream sequences.\nEXAMPLE: --left=0 (use left bound of input sequence as position 0)", metavar="POSITION", default=None)

parser.add_option("--location", dest="location", action="store", type="string", default=None, help=optparse.SUPPRESS_HELP)
parser.add_option("-r", "--ratio", action="store", dest="ratio", type="float", default=2.0, help=optparse.SUPPRESS_HELP)
parser.add_option("--heuristic", action="store", dest="heuristic", choices=['all', 'slices', 'score'], default='slices', help='heutistic for extracting windows. Supported heutistics: all, slices, score\nDEFAULT: %default\nEXAMPLE: --heuristic=all')
parser.add_option("--count", action="store", dest="count", type="string", default='hash', help=optparse.SUPPRESS_HELP)
parser.add_option("--spacing", dest="spacing", action="store", type="string", default='1:1', help=optparse.SUPPRESS_HELP)
parser.add_option("-e", "--error", action="store", dest="error", type="int", default=0, help=optparse.SUPPRESS_HELP)
parser.add_option("--debug", action="store_true", dest="debug", help=optparse.SUPPRESS_HELP)

parser.add_option("-v", "--verbosity", action="store", dest="verbosity", type="int", help="Set verbosity to level LEVEL\nEXAMPLE: --verbosity=2", metavar="LEVEL", default=0)
(options, args) = parser.parse_args()

cli.VERBOSITY = options.verbosity
if options.debug:
    cli.VERBOSITY = 10

########################################
#
# RUN FUNCTION
#
########################################

def run(args, options):
    #if options.help:
    #    print globals()['__doc__'] % {'usage' : USAGE, 'version' : VERSION}
    #    return

    if options.l == None:
        parser.error("options -l is required. You must provide at least an oligomer length")


    timer = Core.timer.Timer()
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
    bgfile = ''
    options.bgwindow = options.bgwindow or 1000000
    if options.bgwindow > s.location[1] - s.location[0] + 1:
        options.bgwindow = s.location[1] - s.location[0] + 1

    if options.bgfile:
        bg = op.oload(options.bgfile)
        bgfile = options.bgfile
    else:
        bgparams = dict(spacing=spacing, count=options.count, error=options.error)
        bg = ormlib.Bg(location=location, W=options.bgwindow, l=options.l, strand=options.strand, overlap=options.overlap, params=bgparams)
        if options.markov >= 0:
            bg.build_markov(s, order=options.markov)
            bgfile = 'Markov order %d from input' % options.markov
        elif options.bgoligo:
            bg.build_from_oligo_file(options.bgoligo, location)
            bgfile = options.bgoligo
        elif options.bgoligomarkov:
            bg.build_markov_from_oligo_file(options.bgoligomarkov, location)            
            bgfile = options.bgoligomarkov
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
              'window' : options.window,
              'window_group' : options.window_group
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
    output = Core.IO.utils.uopen(options.output, mode='w')


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
                                          version = VERSION,
                                          bgfile = bgfile
                                          )

    print >> output, ormlib.format_output(r, COLUMN_HEADER, COLUMN_HEADER_CHAR, FORMAT_ROW)


########################################
#
# 
#
########################################

if __name__ == '__main__':
    try:
        run(args, options)
        #import Core.Devel.devel as devel
        #devel.prof('run(args, options)')
    except KeyboardInterrupt:
        sys.stderr.write('\n')
        sys.stderr.flush()
        sys.exit()
    except SystemExit:
        pass
    except:
        cli.error('Fatal Error')
        if options.debug:
            raise
            #import traceback
            #traceback.print_exc(file=sys.stdout)
