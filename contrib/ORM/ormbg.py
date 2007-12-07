#!/usr/bin/env python
import os, sys, glob

PYTHONPATH = glob.glob1(os.path.abspath(sys.path[0]), '*.zlib')
[sys.path.insert(1, os.path.join(sys.path[0], path)) for path in PYTHONPATH]
import optparse

import Core.Persistence.op as op
import Core.UI.cli as cli
import Core.utils
import Core.types
import Bio.sequence

import ORMLib.ormlib as ormlib

DEBUG = 1
cli.UPDATE = 1

########################################
#
# COMMAND LINE
#
########################################
version = '2007.12'
usage = '''
    usage: %prog -i SEQ -o BG -w WORD_LENGTH
                     [-o outputfile] 
                     [-h | --help]

Example: ormbg.py -w 8 -i hg-2kbp-w8.bg hg17-rnd10000-2kbp.fa -o human.bg -v2'''
description = '''Compute bacground model for ORM'''

parser = optparse.OptionParser(usage, version=version, description=description)
parser.add_option("-s", "--strand", dest="strand", help="strand + or +- (default=+-)", choices=['+', '+-'], default='+-')
parser.add_option("-p", "--overlap", dest="overlap", help="allow overlapping word occurences", action="store_true", default=False)
parser.add_option("-v", "--verbosity", action="store", dest="verbosity", type="int", default=0)
parser.add_option("-w", "--word_length", dest="l", help="word length", action="store", type="int", default=2)
parser.add_option("-W", "--window", dest="W", help="window size", action="store", type="int", default=200)
parser.add_option("--location", dest="location", help="region to scan (example -2000:-1)", action="store", type="string", default=None)
parser.add_option("-o", "--output", action="store", dest="output", default='out.bg')
parser.add_option("-i", "--input", action="store", dest="input", default=sys.stdin)
parser.add_option("-m", "--markovorder", dest="markovOrder", help="use markov model", action="store", type="int", default=-1)
parser.add_option("--count", action="store", dest="count", type="string", default='hash', help='method for counting oligos (hash or tree)')
parser.add_option("--spacing", dest="spacing", help="spacing range (example 0:10)", action="store", type="string", default='1:1')
parser.add_option("-e", "--error", action="store", dest="error", type="int", help="allow e error in motif AAnAA (default=0)", default=0)
parser.add_option("--right", dest="right", action="store", type="int", help="use sequence right bound as reference (ok for upstream)", default=None)
parser.add_option("--left", dest="left", action="store", type="int", help="use sequence left bound as reference (ok for downstream)", default=None)

(options, args) = parser.parse_args()

cli.VERBOSITY = options.verbosity
	

########################################
#
# RUN FUNCTION
#
########################################
def run(args, options):
    if options.left == None and options.right == None:
        options.right = -1
    s = Bio.sequence.fasta2sequences(options.input, rightPosition=options.right, leftPosition=options.left)

    if options.location:
        location = ( int(options.location.split(':')[0]), int(options.location.split(':')[1]) )
    else:
        location = s.location

    if options.W > s.location[1] - s.location[0] + 1:
        options.W = s.location[1] - s.location[0] + 1

    params = {
        'spacing' : ( int(options.spacing.split(':')[0]), int(options.spacing.split(':')[1]) ),
        'count'     : options.count,
        'error' : options.error,
        }

    bg = ormlib.Bg(location=location, W=options.W, l=options.l, strand=options.strand, overlap=options.overlap, params=params)

    if options.markovOrder >= 0:
        bg.build_markov(s, order=options.markovOrder)
    else:
        bg.build(s)

    op.osave(bg, options.output)


if __name__ == '__main__':
    try:
        run(args, options)
    except:
        #import traceback
        #traceback.print_exc(file=sys.stdout)
        cli.error('Error while running ormbg')
        if DEBUG:
            raise


