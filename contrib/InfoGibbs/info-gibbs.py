#!/usr/bin/env python

'''NAME
        info-gibbs.py

VERSION
        %(version)s

AUTHOR
        2006-2006 by Matthieu Defrance (defrance@scmbb.ulb.ac.be)

DESCRIPTION
        Gibbs samling algorithm for motifs discovery

CATEGORY
        sequences
        pattern discovery'''
USAGE = '''
        info-gibbs.py   -l MOTIF_LENGTH
                [-i inputfile]
                [-o outputfile] 
                [-h | -help]
'''
VERSION='2007.12.10'

import os, sys, glob
PYTHONPATH = glob.glob1(os.path.abspath(sys.path[0]), '*.zlib')
[sys.path.insert(1, os.path.join(sys.path[0], path)) for path in PYTHONPATH]
from Bio.scripthelper import *

import Bio.markovsimple as markov
import Lib.rand as rand
import Lib.motif as motif
import Lib.matrix as matrix
import Lib.gibbs93 as gibbs93
import Lib.gibbs95 as gibbs95
import Lib.finalcycle as finalcycle
import Core.cache

  
parser = optparse.OptionParser(usage=USAGE, description=globals()['__doc__'] % {'version': VERSION}, version=VERSION, add_help_option=1, formatter=Formatter())
parser.add_option("-l", "--length", action="store", dest="l", type="int", default=None, help="set oligomer length to L. REQUIRED ARGUMENTS!\nEXAMPLE: --length=7",)
parser.add_option("-i", "--input", action="store", dest="input", default=sys.stdin, metavar="FILE", help="read sequence from FILE (must be in FASTA format)\nif not specified, the standard input is used.",)
parser.add_option("-o", "--output", action="store", dest="output", default=sys.stdout, metavar="FILE", help="output results to FILE\nif not specified, the standard output is used.")
parser.add_option("-s", "--strand", dest="strand", choices=['+', '+-'], default='+-', help="search in foward strand (+) or in both strands (+-)\nDEFAULT:%default\nEXAMPLE: --strand=+",)

parser.add_option("--iter", action="store", dest="iter", type="int", default=100, help="maximum number of Gibbs sampling iterations\nDEFAULT: %default")
parser.add_option("--EM", action="store", dest="EM", type="int", default=0, help="number of Expecation Maximization cycles\nDEFAULT: %default")
parser.add_option("--words", action="store", dest="words", type="int", default=10, help="number of motif occurrences that are expected to be found\nDEFAULT: %default")
parser.add_option("--motifs", action="store", dest="motifs", type="int", default=1, help="number of motif to extract\nDEFAULT: %default ")
parser.add_option("--nrun", action="store", dest="nrun", type="int", default=1, metavar="N", help="try to run the Gibbs sampling N times")
parser.add_option("--finalcycle", action="store_true", dest="finalcycle", help="Try to optimise number of words in a final cycle")
parser.add_option("--shift", action="store", type="int", dest="shift", default=10, metavar="N", help="try to shift the motif evey N iteration\nDEFAULT: %default")

parser.add_option("-m", "--markov", dest="markov", action="store", type="int", default=0, metavar="ORDER", 
    help="use Markov model of order ORDER as backgroundarkov model parameters are computed from the input sequence\nEXAMPLE: --markov=2 (Markov chain of order 2) generated with input sequences)\nDEFAULT: %default")
parser.add_option("--bgfile", action="store", dest="bg", type="str", default=None, metavar="FILE", help="use oligo-analysis background model")
parser.add_option("--type", action="store", dest="type", default="93", choices=['93', '95'], help="use the TYPE method. supported TYPE are 93,95")
parser.add_option("--seed", action="store", dest="seed", type="int", default=None, help="use SEED for random number intialization")
parser.add_option("--startwords", action="store", dest="startwords", type="str", default=None, metavar="FILE", help="read words in FILE and use them as a starting point")
parser.add_option("--startmatrix", action="store", dest="startmatrix", type="str", default=None, metavar="FILE", help="read matrix form FILE and use it as a starting point")
parser.add_option("-v", "--verbosity", action="store", dest="verbosity", type="int", help="Set verbosity to level LEVEL\nEXAMPLE: --verbosity=2", metavar="LEVEL", default=0)


parser.add_option("--fulloutput", action="store_true", dest="fulloutput", default=False, help=optparse.SUPPRESS_HELP)
parser.add_option("--top", action="store", dest="top", type="int", default=1, help=optparse.SUPPRESS_HELP)
parser.add_option("--dir", action="store", dest="dir", type="str", default="", help=optparse.SUPPRESS_HELP)
parser.add_option("--minwords", action="store", dest="minwords", type="int", default=0, help=optparse.SUPPRESS_HELP)
parser.add_option("--maxwords", action="store", dest="maxwords", type="int", default=0, help=optparse.SUPPRESS_HELP)
parser.add_option("--alpha", action="store", dest="alpha", type="float", default=1, help=optparse.SUPPRESS_HELP)
parser.add_option("--sampling", action="store", dest="sampling", type="str", default="llr", help=optparse.SUPPRESS_HELP)
parser.add_option("--normalize", action="store", type="int", dest="normalize", default=0, help=optparse.SUPPRESS_HELP)
parser.add_option("--percent", action="store", dest="percent", type="float", default=1.0, help=optparse.SUPPRESS_HELP)

parser.add_option("--debug", action="store_true", dest="debug", help=optparse.SUPPRESS_HELP)

(options, args) = parser.parse_args()

cli.VERBOSITY = options.verbosity
if options.debug:
    cli.VERBOSITY = 10

cli.UPDATE = 1

def run(args, options):
    #if options.help:
    #    print globals()['__doc__'] % {'usage' : USAGE, 'version' : VERSION}
    #    return

    if options.l == None:
        parser.error("options -l is required. You must provide at least a motif length")
        
    seed = random_seed(options.seed)

    # INPUT
    sequences = read_fasta(options.input, options.strand)

    # PARAMETERS
    parameters = {
        'command' : ' '.join(sys.argv),
        'input' : options.input,
        'iter' :  options.iter,
        'words' : options.words,
        'wordMinLength' : options.l,
        'wordMaxLength' : options.l,
        'l' : options.l,
        'EM' : options.EM,
        'p' : options.words /  float(sum([len(s) for s in sequences])) ,
        'N' : sum([len(s) - options.l +1 for s in sequences]),
        'e' : options.words,
        'type' : options.type,
        'alpha' : options.alpha,
        'seed' : seed,
        'minwords' : options.minwords or options.words,
        'maxwords' : options.maxwords or options.words,
        'sampling' : options.sampling,
        'normalize' : options.normalize,
        'shift' : options.shift,
        'percent' : options.percent
        }

    # OUTPUT
    p = dict(l=parameters['wordMinLength'], alpha=parameters['alpha'])
    
    output = new_output(options.input, options.output, options.dir, '', p)

    # BG MODEL
    if options.bg:
        mm = markov.oligo2MM(options.bg)
    else:
        mm = markov.MM(options.markov)
        mm.learn(sequences)
        #print mm.priori

    mmcache = Core.cache.MemoryCache(mm.logP)
    mm.logP_cached = mmcache.__call__

    # STARTING POINT
    if options.startwords != None:
        sites = open(options.startwords).read().split()
        start = motif.words2motif(sites, sequences, mm, parameters)
    if options.startmatrix != None:
        m = matrix.tab2matrix(options.startmatrix)
        start = finalcycle.matrix2motif(m, sequences, mm, parameters)


    # MAIN
    sites = motif.all_sites(sequences, mm, parameters)
    results = []

    for m in range(options.motifs):
        currentMotifs = []
        parameters['best'] = None
        parameters['info'] = cli.Info(options.nrun * 1 * (options.iter + options.EM)  * (options.maxwords - options.minwords +1), 1)
        
        
        for n in range(options.nrun):
            parameters['tab'] = Bio.tabular.Tab(['i', 'bestIC', 'bestWord', 'bestN', 'IC', 'word', 'N', 'time'])
            parameters['clock'] = Core.Sys.timer.Timer() 
            parameters['nrun'] = n+1

            if options.startwords == None and options.startmatrix == None:
                start = motif.init(sequences, mm, sites, parameters)

            # Gibbs sampling
            if options.type.startswith('95'):
                best = gibbs95.gibbs(sequences, start, mm, parameters)
                results.append(best)
            else:
                for nwords in range(parameters['minwords'], parameters['maxwords']+1):
                    parameters['words'] = nwords
                    best = gibbs93.gibbs(sequences, start, mm, sites, parameters)
        
            if options.finalcycle:
                best = finalcycle.final(best, parameters)

            if best.fitness() > parameters['best'].fitness():
                 parameters['best'] = best

            if options.fulloutput:
                output1 = new_output(options.input, options.output, options.dir, '-%d.out' % n, parameters)
                print >> output1, parameters['tab'].to_txt(parameters, 0, 1)
                output1.flush()
            
            currentMotifs.append(best)


        winner = sorted(currentMotifs, reverse=True)[0]
        results.append(winner)
        sites.difference_update(winner.sp())

        print >> output, winner.full_str()
        if m < options.motifs - 1:
            print >> output, '//'
        output.flush()

    #results = sorted(results)
    #top = options.motifs
    #print >> output, '\n//\n'.join([results[-i].full_str() for i in range(1, top+1)])


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



