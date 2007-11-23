#!/usr/bin/env python

'''
Gibbs sampling for Matrix (PSSM) Discovery

'''
import os, sys, glob
#print sys.version
PYTHONPATH = glob.glob1(os.path.abspath(sys.path[0]), '*.zlib')
[sys.path.insert(1, os.path.join(sys.path[0], path)) for path in PYTHONPATH]
from Bio.scripthelper import *

import Lib.markov as markov
import Lib.rand as rand
import Lib.motif as motif
import Lib.matrix as matrix
import Lib.gibbs93 as gibbs93
import Lib.gibbs95 as gibbs95
import Lib.finalcycle as finalcycle
cli.UPDATE = 1
DEBUG = 1

usage = """
        %prog   [-i inputfile]
                [-o outputfile] 
                [-h | -help]

        %prog ...
"""

parser = optparse.OptionParser(usage)
parser.add_option("-v", "--verbosity", action="store", dest="verbosity", type="int", default=1)
parser.add_option("-o", "--output", action="store", dest="output", default=sys.stdout)
parser.add_option("-i", "--input", action="store", dest="input", default=sys.stdin)
parser.add_option("-l", "--motifLength", action="store", dest="motifLength", type="int", help="motif length", default=6)
parser.add_option("--iter", action="store", dest="iter", type="int", help="", default=100)
parser.add_option("--nrun", action="store", dest="nrun", type="int", help="number of run", default=1)
parser.add_option("--words", action="store", dest="words", type="int", help="", default=10)
parser.add_option("--motifs", action="store", dest="motifs", type="int", help="number of motifs to extract", default=1)
parser.add_option("--minwords", action="store", dest="minwords", type="int", help="", default=0)
parser.add_option("--maxwords", action="store", dest="maxwords", type="int", help="", default=0)
parser.add_option("--dir", action="store", dest="dir", type="str", default="")
parser.add_option("--type", action="store", dest="type", type="str", default="gibbs")
parser.add_option("--top", action="store", dest="top", type="int", default=1)
parser.add_option("--EM", action="store", dest="EM", type="int", help="number of EM iteration", default=10)
parser.add_option("--alpha", action="store", dest="alpha", type="float", default=1)
parser.add_option("--seed", action="store", dest="seed", type="int", help="set random seed", default=None)
parser.add_option("--fulloutput", action="store_true", dest="fulloutput", default=False)
parser.add_option("--bg", action="store", dest="bg", type="str", default=None)
parser.add_option("--startfrom", action="store", dest="startfrom", type="str", default=None)
parser.add_option("--sampling", action="store", dest="sampling", type="str", default="llr")
parser.add_option("--normalize", action="store", type="int", dest="normalize", default=0)
parser.add_option("--strand", action="store", type="str", dest="strand", default='+-')
parser.add_option("--shift", action="store", type="int", dest="shift", default=10, help="try to shift everyn n iteration")
(options, args) = parser.parse_args()

cli.VERBOSITY = options.verbosity


def run(args, options):
    seed = random_seed(options.seed)

    # INPUT
    sequences = read_fasta(options.input, options.strand)

    # PARAMETERS
    parameters = {
        'command' : ' '.join(sys.argv),
        'input' : options.input,
        'iter' :  options.iter,
        'words' : options.words,
        'wordMinLength' : options.motifLength,
        'wordMaxLength' : options.motifLength,
        'l' : options.motifLength,
        'EM' : options.EM,
        'p' : options.words /  float(sum([len(s) for s in sequences])) ,
        'N' : sum([len(s) - options.motifLength +1 for s in sequences]),
        'e' : options.words,
        'type' : options.type,
        'alpha' : options.alpha,
        'seed' : seed,
        'minwords' : options.minwords or options.words,
        'maxwords' : options.maxwords or options.words,
        'sampling' : options.sampling,
        'normalize' : options.normalize,
        'shift' : options.shift,
        }

    # OUTPUT
    p = dict(l=parameters['wordMinLength'], alpha=parameters['alpha'])
    
    output = new_output(options.input, options.output, options.dir, '.mat', p)

    # BG MODEL
    if options.bg:
        mm = markov.oligo2MM(options.bg)

    else:
        mm = markov.MM(0)
        mm.learn(sequences)
        #print mm.priori

    # MAIN
    results = []


    # STARTING POINT
    if options.startfrom != None:
        sites = open(options.startfrom).read().split()
        #print sites
        start = motif.words2motif(sites, sequences, mm, parameters)
        #print start.full_str()


    mask = motif.Mask()


    for m in range(options.motifs):
        currentMotifs = []
        parameters['best'] = None
        parameters['info'] = cli.Info(options.motifs * options.nrun * (options.iter + options.EM) * (options.maxwords - options.minwords +1), 1)
        
        
        for n in range(options.nrun):
            parameters['tab'] = Bio.tabular.Tab(['i', 'bestIC', 'bestWord', 'bestN', 'IC', 'word', 'N', 'time'])
            parameters['clock'] = Core.Sys.timer.Timer() 
            parameters['nrun'] = n+1

            if options.startfrom == None:
                start = motif.init(sequences, mm, mask, parameters)


            if options.type.startswith('95'):
                best = gibbs95.gibbs(sequences, start, mm, parameters)
                results.append(best)
            else:
                for nwords in range(parameters['minwords'], parameters['maxwords']+1):
                    parameters['words'] = nwords
                    best = gibbs93.gibbs(sequences, start, mm, mask, parameters)
        
            #final cycle
            a = best.fitness()
            best = finalcycle.final(best, parameters)
            b = best.fitness()

            if best.fitness() > parameters['best'].fitness():
                 parameters['best'] = best

            if options.fulloutput:
                output1 = new_output(options.input, options.output, options.dir, '-%d.out' % n, parameters)
                print >> output1, parameters['tab'].to_txt(parameters, 0, 1)
                output1.flush()
            
            currentMotifs.append(best)
        sys.stderr.write('\n')

        winner = sorted(currentMotifs, reverse=True)[0]
        results.append(winner)
        #mask
        mask.update(winner.spl)


    results = sorted(results)
    #top = min(len(results), options.top)
    top = options.motifs
    #print >> output, '\n//\n'.join([results[-i].str() for i in range(1, min(len(results)+1, options.top+1))])
    print >> output, '\n//\n'.join([results[-i].full_str() for i in range(1, top+1)])

    #print '\n//\n'.join([results[-i].full_str() for i in range(1, top+1)])





if __name__ == '__main__':
    try:
        run(args, options)
        #import Core.Devel.devel as devel
        #devel.prof('run(args, options)')
    except KeyboardInterrupt:
        sys.stderr.write('\n')
        sys.stderr.flush()
        sys.exit()
    except:
        cli.error('Error while running gibbs')
        if DEBUG:
            raise



