#!/usr/bin/env python

"""A gibbs sampler that optimize the information content.
Usage: infoGibbs [-w] -i [options]

Required :
    -w / --width
    	Give the width (size) of the ma
    -i / --input
	Give the input sequences file (FASTA format is required)
Options:
    -h / --help
        Print this message and exit.
    -c / --create [your_frequence_file_name]
    	Precise that you create the frequence file related to your input file
    -f / --frequence [your_frequence_file_name]
    	Give the frequence file if it exists
    -g / --gibbs
    	Give the gibbs version you want to use :
    		1 : Gibbs which optimize the log-likelihood
    		2 : Gibbs which optimize the information content	
    -n / --number
    	Give the the number of iteration of the gibbs algorithm
    -e / --exp
    	Give the value of the exposant you want to use for discriming the infoContent Value during the update    	
    -s / --seed
    	Give the seed for the initialisation of the matrix
    -o / --output
    	Give the output fileName to get the result in a file
    -v / --verbosity
    	Give the verbosity level you want :
    		3 : Normal verbosity
    		4 : ...
    		Todo : Complete it

Defaut Options :
    - Create a frequence file of order 0
    - GibbsVersion : 2
    - Number iterations : 3000
    - Verbosity : 1
    - Width : 8
    - Seed : a random number between [0 and 1000000] (generate by random.unifrom([range]))
    - Exposant : 20 for inforMationContent and 1 for Weight
    
Remark :
    The order of the arguments is important
"""


import orig_sequence
import markov
import freq2background
import motif
import gibbsInform
import gibbs93
import time
import sys

import getopt

def usage(code, msg=''):
    print >> sys.stderr, __doc__
    if msg:
        print >> sys.stderr, msg
    sys.exit(code)
    
def main():

#****************************************
#		Read the options	*
#****************************************
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:co:f:g:v:n:w:s:e:r:", ["help", "input=","output=", "create", "frequence=", "gibbs=", "verbosity=", "number=","width=","seed=","exp="])
    
	except getopt.error, msg:
		usage(1, msg)
	
	alph = ['a','c','g','t']
	
	width = 8
	nbIterations = 3000
	ordre = 1 
	freqFile = "resultat.freq"
	inputFile = "input.fasta"
	create = 0
	verbose = 1
	endProgram = 0
	gibbsVersion = 2
	matrixSeed = -1
	exp = -1
	res = ''
	
	s = ''
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage(0)
			sys.exit()
			
		if opt in ("-w", "--width"):
			width = int(arg)
			
		if opt in ("-i", "--input"):
			inputFile = arg
			s = orig_sequence.OrigSequences(alph, width, inputFile)
			
		if opt in ("-c", "--create"):
			tmp = inputFile.split(".fasta")
			freqFile = tmp[0] + "_markov_degre_"+ str(ordre) +".freq"
			mar = markov.Markov(ordre,alph, freqFile)	
			mar.createFrequences(s)
			mar.writeFrequences()
		
		if opt in ("-f", "--frequence"):
			freqFile = arg
			
		if opt in ("-g", "--gibbs"):
			gibbsVersion = int(arg)
			
		if opt in ("-n", "--number"):
			nbIterations = int(arg)
		
		if opt in ("-v", "--verbosity"):
			verbose = int(arg)
		
		if opt in ("-s", "--seed"):
			matrixSeed = int(arg)
		if opt in ("-e", "--exp"):
			exp = int(arg)
		if opt in ("-o", "--output"):
			res = arg
	if len(sys.argv) == 1 :
		usage(0)
		sys.exit()
	
	back = freq2background.Background(ordre,alph, freqFile)
	back.freq2Background()
		
    	mot = motif.motif(width,len(s.sequences),s, matrixSeed, back)
		
	
	fileResult = ''
	
	if gibbsVersion == 1 :
		if len(res) > 0 :
			fileResult = open(res,'w')
		if exp == -1 :
			exp = 1
		p = gibbs93.Gibbs93(width,len(s.sequences),alph,s,mot, back, verbose, nbIterations, exp)
			
		if len(res) > 0 :
			completeCommand = ' '.join(sys.argv[0:len(sys.argv)])
			
			fileResult.write('CM\t%s\nXX\n'%(completeCommand))
			fileResult.write('SD\t%s\nXX\n'%(mot.matrixSeed))
			fileResult.write('OP\tweight\nXX\n')
			fileResult.write('AM\t%s\nXX\n'%(exp))
			fileResult.write('IT\t%s\nXX\n'%(nbIterations))
			fileResult.write('IB\t%s\nXX\n'%(100)) #TODO: parametrer !!
			
			cpt = 0
			while cpt < len(s.sequences):
				seqName = s.seqName[cpt]
				size = s.size[cpt]
				fileResult.write('SQ\t%s; %sbp\n'%(seqName,size))
				cpt += 1
			fileResult.write('XX\n')
			fileResult.write('PR \tA:%s; C:%s; G:%s; T:%s\nXX\n\n'%(back.p['a'],back.p['c'],back.p['g'],back.p['t']))
			fileResult.write('//\n')
		
		p.gibbs(s, mot, back, fileResult)
		
		if len(res) > 0 :
			fileResult.close()
		#print 'Pas encore de support pour le Gibbs 93'
	elif gibbsVersion == 2 :
		if len(res) > 0 :
			fileResult = open(res,'w')
		if exp == -1 :
			exp = 20
		p = gibbsInform.GibbsInform(width,len(s.sequences),alph,s,mot,back, verbose, nbIterations, exp)
		
		if len(res) > 0 :
			completeCommand = ' '.join(sys.argv[0:len(sys.argv)])
			
			fileResult.write('CM\t%s\nXX\n'%(completeCommand))
			fileResult.write('SD\t%s\nXX\n'%(mot.matrixSeed))
			fileResult.write('OP\tInfo\nXX\n')
			fileResult.write('AM\t%s\nXX\n'%(exp))
			fileResult.write('IT\t%s\nXX\n'%(nbIterations))
			fileResult.write('IB\t%s\nXX\n'%(100)) #TODO: parametrer !!
		
			cpt = 0
			while cpt < len(s.sequences):
				seqName = s.seqName[cpt]
				size = s.size[cpt]
				fileResult.write('SQ\t%s; %sbp\n'%(seqName,size))
				cpt += 1
			fileResult.write('XX\n')
			fileResult.write('PR \tA:%s; C:%s; G:%s; T:%s\nXX\n\n'%(back.p['a'],back.p['c'],back.p['g'],back.p['t']))
			fileResult.write('//\n')
			
		p.gibbs(s, mot, back, fileResult)
		
		if len(res) > 0 :
			fileResult.close()
	
	if verbose != 1 :
		print 'Seed of the initial matrix :%s'%(mot.matrixSeed)

if __name__ == '__main__':
	t = time.time()
	
	main()
		
	#print time.time()-t
