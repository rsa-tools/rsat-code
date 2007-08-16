import sys
import markov

import time
#################################################################
#								#
#			Class OrigSequences			#
#								#
#################################################################
	
class OrigSequences :

	def __init__(self,alphabet,width, inputFile):
	#sequences is a list of sequences (input sequences)
		self.sequences = []
	#size is a list of the sequences's sizes
		self.size = []
	#indice is the indice of the beginning of motif sequences
	#	self.indice = []
		self.alphabet = alphabet
		self.width = width
		self.pseudocount = []
		self.inputFile = inputFile
		self.seqName = []
	# function to readA fasta file and then create the structure which will contain the sequences		
		self.readFasta()
	
	# function to count the number of different sequences of a fasta file
	def countSeq (self):
		cpt = 0
		f = open(self.inputFile,'r')
		for line in f:
			if line.find('>') >= 0:
				cpt += 1
		print cpt
		
	# function to read a fasta file and fill the structure
	def readFasta (self):
		try:
			# the input file (fasta) is on the standard input
			f = open(self.inputFile,'r')
			try :
				l = f.readline()
				#If the file don't begin with a sequence
				while l.find('>') < 0 :
					l = f.readline ()
					
				eof = 1
				
				l = l.lstrip('>')
				l = l.rstrip('\n')
				self.seqName.append(l)
				
				
				while eof > 0 :
					line = f.readline()
					sequence = ''
					#for a sequence
					while line.find('>') < 0 and len(line) > 0 :
						#we remove the '\n'
						line = line.rstrip('\n')
						#concatenate sequence with the line
						sequence += line
						line = f.readline()
					line = line.lstrip('>')
					line = line.rstrip('\n')
					self.seqName.append(line)
					
					#if we have a sequence we add it to our sequences list
					#we note the size of the sequence too
					sequence = sequence.lower()
					splitList = sequence.split('n')
					acceptSeq = 0
					for elem in splitList :
						if len(elem) >= self.width :
							acceptSeq = 1
					if acceptSeq == 1:
						self.sequences.append(sequence)
						self.size.append(len(sequence))
					#we are at the end of file
					if len(line) == 0:
						eof = 0
			#File empty			
			except IOError, (errno, strerror):
				print "I/O error(%s): %s" % (errno, strerror)
		#The path to the file is bad		
		except IOError, (errno, strerror):
			print "I/O error(%s): %s" % (errno, strerror)
	
	def pseudoCount (self):
		for i in range(0,len(self.alphabet)) :
			self.pseudocount.append(0)
			for seq in self.sequences :
				self.pseudocount[i] += seq.count(self.alphabet[i])
	
#if __name__ == '__main__':
#	alph = ['a','c','g','t']
#	t = time.time()
#	s = OrigSequences(alph)
#	s.pseudoCount()
    	#m = motif.motif(12,len(s.size),s)
#    	print ('Temps calculs : ',time.time()-t)
