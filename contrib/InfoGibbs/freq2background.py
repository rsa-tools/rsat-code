import itertools
import sys
import time

class Elem :
	def __init__(self):
		self.occurences = 0
		self.frequence = 0
		
class Background :

	def __init__(self,degre,alphabet, freqFile):

		self.degre = degre
		self.alphabet = alphabet
		
		self.subPrefixes = []
		self.subMapping = dict()
		
		self.frequences = []
		self.prefixes = []
		self.mapping = dict()
		
		self.p = dict()
		self.freqFile = freqFile
		
	
	def freq2Background (self) :

		if self.readFileFreq() == 1:
			self.initSubPrefixes ()
			self.createMapping ()
			self.createSubMapping ()
			self.createP ()
		
	def readFileFreq (self) :

		try :
			#print self.freqFile,'********'
			f = open(self.freqFile,'r')
			#f = open('resultat.freq', 'r')
			for lines in f:
				#print lines
				lines = lines.rstrip('\n')
				tmp = lines.split('\t')
				if len(tmp) == 3 :
					self.prefixes.append(tmp[0])
					self.frequences.append(Elem())
					self.frequences[len(self.prefixes)-1].occurences = float(tmp[1])
					self.frequences[len(self.prefixes)-1].frequence = float(tmp[2])
				else :
					print "Error : bad file format"
					return 0
					
			return 1		
		#Wrong path to the file		
		except IOError, (errno, strerror):
			print "I/O error(%s): %s" % (errno, strerror)	
	
	def initSubPrefixes (self) :
		s = []
		for i in range(0,self.degre-1):
			s.append(self.alphabet[0])	
		self.createSubPrefixes(s,0)
		
	def createSubPrefixes (self,s,c) :
		if c == self.degre-1:
			prefixe = ''
			for i in range(0,len(s)):
				prefixe += s[i]
				
			self.subPrefixes.append(prefixe)
		else :
			i= 0
			while i < len(self.alphabet) :
				s[c] = self.alphabet[i]
				self.createSubPrefixes(s,c+1)
				i += 1
				
	def createMapping (self):
		nbPrefixes = len(self.prefixes)
		cpt = 0
		indice = []
		while cpt < nbPrefixes :
			indice.append(cpt)
			cpt += 1
		self.mapping = dict(itertools.izip(self.prefixes,indice))

	
	def createSubMapping (self):
		nbPrefixes = len(self.prefixes)
		frequences = []
		i = 0
		while i < nbPrefixes:
			j = 0
			tmp = 0
			while j<len(self.alphabet):
				tmp += self.frequences[i+j].occurences
				j += 1
			frequences.append(tmp)
			i += len(self.alphabet)
			
		self.subMapping = dict(itertools.izip(self.subPrefixes,frequences))
	
	def createP (self):
		frequences = []
		for elem in self.prefixes :
			frequences.append(self.getProbSub(elem))
		self.p = dict(itertools.izip(self.prefixes,frequences))
			
	def getProbSub (self,sub):
		occurence = self.frequences[self.mapping[sub]].occurences
		probA = self.subMapping[sub[0:len(sub)-1]]
		return float(occurence)/float(probA)
		
		
	def printMatTrans(self):
		s = '\t' + 'a' + '\t\t' + 'c' + '\t\t' + 'g' + '\t\t' +'t'
		print s
		for elem in self.subPrefixes :
			s = ''
			s += elem +'\t'
			for lettre in self.alphabet:
				sequence = elem+lettre
				s += str(self.getProbSub(sequence)) + '\t'
			print s
			
#if __name__ == '__main__':
#    	alph = ['a','c','g','t']
    	
#	t = time.time()
	
#	back = Background(2,alph)
#	back.freq2Background()
	#back.printMatTrans()
#	print ('Temps calculs : ',time.time()-t)
