import itertools
import sys
import time
import orig_sequence

#################################################################
#								#
#			Class ElemMarkov			#
#								#
#################################################################

#this class constitute elem of the markov class
class ElemMarkov :
	def __init__(self):
		self.occurences = 0
		self.frequence = 0


#################################################################
#								#
#			Class Markov				#
#								#
#################################################################

#this class will calculate for a certain order and an alphabet, the markovian background model
class Markov :
	def __init__(self,degre,alphabet, freqFile):
		#degre is the order of the markov chain
		self.degre = degre
		self.alphabet = alphabet
		
		#frenquences is a list of ElemMarkov
		self.frequences = []
		#prefixes will contained all the prefixes for the markov chain
		self.prefixes = []
		
		#mapping is a dictionnary which will do the mapping between the prefixes in 
		#the sequences and the indices in frequences
		self.mapping = dict()
		
		self.freqFile = freqFile
		#createPrefixes is a function which will calculate all the prefixes for the markov chain
		self.createPrefixes()
		#createMapping is a function which will create the mapping
		self.createMapping()
		#initFrequences is a function which will create the list frequences with the correct number of elements
		self.initFrequences()
		
		
	def createPrefixes (self):
		prefixes = []
			
		#s is the prefixe and we constitute the first one, the prefixe only constituted by the first
		#alphabet letter
		s = []
		for i in range(0,self.degre):
			s.append(self.alphabet[0])
			
		self.recurPrefixes(s,0)
	
	#this recursive function will calculate all the pefixes	
	def recurPrefixes (self, s, c):
		
		#if this test is successfull, we have a new refixe
		if c == self.degre:
			prefixe = ''
			#we concatenate the different letter of the prefixe
			for i in range(0,len(s)):
				prefixe += s[i]
			#we add the prefixe
			self.prefixes.append(prefixe)
		else :
			i= 0
			while i < len(self.alphabet) :
				#we generate the new letters for the position c
				s[c] = self.alphabet[i]
				#this is the recursive call
				self.recurPrefixes(s,c+1)
				i += 1 
	
	def createMapping (self):
		nbPrefixes = len(self.prefixes)
		indices = []
		
		for i in range(0,nbPrefixes) :
			indices.append(i)
		#we are giving an arbitrary order too each prefixes
		self.mapping = dict(itertools.izip(self.prefixes,indices))
		
	def initFrequences (self) :
		nbPrefixes = len(self.prefixes)
		for i in range(0,nbPrefixes):
			self.frequences.append(ElemMarkov())
		
	#this function will calculate the frequences and occurences of each prefixe	
	#it needs the original sequences
	def createFrequences (self, origSequences):

		nbSeq = len(origSequences.size)
		nbPrefixes = len(self.prefixes)
		sumOccurences = 0.0
		
		for i in range(0,nbSeq):
			#the part of the sequence that we must considered (the "degre-1" last elements will
			#never be considered
			size = origSequences.size[i]-self.degre+1
			s = origSequences.sequences[i]
			#we only use the lower case (important for the mapping)
			s = s.lower()
			for j in range(0,size):
				#sub is a prefixe of size degre
				sub = s[j:j+self.degre]
				#passed the prefixes which contained 'n'
				#if this test is successfull, the prefixes don't containted 'n' character
				if sub.find('n') < 0 :
					self.frequences[self.mapping[sub]].occurences += 1
				#if the prefixe contained 'n', we passed the 'degre' prefixes which will contained 'n' too
				else:
					j += self.degre	
					
		#we calculate the total number of prefixes in the sequences
		for i in range(0,nbPrefixes):
			sumOccurences += self.frequences[i].occurences
			
		#for i in range(0,len(self.prefixes)):
		#	print (self.prefixes[i], self.frequences[i].occurences)
		#print sumOccurences
		
		#we calculate the frequence of each prefixe
		for i in range(0,nbPrefixes):
			self.frequences[i].frequence = float(self.frequences[i].occurences) / sumOccurences
			
	#this function write the results in a file	
	def writeFrequences (self) :
		
		try :
			o = open(self.freqFile,'w') 
			cpt = 0
			#the elements are tab-delimited
			for elem in self.prefixes :
				s = elem + '\t' + str(self.frequences[cpt].occurences) + '\t' + str(self.frequences[cpt].frequence) + '\n'
				o.write(s)
				cpt += 1
			o.close()
		except IOError, (errno, strerror):
			print "I/O error(%s): %s" % (errno, strerror)
			
#if __name__ == '__main__':
#	t = time.time()
#	alph = ['a','c','g','t']
#	s = orig_sequence.OrigSequences(alph,12)
    	#m = motif.motif(12,len(s.size),s)
#    	print ('Temps calculs : ',time.time()-t)
    	
#    	t = time.time()
#	mar = Markov(2,alph)
#	mar.createFrequences(s)
#	mar.writeFrequences()
#	print ('Temps calculs : ',time.time()-t)
