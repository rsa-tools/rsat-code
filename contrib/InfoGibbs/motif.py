import random
import orig_sequence

#################################################################
#								#
#			Class Motif				#
#								#
#################################################################

class motif :
	def __init__(self,width,height,origSeq,matrixSeed):
		#width will contain the width of the motif
		self.width = width
		#height will contain the number of sequences
		self.height = height
		#seed for the initialisation of the matrix
		self.matrixSeed = matrixSeed
		#indice will contain the beginning of the alignement sequence
		self.indice = []
		
		self.fillIndiceRandom(origSeq)
		
	def fillIndiceRandom (self, origSeq):
		if self.matrixSeed != -1 :
			random.seed(self.matrixSeed)
		else :
			self.matrixSeed = int(random.randrange(0,1000000))
			random.seed(self.matrixSeed)
		
		for i in range(0,self.height) :
			if origSeq.size[i] != self.width :
				tmpSeq = origSeq.sequences[i]
				find = 0
				
				while find == 0:
					number = random.randrange(0,origSeq.size[i]-self.width,1)
					if tmpSeq[number:number+self.width].find('n') < 0 :
						find = 1 
						
				self.indice.append(number)
			else :
				self.indice.append(0)
