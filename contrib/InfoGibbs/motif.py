import random
import orig_sequence

#################################################################
#								#
#			Class Motif				#
#								#
#################################################################

class motif :
	def __init__(self,width,height,origSeq,matrixSeed, back):
		#width will contain the width of the motif
		self.width = width
		#height will contain the number of sequences
		self.height = height
		#seed for the initialisation of the matrix
		self.matrixSeed = matrixSeed
		#indice will contain the beginning of the alignement sequence
		self.indice = []
		
		#self.fillIndiceRandom(origSeq)
		self.fillPropProb(origSeq, back)
		
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
				
	def fillPropProb (self, origSeq, back) :
		if self.matrixSeed != -1 :
			random.seed(self.matrixSeed)
		else :
			self.matrixSeed = int(random.randrange(0,1000000))
			random.seed(self.matrixSeed)

		for i in xrange (self.height) :
			listeProba = []
			listeIndice = []
			#listeProbaNormal = []
			if origSeq.size[i] != self.width :
				tmpSeq = origSeq.sequences[i]
				cpt = 0
				size = len(tmpSeq)-self.width+1
				while cpt < size :
					if tmpSeq[cpt:cpt+self.width].find('n') < 0 :
						prob = 1 
						for j in xrange(self.width) :
							prob = prob * back.p[tmpSeq[cpt+j]]
						listeIndice.append(cpt)
						#listeProbaNormal.append(prob)
						listeProba.append(1.0/prob)
					cpt += 1
					
				if len(listeProba) > 1 :
					maxVal = max(listeProba)
					minVal = min(listeProba)
					if minVal != maxVal :
						saveVal = 0
						listeProbaNorm = []
						listeProbaNormAmpl = []
						listeProbaAmplCum = []
						
						for j in xrange(len(listeProba)) :
							normVal = ((listeProba[j] - minVal)/(maxVal - minVal))
							listeProbaNorm.append(normVal)
							normValAmpl = normVal**3
							listeProbaNormAmpl.append(normValAmpl)
							saveVal += normValAmpl
							listeProbaAmplCum.append(saveVal)
							
						#if i == 0 :
						#	res = open('biaisInit.tab','w')
						#	res.write('Indice\tProb\n')
						#	for k in xrange(len(listeProba)) :
						#		res.write('%s\t%s\n'%(listeIndice[k],listeProbaNormal[k]))
						#	res.close()
							
						randNumber = random.uniform(listeProbaAmplCum[0],listeProbaAmplCum[len(listeProba)-1])
						
						find = 0
						cpt = 0
						while not find :
							if listeProbaAmplCum[cpt] > randNumber :
								find = 1
								self.indice.append(listeIndice[cpt])
							cpt += 1
						
					else :
						indice = listeProba.index(maxVal)
						self.indice.append(listeIndice[indice])
				else :
					self.indice.append(listeIndice[0])
			else :
				self.indice.append(0)
		print 'End Init Positions'
