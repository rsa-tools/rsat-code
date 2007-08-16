import itertools
import time
import math
import random
import heapq

#************************#
#	Class of Gibbs	 #
#************************#
class Gibbs93 :
	def __init__(self,width,height, alphabet, seq, m, back, verbose, nbIterations, exp):
		#width of the motif
		self.width = width
		#height of the matrix (number of sequences)
		self.height = height
		#Totla number of letter
		self.nbTotLetter = 0
		#Sum of the peusdocounts
		self.B = 1
		self.exp = exp
		#self.exp = 1
		#Verbosity
		self.verbose = verbose
		#NbIterations
		self.nbIterations = nbIterations
		#Counts
		self.c = []
		#Probability q of the matrix
		self.q = []
		#Prior probability
		self.p = []
		#Pseudocounts
		self.b = []
		#Alphabet
		self.alphabet = alphabet
		#Mapping between : alphabet<>indice
		self.mapping = dict()
		#Mapping between : indice<>alphabet
		self.inverseMap = dict()
		
		#initializing the c and q matrixes
		
		self.c = [[0]*len(self.alphabet) for i in xrange(self.width)]
		self.q = [[0]*len(self.alphabet) for i in xrange(self.width)]
		self.b = [0]*len(self.alphabet)
		
		for letter in self.alphabet :
			self.p.append(back.p[letter])
			
		self.initNbTotLetter(seq)
		self.initB()
		self.createMapping()
		self.initC(seq,m)
			
	def createMapping (self):
		#fill indice with the indices between 0 and length of the alphabet
		indice = [i for i in xrange(len(self.alphabet))]

		self.mapping = dict(itertools.izip(self.alphabet,indice))
		self.inverseMap = dict(itertools.izip(indice,self.alphabet))
	
	def initNbTotLetter (self,seq):
		for sequence in seq.sequences :
			for i in range(0,len(sequence)):
				if sequence[i] != 'n':
					self.nbTotLetter += 1
	
	def initB (self) :
		#we are calculating the pseudo counts for each letter proportionally to the prior probability
		for i in xrange(len(self.alphabet)) :
			self.b[i] = self.B * self.p[i]
			
	#We initialize the total number of count for each letter in the initial matrix
	def initC (self,seq,m):
		numSeq = 0
		self.c = [[0]*len(self.alphabet) for i in xrange(self.width)]
		#for each indice of the beginning of the motif in each sequence
		for indice in m.indice :
			#we stock the sequence in tmpSeq
			tmpSeq = seq.sequences[numSeq]
			for i in xrange(self.width):
				#for each column i, we are countig the number of letter
				#we use the mapping to retrieve the inidice corresponding to the letter
				self.c[i][self.mapping[tmpSeq[indice+i]]] += 1
			numSeq +=1
	
	#We must suppress the counts of the letters of subsequence sequence subSeq when we are scanning originial sequence of subSeq	
	def subC (self,subSeq):
		#for all the indices of the motif
		for i in range(0,self.width):
			indice = self.mapping[subSeq[i]]
			self.c[i][indice] -= 1
	
	#We must add the counts of the letters of subsequence subSeq when we choose a new subsequence for the matrix
	def addC (self,subSeq):
		#for all the indices of the motif
		for i in range(0,self.width):
			indice = self.mapping[subSeq[i]]
			self.c[i][indice] += 1
		
	#We must calculate the Q_{ij}
	#Formule : Q_{ij} = \frac{c_{ij} + b_j}{N-1 + B}	
	def calculateQ (self):
		for i in range(0,self.width):
			for j in range(0,len(self.alphabet)):
				self.q[i][j] = float(self.c[i][j] + self.b[j])/float(self.height-1 + self.B)		
	
	def chooseSequence (self):
		return random.randrange(0,self.height,1)
	
	def chooseSequencePropToWeight (self, back, seq, m):
		scoreList = self.calculRowScore(back,seq,m)
		mappingScoreIndice = dict(itertools.izip(scoreList,xrange(len(scoreList))))

		scoreList.sort()
		
		mappingScoreOrderedIndice = dict(itertools.izip(xrange(len(scoreList)),scoreList))
		
		if len(scoreList) > 1 :
			maxVal = max(scoreList)
			minVal = min(scoreList)
			
			
			listeScoreNorm = []
			listeScoreNormAmpl = []
			listeScoreAmplCum = []
			saveVal = 0
			for i in xrange(len(scoreList)) :
				normVal = (scoreList[i] - minVal)/(maxVal - minVal)
				listeScoreNorm.append(1-normVal)
				normValAmpl = normVal**self.exp
				listeScoreNormAmpl.append(normValAmpl)
				saveVal += normValAmpl
				listeScoreAmplCum.append(saveVal)
		else :
			listeIseqNorm = []
			listeIseqNormAmpl = []
			listeIseqAmplCum = []
			
			listeIseqNorm.append(1)
			listeIseqNormAmpl.append(1)
			listeIseqAmplCum.append(1)
			
		randNumber = random.uniform(0,listeScoreAmplCum[0])
		
		bestApprox = self.findBestApprox(randNumber, listeScoreAmplCum)
		indice = listeScoreAmplCum.index(bestApprox)
		bestIndice = mappingScoreIndice[mappingScoreOrderedIndice[indice]]
		
		return bestIndice
		
	def calculRowScore (self, back, seq, m):
		weightMat = self.calculWeightMat(self.c, back)
					
		Wmax = 0
		Wmin = 0
		
		for elem in weightMat :
			Wmax += max(elem)
			Wmin += min(elem)
		
		cpt = 0
		scoreList = []
		for sequence in seq.sequences :
			indice = m.indice[cpt]
			tmpScore = 0
			for i in xrange(self.width) :
				tmpScore += weightMat[i][self.mapping[sequence[indice+i]]]
			scoreList.append(tmpScore)
			cpt += 1
		
		for i in xrange(len(scoreList)) :
			scoreList [i] = (scoreList[i] - Wmin)/(Wmax-Wmin)

		return scoreList
	
	def calculWeightMat (self, count, back) :
		weightMat = [[0]*len(self.alphabet) for i in xrange(self.width)]
		
		for i in xrange(self.width) :
			for j in xrange(len(self.alphabet)) :
				if count[i][j] != 0 :
					Q = (count[i][j]+self.b[j])/(self.height+self.B)
					P = back.p[self.inverseMap[j]]
					weightMat[i][j] = Q/P
		return weightMat
	
	def calculTheorScoreDistrib (self, decimals, matrix, prob):
		nbCol = self.width
		nbLetter = len(self.alphabet)
		tronc = 10**(decimals)

		scoreList = []
		probaList = []
		
		scoreList.append(0)
		probaList.append(1.0)
		
		for col in xrange(nbCol) :
			currentScoreList = []
			currentProbaList = []
			for letter in xrange(nbLetter) :
				probLetter = prob[self.inverseMap[letter]]
				scoreResidu = matrix[col][letter]*tronc
				
				for prevScore in scoreList :
					currentScore = prevScore+scoreResidu
					
					currentScore = currentScore
					currentScore = math.floor(currentScore)
					indicePrevious = scoreList.index(prevScore)
					
					try :
						indice = currentScoreList.index(currentScore)
						currentProbaList[indice] += probaList[indicePrevious]*probLetter
					except ValueError :
						currentScoreList.append(currentScore)
						currentProbaList.append(probaList[indicePrevious]*probLetter)
			
			scoreList = currentScoreList[:]
			probaList = currentProbaList[:]
		
		mappingScoreProba = dict(itertools.izip(scoreList,probaList))
		
		maxScore = max(scoreList)
		minScore = min(scoreList)

		Wmin = 10**10
		Wmax = -Wmin
		
		for elem in matrix :
			if max(elem) > Wmax :
				Wmax = max(elem)
			elif min(elem) < Wmin :
				Wmin = min(elem)
		
			
		Wmin = Wmin*tronc
		Wmin = math.floor(Wmin)
		Wmax = Wmax*tronc
		Wmax = math.floor(Wmax)	
			
		distribMin = min(minScore, Wmin)
		distribMax = max(maxScore, Wmax)
		
		breakMin = int(distribMin) -1
		breakMax = int(distribMax) +1
		
		scoreSorted = []
		scoreSortedInv = []
		
		for breakScore in range(breakMin, breakMax) :
			scoreSorted.append(breakScore)
			scoreSortedInv.insert(0,breakScore)
		
		#Not very usefull may be deleted
		scoreProbaCum = []
		tmpVal = 0
		for score in scoreSorted :
			try :
				proba = mappingScoreProba[score]
				tmpVal += proba
				scoreProbaCum.append(tmpVal)
				
			except KeyError :
				pass
		#end unusefull
		
		scoreProbaCumInv = []
		tmpScore = []
		tmpVal = 0
		for score in scoreSortedInv :
			try :
				proba = mappingScoreProba[score]
				tmpScore.append(score)
				tmpVal += proba
				scoreProbaCumInv.append(tmpVal)
			except KeyError :
				pass
		
		return dict(itertools.izip(tmpScore,scoreProbaCumInv))
		
	def sampleSequence2 (self, seq, m, back, iterations):
	
		if iterations < self.nbIterations-100 :
			indiceSeq = self.chooseSequence()
		else :
			indiceSeq = self.chooseSequencePropToWeight(back,seq, m)
	
		indice = 0
		probQ = 1
		probP = 1
		
		save = 0
		jumpSeq = 0
		
		listeIndice = []
		listeScore = []		
		
		indiceSeq = self.chooseSequence()
		
		size = seq.size[indiceSeq] - self.width + 1
		tmpSeq = seq.sequences[indiceSeq]
		
		subSeq = tmpSeq[indice:indice+self.width]
		
		while subSeq.find('n')>= 0 and indice < size:
			indice += 1
			subSeq = tmpSeq[indice:indice+self.width]
		
		#We supress the counts of the subsequence originally selected in our matrix
		self.subC(tmpSeq[m.indice[indiceSeq]:m.indice[indiceSeq]+self.width])
		#We calculate Q
		self.calculateQ()
			
		for i in range(0,self.width):
			probQ *= self.q[i][self.mapping[subSeq[i]]]
			probP *= self.p[self.mapping[subSeq[i]]]
			
		maxQP = probQ/probP

		listeIndice.append(indice)
		listeScore.append(maxQP)
		
		oldIndice = indice
		indice += 1
		
		while indice < size :
			probQ = 1
			subSeq = tmpSeq[indice:indice+self.width]
			
			if subSeq.find('n') < 0:
				for i in range(0,self.width):
					probQ *= self.q[i][self.mapping[subSeq[i]]]
					
				#If we jump some sequence letters because of 'n', we can't apply the unique modifications letter
				#to the probability P of the sub sequence
				if jumpSeq :
					probP = 1
					for i in range(0,self.width):
						probP *= self.p[self.mapping[subSeq[i]]]
					jumpSeq = 0
				else :
					probP = probP * (self.p[self.mapping[tmpSeq[indice+self.width-1]]]/self.p[self.mapping[tmpSeq[oldIndice]]])
				
				#tmpQP is the score of the subsequence corresponding to the matrix
				tmpQP = probQ/probP
				
				listeScore.append(tmpQP)
				listeIndice.append(indice)
				
				oldIndice = indice
			
			else :
			#	while indice < size and tmpSeq[indice:indice+self.width].find('n') >= 0 :
			#		indice = indice + self.width
				jumpSeq = 1
		
			indice += 1
		
		#if iterations == 0 :
		#	fileExport = open('logLikeListe.tab', 'w')
		#	fileExport.write('Indice\tLogLike\n')
		#	cpt = 0
		#	while cpt < len(listeScore) :
		#		fileExport.write('%s\t%s\n'%(listeIndice[cpt],listeScore[cpt]))
		#		cpt += 1
		#	fileExport.close()
		#we choose a random number proportional to a uniform function
		#cumulListe is a list of tuple. Each elem of the list is tuple of two elem : firstly the indice
		# and secondly the cumul score		
			
							
		mappingScoreIndice = dict(itertools.izip(listeScore,listeIndice))

		listeScore.sort()
		
		mappingScoreOrderedIndice = dict(itertools.izip(xrange(len(listeScore)),listeScore))
		
		listeScoreNorm = []
		listeScoreNormAmpl = []
		listeScoreAmplCum = []
		
		if len(listeScore) > 1 :
			maxVal = max(listeScore)
			minVal = min(listeScore)
			
			saveVal = 0
			for i in xrange(len(listeScore)) :
				normVal = (listeScore[i] - minVal)/(maxVal - minVal)
				listeScoreNorm.append(normVal)
				normValAmpl = normVal**self.exp
				listeScoreNormAmpl.append(normValAmpl)
				saveVal += normValAmpl
				listeScoreAmplCum.append(saveVal)
			
			maxVal = max(listeScoreAmplCum)
			minVal = min(listeScoreAmplCum)

		else :
			
			listeIseqNorm.append(1)
			listeIseqNormAmpl.append(1)
			listeIseqAmplCum.append(1)


		randNumber = random.uniform(0,listeScoreAmplCum[len(listeScore)-1])
		#val = random.uniform(0,save)
		
		#We fing the best approximation of val with all the scores that we obtained previously. Our choice
		#is then proportional to the scores value (big score have more chance to be pick than little one)
		#The value returned by findBestApprox is a tuple identical to the tuple in cumulListe	
		
		bestApproxVal = self.findBestApprox(randNumber, listeScoreAmplCum)
		indice = listeScoreAmplCum.index(bestApproxVal)
		bestIndice = mappingScoreIndice[mappingScoreOrderedIndice[indice]]
		bestApproxVal = mappingScoreOrderedIndice[indice]
				
		
		#We update the informations for our new sequence
		m.indice[indiceSeq] = bestIndice
		
		#We update the count matrix with our new sub sequence
		self.addC(tmpSeq[bestIndice:bestIndice+self.width])
		
		if self.verbose == 3 :
			print "Pos \t %s \t Seq \t %s \t Q \t %s \t P \t %s \t A \t %s \t InfoContent \t %s"%(bestIndice,subSeq,probQ, probP,probQ/probP, self.calculateInformationContent())
		
		if self.verbose == 4 :
			self.printConsensus()
			
		return bestApproxVal
		
	#Given a value, we must find the nearest cumul value to this value	
	def findBestApprox (self, val, listeIseq):	
		end = 0
		borneInf = 0
		borneSup = len(listeIseq)
		indice = (borneSup-borneInf)/2
		bestApprox = listeIseq[indice]
		
		if val == listeIseq[0] or val < listeIseq[0] :
			bestApprox = listeIseq[0]
			end = 1
			#if val < listeIseq[0] :
			#	print 'inf%s'%(val)
		
		while bestApprox != val and not end :
			if bestApprox > val :
				if indice > 0 and listeIseq[indice-1] >= val :
					if (bestApprox-val) > (listeIseq[indice-1]-val) :
						borneSup = indice
						indice -= (borneSup - borneInf)/2
					else :
						end = 1
				else :
					end = 1
					
			elif indice < len(listeIseq)-1 :
				borneInf = indice
				indice += (borneSup - borneInf)/2
			else :
				end = 1 
			
			bestApprox = listeIseq[indice]

		return bestApprox
	
	def calculeMatrixInfoContent(self,seq, ak, back):
		tmpCount =[[0]*len(self.alphabet) for i in xrange(self.width)]
		
		numSeq = 0
		for indice in ak:
			tmpSeq = seq.sequences[numSeq]
			for i in xrange(self.width):
				tmpCount[i][self.mapping[tmpSeq[indice+i]]] += 1
			numSeq +=1
			
		#print 'Nombre de sequences : ',self.height
		#for i in xrange(self.width) :
		#	print tmpCount[i]
		infoContent = 0
		
		
		for j in xrange(self.width) :
			for i in self.alphabet :
				freq= float(tmpCount[j][self.mapping[i]]+self.b[self.mapping[i]])/(self.height+self.B)
				if freq > 0 :
					infoContent += freq * math.log(freq/back.p[i])
					
		#print 'Background : ',back.p
		#print 'infoContent : %s'%(infoContent)	
		#print 'infoContent / nbSeq : %s'%(infoContent/self.height)	
		#print '---------------------------------'		
		return infoContent
		
	def calculMatrixLogLike(self,seq, ak, back):
		tmpCount =[[0]*len(self.alphabet) for i in xrange(self.width)]
		
		numSeq = 0
		for indice in ak:
			tmpSeq = seq.sequences[numSeq]
			for i in xrange(self.width):
				tmpCount[i][self.mapping[tmpSeq[indice+i]]] += 1
			numSeq +=1
			
		logLike = 0
		for j in xrange(self.width) :
			for i in self.alphabet :
				freq2 = float(tmpCount[j][self.mapping[i]]+self.b[self.mapping[i]])
				freq= float(freq2)/(self.height+self.B)
				logLike += freq2 * math.log(freq/back.p[i])
						
		return logLike
	
	def phaseShift(self, seq, m, back):
		ak = m.indice[:]
		
		maxShift = 3
		
		bestAk = m.indice[:]
		bestInfoContent = self.calculeMatrixInfoContent(seq, m.indice, back)	
		
		for shift in range(-maxShift,maxShift+1) :
			if shift != 0 :
				cpt = 0
				for indice in m.indice :
					if shift < 0:
						if indice + shift < 0 :
							if indice != 0 and seq.sequences[cpt][0:indice].find('n') < 0 :
								ak[cpt] = 0
							else :
								ak[cpt] = m.indice[cpt]
						else :
							if seq.sequences[cpt][indice+shift:indice].find('n') < 0 :
								ak[cpt] = indice+shift
							else :
								ak[cpt] = m.indice[cpt]
					else :
						if indice+shift >= seq.size[cpt]-self.width :
							if seq.sequences[cpt][seq.size[cpt]-self.width : seq.size[cpt]].find('n') < 0  :
								ak[cpt] = seq.size[cpt]-self.width
							else :
								ak[cpt] = m.indice[cpt]
						else :
							if seq.sequences[cpt][indice+shift:indice+shift+self.width].find('n') < 0 :
								ak[cpt] = indice+shift
							else :
								ak[cpt] = m.indice[cpt]
					cpt += 1		
					#if shift < 0 and indice + shift < 0 :
					#	if seq.sequences[cpt][0:indice].find('n') < 0 :
					#		ak[cpt] = 0
					#	else :
					#		ak[cpt] = m.indice[cpt]
					#elif shift > 0 and indice + shift >= seq.size[cpt]-self.width :
					#	if seq.sequences[cpt][indice:seq.size[cpt]-self.width].find('n') < 0  :
					#		ak[cpt] = seq.size[cpt]-self.width-1
					#	else :
					#		ak[cpt] = m.indice[cpt]
					#elif shift <0 and seq.sequences[cpt][indice+shift:indice].find('n') < 0  :
					#	ak[cpt] = indice+shift
					#elif seq.sequences[cpt][indice+shift:indice].find('n') <= 0 :
					#	ak[cpt] = indice+shift
					#else :
					#	ak[cpt] = m.indice[cpt]
					
				infoContent = self.calculeMatrixInfoContent(seq, ak, back)
				
				if bestInfoContent < infoContent :
					bestInfoContent = infoContent
					bestAk = ak[:]
				
		m.indice = bestAk
		self.initC (seq,m)
		
		return bestInfoContent		
	
	def finalCycle (self, seq, back, m, fileResult):
		output = 0
		if str(type(fileResult)).count('file') :
			output = 1
			
		weightMat = self.calculWeightMat(self.c, back)
		
		#Get the P-Values of the sites
		decimals = 2
		mappingScoreProba = self.calculTheorScoreDistrib (decimals, weightMat, back.p)
		
		totalNumberOfPosition = 0
		#Calcul the number of possible positions
		for size in seq.size :
			totalNumberOfPosition += size-self.width+1
		
		infoContentFinal = 0
		
		listeSeq = seq.sequences [:]
		listeIndiceSeq = range(0, len(seq.sequences))
		oldVal = 0
		epsilon = 1000
		iteration = 0
		breakVal = 0.5
		finalMatrix = []
		while iteration < 10 and epsilon > 0.3 :
			finalMatrix = []
			
			maxSeq = self.height*2
			for sequence in listeSeq :
				size = len(sequence)-self.width+1
				indice = 0
				maxEval = 10**10
				while indice < size :
					tmpScore = 0
					if sequence[indice : indice + self.width].find('n') < 0 :
						cpt = 0
						for letter in sequence[indice : indice + self.width] :
							tmpScore += math.floor(weightMat[cpt][self.mapping[letter]]*(10**decimals))
							cpt += 1
							
						reussi = 1
						save = tmpScore
						pVal = 0
						eVal = 0
						while reussi :
							try :
								pVal = mappingScoreProba[tmpScore]
								eVal = pVal * totalNumberOfPosition
								reussi = 0
							except KeyError :
								if reussi == 1 :
									if tmpScore > max(mappingScoreProba) :
										reussi = 2
										tmpScore = save
									else :
										tmpScore += 1
								elif reussi == 2 :
									if tmpScore < min(mappingScoreProba) :
										reussi = 0
									else :
										tmpScore -= 1
						
						
						if maxEval > eVal :
							numSeq = listeSeq.index(sequence)
							numSeq = listeIndiceSeq[numSeq]
							finalMatrix.append((eVal, sequence[indice:indice+self.width], numSeq, indice, indice+self.width-1, pVal, tmpScore))
						
						if len(finalMatrix) > maxSeq :
							maxEval = max (finalMatrix)
							elem = finalMatrix[finalMatrix.index(maxEval)]
							finalMatrix.remove(elem)
							maxEval = max (finalMatrix)	
						
					indice += 1
			
			cpt = 0
			while cpt < len(finalMatrix) :
				elem = finalMatrix[cpt]
				if elem[0] > breakVal :
					finalMatrix.remove(elem)
				else :
					cpt += 1
				
			#for elem in finalMatrix :
			#	print elem[0:7]
					
					
			if len(finalMatrix) > 0 :
				tmpCount =[[0]*len(self.alphabet) for i in xrange(self.width)]
				
				for i in xrange(self.width) :
					for j in xrange(len(self.alphabet)) :
						tmpCount[i][j] = self.b[j]
				
				for elem in finalMatrix :
					sequenceTmp = elem[1]
					cpt = 0
					for letter in sequenceTmp :
						tmpCount[cpt][self.mapping[letter]] += 1
						cpt += 1
			
				infoContent = 0
				for j in xrange(self.width) :
					for i in self.alphabet :
						freq= float(tmpCount[j][self.mapping[i]]+self.b[self.mapping[i]])/(len(finalMatrix)+self.B)
						if freq > 0 :
							infoContent += freq * math.log(freq/back.p[i])
						else :
							print 'Bizarre'
				
					
				#print '----------------'
				#print ' Information Content Value: ',infoContent
				
				infoContentFinal = infoContent
								
				weightMat = [[0]*len(self.alphabet) for i in xrange(self.width)]						
						
				for i in xrange(self.width) :
					for j in xrange(len(self.alphabet)) :
						Q = (tmpCount[i][j]+self.b[j])/(self.height+self.B)
						P = back.p[self.inverseMap[j]]
						weightMat[i][j] = Q/P
				
				
				listeIndiceSeqTmp = []
				for elem in finalMatrix :
					listeIndiceSeqTmp.append(elem[2])
				
				for elem in listeIndiceSeq :
					if listeIndiceSeqTmp.count(elem) == 0 :
						try :
							indiceSeq = listeIndiceSeq.index(elem)
							listeIndiceSeq.remove(indiceSeq)
							tmpSeq = seq.sequences[indiceSeq]
							listeSeq.remove(tmpSeq)
						except ValueError :
							pass
				
				epsilon = abs(infoContent - oldVal)
				oldVal = infoContent
				iteration += 1	
				breakVal = max(finalMatrix)[0]
				
				mappingScoreProba = self.calculTheorScoreDistrib (decimals, weightMat, back.p)
			else :
				epsilon = 0
				
			consensus = ''
			cpt = 0
			for elem in weightMat :
				indiceMax = weightMat[cpt].index(max(elem))
				consensus += self.inverseMap[indiceMax]
				cpt += 1
				
			#print '----------------'
			#print 'Motif Final :', consensus
		
		if output :
			datum = time.localtime()[0:6]
			AC = ''
			for elem in datum :
				AC += str(elem)
			fileResult.write('AC\t%s\nXX\n'%(AC))
			fileResult.write('TY\tfinal\nXX\n')
			
			nbZeros = 1
			rest = self.width+1
			while rest > 1 :
				rest = rest/10
				nbZeros += 1
				
			
			tmpCount =[[0]*len(self.alphabet) for i in xrange(self.width)]
			for i in xrange(self.width) :
					for j in xrange(len(self.alphabet)) :
						tmpCount[i][j] = self.b[j]
				
			for elem in finalMatrix :
				sequenceTmp = elem[1]
				cpt = 0
				for letter in sequenceTmp :
					tmpCount[cpt][self.mapping[letter]] += 1
					cpt += 1
			fileResult.write("P0\tA\tC\tG\tT\n")	
			if len(finalMatrix) > 0 :
				# pI =  "%0*.f" % (nbZeros,2) : writing the correct id for P for the correct number of zeros
				for i in xrange(self.width) :
					counts = tmpCount[i][:]
					for j in xrange(len(counts)) :
						counts[j] = int(math.floor(counts[j]))
						letter = self.inverseMap[counts.index(max(counts))]
					fileResult.write('%s\t%s\t%s\t%s\t%s\t%s\n'%("%0*.f" % (nbZeros,i+1), counts[0], counts[1], counts[2], counts[3], letter.upper()))
				fileResult.write('XX\n')

				#finalMatrix.append((eVal, sequence[indice:indice+self.width], numSeq, indice, indice+self.width-1, pVal, tmpScore))
				
				weightMat = [[0]*len(self.alphabet) for i in xrange(self.width)]						
						
				for i in xrange(self.width) :
					for j in xrange(len(self.alphabet)) :
						Q = (tmpCount[i][j]+self.b[j])/(self.height+self.B)
						P = back.p[self.inverseMap[j]]
						weightMat[i][j] = Q/P
				
				cpt = 0
				for elem in finalMatrix :
					tmpSeq = elem[1]
					seqId = seq.seqName[elem[2]]
					start = elem[3]
					end = start+self.width-1
					
					weight = 0
					cpt2 = 0
					for letter in tmpSeq :
						weight += weightMat[cpt2][self.mapping[letter]]
						cpt2 += 1
						 
					fileResult.write('BS\t%s; %s; %s; %s; strand D; %s\n'%(tmpSeq,seqId,start,end,weight))
					cpt += 1
				fileResult.write('XX\n')
				
					
				IC = 0
				LL = 0
				for j in xrange(self.width) :
					for i in self.alphabet :
						freq2 = float(tmpCount[j][self.mapping[i]]+self.b[self.mapping[i]])
						freq = freq2/(len(finalMatrix)+self.B)
						if freq > 0 :
							IC += freq * math.log(freq/back.p[i])
							LL += freq2 * math.log(freq/back.p[i])
				
				ICC = IC/self.width
				LLC = LL/self.width
				fileResult.write('IC\t%s\nXX\n'%(IC))
				fileResult.write('CC\t%s\nXX\n'%(ICC))
				fileResult.write('LL\t%s\nXX\n'%(LL))
				fileResult.write('CL\t%s\nXX\n'%(LLC))
				fileResult.write('//')
			else :
				for i in xrange(self.width) :
					fileResult.write('%s\t%s\t%s\t%s\t%s\t%s\n'%("%0*.f" % (nbZeros,i), 0, 0, 0, 0, '-'))
				fileResult.write('XX\n')
				fileResult.write('IC\t%s\nXX\n'%(0))
				fileResult.write('CC\t%s\nXX\n'%(0))
				fileResult.write('LL\t%s\nXX\n'%(0))
				fileResult.write('CL\t%s\nXX\n'%(0))
				fileResult.write('//')
			

		returnedList = str(infoContentFinal)+'\t'+consensus
		
		print 'Consensus du motif final : ',consensus
		print 'Contenu informationel du motif trouve : ',infoContentFinal

		return returnedList
	#We calculate the information content like in Hertz 1999
	#info = \sum_{i=1}^{L} \sum_{j = 1}{A} f_{ij} ln(\frac{f_{ij}}{p_j})	
	def calculateInformationContent (self):
		infoContent = 0
		for i in range(0,self.width):
			for j in range(0, len(self.alphabet)):
				ratio = (self.c[i][j]+self.b[j])/float(self.height+self.B)
				infoContent += ratio*math.log(ratio/(self.p[j]))
		return infoContent
	
	#We print de strict consensus
	def printConsensus (self) :
		#weightMat is the weight matrix for the motif
		weightMat = [[0]*len(self.alphabet) for i in xrange(self.width)]
		#maxi will be the consensus
		consensus = ['']*self.width
		
		for i in xrange(self.width) :
			for j in xrange(len(self.alphabet)) :
				#Reflechir au cas ou c_{ij} = 0
				#if self.c[i][j] == 0 :
				#	weightMat[i][j] = math.log(0.000000001 / (self.p[j]*self.height))
				if self.c[i][j] != 0 :
					#the weight is defied by : ln(\frac{f_{ij}}{p_i}) 
					weightMat[i][j] = math.log(float(self.c[i][j]) / (self.p[j]*self.height))
					
		for i in xrange(self.width) :
			consensus[i] = self.inverseMap[weightMat[i].index(max(weightMat[i]))] 
			
		print "consensus :", consensus
	
	def printMotif (self, m):
		weightMat = []
		
		for j in range (0, self.width):
			tmp = [0] * len(self.alphabet)
			weightMat.append(tmp)
			
		for i in range (0, self.width):
			for j in range (0, len(self.alphabet)):
				if self.c[i][j] != 0 :
					#weightMat[i][j] = math.log(self.c[i][j]/(self.p[j]*self.height))
					weightMat[i][j] = math.log(float(self.c[i][j])/(float(self.p[j])*self.height))
					#print weightMat[i][j]
				else :
					#print "\t \t",0
					weightMat[i][j] = 0		
			print weightMat[i]
		
		resultList = []	
		for i in range (0, self.width):	
			maxi = max(weightMat[i])
			resultList.append(self.inverseMap[weightMat[i].index(maxi)])
			
		print resultList
	
	def printToFile (self, seq, m):
		inputFile = seq.inputFile
		inputFile = inputFile.split(".fasta")
		
		outputFile = inputFile[0]+"_93.res"
		f = open(outputFile,'w')
		
		f.write("BF \t %s \n"%(inputFile[0]))
		f.write("XX \n")
		f.write("P0 \t A \t C \t G \t T \n")
		for i in xrange(self.width) :
			nbLetter = self.getLettersSubSequence(i, seq, m)
			mostRepresentedLetter = self.inverseMap[nbLetter.index(max(nbLetter))]
			f.write("%s \t %s \t %s \t %s \t %s \t %s \n"%(i,nbLetter[0],nbLetter[1],nbLetter[2],nbLetter[3], mostRepresentedLetter.upper()))
		f.write("XX \n")
		
		for i in xrange(self.height) :
			subSeq = seq.sequences[i][m.indice[i] : m.indice[i]+self.width]
			f.write("BS %s; %s \n"%(subSeq,seq.seqName[i]+"_"+str(m.indice[i])+"_"+str(m.indice[i]+self.width)+"_"+"D"))
		
		#f = open(self.freqFile,'r')
		f.close()
		
	def getLettersSubSequence(self, indice, seq, m) :
		res = [0] * len(self.alphabet)
		
		for i in xrange(self.height) :
			res[self.mapping[seq.sequences[i][m.indice[i] + indice]]] += 1
			
		return res
	
	def gibbs (self, seq, m, back, fileResult):
		nbIter = self.nbIterations
		#resultList = []
		#resultInform = []
		res = 0
		maxi = 0
		mini = 100000
		cpt = 0
		for i in range(0,nbIter):
			#print (self.sampleSequence (seq, m, back))
			#res = self.sampleSequence2 (seq, m, back)
			#res = self.sampleSequence (seq, m, back)
			#res = self.sampleSequenceLog (seq, m, back)
			
			res = self.sampleSequence2 (seq, m, back, i)
			
			if cpt%100 == 0 : #and cpt < nbIter/2 :
				res = self.phaseShift(seq,m,back)
			
			cpt += 1
			#if res > maxi :
			#	maxi = res
			#elif res < mini :
			#	mini = res
			#resultList.append((i,res))
			#resultInform.append((i,self.calculateInformationContent()))
		
		for elem in range(0,self.height) :
			print seq.sequences[elem][m.indice[elem] : m.indice[elem]+self.width]
		
		print ("Le motif trouve :")	
		self.printMotif (m)
		
		if str(type(fileResult)).count('file') :
				self.printToResFile(seq, m, back, fileResult)
	
		print 'Cycle Final :'	
		self.finalCycle (seq, back, m, fileResult)
		
		
		#for elem in range(0,self.height) :
		#	print seq.sequences[elem][m.indice[elem] : m.indice[elem]+self.width]
		#print ("----------------\nLe motif trouve :")	
		#self.printMotif (m)
		
		#self.printToFile(seq, m)
		
		#print 'Contenu informationel : %s'%(self.calculateInformationContent())
		

	def printToResFile(self, seq, m, back, fileResult):
		datum = time.localtime()[0:6]
		AC = ''
		for elem in datum :
			AC += str(elem)
		fileResult.write('AC\t%s\nXX\n'%(AC))
		fileResult.write('TY\tend of sampling\nXX\n')
		
		nbZeros = 1
		rest = self.width+1
		while rest > 1 :
			rest = rest/10
			nbZeros += 1
		
		fileResult.write("P0\tA\tC\tG\tT\n")
		# pI =  "%0*.f" % (nbZeros,2) : writing the correct id for P for the correct number of zeros
		for i in xrange(self.width) :
			counts = self.c[i][:]
			for j in xrange(len(counts)) :
				counts[j] = int(math.floor(counts[j]))

			letter = self.inverseMap[counts.index(max(counts))]
			fileResult.write('%s\t%s\t%s\t%s\t%s\t%s\n'%("%0*.f" % (nbZeros,i+1), counts[0], counts[1], counts[2], counts[3], letter.upper()))
		fileResult.write('XX\n')
		weightMat = self.calculWeightMat(self.c, back)
		
		cpt = 0
		for sequence in seq.sequences :
			tmpSeq = sequence[m.indice[cpt] : m.indice[cpt]+self.width]
			seqId = seq.seqName[cpt]
			start = m.indice[cpt]
			end = start+self.width-1
			
			weight = 0
			cpt2 = 0
			for letter in tmpSeq :
				weight += weightMat[cpt2][self.mapping[letter]]
				cpt2 += 1
				 
			fileResult.write('BS\t%s; %s; %s; %s; strand D; %s\n'%(tmpSeq,seqId,start,end,weight))
			cpt += 1
		fileResult.write('XX\n')
		IC = self.calculeMatrixInfoContent(seq, m.indice, back)
		ICC = IC/self.width
		LL = self.calculMatrixLogLike(seq, m.indice, back)
		LLC = LL/self.width
		fileResult.write('IC\t%s\nXX\n'%(IC))
		fileResult.write('CC\t%s\nXX\n'%(ICC))
		fileResult.write('LL\t%s\nXX\n'%(LL))
		fileResult.write('CL\t%s\nXX\n'%(LLC))
		fileResult.write('//\n')
