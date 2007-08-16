import itertools
import time
import math
import random
import sys

class GibbsInform :
	def __init__(self,width,height, alphabet, seq, m,back, verbose, nbIterations, exp):
		self.width = width
		self.height = height
		self.k = 1
		
		self.verbose = verbose
		self.nbIterations = nbIterations
		self.exp = exp
		
		self.alphabet = alphabet
		self.mapping = dict()
		self.inverseMap = dict()
		
		self.c = [[0]*len(self.alphabet) for i in xrange(self.width)]
		self.informCol = [[0]*len(self.alphabet) for i in xrange(self.width)]
			
		self.createMapping()
		self.initCount(seq,m, back)
			
	def createMapping (self):
		indice = []
		for i in range(0,len(self.alphabet)) :
			indice.append(i)
		self.mapping = dict(itertools.izip(self.alphabet,indice))
		self.inverseMap = dict(itertools.izip(indice,self.alphabet))
	
	def initCount (self,seq,m, back):
		numSeq = 0
		self.c = [[0]*len(self.alphabet) for i in xrange(self.width)]
		
		for i in xrange(self.width) :
			for j in xrange(len(self.alphabet)) :
				self.c[i][j] = self.k * back.p[self.inverseMap[j]]
		
		for indice in m.indice:
			tmpSeq = seq.sequences[numSeq]
			for i in range(0,self.width):
				self.c[i][self.mapping[tmpSeq[indice+i]]] += 1
			numSeq +=1
	
	def subSeqCount (self, subSeq) :
		for i in xrange(self.width):
			self.c[i][self.mapping[subSeq[i]]] -= 1.0
			
	def recalculateCount (self, subSeq):
		for i in xrange(self.width):
			self.c[i][self.mapping[subSeq[i]]] += 1.0
	
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
	
	def calculWeightMat (self, count, back) :
		weightMat = [[0]*len(self.alphabet) for i in xrange(self.width)]
		
		for i in xrange(self.width) :
			for j in xrange(len(self.alphabet)) :
				if count[i][j] != 0 :
					weightMat[i][j] = math.log((float(count[i][j])/(self.height+self.k))/back.p[self.inverseMap[j]])
		return weightMat
	
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
	
	#Special function which calculate the informContent of a column of the matrix but whithout one sequence. There is
	#4 choices of inforcontent per column if the alphabet is A,C,G,T.
	def calculInformCol (self, back) :
		for col in xrange(self.width) :
			for letter in xrange(len(self.alphabet)) :
				self.c[col][letter] += 1
				self.informCol[col][letter] = 0
				
				for j in self.alphabet :
					freq = self.c[col][self.mapping[j]]/float(self.height+self.k)
					#freq2 = self.c[col][self.mapping[j]]
					if freq > 0 :
						self.informCol[col][letter] += freq * math.log(freq/back.p[j])
					else :
						print 'Bizarre'
						
				self.c[col][letter] -= 1
		
	def sampleSequenceFast (self, seq, m, back, iterations):
	#Add log likelyhood score == infocontent * -N (cf article hertz) << verif, if true, only one class :D [Urgent]
		indice = 0
		
		if iterations < self.nbIterations-100 :
			indiceSeq = self.chooseSequence()
		else :
			indiceSeq = self.chooseSequencePropToWeight(back,seq, m)

		bestIndice = m.indice[indiceSeq]
		
		size = seq.size[indiceSeq]
		size = size - self.width + 1
		tmpSeq = seq.sequences[indiceSeq]
		
		self.subSeqCount(tmpSeq[bestIndice:bestIndice+self.width])
		self.calculInformCol (back)
		
		listeIndice = []
		listeIseq = []
		
		maxInform = 0 
		
		while indice < size :
			subSeq = tmpSeq[indice:indice+self.width]
			if subSeq.find('n') < 0 :
				tmpIseq = 0
				cpt = 0
				for letter in subSeq :
					tmpIseq += self.informCol[cpt][self.mapping[letter]]
					cpt += 1
					
				listeIndice.append(indice)
				listeIseq.append(tmpIseq)
				
			indice += 1
		
		#if iterations == 0 :
		#	fileExport = open('logLike_0.tab', 'w')
		#	fileExport.write('Indice\tContenuInformationel\n')
		#	cpt = 0
		#	while cpt < len(listeIseq) :
		#		fileExport.write('%s\t%s\n'%(listeIndice[cpt],listeIseq[cpt]))
		#		cpt += 1
		#	fileExport.close()
		
		#if iterations == self.nbIterations - 1 :
		#	fileExport = open('logLike_3000.tab', 'w')
		#	fileExport.write('Indice\tContenuInformationel\n')
		#	cpt = 0
		#	while cpt < len(listeIseq) :
		#		fileExport.write('%s\t%s\n'%(listeIndice[cpt],listeIseq[cpt]))
		#		cpt += 1
		#	fileExport.close()
		
		mappingIseqIndice = dict(itertools.izip(listeIseq,listeIndice))

		listeIseq.sort()
		
		mappingIseqOrderedIndice = dict(itertools.izip(xrange(len(listeIseq)),listeIseq))
		
		listeIseqNorm = []
		listeIseqNormAmpl = []
		listeIseqAmplCum = []
		
		if len(listeIseq) > 1 :
			maxVal = max(listeIseq)
			minVal = min(listeIseq)
			
			saveVal = 0
			for i in xrange(len(listeIseq)) :
				normVal = (listeIseq[i] - minVal)/(maxVal - minVal)
				listeIseqNorm.append(normVal)
				normValAmpl = normVal**self.exp
				listeIseqNormAmpl.append(normValAmpl)
				saveVal += normValAmpl
				listeIseqAmplCum.append(saveVal)
			
			maxVal = max(listeIseqAmplCum)
			minVal = min(listeIseqAmplCum)
			listeIseqAmplCumNorm = []
			for i in xrange(len(listeIseq)) :
				listeIseqAmplCumNorm.append((listeIseqAmplCum[i] - minVal)/(maxVal - minVal))
		else :
			
			listeIseqNorm.append(1)
			listeIseqNormAmpl.append(1)
			listeIseqAmplCum.append(1)
		
		#coeff = 1.0/listeIseq[len(listeIseq)-1]
		#save = 0

		#for i in xrange(len(listeIseq)) :
		#	listeIseq[i] = (coeff*listeIseq[i])**self.exp
			#val += listeIseq[i]
		#	sumCumulIseq.append((i, listeIseq[i]))
		
		
		randNumber = random.uniform(0,listeIseqAmplCum[len(listeIseq)-1])
		
		
		#maxVal = listeIseq[len(listeIseq)-1]
		#Old Trick to choose better segments
		#maxVal = listeIseqNormAmpl[len(listeIseq)-1]
		#mu = maxVal*0.8
		#sigma = maxVal*0.2
		#randNumber = random.gauss(mu,sigma)
		#if randNumber > maxVal :
		#	randNumber = maxVal
		
		bestApprox = self.findBestApprox(randNumber, listeIseqAmplCum)

		
		indice = listeIseqAmplCum.index(bestApprox)
		bestIndice = mappingIseqIndice[mappingIseqOrderedIndice[indice]]
		chooseInform = mappingIseqOrderedIndice[indice]
				
		self.recalculateCount(tmpSeq[bestIndice:bestIndice+self.width])
		m.indice[indiceSeq] = bestIndice
			
		if self.verbose == 3 :
			print "Pos \t %s \t Seq \t %s \t InfoContentMat \t %s"%(indice,subSeq,tmpIseq)
		
		#affiche le concensus de chaque iter	
		if self.verbose == 4 :
			self.printConsensus(back)
			
		return chooseInform		
	
	def findBestApprox (self, val, listeIseq):	
		end = 0
		borneInf = 0
		borneSup = len(listeIseq)
		indice = (borneSup-borneInf)/2
		bestApprox = listeIseq[indice]
		
		if val == listeIseq[0] or val < listeIseq[0] :
			bestApprox = listeIseq[0]
			end = 1
		
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
		
		for i in xrange(self.width) :
			for j in xrange(len(self.alphabet)) :
				tmpCount[i][j] = self.k * back.p[self.inverseMap[j]]
		
		numSeq = 0
		for indice in ak:
			tmpSeq = seq.sequences[numSeq]
			for i in xrange(self.width):
				tmpCount[i][self.mapping[tmpSeq[indice+i]]] += 1
			numSeq +=1
			
		infoContent = 0
		for j in xrange(self.width) :
			for i in self.alphabet :
				freq= float(tmpCount[j][self.mapping[i]])/(self.height+self.k)
				if freq > 0 :
					infoContent += freq * math.log(freq/back.p[i])
				else :
					print 'Bizarre'
						
		return infoContent
		
	def calculMatrixLogLike(self,seq, ak, back):
		tmpCount =[[0]*len(self.alphabet) for i in xrange(self.width)]
		
		for i in xrange(self.width) :
			for j in xrange(len(self.alphabet)) :
				tmpCount[i][j] = self.k * back.p[self.inverseMap[j]]
		
		numSeq = 0
		for indice in ak:
			tmpSeq = seq.sequences[numSeq]
			for i in xrange(self.width):
				tmpCount[i][self.mapping[tmpSeq[indice+i]]] += 1
			numSeq +=1
			
		logLike = 0
		for j in xrange(self.width) :
			for i in self.alphabet :
				freq2 = tmpCount[j][self.mapping[i]]
				freq= float(freq2)/(self.height+self.k)
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
					
				infoContent = self.calculeMatrixInfoContent(seq, ak, back)
				
				if bestInfoContent < infoContent :
					bestInfoContent = infoContent
					bestAk = ak[:]
				
		m.indice = bestAk
		self.initCount (seq,m, back)
		
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
						tmpCount[i][j] = self.k * back.p[self.inverseMap[j]]
				
				for elem in finalMatrix :
					sequenceTmp = elem[1]
					cpt = 0
					for letter in sequenceTmp :
						tmpCount[cpt][self.mapping[letter]] += 1
						cpt += 1
			
				infoContent = 0
				for j in xrange(self.width) :
					for i in self.alphabet :
						freq= float(tmpCount[j][self.mapping[i]])/(len(finalMatrix)+self.k)
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
						if tmpCount[i][j] != 0 :
							weightMat[i][j] = math.log((float(tmpCount[i][j])/len(finalMatrix))/back.p[self.inverseMap[j]])
				
				
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
					tmpCount[i][j] = self.k * back.p[self.inverseMap[j]]
				
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
						if tmpCount[i][j] != 0 :
							weightMat[i][j] = math.log((float(tmpCount[i][j])/len(finalMatrix))/back.p[self.inverseMap[j]])
				
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
						freq2 = float(tmpCount[j][self.mapping[i]])
						freq = freq2/(len(finalMatrix)+self.k)
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
		return returnedList
		
	def gibbs (self, seq, m, back, fileResult):
		
		nbIter = self.nbIterations
		resultList = []
		res = 0
		cpt = 0
		
		for i in range(0,nbIter):
			self.sampleSequenceFast (seq, m, back, cpt)
			
			res = self.calculeMatrixInfoContent(seq, m.indice, back)
				
			if cpt%100 == 0 : 
				self.phaseShift(seq,m,back)
				
			resultList.append((i,res))
			cpt += 1
			
					
		#if self.verbose == 1 :
		#	print >> sys.stdout, self.finalCycle(seq, back, m, resFile) 
			
		if self.verbose == 5 :
			for elem in range(0,self.height) :
				print seq.sequences[elem][m.indice[elem] : m.indice[elem]+self.width]
			
			print ("----------------\nLe motif trouve :")	
			self.printMotif (m, back)
			print ("----------------\nLa matrice de score du motif trouve:")
			self.printScore (back)
			print ("----------------\nInformation Content Value : %s")%(self.calculeMatrixInfoContent(seq, m.indice, back))
			self.printToFile(seq, m)
			
			if str(type(fileResult)).count('file') :
				self.printToResFile(seq, m, back, fileResult)
			
			print ("----------------\nFinal Cycle :")
			tmp = self.finalCycle(seq, back, m, fileResult)

#Usage Function From Class GibbsInform
	def printMotif (self, m, back):
		weightMat = [[0]*len(self.alphabet) for i in xrange(self.width)]
			
		for i in xrange(self.width):
			for j in xrange(len(self.alphabet)):
				if self.c[i][j] != 0 :
					weightMat[i][j] = math.log(float(self.c[i][j])/(back.p[self.inverseMap[j]]*(self.height+self.k)))
				else :
					weightMat[i][j] = 0		
		listeMaxi = []		
		for i in range (0, self.width):	
			maxi = max(weightMat[i])
			listeMaxi.append(self.inverseMap[weightMat[i].index(maxi)])
		print listeMaxi
		
	def printConsensus (self, back) :
		weightMat = [[0]*len(self.alphabet) for i in xrange(self.width)]
		maxi = ['']*self.width
		
		for i in xrange(self.width) :
			for j in xrange(len(self.alphabet)) :
				if self.c[i][j] != 0 :
					weightMat[i][j] = math.log10((float(self.c[i][j])/(self.height+self.k))/back.p[self.inverseMap[j]])
					
		for i in xrange(self.width) :
			maxi[i] = self.inverseMap[weightMat[i].index(max(weightMat[i]))] 
			
		print "consensus :", maxi	
	def printScore (self, back):
		weightMat = [[0]*len(self.alphabet) for i in xrange(self.width)]
		maxi = ['']*self.width
		
		for i in xrange(self.width) :
			for j in xrange(len(self.alphabet)) :
				if self.c[i][j] != 0 :
					weightMat[i][j] = math.log((float(self.c[i][j])/(self.height+self.k))/back.p[self.inverseMap[j]])
		
		for i in xrange(self.width) :
			print weightMat[i]
			
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
		
		
		
		
		
	def printToFile (self, seq, m):
		inputFile = seq.inputFile
		inputFile = inputFile.split(".fasta")
		
		outputFile = inputFile[0]+"_inform.res"
		f = open(outputFile,'w')
		f.write("// \nXX \n")
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
		f.write("//")
		#f = open(self.freqFile,'r')
		f.close()
		
	def getLettersSubSequence(self, indice, seq, m) :
		res = [0] * len(self.alphabet)
		
		for i in xrange(self.height) :
			res[self.mapping[seq.sequences[i][m.indice[i] + indice]]] += 1
			
		return res
