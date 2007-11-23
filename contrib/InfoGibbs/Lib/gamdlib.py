
import random
from copy import copy

import ga

from matrix import *


class Individual(object):
    def __init__(self, sequences, bg, spl=None):

        self.sequences = sequences
        self.bg = bg # Background
        self.spl = spl or [] #sequence, position, length
        self._fitness = None
        self._matrix = None
        #self.priori = sites2bernoulli(sequences)

    def get_words(self):
        words = []
        for s,p,l in self.spl:
            words += [ self.sequences[s][p:p+l] ]
        return words

    #def __repr__(self):
    #    str = [ 'IC = %.3f %d' % (self.fitness(), len(self.spl)) ]
    #    words = self.get_words()
    #    matrix = words2matrix(words, self.bg)
    #    for s,p,l in self.spl:
    #        str += [ '%i,%i,%i' % (s,p,l) + ' %s' % self.sequences[s][p:p+l] ]
    #    return '\n'.join(str)

    def consensus(self):
        return consensus(self.matrix(), self.bg.priori)

    def fitness(self):
        if self._fitness == None:
            #self._fitness = llr(self.get_words(), self.matrix(), self.bg)
            self._fitness = Iseq_markov(self.get_words(), self.matrix(), self.bg)
        return self._fitness

    def matrix(self):
        return words2matrix(self.get_words(), self.bg.priori)        

    def update(self):
        self.spl.sort()
        self._fitness = None
        self.fitness()


    def __eq__(self, other):
        return self.spl == other.spl

    def __cmp__(self, other):
        return cmp(self.fitness(), other.fitness())


    def copy(self):
        return Individual(self.sequences, self.bg, copy(self.spl))        

    def str(self):
        title = [ 'IC=%.3f' % (self.fitness()) ]
        title += ['N=%d' % len(self.spl) ]
        title += ['Consensus=%s' % self.consensus() ]    
        return matrix2tab(self.matrix(), '\t'.join(title))


    def full_output(self):
        str = [ 'Information Content = %.3f' % (self.fitness()) ]
        str += ['Number of words = %d' % len(self.spl) ]
        str += ['Consensus = %s' % self.consensus() ]

        str += [ '' ]
        #str += [ 'ALIGNED SEQUENCES' ]
        words = self.get_words()
        matrix = self.matrix()
        for s,p,l in self.spl:
            str += [ '%i,%i,%i' % (s,p,l) + ' %s' % self.sequences[s][p:p+l] ]
        str += [ '' ]
        #str += [ 'ALIGNED SEQUENCES' ]
        str += [ '\n'.join(self.get_words()) ]
        str += [ '' ]
        #str += [ 'MATRIX' ]
        str += [ matrix2atb(self.matrix()) ]
        str += [ '' ]
        str += [ '' ]
        return '\n'.join(str)



def fitness(i):
    return i.fitness()


def shifted(i, x):
    newI = Individual(i.sequences, i.bg)
    for k in range(len(i.spl)):
        s,p,l = i.spl[k]
        if p+x <= len(i.sequences[s])-l and p+x >= 0:
            newI.spl.append( (s,p+x,l) )
    newI.update()
    return newI


###
#
# Crossover
#
###
def crossover(i1, i2):
    li1 = len(i1.spl)
    li2 = len(i2.spl)
    c = random.randint(1, min(li1, li2) - 1)
    
    p1 = random.sample(range(li1), c)
    p2 = random.sample(range(li2), c)

    newI1 = i1.copy()
    newI2 = i2.copy()

    for k in range(len(p1)):
        spl1 = i1.spl[p1[k]]
        spl2 = i2.spl[p2[k]]        

        #interchange if ok
        if not spl1 in i2.spl and not spl2 in i1.spl:
            newI1.spl[p1[k]] = spl2
            newI2.spl[p2[k]] = spl1

    newI1.update()
    newI2.update()

    return newI1, newI2


def crossover_full(i1, i2):
    li1 = len(i1.spl)
    li2 = len(i2.spl)
    c = random.randint(1, min(li1, li2) - 1 )
    
    newI1 = i1.copy() #Individual(i1.sequences, i1.bg)
    newI2 = i2.copy()

    if li1 <= 2 or li2 <=2:
        return i1, i2


    for k in range(min(li1,li2)):
        spl1 = i1.spl[k]
        spl2 = i2.spl[k]        

        alpha = random.random()

        #crossover
        if alpha < 0.8:
            if not spl1 in newI2.spl and not spl2 in newI1.spl:
                newI1.spl.remove(spl1)
                newI2.spl.remove(spl2)

                newI1.spl.append(spl2)
                newI2.spl.append(spl1)


        #i1->i2
        elif alpha < 0.9 and not spl1 in newI2.spl:
                newI1.spl.remove(spl1)
                newI2.spl.remove(spl2)

                newI2.spl.append(spl1)

        #i2->i1
        #elif alpha < 1.0 and not spl2 in newI1.spl:
        #    newI1.spl.append(spl2)


    newI1.update()
    newI2.update()
    return newI1, newI2


###
#
# Mutation
#
###




def mutation_soft(i1):
    # choose a new word
    spl = choose_word(i1, i1.spl[0][-1], i1.spl[0][-1])

    i = i1.copy()
    
    #replace an existing word with that
    if not spl in i.spl: 
        x = random.randint(0, len(i.spl)-1)
        i.spl[x] = spl
        i.update()
    return i


def mutation_hard(i1):
    i = i1.copy()
    #choose a number of words
    l = len(i.spl)
    for k in random.sample(range(l), random.randint(1, l)):

        # choose a new word
        spl = choose_word(i, i.spl[0][-1], i.spl[0][-1])

        #replace an existing word with that
        if not spl in i.spl: 
            i.spl[k] = spl
            i.update()
    return i


#def mutation_insert(i1):
#    # choose a new word
#    spl = choose_word(i1, i1.spl[0][-1], i1.spl[0][-1])
#    i = i1.copy()
#    if not spl in i.spl: 
#        i1.spl.append(spl)
#        i.update()
#    return i



def mutation(i1):
    alpha = random.random()

    if alpha < 0.8:
        return mutation_soft(i1)

    elif alpha < 0.9:
        return mutation_shift(i1)

    #elif alpha < 0.95:
    #    return mutation_indel(i1)

    elif alpha <= 1.0:
        return mutation_hard(i1)
    
    return i1



# mutation 2

def mutation_shift(i1):
    shift = random.randint(-2,2)
    newI = shifted(i1, shift)
    
    if len(i1.spl) != len(newI.spl):
        return i1
    else:
        return newI


def mutation_indel(i1):
    if len(i1.spl) <= 2:
        return i1

    # choose a new word
    spl = choose_word(i1, i1.spl[0][-1], i1.spl[0][-1])

    i = i1.copy()

    alpha = random.random()

    #insert
    if alpha < 0.01:
        if not spl in i.spl: 
            i1.spl.append(spl)
            i.update()

    #delete
    elif alpha < 1.0:
        x = random.randint(0, len(i.spl)-1)
        del i.spl[x]
        i.update()


    return i


def mutation_final(i):
    alpha = random.random()
    if alpha < 0.4:
        return mutation_shift(i)
    else:
        return mutation_soft(i)
    

###
#
# Initial population
#
###
def choose_word(i, MIN_WORD_LENGTH=4, MAX_WORD_LENGTH=4):
    s = random.randint(0, len(i.sequences)-1)
    #l = random.randint( MIN_WORD_LENGTH, min(MAX_WORD_LENGTH, len(i.sequences[s])) )
    l = MIN_WORD_LENGTH
    p = random.randint(0, len(i.sequences[s])-l)
    return s,p,l


def initial_population_boosted(sequences, bg, parameters):
    max = -10000
    for n in range(100):
        population = initial_population(sequences, bg, parameters)
        f = population[-1].fitness()
        if f > max:
            best = population
            max = f

    return best


def initial_population(sequences, bg, parameters):
    """

    """

    size = parameters['populationSize']
    wordMinLength = parameters['wordMinLength']
    wordMaxLength = parameters['wordMaxLength']
    wordsPerIndividual = parameters['wordsPerIndividual']


    population = []
    while (len(population) < size):

        #create individual
        i = Individual(sequences, bg)
        while(len(i.spl)) < wordsPerIndividual:
            w = choose_word(i, wordMinLength, wordMaxLength)
            if w not in i.spl:
                i.spl += [ w ]
        #

        if i not in population:
            i.update()
            population += [ i ]

    population.sort()

    return population


def initial_population_bias(sequences, bg, parameters):
    """

    """

    size = parameters['populationSize']
    wordMinLength = parameters['wordMinLength']
    wordMaxLength = parameters['wordMaxLength']
    wordsPerIndividual = parameters['wordsPerIndividual']

    l = wordMaxLength

    priori = bg

    f = []
    spl = []

    
    #build choice function
    for s in range(len(sequences)):
        for p in range(0, len(sequences[s])-l+1):
            word = sequences[s][p:p+l]
            Score = 1.0 - word_prob(word, priori)
            f += [ Score ]
            spl += [ (s,p,l) ]
    
    choose_word_bias = wchoice(spl, f)




    population = []
    while (len(population) < size):

        #create individual
        i = Individual(sequences, bg)

        while(len(i.spl)) < wordsPerIndividual:
            w = choose_word_bias()
            if w not in i.spl:
                i.spl += [ w ]
            
        #

        if i not in population:
            i.update()
            population += [ i ]

    population.sort()

    return population




def print_population(p, info='', printIndividuals=True):
    print '<<<<<<<< %s >>>>>>>>' % info
    f = [i.fitness() for i in p]
    f.sort()

    print 'Best = %.3f' % f[-1]

    #print individuals
    if not printIndividuals:
        return
    count = 0
    for i in p:
        count +=1
        print '[%d]' % count
        print i.full_output()
        print ''


###
#
# Tests
#
###

def run():
    sequences = ['ATTAAAATGGTT', 'GAAAAGCCGTGT', 'ATGAAAAGGCTGGATGC', 'TTTGAAAAT']
    begin = create_initial_population(6, sequences)
    ga.mutation = mutation2
    ga.crossover = crossover
    ga.fitness = fitness
    final = ga.ga(begin)
    #print_population(begin, 'starting')
    print_population(final, 'final')

    print 'Top 3'
    print final[-1]
    print final[-2]
    print final[-3]

if __name__ == '__main__':
    run()


    

