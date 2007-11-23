"""
   1. Choose initial population
   2. Evaluate the fitness of each individual in the population
   3. Repeat
         1. Select best-ranking individuals to reproduce
         2. Breed new generation through crossover and mutation (genetic operations) and give birth to offspring
         3. Evaluate the individual fitnesses of the offspring
         4. Replace worst ranked part of population with offspring
   4. Until <terminating condition>
"""

import random

from Stats.utils import mean, var


import Core.UI.cli as cli


def fitness(i1):
    raise NotImplementedError


def crossover(i1, i2):
    raise NotImplementedError

def create_initial_population(size):
    raise NotImplementedError


def mutate(i1):
    raise NotImplementedError


def roll(probs):
    """
    Pipped dice roll
    probs -- list of probalitity for each value (example [0.25,0.25,0.25,0.25]) value are in [0..len(l)[
    
    return: value (int)
    """
    a = 0.0
    bounds = []
    for i in range(len(probs)):
        bounds.append( [a, a+probs[i]])
        a += probs[i]

    r = random.random()
    i = 0
    for a,b in bounds:
        if r >= a and  r < b:
            return i
        i += 1
    raise Error


#def sort(population):
#    population.sort( lambda i1, i2 : cmp(fitness(i1), fitness(i2)) )
    
def tournament_selection(population, size=3):
    assert(len(population) - size - 1 > 0) 
    
    r = range(len(population))

    x1 = random.sample(r, 3)[-1]
    r.remove(x1)
    x2 = random.sample(r, 3)[-1]

    return population[x1], population[x2]


def ga(population, parameters):
    """
    parameters dict
    Pm : mutation prob
    Pc : crossover prob
    
    maxGeneration : number of generation
    tournamentSize : tournament size
    """

    Pm = parameters['Pm']
    Pc = parameters['Pc']
    maxGeneration = parameters['maxGeneration']
    tournamentSize = parameters['tournamentSize']
    populationSize = parameters['populationSize']
    elitism = parameters['elitism']

    #1. Choose initial population
    #population = create_initial_population(populationSize)

    info = cli.Info(maxGeneration+1, 1)
    #2. Evaluate the fitness of each individual in the population
   

    #3. Repeat
    generationCount = 0
    while True:

        population.sort()

        #
        # STATS
        #
        best = population[-1].fitness() 
        bestConsensus = population[-1].consensus()
        #worst = population[0].fitness()
        l = [i.fitness() for i in population]
        mu, sigma = mean(l), var(l)

        avgWords = mean([len(i.spl) for i in population])
        nw = len(i.spl)

        info('gen=%-4d best=%.3f %s %d mean=%.3f var=%.4f avgW=%.1f' % (generationCount, best, bestConsensus, nw, mu, sigma, avgWords))
        parameters['tab'].append({'gen' : generationCount,'IC' : best, 'word' : bestConsensus, 'nw' : nw, 'mean' : mu, 'var' : sigma, 'time' : parameters['clock'].get_value()}) 

        #
        # BREAK
        #
        if generationCount >= maxGeneration:
            break

        #
        # NEW GENERATION
        #
        generationCount += 1
        offspring = []

        for i in range(0, len(population), 2):
            #1. Select best-ranking individuals to reproduce
            
            cli.info('====tournament', tournamentSize)

            j1, j2 = tournament_selection(population, 3)

            cli.info('%s\n*%s' % (j1,j2), 5); cli.info('tournament======', 5)

            #2. Breed new generation through crossover and mutation (genetic operations) and give birth to offspring
            if random.random() <= Pc:
                cli.info('=====crossover', 4); cli.info('(input)%s\n*%s' % (j1,j2), 4)

                j1, j2 = crossover(j1, j2)

                cli.info('(output)%s\n*%s' % (j1,j2), 4); cli.info('crossover=====', 4)
                
            if random.random() <= Pm:

                cli.info('=====mutation', 4); cli.info('(input)%s\n*%s' % (j1,j2), 4)

                j1, j2 = mutation(j1), mutation(j2)

                cli.info('(outplut)%s\n*%s' % (j1,j2), 4); cli.info('mutation======', 4)
            
                
            if j1 not in population and j1 not in offspring:
                offspring += [j1]
            if j2 not in population and j2 not in offspring:
                offspring += [j2]


        #3. Evaluate the individual fitnesses of the offspring
        population = population + offspring
        population.sort()

        #4. Replace worst ranked part of population with offspring

        #hard selection
        #population = population[-populationSize:]
        
        #soft selection
        top = population[min(-1,-int(populationSize*elitism)):]
        random.shuffle(population)
        population = population[:populationSize-len(top)] + top


    return population
