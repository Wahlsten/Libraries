import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import animation
import math

def GeneticAlgorithm(f, populationSize = 10,  numGenerations = 10, geneSize = 8):

    num_generations = numGenerations
    population_size = populationSize
    gene_size = geneSize
    func = f

def InitializeChromosomes(population_size, chromosome_size):
    ''' Initializing Chromosomes'''
    
    chromosome_list = []

    for k in range(population_size):

        chromosome_list.append(InitChromosome(chromosome_size))

    return chromosome_list

def InitChromosome( chromosome_size ):
    ''' Initializing chromosome '''
    
    chromosome = np.random.randint(0, 2, chromosome_size)

    return chromosome

def DecodeChromosome( chromosome ):
    ''' Decoding chromosome '''

    exp_upper = np.floor(len(chromosome) / 2.0)
    exp_lower = np.ceil(len(chromosome) / 2.0) - 1
    exponents = np.linspace(exp_lower, exp_upper)

    chromosome_decoded = 0.0

    for k in range(len(chromosome)):

        chromosome_decoded = chromosome_decoded + chromosome[k]

    return chromosome_decoded
    
def Crossover( chromosome1, chromosome2, cross_left, cross_right ):
    ''' Crossover '''
    
    crossover_new = np.concatenate((chromosome1[:cross_left], chromosome2[cross_left:cross_right], chromosome1[cross_right:]), axis=None)

    return crossover_new

def Mutate(chromosome, mutate_prob):
    ''' Mutate '''
    
    mutated_chromosome = np.zeros(len(chromosome))
    randomVec = np.random.rand(len(chromosome))
    
    for k in range(len(chromosome)):

        r = randomVec[k]

        if r < mutate_prob:

            if chromosome[k] == 1:

                mutated_chromosome[k] = 0
            else:

                mutated_chromosome[k] = 1
        else:

            mutated_chromosome[k] = chromosome[k]

    return mutated_chromosome

def Tournament(chromosome_list, performance_list):
    ''' Tournament '''
    
    cross_left = 3
    cross_right = 6

    crossover_chromosome_list = []

    for k in range(len(performance_list)):
        
        # Find parents
        parent1 = FindParent(performance_list)
        parent2 = FindParent(performance_list)

        # Crossover
        crossover_chromosome = Crossover(chromosome_list[parent1], chromosome_list[parent2], cross_left, cross_right)
        crossover_chromosome_list.append(crossover_chromosome)

    
    return crossover_chromosome_list

def FindParent(performance_list):
    ''' Find parent '''

    # Normalize performance_list

    # Cumulative sum
    
    r = np.random.rand()

    for k in range(len(performance_list)):
        if r < performance_list[k]:
            return k

    return len(performance_list)-1


