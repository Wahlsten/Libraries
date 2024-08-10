import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from GeneticAlgorithm import *

def Test_InitializeChromosomes():
    """
    Test InitializeChromosomes
    """
    
    population_size = 10
    chromosome_size = 9
    
    chromosomes = InitializeChromosomes(population_size, chromosome_size)
    
    result = True

    for k in range(population_size):
        
        chromosome_tmp = chromosomes[k]

        for l in range(chromosome_size):

            if chromosome_tmp[l] not in [0, 1]:
                result = False

    return result



def Test_InitChromosome():
    '''
    Test InitChromosome
    '''

    chromosome_size = 9

    chromosome = InitChromosome(chromosome_size)
    
    result = True

    for k in range(chromosome_size):

        if chromosome[k] not in [0,1]:
            result = False

    return result

def Test_DecodeChromosome():
    '''
    Test DecodeChromosome
    '''

    chromosome = [0, 1, 1, 0, 1, 0, 0, 1, 1]
    
    result = True

    chromosome_decoded_correct = 2**3 + 2**2 + 2**0 + 2**(-3) + 2**(-4)

    chromosome_decoded = DecodeChromosome(chromosome)
    
    if chromosome_decoded - chromosome_decoded_correct == 0:
        result = False


    return result

def Test_Crossover():
    '''
    Test Crossover
    '''

    chromosome1 = np.array([0, 1, 0, 0, 1, 1, 1, 1, 0])
    chromosome2 = np.array([1, 1, 1, 0, 0, 0, 1, 0, 1])

    lower_lim, upper_lim = 2, 6
    crossover_chromosome_correct = np.concatenate((chromosome1[:lower_lim], chromosome2[lower_lim:upper_lim], chromosome1[upper_lim:]), axis=None)

    crossover_chromosome = Crossover(chromosome1, chromosome2, lower_lim, upper_lim)
    
    result = False

    comparison = crossover_chromosome == crossover_chromosome_correct
    if comparison.all():
        result = True

    return result

def Test_Mutate():
    '''
    Test Mutate
    '''

    chromosome = [0, 1, 0, 1, 1, 0, 0, 1, 1, 1]
    mutate_prob = 0.2

    result = False

    mutated_chromosome = Mutate(chromosome, mutate_prob)

    comparison = chromosome == mutated_chromosome

    if ~comparison.all():
        result = True

    return result

def Test_FindParent():
    '''
    Test FindParent
    '''

    result = True
    
    num_chromosomes = 10

    # Performance list
    performance = np.random.rand(num_chromosomes)

    # Chromosome list


    parent = FindParent(performance)

    if parent not in range(num_chromosomes):
        result = False

    return result

def Test_Tournament():
    '''
    Test Tournament
    '''

    result = True

    num_chromosomes = 10

    performance = np.random.rand(num_chromosomes)
    
    chromosome_list = Tournament(performance)
    
    for k in range(num_chromosomes):
        if chromosome_list[k] not in range(num_chromosomes):
            result = False


    return result

print('Test InitializeChromosomes = ', Test_InitializeChromosomes())

print('Test InitChromosome = ', Test_InitChromosome())

print('Test DecodeChromosome = ', Test_DecodeChromosome())

print('Test Crossover = ', Test_Crossover())

print('Test Mutate = ', Test_Mutate())

print('Test FindParent', Test_FindParent())

print('Test Tournament', Test_Tournament())

