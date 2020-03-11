## Ford Fishman

import numpy as np
import PhageReceptor; import Population

NUCLEOTIDES = ("A","C","G","T")

class Phage():

    """
    Phage
    Infects Strains with the right PhageReceptors
    --------------------
    name (str): the name of the phage
    genomeLength (int): number of genes in genome
    population (int): population size 
    targetPop (Population): the population this phage infects
    """

    def __init__(self, name:str, phageReceptor, targetPop:Population.Population,genome:str = None, genomeLength:int = 100, pop = 100):

        self.__name = name
        self.__pop = pop
        self.__phageReceptor = phageReceptor
        self.__targetPop = targetPop

        if not genome is None: 
            self.__genome = genome

        else: # Generates a pseudo-genome for the phage
            # Essentially provides it with pseudo-spacers
            self.__genome = np.random.choice(NUCLEOTIDES, size=genomeLength, replace=True)
    
    """
    Attribute functions
    """
    def genome(self):
        return self.__genome

    def phageReceptor(self):
        return self.__phageReceptor
    
    def name(self):
        return self.__name
    
    def pop(self):
        return(self.__pop)

    def targetPop(self):
        return(self.__targetPop)

    """
    Other functions
    """

