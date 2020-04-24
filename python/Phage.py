## Ford Fishman

import numpy as np
import PhageReceptor; import Population
from Enums import Mutation
import general as gen

NUCLEOTIDES = ("A","C","G","T") # tuple of nucleotides

##########################################################################################################

class Phage():

    """
    Phage
    Infects Strains with the right PhageReceptors
    --------------------
    name (str): the name of the phage
    genomeLength (int): number of genes in genome
    population (int): population size 
    # strains (set(str)): set of strains this phage can infect
    """

    def __init__(self, name:str, receptor:PhageReceptor, genome:str = None, genomeLength:int = 100, pop:float = 1, fitness:float = 1):

        self.__name = name
        self.__pop = pop
        self.__receptor = receptor
        self.__fitness = fitness
        # self.__strains = set()
        # self.__targetPop = targetPop

        if not genome is None: 
            self.__genome = genome

        else: # Generates a pseudo-genome for the phage
            # Essentially provides it with pseudo-spacers
            self.__genome = "".join(np.random.choice(NUCLEOTIDES, size=genomeLength, replace=True))

##########################################################################################################

    """
    Main timestep function
    """

    def timestep(self, Ns:float, bP:float, absP:float, dP:float):
        """
        Ns (float): number of susceptible hosts
        bP (float): burst size of phage
        absP (float): absorption rate of phage
        dP (float): decay rate of phage  
        """

        Np = self.__pop # phage pop

        self.__pop += absP*(bP-1)*Ns*Np*self.__fitness - dP*Np

        if self.__pop < 1: self.__pop = 0
        return None

##########################################################################################################
    
    """
    Attribute functions
    """
    def genome(self):
        return self.__genome

    def receptor(self):
        return self.__receptor
    
    def name(self):
        return self.__name
    
    def pop(self):
        return(self.__pop)

    # def targetPop(self):
    #     return(self.__targetPop)

##########################################################################################################

    """
    Other functions
    """

    @gen.dispatch_on_value
    def mutate(self, mutation) -> str:
        pass # only runs if a mutation occurs that's not in thhe class
        

    @mutate.register(Mutation.SNP)
    def _(self, mutation) -> str:
        # print(self.__genome)
        nt = np.random.choice(NUCLEOTIDES) # the new nucleotide

        gLength = len(self.__genome) # genome length
        i = np.random.choice( range(0, gLength) )
        # make genome into a list to change position
        genomeList = list(self.__genome) 
        genomeList[i] = nt
        newGenome = "".join(genomeList)
        # print(newGenome)
        return newGenome

    @mutate.register(Mutation.DELETION)
    def _(self, mutation) -> str:
        
        gLength = len(self.__genome) # genome length
        i = np.random.choice( range(0, gLength) )
        # make genome into a list to delete position
        genomeList = list(self.__genome)
        genomeList.pop(i)
        newGenome = "".join(genomeList)

        return newGenome



##########################################################################################################

    """
    think about making spacers of variable cost
    """
"""
print testing
"""
a = "".join(np.random.choice(NUCLEOTIDES, size=20, replace=True))
# print(len(a))