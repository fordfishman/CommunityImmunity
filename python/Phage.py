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

    def __init__(self, name:str, absp:float, beta:float, d:float, receptor:PhageReceptor, genome:str = None, genomeLength:int = 100, pop:float = 1, fitness:float = 1):

        self.name = name
        self.pop = pop
        self.receptor = receptor
        self.fitness = fitness
        # self.__strains = set()
        # self.__targetPop = targetPop

        self.absp = absp
        self.beta = beta
        self.d = d
        self.newInfections = 0
        self.lysisEvents = 0

        if not genome is None: 
            self.genome = genome

        else: # Generates a pseudo-genome for the phage
            # Essentially provides it with pseudo-spacers
            self.genome = "".join(np.random.choice(NUCLEOTIDES, size=genomeLength, replace=True))

##########################################################################################################

    """
    Main timestep function
    """

    def timestep(self):
        """
        Ns (float): number of susceptible hosts
        inf:
        l:
        """

        # absp = self.absp
        beta = self.beta
        d = self.d

        Np = self.pop # phage pop
        newInf = self.newInfections
        lysisEvents = self.lysisEvents
        
        # self.__pop += absp*(beta-1)*Ns*Np*self.__fitness - d*Np
        # self.pop += l*beta*inf - absp*Np*Ns - d*Np
        self.pop += beta*lysisEvents - newInf - d*Np
        # if beta*inf < 1 and self.pop < 1: self.pop = 0
        if self.pop < 1: self.pop = 0

        self.lysisEvents = 0
        self.newInfections = 0
        return None

##########################################################################################################
    
    """
    Attribute functions
    """
    # def genome(self):
    #     return self.__genome

    # def receptor(self):
    #     return self.__receptor
    
    # def name(self):
    #     return self.__name
    
    # def pop(self):
    #     return(self.__pop)

    # def fitness(self):
    #     return(self.__fitness)

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

        gLength = len(self.genome) # genome length
        i = np.random.choice( range(0, gLength) )
        # make genome into a list to change position
        genomeList = list(self.genome) 
        genomeList[i] = nt
        newGenome = "".join(genomeList)
        # print(newGenome)
        return newGenome

    @mutate.register(Mutation.DELETION)
    def _(self, mutation) -> str:
        
        gLength = len(self.genome) # genome length
        i = np.random.choice( range(0, gLength) )
        # make genome into a list to delete position
        genomeList = list(self.genome)
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
# a = {"1":{"a":3,"b":4},"2":{"a":100}}

# print(*a)
