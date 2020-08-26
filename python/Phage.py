## Ford Fishman

import numpy as np; import pandas as pd
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

    def __init__(self, name:str, adsp:float, beta:float, d:float, receptor:PhageReceptor, genome:str = None, genomeLength:int = 100, pop:float = 1, fitness:float = 1):

        self.name = name
        self.pop = pop
        self.receptor = receptor
        self.fitness = fitness
        # self.__strains = set()
        # self.__targetPop = targetPop
        self.record = None
        self.type = 'phage'

        self.adsp = adsp
        self.beta = beta
        self.d = d
        self.adsorbed = 0
        self.lysisEvents = 0

        self.record = gen.initRecord()

        if not genome is None: 
            self.genome = genome

        else: # Generates a pseudo-genome for the phage
            # Essentially provides it with pseudo-spacers
            self.genome = "".join(np.random.choice(NUCLEOTIDES, size=genomeLength, replace=True))

##########################################################################################################

    """
    Main timestep function
    """

    def timestep(self, step:int):
        """
        step (int): current timestep
        """

        i = len(self.record) # how long this phage has been around
        # adsp = self.adsp
        beta = self.beta
        d = self.d

        Np = self.pop # phage pop
        adsorbed = self.adsorbed
        lysisEvents = self.lysisEvents
        
        # self.__pop += adsp*(beta-1)*Ns*Np*self.__fitness - d*Np
        # self.pop += l*beta*inf - adsp*Np*Ns - d*Np
        self.pop += beta*lysisEvents - adsorbed - d*Np
        # if beta*inf < 1 and self.pop < 1: self.pop = 0
        # if self.pop < 1: self.pop = 0
        
        # record data from this timestep
        # columns: "timestep", "name", "pop", "dpop", "dpop_pop","type", "spacers"
        self.record.loc[i] = [step, self.name, self.pop, self.pop-Np, (self.pop-Np)/Np, self.type, None]

        self.lysisEvents = 0
        self.adsorbed = 0
        return None


##########################################################################################################

    """
    Other functions
    """

    @gen.dispatch_on_value
    def mutate(self, mutation) -> str:
        pass # only runs if a mutation occurs that's not in the class
        

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
