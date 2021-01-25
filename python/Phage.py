## Ford Fishman

import numpy as np; import pandas as pd
import PhageReceptor
from Enums import Mutation, Type
from copy import deepcopy
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

    def __init__(self, name:str, adsp:float, beta:float, d:float, receptor:PhageReceptor, protospacers:set, genome:str = None,genomeLength:int = 100, pop:float = 1, fitness:float = 1):

        self.name = name
        self.pop = pop
        self.receptor = receptor
        self.fitness = fitness
        self.record = None
        self.type = 'phage'
        self.adsp = adsp
        self.beta = beta
        self.d = d
        self.adsorbed = 0
        self.lysisEvents = 0
        self.descendents = 0 # the number of descendent phage types from this phage

        self.record = gen.initRecord()

        # if not protospacers is None:
        self.protospacers = protospacers

        # else:
        #     self.protospacers = set()

        #     for i in range(numProto):

        #         self.protospacers.add( gen.generateName(Type.PROTO) )
            

##########################################################################################################

    """
    Main timestep function
    """

    def timestep(self, step:int):
        """
        step (int): current timestep
        """

        i = len(self.record) # how long this phage has been around
        beta = self.beta
        d = self.d

        Np = self.pop # phage pop
        adsorbed = self.adsorbed
        lysisEvents = self.lysisEvents
        
        self.pop += beta*lysisEvents - adsorbed - d*Np
        
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
    def numProto(self):
        return len(self.protospacers)

    @gen.dispatch_on_value
    def mutate(self, mutation) -> str:
        pass # only runs if a mutation occurs that's not in the class
        

    @mutate.register(Mutation.PROTOCHANGE)
    def _(self, mutation, protoName=None) -> set:
        """
        Alter a protospacer currently in the strain
        """
        protospacers = deepcopy(self.protospacers)
        protospacer = np.random.choice(list(protospacers))
        protospacers.discard(protospacer)
        protospacers.add(protoName)

        return protospacers


    @mutate.register(Mutation.PROTOADD)
    def _(self, mutation, protoName=None) -> set:
        """
        Add a new protospacer to genome
        """
        protospacers = deepcopy(self.protospacers)
        protospacers.add(protoName)

        return protospacers


    @mutate.register(Mutation.PROTODELETE)
    def _(self, mutation) -> set:
        """
        Delete protospacer from genome
        """
        protospacers = deepcopy(self.protospacers)
        protospacer = np.random.choice(list(protospacers))
        protospacers.discard(protospacer)

        return protospacers


##########################################################################################################

    """
    think about making spacers of variable cost
    """
"""
print testing
"""
# a = {"1":{"a":3,"b":4},"2":{"a":100}}

# print(*a)
