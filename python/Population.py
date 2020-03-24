## Ford Fishman

import numpy as np
import Strain; import Phage; import PhageReceptor; import Crispr

class Population():
    """
    Population
    A set of strains located together
    -----------------
    strains (dict): set of strains making up a single population
    """
    def __init__(self, name:str, strains:dict):

        self.__name = name
        self.__strains = strains
        self.__updatePopSize()

    
    """
    Attribute functions
    """

    def name(self):
        return self.__name

    def getStrain(self, strainName:str):
        return self.__strains[strainName]

    def popSize(self):
        return self.__popSize

    """
    Main timestep function
    """
    def timestep(self, N:int, a:float, b:float, c:float):
        """
        N (int): total community size
        a (float): competition coefficient
        b (float): intrinsic growth rate
        c (float): cost of crispr
        """
        # pop = {strainName:strain.timestep(N,a,b,c) for (strainName,strain) in self.__strains.items()} # perform timestep on all strains
        # self.__strains = pop
        for strainName in self.__strains:
            self.__strains[strainName].timestep(N,a,b,c)
        self.__updatePopSize() # re-calculate population size

        return None

    """
    Other functions
    """

    def __updatePopSize(self):
        """Re-total population size across all strains"""
        strains = self.__strains
        popSize = 0

        for i in strains:
            popSize += strains[i].pop()
        
        self.__popSize = popSize

        return None

    def newSpacer(self, strainName:str, newStrainName:str, phage:Phage):

        """Adds a new strain based on another but with a new spacer"""

        strain = self.getStrain(strainName)
        crispr = strain.crispr() 
        crispr.newSpacer( genome = phage.genome() ) # generate a new spacer based on phage genome

        newStrain = Strain.Strain( # organism with new spacer is a new strain, one starting cell
            name=newStrainName,
            crispr=crispr,
            phReceptors=strain.phReceptors()
        )

        self.__strains[newStrainName] = newStrain

        self.__updatePopSize() # re-calculate population size

        return None
    # mutations of receptors, new spacers?

    def receptorMod(self, strainName:str, newStrainName:str, newReceptorName:str):
        """At random, a strain in population will have a random receptor modified"""
        strain = np.random.choice( # random strain chosen
            a = self.__strains,
            replace = True
        )

        receptor = np.random.choice( # random receptor chosen
            a = strain.phReceptors(),
            replace = True
        )

        newStrain = Strain.Strain( # organism with new spacer is a new strain, one starting cell
            name=newStrainName,
            crispr=strain.crispr(),
            phReceptors=strain.phReceptors()
        )
        
        fitness = np.random.random()

        newReceptor = PhageReceptor.PhageReceptor(
            name = newReceptorName,
            fitness = fitness
        )

        newStrain.addReceptor(newReceptor)
        newStrain.removeReceptor(receptor.name())

        self.__strains[newReceptorName] = newReceptor

        self.__updatePopSize() # re-calculate pop size

        return None
        



