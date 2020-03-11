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
    def __init__(self, strains:dict):
        self.__strains = strains
    
    """
    Attribute functions
    """

    def getStrain(self, strainName:str):
        return self.__strains[strainName]

    """
    Other functions
    """

    def newSpacer(self, strainName:str, newStrainName:str, phage:Phage.Phage):

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

        return None
        



