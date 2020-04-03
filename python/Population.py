## Ford Fishman

import numpy as np
import Strain; import Phage; import PhageReceptor; import Crispr

##########################################################################################################

class Population():

    """
    Population
    A set of strains located together
    -----------------
    name (str)
    strains (dict(str:Strain)): dict of strains making up a single population
    phageSet (dict(str:str)): dict of strain name to infectable phages
    """

    def __init__(self, name:str, strains:dict):

        self.__name = name
        self.__strains = strains
        self.__updatePopSize()
        self.__phageSet = dict()

##########################################################################################################
     
    """
    Attribute functions
    """

    def name(self):
        return self.__name

    def getStrain(self, strainName:str):
        return self.__strains[strainName]
    
    def strains(self):
        return self.__strains

    def popSize(self):
        return self.__popSize

    def phageSet(self):
        return self.__phageSet

##########################################################################################################

    """
    Main timestep function
    """

    def timestep(self, N:int, a:float, b:float, c:float,  y:float, absP:float, p:dict):
        """
        N (int): total community size
        a (float): competition coefficient
        b (float): intrinsic growth rate
        c (float): cost of crispr
        """
        
        for strainName in self.__strains:

        #     p = self.__totalPhages( 
        #         phageSet = self.__phageSet[strainName],
        #         phageDict = p
        # )
            phageTotal = p[strainName] # the number of infecting phages per strain

            self.__strains[strainName].timestep(N=N,a=a,b=b,c=c,y=y,absP=absP,p=phageTotal)

        self.__updatePopSize() # re-calculate population size

        return None

##########################################################################################################

    """
    Other functions
    """

    def vulnerableStrains(self, phageGenome:str, receptor:str):

        strains = self.__strains
        # vStrains = dict() # dictionary for storing vulnerable strains
        vStrains = set()

        """optimize when I add in event types"""

        for strain in strains.values():

            if strain.isVulnerable(receptor) and not strain.isImmune(phageGenome):

                vStrains.add(strain)

        return vStrains

    def subTotal(self, phageGenome:str, receptor:str):
        """Returns the total number of hosts susceptible to a phage"""
        total = 0

        for strain in self.vulnerableStrains(phageGenome, receptor):

            total += strain.pop()

        return total

    def updatePhageSet(self, receptor:str, phageGenome:str, phageName:str):
        """
        
        """
        vStrains = self.vulnerableStrains(phageGenome, receptor)

        for strain in vStrains:

            self.__phageSet[strain.name()] = strain.phages()

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
        
##########################################################################################################

    """
    Private methods
    """
    def __updatePopSize(self):
        """Re-total population size across all strains"""
        strains = self.__strains
        popSize = 0

        for i in strains:
            popSize += strains[i].pop()
        
        self.__popSize = popSize

        return None

    def __totalPhages(self, phageSet, phageDict:dict):
        """
        Adds phage pops that a strain is vulnerable to
        """
        total = 0
        
        for phage in phageSet:

            total += phageDict[phage]

        return total