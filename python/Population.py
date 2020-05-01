## Ford Fishman

import numpy as np
from copy import copy, deepcopy
import Strain; import Phage; import PhageReceptor; import Crispr
import general as gen; from Enums import Type

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
        self.__resSize = 0
        self.__vulnSize = 0

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
    
    def resSize(self):
        return self.__resSize
    
    def vulnSize(self):
        return self.__vulnSize

    def phageSet(self):
        return self.__phageSet

    def richness(self):
        return len(self.__strains)
    
    def spacerRichness(self)->list:
        """each strain has a value"""
        strains = self.__strains.values()
        return [ strain.crisprLength() for strain in strains ]

##########################################################################################################

    """
    Main timestep function
    """

    def timestep(self, N:int, a:float, b:float, c:float,  y:float, absP:float, p:dict, pS:float):
        """
        N (int): total community size
        a (float): competition coefficient
        b (float): intrinsic growth rate
        c (float): cost of crispr
        """
        strains = deepcopy(self.__strains)
        newStrains = list()
        extinctStrainNames = set()
        resSize = 0 
        vulnSize = 0 # number of vulnerable hosts at this timestep

        for strainName,strain in strains.items():

            if strain.pop() == 0:
                extinctStrainNames.add(strainName)
                print("1 host strain dead")
                continue

            phages = p[strainName] # the phages infecting this strain
            phageTotalList = [phage.pop() for phage in phages.values()]
            phageTotal = sum(phageTotalList) # the number of infecting phages per strain

            strains[strainName].timestep(N=N,a=a,b=b,c=c,y=y,absP=absP,p=phageTotal)

            if strain.hasCost():
                resSize += strains[strainName].pop()
            else:
                vulnSize += strains[strainName].pop()

            for phage in phages.values():

                n = absP*phage.pop()*strain.pop() # number of infections

                if n == 0:
                    continue
                
                newStrains = [*newStrains, *self.newSpacer(n,p=pS,strain=strain, phage=phage) ]

        newStrains_dict = {newStrain.name():newStrain for newStrain in newStrains}
        strains.update(newStrains_dict)
        finalStrains = {strainName:strain for strainName,strain in strains.items() if not strainName in extinctStrainNames}
        self.__strains = finalStrains
        
        self.__updatePopSize() # re-calculate population size
        self.__resSize = resSize
        self.__vulnSize = vulnSize

        return None

##########################################################################################################

    """
    Other functions
    """

    def vulnerableStrains(self, phageGenome:str, receptor:str):

        # strains = deepcopy( self.__strains )
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
        # vStrains = deepcopy(self.vulnerableStrains(phageGenome, receptor))
        vStrains = self.vulnerableStrains(phageGenome, receptor)

        for strain in vStrains:

            self.__phageSet[strain.name()] = strain.phages()

        return None

    @gen.runProcess
    def newSpacer(self, *args, p:float=0, strain:Strain=None, phage:Phage=None, num:int=0):

        """Adds a new strain based on another but with a new spacer"""

        if strain is None or phage is None:
            raise KeyError("Need to specify strain name and phage")

        oldSpacers = deepcopy(strain.crispr().spacers())
        # oldSpacers = strain.crispr().spacers()


        crispr = Crispr.Crispr( oldSpacers
        ) 

        spacer = crispr.makeSpacer( genome = phage.genome() ) # generate a new spacer based on phage genome
        crispr.addSpacer(spacer)
            
        newName = gen.generateName(Type.STRAIN, len(self.__strains)+1)

        newStrain = Strain.Strain( # organism with new spacer is a new strain, one starting cell
            name=newName,
            crispr=crispr,
            phReceptors=strain.phReceptors(),
            pop = 1.0
        )
        # newStrain.addSpacer(spacer)

        return newStrain
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