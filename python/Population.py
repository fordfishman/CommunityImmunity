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

        self.name = name
        self.strains = strains
        self.__updatePopSize()
        self.phageSet = dict()
        self.resSize = 0
        self.vulnSize = 0
        self.infected = dict() # current number of infections due to phage x

##########################################################################################################
     
    """
    Attribute functions
    """

    # def name(self):
    #     return self.__name

    # def getStrain(self, strainName:str):
    #     return self.strains[strainName]

    # def popSize(self):
    #     return self.__popSize

    # def ipopSize(self):
    #     return self.__ipopSize

    # def infected(self):
    #     return self.__infected
    
    # def resSize(self):
    #     return self.__resSize
    
    # def vulnSize(self):
    #     return self.__vulnSize

    # def phageSet(self):
    #     return self.__phageSet

    def richness(self):
        return len(self.strains)
    
    def spacerRichness(self)->list:
        """each strain has a value"""
        strains = self.strains.values()
        return [ len(strain.crispr) for strain in strains ]

##########################################################################################################

    """
    Main timestep function
    """

    def timestep(self, N:int, p:dict, pS:float, l:float):
        """
        N (int): total community size
        a (float): competition coefficient
        b (float): intrinsic growth rate
        c (float): cost of crispr
        """
        strains = deepcopy(self.strains)
        newStrains = list()
        extinctStrainNames = set()
        # allPhages = set() # all phage names that are infecting
        # infected = dict()
        resSize = 0 
        vulnSize = 0 # number of vulnerable hosts at this timestep


        for strainName,strain0 in strains.items():

            strain = deepcopy(strain0)

            if strain.pop == 0 and strain.ipop == 0:
                extinctStrainNames.add(strainName)
                continue
            
            phages = p[strainName] # the phages infecting this strain
            
            absorbedPhages = 0

            for phage in phages.values():

                absorbedPhages += phage.pop*phage.absp

                n = self.__activeInfections(phage.absp, strain, phage, l)  # number of infections

                self.infected[phage.name] = n

                if n == 0:
                    continue
                
                newStrains += self.newSpacer(n,p=pS,strain=strain, phage=phage) 

            strains[strainName].timestep(N=N,absP=absorbedPhages, l=l)
            if strain.hasCost():
                resSize += strains[strainName].pop
            else:
                vulnSize += strains[strainName].pop

        newStrains_dict = {newStrain.name:newStrain for newStrain in newStrains}
        strains.update(newStrains_dict)
        finalStrains = {strainName:strain for strainName,strain in strains.items() if not strainName in extinctStrainNames}
        self.strains = finalStrains
        
        self.__updatePopSize() # re-calculate population size
        self.resSize = resSize
        self.vulnSize = vulnSize
        # self.__infected = infected

        return None

##########################################################################################################

    """
    Other functions
    """

    # def vulnerableStrains(self, phageGenome:str, receptor:str):

    #     # strains = deepcopy( self.strains )
    #     strains = self.strains 
    #     # vStrains = dict() # dictionary for storing vulnerable strains
    #     vStrains = {strain for strain in strains.values() if strain.isVulnerable(receptor) and not strain.isImmune(phageGenome)} 

    #     """optimize when I add in event types"""

    #     # for strain in strains.values():

    #     #     if strain.isVulnerable(receptor) and not strain.isImmune(phageGenome):

    #     #         vStrains.add(strain)

    #     return vStrains



    @gen.runProcess
    def newSpacer(self, *args, p:float=0, strain:Strain=None, phage:Phage=None):

        """Adds a new strain based on another but with a new spacer"""

        if strain is None or phage is None:
            raise KeyError("Need to specify strain name and phage")
        oldSpacers = strain.crispr.spacers
        
        a = strain.a
        b = strain.b 
        c = strain.c 
        y = strain.y

        crispr = Crispr.Crispr( oldSpacers) 

        spacer = crispr.makeSpacer( genome = phage.genome ) # generate a new spacer based on phage genome
        crispr.addSpacer(spacer)
            
        newName = gen.generateName(Type.STRAIN)

        newStrain = Strain.Strain( # organism with new spacer is a new strain, one starting cell
            name=newName,
            a=a,
            b=b,
            c=c,
            y=y,
            crispr=crispr,
            phReceptors=strain.phReceptors,
            pop = 1.0
        )
        
        return newStrain
    # mutations of receptors, new spacers?

    def receptorMod(self, strainName:str, newStrainName:str, newReceptorName:str):
        """At random, a strain in population will have a random receptor modified"""
        # strain = np.random.choice( # random strain chosen
        #     a = self.strains,
        #     replace = True
        # )

        # receptor = np.random.choice( # random receptor chosen
        #     a = strain.phReceptors(),
        #     replace = True
        # )

        # newStrain = Strain.Strain( # organism with new spacer is a new strain, one starting cell
        #     name=newStrainName,
        #     crispr=strain.crispr(),
        #     phReceptors=strain.phReceptors()
        # )
        
        # fitness = np.random.random()

        # newReceptor = PhageReceptor.PhageReceptor(
        #     name = newReceptorName,
        #     fitness = fitness
        # )

        # newStrain.addReceptor(newReceptor)
        # newStrain.removeReceptor(receptor.name())

        # self.strains[newReceptorName] = newReceptor

        # self.__updatePopSize() # re-calculate pop size

        return None
        
##########################################################################################################

    """
    Private methods
    """
    def __updatePopSize(self):
        """Re-total population size across all strains"""
        strains = self.strains
        popSize = 0
        ipopSize = 0

        for strain in strains.values():
            popSize += strain.pop
            ipopSize += strain.ipop
        
        self.popSize = popSize
        self.ipopSize = ipopSize

        return None

    def __totalPhages(self, phageSet, phageDict:dict):
        """
        Adds phage pops that a strain is vulnerable to
        """
        total = 0
        
        for phage in phageSet:

            total += phageDict[phage]

        return total

    def __activeInfections(self, absP:float, strain:Strain, phage:Phage,l:float):
        """
        Returns current number of infections due to this phage
        """
        newInfections = strain.pop*absP*phage.pop
        phageName = phage.name
        totalInfections = 0
        oldInfections = self.infected.get(phageName, 0)
        # try:
        #     self.__infected[phageName]
        # except KeyError:
        #     totalInfections = newInfections
        # else:
        # oldInfections = self.__infected[phageName]
        lysisEvents = l*oldInfections
        totalInfections = newInfections + oldInfections - lysisEvents

        return totalInfections

