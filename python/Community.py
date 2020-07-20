## Ford Fishman

from copy import copy, deepcopy
import timeit
from Phage import Phage; from Population import Population
from Strain import Strain; from Crispr import Crispr
import general as gen; from Enums import Mutation, Type
import pandas as pd; import numpy as np
##########################################################################################################

class Community():

    """
    population (dict): dictionary of populations by name
    phages (dict): dictionary of phages by name
    k (int): carrying capacity
    """

    def __init__(self, pS, m, l, strains:dict, phages:dict=None):
        self.pS = pS
        self.m = m
        self.l = l
        self.populations = dict()
        self.strains = strains
        self.phages = phages
        self.comSizeOverTime = list()
        self.phagePopOverTime = list()
        self.imOverTime = list()
        self.susOverTime = list()
        self.__updateComSize()
        # self.__updateStrainPhageDF()
        self.strainTimes = list()
        self.phageTimes = list()
        self.dfTimes = list()
        self.otherTimes = list()

##########################################################################################################


    """
    Attribute functions
    """

    def richness(self)->int:
        """Returns strain richness for entire community"""
        return len(self.strains)

    def spacerRichness(self)->list:
        """per species spacer richness"""
        pops = self.populations.values()

        return [ pop.spacerRichness() for pop in pops ] 





##########################################################################################################

    """
    Main timestep function
    """
    def timestep(self, step:int):
        """
        N (int): total community size
        a (float): competition coefficient
        b (float): intrinsic growth rate
        c (float): cost of crispr
        y (float):
        bP (float): phage birth rate
        adsP (float): aborption rate of phage
        dP (float): decay rate of phage
        pS (float): probability of spacer forming during infection per host
        """

        pops = deepcopy(self.populations)
        phages = deepcopy(self.phages)
        strains = deepcopy(self.strains)
        pS = self.pS
        m = self.m
        l = self.l
        N = self.totalComSize
        totalImmune = 0
        totalSusceptible = 0
        # infected = dict()
        newPhages = list()
        newStrains = list()
        extinctPhageNames = set()
        extinctStrainNames = set()


        for strainName, strain in strains.items():

            if strain.pop == 0 and strain.ipop == 0:
                extinctStrainNames.add(strainName)
                continue

            if strain.hasCost():
                totalImmune += strains[strainName].pop
            else:
                totalSusceptible += strains[strainName].pop

            for phageName,phage in phages.items():
                # consider changing exinction parameters

                if strain.isInfectable(phage): # can phage infect this host?
                    # newInfections = np.random.binomial(n=phage.pop*strain.pop,p=phage.adsp)
                    newInfections = phage.pop*phage.adsp*strain.pop

                else:
                    newInfections = 0
                    # newInfections = np.random.binomial(n=phage.pop*strain.pop,p=1e-5*phage.adsp)

                currentInfections = strain.infections.get(phageName, 0)
                phage.newInfections += newInfections # phage keeps track of its infections
                lysisEvents = currentInfections * l # how many infections are now lysing?
                # lysisEvents = np.random.binomial(n=currentInfections,p=l) # how many infections are now lysing?
                phage.lysisEvents += lysisEvents
                newVirions = lysisEvents*phage.beta
                strain.infections[phageName] = currentInfections + newInfections - lysisEvents

                newPhages += self.phageMutation(newVirions, p=m, phage=deepcopy(phage))
                newStrains += self.newSpacer(newInfections,p=pS, strain=deepcopy(strain), phage=phage)

            strains[strainName].timestep(N=N,l=l)

        for phageName, phage in phages.items():

            if phage.lysisEvents*phage.beta + phage.pop < 0.1: # remove extinct phages from community
                extinctPhageNames.add(phageName)
                continue

            phages[phageName].timestep()

        newPhages_dict = {phage.name:phage for phage in newPhages}

        phages.update(newPhages_dict)
        finalPhages = {phageName:phage for phageName,phage in phages.items() if not phageName in extinctPhageNames}
        newStrains_dict = {newStrain.name:newStrain for newStrain in newStrains}
        strains.update(newStrains_dict)
        finalStrains = {strainName:strain for strainName,strain in strains.items() if not strainName in extinctStrainNames}
            
        self.phages = finalPhages
        self.strains = finalStrains    
        
        self.__updateComSize()

        self.imOverTime.append(totalImmune)
        self.susOverTime.append(totalSusceptible)

        return None

##########################################################################################################

    """
    Other functions
    """

    def phagesPopDict(self):
        """
        returns dictionary of phage name to phage pop (str:float)
        """
        phages = self.phages
        popDict = {phageName:phage.pop for phageName, phage in phages.items()}

        return popDict

    @gen.runProcess
    def phageMutation(self, *args, p:float=0, phage:Phage=None) -> Phage:
        """
        """
        if phage is None:

            raise KeyError("Need to specify strain name and phage")

        newGenome = phage.mutate(Mutation.SNP)
        newName = gen.generateName(Type.PHAGE)
        mod = phage.fitness
        # mod = 1 
        fitness = np.random.uniform(0,1.5) * mod # some cost to mutating genome
        # fitness = 1

        newPhage = Phage(
            name=newName, 
            adsp = phage.adsp, 
            beta = phage.beta, 
            d = phage.d, 
            receptor=phage.receptor,
            genome=newGenome, 
            fitness = fitness,
            pop = 1.0
        )
        
        return newPhage

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

        crispr = Crispr( oldSpacers ) 

        spacer = crispr.makeSpacer( genome = phage.genome ) # generate a new spacer based on phage genome
        crispr.addSpacer(spacer)
            
        newName = gen.generateName(Type.STRAIN)

        newStrain = Strain( # organism with new spacer is a new strain, one starting cell
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

##########################################################################################################

    """
    Private methods
    """
    
    def __updateComSize(self):
        """Re-total community size across all populations"""
        # populations = self.populations
        comSize = 0

        # for pop in populations.values():
        #     comSize += pop.popSize
        #     self.strains.update( pop.strains )
        for strain in self.strains.values():
            comSize += strain.pop

        self.totalComSize = comSize
        self.comSizeOverTime.append(comSize)

        phagePop = 0

        for pop in self.phagesPopDict().values():
            phagePop += pop

        self.phagePopOverTime.append(phagePop)

        return None

        
    

##########################################################################################################
    
"""
Print testing
"""


