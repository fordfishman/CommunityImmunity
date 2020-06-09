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
        self.vulnOverTime = list()
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
        absP (float): aborption rate of phage
        dP (float): decay rate of phage
        pS (float): probability of spacer forming during infection per host
        """
        # tHost=0 # for testing
        # time1 = timeit.default_timer()
        pops = deepcopy(self.populations)
        phages = deepcopy(self.phages)
        strains = deepcopy(self.strains)
        pS = self.pS
        m = self.m
        l = self.l
        N = self.totalComSize
        totalImmune = 0
        totalVulnerable = 0
        # infected = dict()
        newPhages = list()
        newStrains = list()
        extinctPhageNames = set()
        extinctStrainNames = set()
        # structuralChange = False # d
        # df = pd.DataFrame(None, index = self.strains.keys(), columns=self.phages.keys())

        # time2 = timeit.default_timer()
        
        # t0 = timeit.default_timer()

        for strainName, strain in strains.items():

            if strain.pop == 0 and strain.ipop == 0:
                extinctStrainNames.add(strainName)
                continue

            if strain.hasCost():
                totalImmune += strains[strainName].pop
            else:
                totalVulnerable += strains[strainName].pop

            for phageName,phage in phages.items():
                # consider changing exinction parameters

                if strain.isInfectable(phage): # can phage infect this host?
                    newInfections = phage.pop*phage.absp*strain.pop

                else:
                    newInfections = 0

                currentInfections = strain.infections.get(phageName, 0)
                phage.newInfections += newInfections # phage keeps track of its infections
                lysisEvents = currentInfections * l # how many infections are now lysing?
                phage.lysisEvents += lysisEvents
                newVirions = lysisEvents*phage.beta
                strain.infections[phageName] = currentInfections + newInfections - lysisEvents

                newPhages += self.phageMutation(newVirions, p=m, phage=deepcopy(phage))
                newStrains += self.newSpacer(newInfections,p=pS, strain=deepcopy(strain), phage=phage)

            strains[strainName].timestep(N=N,l=l)

        for phageName, phage in phages.items():

            if phage.lysisEvents*phage.beta < 0.0001 and phage.pop < 1: # remove extinct phages from community
                extinctPhageNames.add(phageName)
                continue

            phages[phageName].timestep()

                


        # for popName, pop in pops.items(): # iterate over pops
            
        #     popInfected = pop.infected # the phage infections organized by phage for this population
        #     # iterate through all phages
        #     infected0 = {phageName:self.__updateInfections(phageName,infected,newInfections) for phageName, newInfections in popInfected.items()}
        #     infected = infected0

        #     p = self.infectingPhages(pop) # iterate through strainphagedf

        #     # iterate through all strains, and all infecting phages for each strain
        #     pops[popName].timestep(N=self.totalComSize,p=p,pS=pS, l=l) 

        #     totalImmune += pop.imSize
        #     totalVulnerable += pop.vulnSize

        # t1 = timeit.default_timer()

        # self.strainTimes.append(t1-t0)

        # newPhages = list()
        # extinctPhageNames = set()

        # t0 = timeit.default_timer()

        # for phageName, phage in phages.items(): # iterate over all phages

        #     inf = infected.get(phageName, 0.0)
        #     n = l * phage.beta * inf  # number of new phages

            # if n < 0.0001 and phage.pop < 1: # remove extinct phages from community
            #     extinctPhageNames.add(phageName)
            #     continue
        #     #iterate
        #     Ns = self.totalVulnerable(self.phages[phageName]) # hosts vulnerable and infected to this phage

        #     phages[phageName].timestep(Ns=Ns,inf=inf, l=l)

        #     if n == 0:
        #         continue

        #     newPhages += self.phageMutation(phage.beta, inf*l, p=m, phage=deepcopy(phage))

        newPhages_dict = {phage.name:phage for phage in newPhages}

        phages.update(newPhages_dict)
        finalPhages = {phageName:phage for phageName,phage in phages.items() if not phageName in extinctPhageNames}
        newStrains_dict = {newStrain.name:newStrain for newStrain in newStrains}
        strains.update(newStrains_dict)
        finalStrains = {strainName:strain for strainName,strain in strains.items() if not strainName in extinctStrainNames}
            
        # t1 = timeit.default_timer()

        # self.phageTimes.append(t1-t0)
        # time3 = timeit.default_timer()
        # self.populations = pops
        self.phages = finalPhages
        self.strains = finalStrains    
        
        self.__updateComSize()
        # time4 = timeit.default_timer()
        # time5 = timeit.default_timer()
        # self.__updateStrainPhageDF()
        # time6 = timeit.default_timer()
        self.imOverTime.append(totalImmune)
        self.vulnOverTime.append(totalVulnerable)
        # self.dfTimes.append(time6-time5)
        # self.otherTimes.append(time2-time1 + time4-time3)
        return None

##########################################################################################################

    """
    Other functions
    """

    # def vulnerableStrains(self, phageGenome:str, receptorName:str):
    #     """returns a dict of vulnerable strains to the phage"""
    #     strains = self.strains
    #     # vStrains = set()
    #     vStrains = { 
    #         strain:
    #             True if strain.isVulnerable(receptorName) 
    #                 and not strain.isImmune(phageGenome) 
    #             else False 
    #             for strain in strains.values() 
    #         }
        

    #     return vStrains
    # def vulnStrains(self, col:pd.Series):
    #     """returns a dict of vulnerable strains to the phage"""
    #     strains = col.index
    #     phage = col.name
    #     genome = phage.genome
    #     receptor = phage.receptor.name

    #     x = pd.Series(True, name = phage, index = strains)

    #     for strain, isVulnerable in col.iteritems():

    #         x.loc[strain] = strain.isVulnerable(receptor) and not strain.isImmune(genome)

    #     # strains.where(strains.isImmune(phage.genome) & strains.isVulnerable(phage.receptor.name), False)

    #     return x

    # def totalVulnerable(self, phage:str) -> float:
    #     """returns the total number of vulnerable and infected hosts for a specific phage"""
    #     vuln = 0

    #     areStrainsVulnerable = self.StrainPhageDF[phage]

    #     for strain, isVulnerable in areStrainsVulnerable.items():

    #         if isVulnerable:

    #             vuln += self.strains[strain.name].pop
        
    #     return vuln
        
    # def infectingPhages(self, pop:Population):
    #     """
    #     Return the phages that infect each strain in pop (dict) 
    #     """
    #     phages = dict() # dict of phageName:phage
    #     strainNames = set(pop.strains.keys())

    #     for strain, row in self.StrainPhageDF.iterrows():

    #         if strain.name in strainNames:

    #             phages[strain.name] = {phage.name:self.phages[phage.name] for phage, infects in row.items() if infects}

    #     return phages


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
            absp = phage.absp, 
            beta = phage.beta, 
            d = phage.d, 
            receptor=phage.receptor,
            genome=newGenome, 
            fitness = fitness,
            pop = 1
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

        
    # def __updateStrainPhageDF(self):
    #     """
    #     Creates dataframe for identifying which strains are vulnerable to which phage.
    #     Rows are strain names, columns are phage names.
    #     Each entry is a boolean. 
    #     """

    #     strains = self.strains
    #     phages = self.phages

    #     df = pd.DataFrame(True, index = strains.values(), columns=self.phages.values())
    #     # t1 = timeit.default_timer()
    #     df = df.apply(
    #         lambda col: 
    #             pd.Series(
    #                 self.vulnerableStrains( 
    #                     col.name.genome,
    #                     col.name.receptor.name 
    #                     )
    #                 ),
    #         axis=0
    #     )
    #     # t2 = timeit.default_timer()
    #     # t3 = timeit.default_timer()
    #     # df = df.apply(
    #     #     lambda col: 
    #     #         pd.Series(self.vulnStrains(col)),
    #     #     axis=0
    #     # )
    #     # # if self.StrainPhageDF != df: raise NameError("stupid")
    #     # t4 = timeit.default_timer()
    #     self.StrainPhageDF = df

    #     # print((t4-t3)/(t2-t1))
    #     return None
    
    # def __updateInfections(self, phageName:str, oldInfections:dict, newInfections:float)->float:
    #     """
    #     updates active infections for particular phage
    #     """
    #     # phageName = phage.name
    #     activeInfections = 0
        
    #     # try:
    #     #     oldInfections[phageName]
    #     # except KeyError:
    #     #     activeInfections = newInfections
    #     # else:
    #     activeInfections = newInfections + oldInfections.get(phageName, 0.0)

    #     return activeInfections


##########################################################################################################
    
"""
Print testing
"""


