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
    pS (float): 
    strains (dict): dictionary of strains by name
    phages (dict): dictionary of phages by name
    k (int): carrying capacity
    """

    def __init__(self, pS:float, m:float, l:float, strains:dict, phages:dict=None):

        self.pS = pS
        self.m = m
        self.l = l
        # self.populations = dict()
        self.strains = strains
        self.phages = phages

        # will keep history of all strains and phages
        self.allStrains = strains
        self.allPhages = phages

        # total population at each timestep
        self.NList = list() # hosts
        self.PList = list() # phage
        self.IList = list() # spacer variants
        self.SList = list() # initial host variant

        self.__updateComSize()
        self.strainTimes = list()
        self.phageTimes = list()
        self.dfTimes = list()
        self.otherTimes = list()
        self.record = gen.initRecord()

##########################################################################################################


    """
    Attribute functions
    """

    def richness(self)->int:
        """Returns strain richness for entire community"""
        return len(self.strains)

    def spacerRichness(self)->list:
        """how many unique spacers are in this community?"""

        uniqueSpacers = set()
        

        for strain in self.strains.values():

            spacers = strain.crispr.spacers # the set of spacers this strain has
            uniqueSpacers.update(spacers)

        return len(uniqueSpacers)

    def fullRecord(self):
        """Record of all strain and phage data from the beginning to the current step"""
        record = self.record

        for phage in self.phages.values():
            
            record = record.append( phage.record, ignore_index = True)  

        for strain in self.strains.values():

            record = record.append( strain.record, ignore_index = True )  
        
        return record

##########################################################################################################

    """
    Main timestep function
    """
    def timestep(self, step:int, maxtimestep=1000):
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

        # pops = deepcopy(self.populations)
        phages = deepcopy(self.phages)
        strains = deepcopy(self.strains)
        pS = self.pS
        m = self.m
        l = self.l
        N = self.N # total cells in system
        totalImmune = 0
        totalSusceptible = 0
        newPhages = list()
        newStrains = list()
        # keep track of each member going extinct
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

                if strain.isInfectable(phage): # can phage infect this host?
                    newInfections = phage.pop*phage.adsp*strain.pop
                    failedInfections = 0

                else:
                    # new infections due to CRISPR failure rate
                    newInfections = phage.pop*phage.adsp*strain.pop*strain.f
                    # failed infections remove phage from system
                    failedInfections = phage.pop*phage.adsp*strain.pop*(1-strain.f)

                currentInfections = strain.infections.get(phageName, 0)
                phage.adsorbed += newInfections + failedInfections # phage keeps tracked of adsorbed members
                lysisEvents = currentInfections * l # how many infections are now lysing?
                # lysisEvents = np.random.binomial(n=currentInfections,p=l) 
                phage.lysisEvents += lysisEvents
                newVirions = lysisEvents*phage.beta
                strain.infections[phageName] = currentInfections + newInfections - lysisEvents

                newPhages += self.phageMutation(newVirions, p=m, phage=deepcopy(phage))

                if not strain.crispr is None: # if the strain has a crispr locus
                    newStrains += self.newSpacer(newInfections,p=pS, strain=deepcopy(strain), phage=phage)

            strains[strainName].timestep(N=N,l=l,step=step)

        
        # phage timestep 
        for phageName, phage in phages.items():
            proxyPop = phage.lysisEvents*phage.beta + phage.pop
            if proxyPop < 1: # remove extinct phages from community
                n = np.random.uniform()

                if n<0.01 or proxyPop<0.1: extinctPhageNames.add(phageName)

            # if proxyPop < 0.8: # remove extinct phages from community
            #     extinctPhageNames.add(phageName)
            #     continue

            phages[phageName].timestep(step)

        # update dictionaries to account for new and extinct strains and phages
        newPhages_dict = {phage.name:phage for phage in newPhages}
        phages.update(newPhages_dict)
        self.allPhages.update(newPhages_dict)
        finalPhages = {phageName:phage for phageName,phage in phages.items() if not phageName in extinctPhageNames}

        newStrains_dict = {newStrain.name:newStrain for newStrain in newStrains}
        strains.update(newStrains_dict)
        self.allStrains.update(newStrains_dict)
        finalStrains = {strainName:strain for strainName,strain in strains.items() if not strainName in extinctStrainNames}
        
        # Add records of extinct members to master record

        for phageName in extinctPhageNames:
            
            self.record = self.record.append( self.phages[phageName].record, ignore_index = True )  

        for strainName in extinctStrainNames:

            self.record = self.record.append( self.strains[strainName].record, ignore_index = True ) 

        # update community attributes 

        self.phages = finalPhages
        self.strains = finalStrains    
        
        self.__updateComSize()

        self.IList.append(totalImmune)
        self.SList.append(totalSusceptible)

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
        Generates a phage mutant with qualities based upon the original strain.
        This function is an argument of gen.runProcess(), and is called as many times
        as new phages are formed.
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
        f = strain.f 
        y = strain.y

        crispr = Crispr( oldSpacers ) 

        spacer = crispr.makeSpacer( genome = phage.genome ) # generate a new spacer based on phage genome
        crispr.addSpacer(spacer)
            
        newName = gen.generateName(Type.STRAIN)

        newStrain = Strain( # organism with new spacer is a new strain, one starting cell
            name=newName,
            a=a, b=b, c=c, f=f, y=y,
            crispr=crispr,
            phReceptors=strain.phReceptors,
            pop = 1.0
        )
        
        return newStrain
    
    def globalInfectionEdges(self, trim:bool=True):
        """
        Finds edges in global infection network
        return: list of tuples
        """

        if trim:

            self.trimNetwork()
            strains, phages = self.trimmedStrains, self.trimmedPhages

        else:

            strains, phages = self.allStrains, self.allPhages

        strainNames, phageNames = list(strains.keys()), list(phages.keys())

        # df = pd.DataFrame(None, columns=phageNames, index=strainNames) # initialize blank dataframe

        edges = list()

        for strainName in strainNames:

            strain = strains[strainName]

            # df.loc[strainName] = [ strain.isInfectable( phages[phageName] ) for phageName in phageNames ]

            edges += [ (strainName, phageName) for phageName in phageNames if strain.isInfectable( phages[phageName] )]

        return edges

    def trimNetwork(self)->None:

        strains, phages = self.allStrains, self.allPhages

        record = self.fullRecord()

        self.trimmedStrains = {strainName:strain for strainName,strain in strains.items() if record[record['name']==strainName]['pop'].max()>=100}
        self.trimmedPhages = {phageName:phage for phageName,phage in phages.items() if record[record['name']==phageName]['pop'].max()>=100}
        print('Trimmed Strains:\t%s'% len(self.trimmedStrains))
        print('Trimmed Phages:\t%s'% len(self.trimmedPhages))

        return None





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

        self.N = comSize
        self.NList.append(comSize)

        phagePop = 0

        for pop in self.phagesPopDict().values():
            phagePop += pop

        self.PList.append(phagePop)

        return None

        
    

##########################################################################################################
    
"""
Print testing
"""


