## Ford Fishman

from copy import copy, deepcopy
import timeit
import Phage; import Population; import Strain
import general as gen; from Enums import Mutation, Type
import pandas as pd; import numpy as np
##########################################################################################################

class Community():

    """
    population (dict): dictionary of populations by name
    phages (dict): dictionary of phages by name
    k (int): carrying capacity
    """

    def __init__(self, populations:dict, phages:dict=None):
        self.populations = populations
        self.strains = dict()
        self.phages = phages
        self.comSizeOverTime = list()
        self.phagePopOverTime = list()
        self.resOverTime = list()
        self.vulnOverTime = list()
        self.__updateComSize()
        self.__initStrainPhageDF()
        self.strainTimes = list()
        self.phageTimes = list()
        self.otherTimes = list()

##########################################################################################################


    """
    Attribute functions
    """

    # def getPopulation(self, popName):

    #     try: 
    #         self.populations[popName]

    #     except KeyError:
    #         print("Population does not exist")

    #     else: 
    #         return self.populations[popName]

    # def getPhage(self, phageName):

    #     try: 
    #         self.phages[phageName]

    #     except KeyError:
    #         print("Phage does not exist")

    #     else: 
    #         return self.phages[phageName]
    

    # def totalComSize(self):
    #     return self.__totalComSize

    # def comSizeOverTime(self):
    #     return self.__comSizeOverTime
    
    # def resOverTime(self):
    #     return self.__resOverTime
    
    # def vulnOverTime(self):
    #     return self.__vulnOverTime

    # def phagePopOverTime(self):
    #     return self.__phagePopOverTime
    
    # def phages(self):
    #     """Remove, only using for testing"""
    #     return self.__phages

    def richness(self)->int:
        """Returns strain richness for entire community"""
        richness = 0

        for pop in self.populations.values():
            richness += pop.richness()

        return richness

    def spacerRichness(self)->list:
        """per species spacer richness"""
        pops = self.populations.values()

        return [ pop.spacerRichness() for pop in pops ] 





##########################################################################################################

    """
    Main timestep function
    """
    def timestep(self, pS:float, m:float):
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
        time1 = timeit.default_timer()
        pops = deepcopy(self.populations)
        phages = deepcopy(self.phages)
        totalResistant = 0
        totalVulnerable = 0
        infected = dict()
        time2 = timeit.default_timer()
        

        t0 = timeit.default_timer()

        for popName, pop in pops.items():
            
            popInfected = pop.infected # the phage infections organized by phage for this population

            infected0 = {phageName:self.__updateInfections(phageName,infected,newInfections) for phageName, newInfections in popInfected.items()}
            infected = infected0

            p = self.infectingPhages(pop)

            pops[popName].timestep(N=self.totalComSize,p=p,pS=pS)

            totalResistant += pop.resSize
            totalVulnerable += pop.vulnSize
        t1 = timeit.default_timer()

        self.strainTimes.append(t1-t0)

        newPhages = list()
        extinctPhageNames = set()

        t0 = timeit.default_timer()

        for phageName, phage in phages.items(): 

            inf = infected.get(phageName, 0.0)
            n = 0.5 * phage.beta * inf  # number of new phages

            if n < 0.0001 and phage.pop < 1: # remove extinct phages from community
                extinctPhageNames.add(phageName)
                continue

            Ns = self.totalVulnerable(phageName) # hosts vulnerable and infected to this phage

            phages[phageName].timestep(Ns=Ns,inf=inf)

            

            if n == 0:
                continue

            newPhages += self.phageMutation(phage.beta, inf*0.5, p=m, phage=deepcopy(phage))

        newPhages_dict = {phage.name:phage for phage in newPhages}

        phages.update(newPhages_dict)
        finalPhages = {phageName:phage for phageName,phage in phages.items() if not phageName in extinctPhageNames}
        t1 = timeit.default_timer()

        self.phageTimes.append(t1-t0)
        time3 = timeit.default_timer()
        self.populations = pops
        self.phages = finalPhages
        self.__updateComSize()
        self.__updateStrainPhageDF()
        self.resOverTime.append(totalResistant)
        self.vulnOverTime.append(totalVulnerable)
        time4 = timeit.default_timer()
        self.otherTimes.append(time2-time1 + time4-time3)
        return None

##########################################################################################################

    """
    Other functions
    """

    def infection(self):
        """might add this in with IBM"""
        pass

    def vulnerableStrains(self, phageGenome:str, receptorName:str):
        """returns a dict of vulnerable strains to the phage"""
        pops = self.populations
        vStrains = set()

        for pop in pops.values():

            vStrains.update( pop.vulnerableStrains(phageGenome,receptorName) )#  all vulnerable strains from each population

        return vStrains

    def totalVulnerable(self, phageName:str) -> float:
        """returns the total number of vulnerable and infected hosts for a specific phage"""
        vuln = 0

        areStrainsVulnerable = self.StrainPhageDF[phageName]

        for strainName, isVulnerable in areStrainsVulnerable.items():

            if isVulnerable:

                vuln += self.strains[strainName].pop
        
        return vuln
        
    def infectingPhages(self, pop:Population):
        """
        Return the phages that infect each strain in pop (dict) 
        """
        phages = dict() # dict of phageName:phage
        strainNames = set(pop.strains.keys())

        for strainName, row in self.StrainPhageDF.iterrows():

            if strainName in strainNames:

                phages[strainName] = {phageName:self.phages[phageName] for phageName, infects in row.items() if infects}

        return phages

    def infectedHosts(self):
        """
        """
        pass

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

        newPhage = Phage.Phage(
            name=newName, 
            absp = phage.absp, 
            beta = phage.beta, 
            d = phage.d, 
            receptor=phage.receptor,
            genome=newGenome, 
            fitness = fitness,
            pop = 100
        )
        
        return newPhage



    

##########################################################################################################

    """
    Private methods
    """
    
    def __updateComSize(self):
        """Re-total community size across all populations"""
        populations = self.populations
        comSize = 0

        for pop in populations.values():
            comSize += pop.popSize
            self.strains.update( pop.strains )

        self.totalComSize = comSize
        self.comSizeOverTime.append(comSize)

        phagePop = 0

        for pop in self.phagesPopDict().values():
            phagePop += pop

        self.phagePopOverTime.append(phagePop)

        return None

    def __initStrainPhageDF(self):
        """
        Creates dataframe for identifying which strains are vulnerable to which phage.
        Rows are strain names, columns are phage names.
        Each entry is a boolean. 
        """

        strains = self.strains

        df = pd.DataFrame(None, index = strains.keys(), columns=self.phages.keys())

        for phageName, phage in self.phages.items():

            vStrains = self.vulnerableStrains(phageGenome=phage.genome, receptorName=phage.receptor.name)
            # print(vStrains) 
            Strains = {strainName:(strain in vStrains) for strainName, strain in strains.items()} # which strains are vulnerable to this phage

            df[phageName] = df.index.map(Strains)

        self.StrainPhageDF = df

        return None
        
    def __updateStrainPhageDF(self):
        """
        Creates dataframe for identifying which strains are vulnerable to which phage.
        Rows are strain names, columns are phage names.
        Each entry is a boolean. 
        """

        strains = self.strains

        df = pd.DataFrame(None, index = strains.keys(), columns=self.phages.keys())

        for phageName, phage in self.phages.items():

            vStrains = self.vulnerableStrains(phageGenome=phage.genome, receptorName=phage.receptor.name)
            # print(vStrains) 
            Strains = {strainName:(strain in vStrains) for strainName, strain in strains.items()} # which strains are vulnerable to this phage

            df[phageName] = df.index.map(Strains)

        self.StrainPhageDF = df

        return None
    
    def __updateInfections(self, phageName:str, oldInfections:dict, newInfections:float)->float:
        """
        updates active infections for particular phage
        """
        # phageName = phage.name
        activeInfections = 0
        
        # try:
        #     oldInfections[phageName]
        # except KeyError:
        #     activeInfections = newInfections
        # else:
        activeInfections = newInfections + oldInfections.get(phageName, 0.0)

        return activeInfections


##########################################################################################################
    
"""
Print testing
"""


