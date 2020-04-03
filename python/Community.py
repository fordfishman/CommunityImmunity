## Ford Fishman

import Phage; import Population; import Strain
import pandas as pd
##########################################################################################################

class Community():

    """
    population (dict): dictionary of populations by name
    phages (dict): dictionary of phages by name
    k (int): carrying capacity
    """

    def __init__(self, populations:dict, phages:dict=None):
        self.__populations = populations
        self.__strains = dict()
        self.__phages = phages
        self.__comSizeOverTime = list()
        self.__phagePopOverTime = list()
        self.__updateComSize()
        self.__updateStrainPhageDF()

##########################################################################################################


    """
    Attribute functions
    """

    def getPopulation(self, popName):

        try: 
            self.__populations[popName]

        except KeyError:
            print("Population does not exist")

        else: 
            return self.__populations[popName]

    def getPhage(self, phageName):

        try: 
            self.__phages[phageName]

        except KeyError:
            print("Phage does not exist")

        else: 
            return self.__phages[phageName]
    
    # def k(self):
    #     return self.__k

    def totalComSize(self):
        return self.__totalComSize

    def comSizeOverTime(self):
        return self.__comSizeOverTime

    def phagePopOverTime(self):
        return self.__phagePopOverTime

##########################################################################################################

    """
    Main timestep function
    """
    def timestep(self, aH:float, bH:float, c:float, y:float, bP:float, absP:float, dP:float, pS:float):
        """
        N (int): total community size
        a (float): competition coefficient
        b (float): intrinsic growth rate
        c (float): cost of crispr
        y (float):
        bP (float):
        absP (float):
        dP (float):
        pS (float): probability of spacer forming during infection
        """
        # tHost=0 # for testing

        pops = self.__populations
        phages = self.__phages

        for popName, pop in pops.items():

            p = self.infectingPhages(pop)
            # p = self.__phagePopOverTime[-1] 
            # print("Phages:\t"+str(p))
            # p = self.phagesPopDict()

            pops[popName].timestep(N=self.totalComSize(),a=aH,b=bH,c=c,y=y,absP=absP,p=p)

        for phageName, phage in phages.items(): 

            totalHosts = self.totalVulnerable(phageName)
            # totalHosts = self.__totalComSize 
            
            # print("Hosts:\t"+str(totalHosts))
            # print(self.__totalComSize)
            # print(self.__populations["pop1"].getStrain("s1").isVulnerable(phage.receptor().name()))
            # print(self.vulnerableStrains(phage.genome(), phage.receptor()))
            phages[phageName].timestep(Ns=totalHosts, bP=bP, absP=absP, dP=dP)

            # self.__populations[popName].timestep(N,aH,bH,c,y,absP,p=self.phagesPopDict())
            
            # print(str(bH-absP*p))
            # vStrains = self.vulnerableStrains( phage.genome(), phage.receptor().name() )

        self.__populations = pops
        self.__phages = phages
        self.__updateComSize()
        self.__updateStrainPhageDF()
        # print(self.__StrainPhageDF)
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
        pops = self.__populations
        vStrains = set()

        for pop in pops.values():

            vStrains.update( pop.vulnerableStrains(phageGenome,receptorName) )#  all vulnerable strains from each population

        return vStrains

    def totalVulnerable(self, phageName:str):
        """returns the total number of vulnerable hosts to a specific phage"""
        # pops = self.__populations
        total = 0

        areStrainsVulnerable = self.__StrainPhageDF[phageName]

        for strainName, isVulnerable in areStrainsVulnerable.items():
            # print(isVulnerable)
            if isVulnerable:

                total += self.__strains[strainName].pop()

        return total

    def infectingPhages(self, pop:Population):
        """
        Return the number of phages that infect each strain in pop (dict) 
        """
        phages = dict() # dict of strainName: Array of bools
        strainNames = set(pop.strains().keys())

        for strainName, row in self.__StrainPhageDF.iterrows():

            if strainName in strainNames:

                phageTotal = [self.phagesPopDict()[phageName] for phageName,infects in row.items() if infects]
                phages[strainName] = sum(phageTotal)

        return phages

    def phagesPopDict(self):
        """
        returns dictionary of phage name to phage pop (str:float)
        """
        phages = self.__phages
        popDict = {phageName:phage.pop() for phageName, phage in phages.items()}

        return popDict


    

##########################################################################################################

    """
    Private methods
    """
    
    def __updateComSize(self):
        """Re-total community size across all populations"""
        populations = self.__populations
        comSize = 0

        for pop in populations.values():
            comSize += pop.popSize()
            self.__strains.update( pop.strains() )

        self.__totalComSize = comSize
        self.__comSizeOverTime.append(comSize)

        phagePop = 0

        for pop in self.phagesPopDict().values():
            phagePop += pop

        self.__phagePopOverTime.append(phagePop)

        return None

    def __updateStrainPhageDF(self):
        """
        """

        strains = self.__strains

        df = pd.DataFrame(None, index = strains.keys(), columns=self.__phages.keys())

        for phageName, phage in self.__phages.items():

            vStrains = self.vulnerableStrains(phageGenome=phage.genome(), receptorName=phage.receptor().name())
            # print(vStrains) 
            Strains = {strainName:(strain in vStrains) for strainName, strain in strains.items()} # which strains are vulnerable to this phage

            df[phageName] = df.index.map(Strains)

        self.__StrainPhageDF = df

        return None
##########################################################################################################
    
"""
Print testing
"""


# testing pandas
import string

a = pd.DataFrame(None, index=range(0,4), columns=list(string.ascii_lowercase[0:5]))
# a.loc[2] = [1,2,3,4,5]

b = {
    1 : 1,
    0 : 2,
    2 : 0,
    3 : 4
}
# print(a.index)
a["a"] = a.index.map(b)
for x,y in a.iterrows():
    for z,d in y.items():
        # print(z,d)
        continue
    continue
# print(a["a"])
# print(a.loc[0]["a"])

c = list(range(0,5))
print(
    # [x for x in c if x!=3]
)
