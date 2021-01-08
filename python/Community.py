## Ford Fishman

from copy import copy, deepcopy
import timeit
from Phage import Phage; from Population import Population
from Strain import Strain; from Crispr import Crispr
import general as gen; from Enums import Mutation, Type
from network import createNetwork, adjacencyMatrix, addToNetwork
import pandas as pd; import numpy as np
##########################################################################################################

class Community():

    """
    pS (float): 
    strains (dict): dictionary of strains by name
    phages (dict): dictionary of phages by name
    k (int): carrying capacity
    """

    def __init__(self, c:float, pS:float, m:float, l:float, strains:dict, phages:dict=None, nameGenerator=None):

        self.base_c = c
        self.pS = pS
        self.m = m
        self.l = l
        self.strains = strains
        self.phages = phages

        self.initStrains(); self.initPhages(); self.initNetwork()

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
        self.nameGenerator = nameGenerator

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
    # def timestep(self, step:int):
    #     """
    #     step (int): current time step

    #     """
    #     phages = deepcopy(self.phages)
    #     strains = deepcopy(self.strains)
    #     pS = self.pS
    #     m = self.m
    #     l = self.l
    #     N = self.N_tot # total cells in system
    #     totalImmune = 0
    #     totalSusceptible = 0
    #     newPhages = list()
    #     newStrains = list()
    #     # keep track of each member going extinct
    #     extinctPhageNames = set()
    #     extinctStrainNames = set()

    #     for strainName, strain in strains.items():

    #         if strain.pop == 0 and strain.ipop == 0:
    #             extinctStrainNames.add(strainName)
    #             continue

    #         if strain.hasCost():
    #             totalImmune += strains[strainName].pop
    #         else:
    #             totalSusceptible += strains[strainName].pop

    #         for phageName,phage in phages.items():

    #             if strain.isInfectable(phage): # can phage infect this host?
    #                 newInfections = phage.pop*phage.adsp*strain.pop
    #                 failedInfections = 0

    #             else:
    #                 # new infections due to CRISPR failure rate
    #                 newInfections = phage.pop*phage.adsp*strain.pop*strain.f
    #                 # failed infections remove phage from system
    #                 failedInfections = phage.pop*phage.adsp*strain.pop*(1-strain.f)

    #             currentInfections = strain.infections.get(phageName, 0)
    #             phage.adsorbed += newInfections + failedInfections # phage keeps tracked of adsorbed members
    #             lysisEvents = currentInfections * l # how many infections are now lysing?
    #             phage.lysisEvents += lysisEvents
    #             newVirions = lysisEvents*phage.beta
    #             strain.infections[phageName] = currentInfections + newInfections - lysisEvents

    #             newPhages += self.phageMutation(newVirions, p=m, phage=deepcopy(phage), nameGenerator=self.nameGenerator)

    #             if not strain.crispr is None: # if the strain has a crispr locus
    #                 newStrains += self.newSpacer(newInfections,p=pS, strain=deepcopy(strain), phage=deepcopy(phage), nameGenerator=self.nameGenerator)


    #         strains[strainName].timestep(N=N,l=l,step=step)

        
    #     # phage timestep 
    #     for phageName, phage in phages.items():
    #         proxyPop = phage.lysisEvents*phage.beta + phage.pop
    #         if proxyPop < 1: # remove extinct phages from community
    #             n = np.random.uniform()

    #             if n<0.01 or proxyPop<0.1: extinctPhageNames.add(phageName)

    #         phages[phageName].timestep(step)

    #     # update dictionaries to account for new and extinct strains and phages
    #     newPhages_dict = {phage.name:phage for phage in newPhages}
    #     phages.update(newPhages_dict)
    #     self.allPhages.update(newPhages_dict)
    #     finalPhages = {phageName:phage for phageName,phage in phages.items() if not phageName in extinctPhageNames}

    #     newStrains_dict = {newStrain.name:newStrain for newStrain in newStrains}
    #     strains.update(newStrains_dict)
    #     self.allStrains.update(newStrains_dict)
    #     finalStrains = {strainName:strain for strainName,strain in strains.items() if not strainName in extinctStrainNames}
        
    #     # Add records of extinct members to master record

    #     for phageName in extinctPhageNames:
            
    #         self.record = self.record.append( self.phages[phageName].record, ignore_index = True )  

    #     for strainName in extinctStrainNames:

    #         self.record = self.record.append( self.strains[strainName].record, ignore_index = True ) 

    #     # update community attributes 
    #     self.phages = finalPhages
    #     self.strains = finalStrains    
        
    #     self.__updateComSize()


    #     return None
    
    def timestep2(self, step:int):

        self.growth()
        self.__updateComSize()
        


        return None


##########################################################################################################

    """
    Other functions
    """

    def initStrains(self):
        """
        Initialize arrays with attributes
        """

        strains = self.strains

        self.strainIDS = list()
        self.a = np.array([])
        self.b = np.array([])
        self.c = np.array([])
        self.f = np.array([])
        self.y = np.array([])
        self.N = np.array([])
        self.I_n = np.array([])

        for strain in strains.values():

            self.strainIDS.append(strain.name)
            self.a = np.append(self.a, strain.a)
            self.b = np.append(self.b, strain.b)
            self.f = np.append(self.f, strain.f)
            self.y = np.append(self.y, strain.y)
            self.N = np.append(self.N, strain.pop)
            self.I_n = np.append(self.I_n, strain.ipop)

            if not strain.hasCost():
                c = 0
            else:
                c = self.base_c

            self.c = np.append(self.c, c)

        self.reshape(phage=False)
        
        return None

    def initPhages(self):
        """
        Initialize arrays with attributes
        """

        phages = self.phages

        self.phageIDS = list()
        self.d = np.array([])
        self.adsp = np.array([])
        self.beta = np.array([])
        self.P = np.array([])
        self.adsorbed = np.array([])
        self.I_p = np.array([])

        for phage in phages.values():

            self.phageIDS.append(phage.name)
            self.d = np.append(self.d, phage.d)
            self.adsp = np.append(self.adsp, phage.adsp)
            self.beta = np.append(self.beta, phage.beta)
            self.adsorbed = np.append(self.adsorbed, phage.adsorbed)
            self.P = np.append(self.P, phage.pop)
            self.I_p = np.append(self.I_p, phage.adsorbed)
        
        self.reshape(host=False)

        return None

    def appendPhages(self, phage):

        self.phageIDS.append(phage.name)
        self.d = np.append(self.d, phage.d)
        self.adsp = np.append(self.adsp, phage.adsp)
        self.beta = np.append(self.beta, phage.beta)
        self.adsorbed = np.append(self.adsorbed, phage.adsorbed)
        self.P = np.append(self.P, phage.pop)
        self.I_p = np.append(self.I_p, phage.adsorbed)
        self.reshape()

        return None
    
    def appendStrains(self, strain):

        self.strainIDS.append(strain.name)
        self.a = np.append(self.a, strain.a)
        self.b = np.append(self.b, strain.b)
        self.f = np.append(self.f, strain.f)
        self.y = np.append(self.y, strain.y)
        self.N = np.append(self.N, strain.pop)
        self.I_n = np.append(self.I_n, strain.ipop)

        if not strain.hasCost():
            c = 0
        else:
            c = self.base_c

        self.c = np.append(self.c, c)  
  
        self.reshape()

        return None
    
    def reshape(self, phage:bool=True, host:bool=True):
        """ 
        turn numpy arrays into vertical rather than dimensionless vectors
        """
        if phage:
            self.d = self.d.reshape(-1,1)
            self.beta = self.beta.reshape(-1,1)
            self.adsorbed = self.adsorbed.reshape(-1,1)
            self.P = self.P.reshape(-1,1)
            self.adsp = self.adsp.reshape(-1,1)
            self.I_p = self.I_p.reshape(-1,1)
        
        if host:
            self.a = self.a.reshape(-1,1)
            self.b = self.b.reshape(-1,1)
            self.c = self.c.reshape(-1,1)
            self.f = self.f.reshape(-1,1)
            self.y = self.y.reshape(-1,1)
            self.N = self.N.reshape(-1,1)

        return None        

    def growth(self):

        self.reshape()

        # infections
        # theoretical number of infections if all hosts are infectible to all phagesS
        I_pot = np.matmul((self.N), (self.P*self.adsp).transpose())
        I = self.A*I_pot #
        I_p = np.sum(I, axis=0).reshape(-1,1) # number of new infections for each phage
        I_n = np.sum(I, axis=1).reshape(-1,1) # number of new infections for each strain
        # self.I_p += I_p
        # self.I_n += I_n

        # hosts
        r = self.b * (1-self.c)
        self.N = ( r*self.N )/(1 + (np.sum(self.N)/self.a)**self.y) - I_n
        self.N = np.where(self.N < 0, 0, self.N)

        # phage
        # self.P += self.beta*self.I_p*self.l - I_p - self.d*self.P
        self.P += self.beta*I_p - I_p - self.d*self.P
        self.P = np.where(self.P < 0, 0, self.P)

        # mutations
        newPhages = self.mutations(I_p)
        
        for phage in newPhages.values():
            self.appendPhages(phage)
        
        # new spacers
        newStrains = self.newSpacers(I)
        
        for strain in newStrains.values():
            self.appendStrains(strain)

        # update 
        self.strains.update(newStrains)
        self.phages.update(newPhages)

        return None

    def phagesPopDict(self):
        """
        returns dictionary of phage name to phage pop (str:float)
        """
        phages = self.phages
        popDict = {phageName:phage.pop for phageName, phage in phages.items()}

        return popDict

    def mutations(self, n)->dict:
        """
        n: list of infection events by phage
        """
        n = list(n.reshape(-1))
        numEvents = np.random.binomial(n,self.m)
        newPhages = dict()
        strains = self.strains
        edges = list()

        if np.count_nonzero(numEvents) < 1:
            return newPhages

        for i, n in enumerate(numEvents):

            if n < 1: # if there are no mutation events for the phage
                    continue 

            phageName = self.phageIDS[i]
            phage = deepcopy(self.phages.get(phageName))

            for ind in range(n):
                newPhage = self.phageMutation(phage)
                newPhages[newPhage.name] = newPhage

                edges += [ (strainName, newPhage.name) for strainName, strain in strains.items() if strain.isInfectable( newPhage )]

        self.addToNetwork(edges, phageNames=list(newPhages.keys()))

        return newPhages

    def phageMutation(self, phage:Phage) -> Phage:
        """
        Generates a phage mutant with qualities based upon the original strain.
        """
        # equal prob. of a change, addition, or deletion of protospacer
        if len(phage.protospacers) >= 1:
            mutation = np.random.choice([Mutation.PROTOCHANGE, Mutation.PROTODELETE, Mutation.PROTOADD])
        else:
            mutation = Mutation.PROTOADD
        # newGenome = phage.mutate(mutation)
        if mutation is Mutation.PROTOADD or mutation is Mutation.PROTOCHANGE:
            protoName = self.nameGenerator.generateName(Type.PROTO)
            protospacers = phage.mutate(mutation, protoName)
        
        else:
            protospacers = phage.mutate(mutation)

        protospacers = phage.mutate(mutation)
        newName = self.nameGenerator.generateName(Type.PHAGE)
        mod = phage.fitness
        # mod = 1 
        fitness = np.random.uniform(0,1.5) * mod # some cost to mutating genome
        fitness = 1

        newPhage = Phage(
            name=newName, 
            adsp = phage.adsp, 
            beta = phage.beta, 
            d = phage.d, 
            receptor=phage.receptor,
            # genome=newGenome, 
            protospacers=protospacers,
            fitness = fitness,
            pop = 1.0
        )
        
        return newPhage

    def newSpacers(self, I)->dict:
        """
        I: numpy array of infection events by host and phage
        """
        numEvents = np.random.binomial(list(I), self.pS)
        newStrains = dict()
        phages = self.phages
        nstrain = I.shape[0]
        nphage = I.shape[1]
        edges = list()

        if np.count_nonzero(numEvents) < 1:
            return newStrains

        for i in range(nstrain):
            for j in range(nphage):

                n = numEvents[i,j]
                if n < 1: # if there are no mutation events for the phag
                    continue 

                strainName = self.strainIDS[i]
                phageName = self.phageIDS[j]
                strain = deepcopy(self.strains.get(strainName))
                phage = deepcopy(self.phages.get(phageName))

                for ind in range(n):
                    newStrain = self.newSpacer(strain, phage)
                    newStrains[newStrain.name] = newStrain

                    edges += [ (newStrain.name, phageName) for phageName, phage in phages.items() if newStrain.isInfectable( phage )]
        
        self.addToNetwork(edges, strainNames=list(newStrains.keys()))
        
        return newStrains

    def newSpacer(self, strain, phage):

        """Adds a new strain based on another but with a new spacer"""

        oldSpacers = strain.crispr.spacers
        
        a = strain.a
        b = strain.b 
        c = strain.c
        f = strain.f 
        y = strain.y

        crispr = Crispr( oldSpacers ) 

        if len(phage.protospacers) !=0:

            spacer = np.random.choice( list(phage.protospacers) )
            crispr.addSpacer(spacer)
            
        newName = self.nameGenerator.generateName(Type.STRAIN)

        newStrain = Strain( # organism with new spacer is a new strain, one starting cell
            name=newName,
            a=a, b=b, c=c, f=f, y=y,
            crispr=crispr,
            phReceptors=strain.phReceptors,
            pop = 1.0
        )
        
        return newStrain
    
    def initNetwork(self):
        """ Create initial adjacency matrix """
        edges = self.globalInfectionEdges()

        B = createNetwork(edges, self.strains.keys(), self.phages.keys())
        self.net = B
        A_df = adjacencyMatrix(B).sort_index().sort_index(axis=1)
        self.A = np.array(A_df)

        return None
    
    def addToNetwork(self, edges:list, strainNames:list=None, phageNames:list=None):

        addToNetwork(self.net, edges, strainNames, phageNames)
        A_df = adjacencyMatrix(self.net).sort_index().sort_index(axis=1)
        self.A = np.array(A_df)

        return None
    
    def globalInfectionEdges(self, trim:bool=True):
        """
        Finds edges in global infection network
        return: list of tuples
        """

        # if trim:

        #     self.trimNetwork()
        #     strains, phages = self.trimmedStrains, self.trimmedPhages

        # else:

        strains, phages = self.strains, self.phages

        strainNames, phageNames = list(strains.keys()), list(phages.keys())


        edges = list()

        for strainName in strainNames:

            strain = strains[strainName]

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

        self.N_tot = np.sum(self.N)
        self.NList.append(self.N_tot)

        self.P_tot = np.sum(self.P)
        self.PList.append(self.P_tot)

        # index of susceptible 
        ind_suc = [0,]
        currentSus = self.N[ind_suc]
        self.SList.append(float(currentSus))
        self.IList.append(float(self.N_tot-currentSus))

        return None

##########################################################################################################
    
"""
Print testing
"""


