## Ford Fishman

from copy import copy, deepcopy
import timeit
from Phage import Phage
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

    def __init__(self, c:float, l:float, strains:dict, phages:dict=None, nameGenerator=None):

        self.base_c = c
        self.l = l
        self.strains = strains
        self.phages = phages

        # keep track of all phages and all strains
        self.allN = list()
        self.allP = list()

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
        # self.record = gen.initRecord()
        self.nameGenerator = nameGenerator

##########################################################################################################


    """
    Attribute functions
    """

    def richness(self)->int:
        """Returns strain richness for entire community"""
        return len(self.N.reshape(-1))
    
    def phageRichness(self)->int:
        return len(self.P.reshape(-1))

    def spacerRichness(self)->list:
        """how many unique spacers are in this community?"""

        uniqueSpacers = set()
        

        for strain in self.strains.values():

            spacers = strain.crispr.spacers # the set of spacers this strain has
            uniqueSpacers.update(spacers)

        return len(uniqueSpacers)

    def fullDF(self):
        """ 
        """
        # strains
        maxStrains = len(self.N.reshape(-1))
        Ns = self.allN

        timesteps = len(Ns)

        strains_mat = np.zeros([maxStrains,timesteps])

        for i,N in enumerate(Ns):
            zeros = np.zeros( maxStrains - len(N) )
            full_N = np.append(N, zeros)
            strains_mat[:,i] = full_N

        df = pd.DataFrame(None)

        for j in range(maxStrains):

            pop = strains_mat[j,:]
            dpop = np.ediff1d(pop, to_begin=0)
            name = self.strainIDS[j]
            strain = self.strains[name]
            type_ = strain.type
            spacers = strain.numSpacers()
            dpop_pop = np.where((dpop==0) & (pop==0), 0, dpop/pop)

            df_j = pd.DataFrame(
                {
                    'timestep': list( range(timesteps) ),
                    'name': self.strainIDS[j],
                    'pop': pop,
                    'dpop': dpop,
                    'dpop_pop': dpop_pop,
                    'type':type_,
                    'spacers': spacers
                }
            )

            df = df.append(df_j, ignore_index=True)

        # phage
        maxPhages = len(self.P.reshape(-1))
        Ps = self.allP

        phages_mat = np.zeros([maxPhages,timesteps])

        for i,P in enumerate(Ps):
            zeros = np.zeros( maxPhages - len(P) )
            full_P = np.append(P, zeros)
            phages_mat[:,i] = full_P

        for j in range(maxPhages):

            pop = phages_mat[j,:]
            dpop = np.ediff1d(pop, to_begin=0)
            name = self.phageIDS[j]
            phage = self.phages[name]
            spacers = phage.numProto()
            dpop_pop = np.where((dpop==0) & (pop==0), 0, dpop/pop)

            df_j = pd.DataFrame(
                {
                    'timestep': list( range(timesteps) ),
                    'name': self.phageIDS[j],
                    'pop': pop,
                    'dpop': dpop,
                    'dpop_pop': dpop_pop,
                    'type':'phage',
                    'spacers': spacers
                }
            )

            df = df.append(df_j, ignore_index=True)

        return df

    # def fullRecord(self):
    #     """Record of all strain and phage data from the beginning to the current step"""
    #     record = self.record

    #     for phage in self.phages.values():
            
    #         record = record.append( phage.record, ignore_index = True)  

    #     for strain in self.strains.values():

    #         record = record.append( strain.record, ignore_index = True )  
        
    #     return record

##########################################################################################################

    """
    Main timestep function
    """
    
    def timestep(self, step:int):

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
        self.pS = np.array([])

        for strain in strains.values():

            self.strainIDS.append(strain.name)
            self.a = np.append(self.a, strain.a)
            self.b = np.append(self.b, strain.b)
            self.f = np.append(self.f, strain.f)
            self.y = np.append(self.y, strain.y)
            self.N = np.append(self.N, strain.pop)
            self.allN.append(self.N)
            self.I_n = np.append(self.I_n, strain.ipop)
            self.pS = np.append(self.pS, strain.pS).reshape(-1,1)

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
        self.m = np.array([])
        self.P = np.array([])
        self.adsorbed = np.array([])
        self.I_p = np.array([])
        

        for phage in phages.values():

            self.phageIDS.append(phage.name)
            self.d = np.append(self.d, phage.d)
            self.adsp = np.append(self.adsp, phage.adsp)
            self.beta = np.append(self.beta, phage.beta)
            self.m = np.append(self.m, phage.m)
            self.adsorbed = np.append(self.adsorbed, phage.adsorbed)
            self.P = np.append(self.P, phage.pop)
            self.allP.append(self.P)
            self.I_p = np.append(self.I_p, phage.adsorbed)
        
        self.reshape(host=False)

        return None

    def appendPhages(self, phage):

        self.phageIDS.append(phage.name)
        self.d = np.append(self.d, phage.d)
        self.adsp = np.append(self.adsp, phage.adsp)
        self.beta = np.append(self.beta, phage.beta)
        self.m = np.append(self.m, phage.m)
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
        self.pS = np.append(self.pS, strain.pS).reshape(-1,1)

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
        newPhages = self.mutations(self.beta*I_p)
        
        for phage in newPhages.values():
            self.appendPhages(phage)
        
        # new spacers
        newStrains = self.newSpacers(I)
        
        for strain in newStrains.values():
            self.appendStrains(strain)

        # update 
        self.strains.update(newStrains)
        self.phages.update(newPhages)
        self.allN.append(self.N)
        self.allP.append(self.P)

        return None


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
            # mutation = np.random.choice([Mutation.PROTOCHANGE, Mutation.PROTODELETE, Mutation.PROTOADD])
            mutation = Mutation.PROTOCHANGE
        else:
            mutation = Mutation.PROTOADD
        if mutation is Mutation.PROTOADD or mutation is Mutation.PROTOCHANGE:
            protoName = self.nameGenerator.generateName(Type.PROTO)
            protospacers = phage.mutate(mutation, protoName)
        
        else:
            protospacers = phage.mutate(mutation)

        protospacers = phage.mutate(mutation)
        newName = self.nameGenerator.generateName(Type.PHAGE)

        p = {
            'adsp': phage.adsp,
            'beta': phage.beta, 
            'd': phage.d, 
            'm': phage.m,
            'fitness':phage.fitness
        }

        if not phage.evoTraits is None:

            for param in phage.evoTraits: # limit scale of change to current magnitude
                p[param] = np.random.normal(0, 0.2) * p[param] + p[param]

        newPhage = Phage(
            name=newName, 
            adsp = p['adsp'], 
            beta = p['beta'], 
            d = p['d'], 
            m = p['m'],
            receptor=phage.receptor,
            # genome=newGenome, 
            protospacers=protospacers,
            fitness = p['fitness'],
            pop = 1.0,
            evoTraits=phage.evoTraits
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

        crispr = Crispr( oldSpacers ) 

        if len(phage.protospacers) !=0:

            spacer = np.random.choice( list(phage.protospacers) )
            crispr.addSpacer(spacer)
            
        newName = self.nameGenerator.generateName(Type.STRAIN)

        p = {
            'a': strain.a,
            'b': strain.b, 
            'c': strain.c, 
            'f': strain.f,
            'y': strain.y,
            'pS': strain.pS
        }

        if not strain.evoTraits is None:

            for param in strain.evoTraits: # limit scale of change to current magnitude
                p[param] = np.random.normal(0, 0.2) * p[param] + p[param]

        newStrain = Strain( # organism with new spacer is a new strain, one starting cell
            name=newName,
            type_ = 'novel',
            a=p['a'], b=p['b'], c=p['c'], f=p['f'], y=p['y'], pS=p['pS'],
            crispr=crispr,
            phReceptors=strain.phReceptors,
            pop = 1.0,
            evoTraits=strain.evoTraits
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

    # def trimNetwork(self)->None:

    #     strains, phages = self.allStrains, self.allPhages

    #     record = self.fullRecord()

    #     self.trimmedStrains = {strainName:strain for strainName,strain in strains.items() if record[record['name']==strainName]['pop'].max()>=100}
    #     self.trimmedPhages = {phageName:phage for phageName,phage in phages.items() if record[record['name']==phageName]['pop'].max()>=100}
    #     print('Trimmed Strains:\t%s'% len(self.trimmedStrains))
    #     print('Trimmed Phages:\t%s'% len(self.trimmedPhages))

    #     return None

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


