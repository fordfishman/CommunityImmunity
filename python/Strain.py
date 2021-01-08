## Ford Fishman

import Spacer; import Crispr; import PhageReceptor
from Phage import Phage
import general as gen
import Error as er

##########################################################################################################

class Strain():

    """
    name (str): tag for this strain
    parent (str): the id of the parent that spawned this strain
    a (float): competition coefficient (hosts)
    b (float): birth rate (hosts/timestep)
    c (float): cost of CRISPR [0,1]
    f (float): failure rate of CRISPR immunity
    crispr (Crispr): CRISPR locus, if none is present, keep as None
    phReceptors (dict): dictionary of this strains phage receptors
    activeReceptors (dict)
    pop (float)
    """
    def __init__(self, name:str, a:float, b:float, c:float, f:float,y:float, crispr:Crispr= None,phReceptors:dict = None, pop:float = 1):
        self.name = name
        self.crispr = crispr
        self.phReceptors = phReceptors
        self.activeReceptors = dict()
        self.pop = pop
        self.ipop = 0 # infected pop
        self.intrinsicFitness = 1 # CHANGE THIS AT SOME POINT WHEN I ADD IN RESOURCES
        self.infections = dict()
        # self.parent = parent
        self.descendents = 0 # the number of descendent strains from this strain

        self.a = a
        self.b = b
        self.c = c
        self.y = y
        self.f = f

        self.record = gen.initRecord()
        
        if not phReceptors is None:
            for receptorName in phReceptors: # for all receptors in a strain
                # if the receptor just became active, add to active dict
                if phReceptors[receptorName].isExpressed:
                    self.activeReceptors[receptorName] = phReceptors[receptorName]

##########################################################################################################
    """
    Main timestep function
    """
    def timestep(self, N:int, l:float, step:int):
        """
        N (int): total host density
        l (float): probability an infection lyses
        step (int): current timestep
        """

        i = len(self.record) # how long this strain has been around

        a = self.a
        b = self.b
        c = self.c 
        y = self.y
        Nh = self.pop

        currentInfections = sum(self.infections.values())
        self.ipop = currentInfections

        # labeling system for record
        strainType = 'novel'

        if not self.hasCost(): # set cost to 0 if strain does not have CRISPR-associated cost
            c = 0
        if self.name == 's1':
            strainType = 'initial'
        
        # fitness 
        r = b * ( 1 - c )
        # self.ipop += currentInfections - l * self.ipop # infected pop 
        
        # self.pop = ( r*Nh - currentInfections)/( 1 + ( N/a )**y ) 
        # based upon Beverton-Holt model
        self.pop = ( r*Nh )/( 1 + ( N/a )**y ) - currentInfections
        # lotka-volterra predation
        # self.pop += r*Nh*(1 - N/a) - currentInfections  
        # Extinction threshold
        if self.pop+self.ipop < 0.1: 
            self.pop = 0
            self.ipop = 0

        if not self.crispr is None:
            spacers = len(self.crispr)
        else:
            spacers = None 

        self.record.loc[i] = [step, self.name, self.pop, self.pop-Nh, (self.pop-Nh)/Nh, strainType, spacers]


        return None

##########################################################################################################

    """
    Other functions
    """
    def addSpacer(self,protospacer:str):
        if not self.crispr is None:

            self.crispr.addSpacer(protospacer)
        return None
        

    def isSusceptible(self, receptor:str): 
        """Does this strain have the phage receptor to be susceptible to this phage"""
        # if the receptor is in strain and is expressed
        return receptor in self.phReceptors and self.phReceptors[receptor].isExpressed  


    def isImmune(self, protospacers:set):
        """Does the strain have CRISPR resistance"""
        isImmune = False

        if self.crispr is None:
            isImmune = False
        else:
            for protospacer in protospacers:

                if  self.crispr.hasSpacer( protospacer ):
                    isImmune = True
                    break

        return isImmune

    def isInfectable(self, phage:Phage):
        """Can this phage infect this strain"""
        return self.isSusceptible(phage.receptor.name) and not self.isImmune(phage.protospacers)

    def addReceptor(self, receptor:PhageReceptor.PhageReceptor):
        """Add a receptor to strain"""
        self.phReceptors[receptor.name] = receptor

        return None

    def removeReceptor(self, receptorName:str):
        """Removes receptor from strain. Use when a receptor is modified (old one is lost, new one is gained)"""
        
        self.phReceptors.pop(receptorName)
        
        return None

    def hasCost(self):
        """Is there a CRISPR-associated cost to this strain?"""
        cost = False # initialize cost to be 0

        if not self.crispr is None: # if strain has a crispr 

            cost = self.crispr.hasCost

        return cost

    
    # def changeReceptorActivity(self, receptorName:str, active:bool):

    #     phReceptors = self.__phReceptors
    #     activeReceptors = self.__activeReceptors

    #     try:
    #         if not receptorName in phReceptors:
    #             raise er.FunctionError
            
    #     except:
    #         print("phage receptor not in strain")
    #         print()

    #     else:

    #         if active:

    #     finally:
    #         return None
    

    
        
        


##########################################################################################################

"""
Print tests
"""

# timesteps = 4
# pS = 0.001 # prob of spacer forming if infection occurs
# b = 1.6 # initial max intrinsic growth
# a = 0.04 # competition coefficient (inter and intraspecific)
# crisprCost = 0.5 # fitness cost of having active CRISPR system

# receptor1 = PhageReceptor.PhageReceptor( name = "r" + "1" ) # change numbering system

# crispr0 = Crispr.Crispr(spacerLength=1)
# print(crispr0.hasCost())
# crispr0.newSpacer("GATG")
# print(crispr0.hasCost())
# strain1 = Strain(
#     name = "s" + "1",
#     crispr = crispr0,
#     phReceptors = {
#         receptor1.name():receptor1
#         },
#     pop = 3
# )
# for i in range(timesteps): 
#     strain1.timestep(strain1.pop(), a, b, crisprCost)

# print(strain1.pop())