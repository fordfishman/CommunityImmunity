## Ford Fishman

import Spacer; import Crispr; import PhageReceptor
import Error as er

##########################################################################################################

class Strain():

    """
    name (str)
    crispr (Crispr)
    phReceptors (dict)
    activeReceptors (dict)
    pop (float)
    intrinsicFitness (float)
    phages (set(str)): names of phages this strain can be infected by
    """
    def __init__(self, name:str, crispr:Crispr= None, phReceptors:dict = None, pop:float = 1):
        self.__name = name
        self.__crispr = crispr
        self.__phReceptors = phReceptors
        self.__activeReceptors = dict()
        self.__pop = pop
        self.__intrinsicFitness = 1 # CHANGE THIS AT SOME POINT WHEN I ADD IN RESOURCES
        self.__phages = set()
        
        if not phReceptors is None:
            for receptorName in phReceptors: # for all receptors in a strain
                # if the receptor just became active, add to active dict
                if phReceptors[receptorName].isExpressed():
                    self.__activeReceptors[receptorName] = phReceptors[receptorName]

##########################################################################################################

    """
    Attribute functions
    """

    def name(self):
        return self.__name

    def phReceptors(self):
        return self.__phReceptors

    def getReceptor(self, receptorName:str):
        return self.__phReceptors[receptorName]

    def crispr(self):
        return self.__crispr
    
    def crisprLength(self):
        return len(self.__crispr)

    def pop(self):
        return self.__pop

    def intrinsicFitness(self):
        return self.__intrinsicFitness

    def phages(self):
        return (self.__phages)

##########################################################################################################
    """
    Main timestep function
    """
    def timestep(self, N:int, a:float, b:float, c:float, y:float, absP:float, p:float ):
        """
        N (int): total host density
        a (float): competition coefficient
        b (float): intrinsic growth rate
        c (float): cost of crispr
        absP (float): absorbance rate of phage to host
        p (float): phage density
        y (float):
        """

        if not self.hasCost(): # set cost to 0 if strain does not have CRISPR-associated cost
            c = 0

        # fitness
        r = b * ( self.__intrinsicFitness - c ) - absP*p
            # self.__pop = ( b*self.__pop )/( 1 + ( N/a )**y ) # Beverton-Holt Model
        # reproduce
        self.__pop = ( r*self.__pop )/( 1 + ( N/a )**y ) 
        if self.__pop < 1: self.__pop = 0
        return None

##########################################################################################################

    """
    Other functions
    """
    def addSpacer(self,spacer:str):
        if not self.__crispr is None:

            self.__crispr.addSpacer(spacer)
        return None
        

    def isVulnerable(self, receptor:str): 
        """Does this strain have the phage receptor to be vulnerable to this phage"""
        # if the receptor is in strain and is expressed
        return receptor in self.__phReceptors and self.__phReceptors[receptor].isExpressed()  


    def isImmune(self, phageGenome:str):
        """Does the strain have CRISPR resistance"""
        crispr = self.__crispr
        
        return crispr.hasSpacer( genome = phageGenome )

    def updatePhages(self, receptor:str, phageGenome:str, phageName:str):
        """Updates set of infectable phages"""

        if self.isVulnerable(receptor) and not self.isImmune(phageGenome):

            self.__phages.add(phageName)
        
        else: 

            self.__phages.remove(phageName)

        return None

    def addReceptor(self, receptor:PhageReceptor.PhageReceptor):
        """Add a receptor to strain"""
        self.__phReceptors[receptor.name()] = receptor

        return None

    def removeReceptor(self, receptorName:str):
        """Removes receptor from strain. Use when a receptor is modified (old one is lost, new one is gained)"""
        
        self.__phReceptors.pop(receptorName)
        
        return None

    def hasCost(self):
        """Is there a CRISPR-associated cost to this strain?"""
        cost = False # initialize cost to be 0

        if not self.__crispr is None: # if strain has a crispr 

            cost = self.__crispr.hasCost()

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
    Private methods
    """

    
        
        


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