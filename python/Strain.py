## Ford Fishman

import Spacer; import Crispr; import PhageReceptor; import Phage
import Error as er

class Strain():
    """
    """
    def __init__(self, name:str, crispr:Crispr.Crispr = None, phReceptors:dict = None, pop:int = 1):
        self.__name = name
        self.__crispr = crispr
        self.__phReceptors = phReceptors
        self.__activeReceptors = dict()
        self.__pop = 1
        self.__intrinsicFitness = 1 # CHANGE THIS AT SOME POINT WHEN I ADD IN RESOURCES
        for receptorName in phReceptors: # for all receptors in a strain
            # if the receptor just became active, add to active dict
            if phReceptors[receptorName].isExpressed():
                self.__activeReceptors[receptorName] = phReceptors[receptorName]


    """
    Attribute functions
    """

    def phReceptors(self):
        return self.__phReceptors

    def getReceptor(self, receptorName:str):
        return self.__phReceptors[receptorName]

    def crispr(self):
        return self.__crispr

    def pop(self):
        return self.__pop

    def intrinsicFitness(self):
        return self.__intrinsicFitness

    """
    Other functions
    """
    # Maybe remove references to Phage in everything but community?
    def isVulnerable(self, phage:Phage.Phage): 
        """Does this strain have the phage receptor to be vulnerable to this phage"""
        receptorTarget = phage.phageReceptor().name()
        # if the receptor is in strain and is expressed
        return receptorTarget in self.__phReceptors and self.__phReceptors[receptorTarget].isExpressed()  


    def isImmune(self, phage:Phage.Phage):
        """Does the strain have CRISPR resistance"""
        crispr = self.__crispr
        
        return crispr.hasSpacer( phageGenome = phage.genome() )

    def growth(self):

        return None

    def addReceptor(self, receptor:PhageReceptor.PhageReceptor):
        """Add a receptor to strain"""
        self.__phReceptors[receptor.name()] = receptor

        return None

    def removeReceptor(self, receptorName:str):
        """Removes receptor from strain. Use when a receptor is modified (old one is lost, new one is gained"""
        
        self.__phReceptors.pop(receptorName)
        
        return None
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

    