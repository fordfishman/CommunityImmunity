## Ford Fishman

import Phage; import Population

class Community():
    """
    population (dict): dictionary of populations by name
    phages (dict): dictionary of phages by name
    k (int): carrying capacity
    """
    def __init__(self, k:int, populations:dict, phages:dict=None):
        self.__populations = populations
        self.__phages = phages

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
    
    """
    Other functions
    """

    def infection(self):
        return None 
    