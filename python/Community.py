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
        self.__k = k
        self.__updateComSize()


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
    
    def k(self):
        return self.__k

    def totalComSize(self):
        return self.__totalComSize

    """
    Main timestep function
    """
    def timestep(self, N:int, a:float, b:float, c:float):
        """
        N (int): total community size
        a (float): competition coefficient
        b (float): intrinsic growth rate
        c (float): cost of crispr
        """
        # pop = {popName:pop.timestep(N,a,b,c) for (popName,pop) in list(self.__populations.items())} # perform timestep on all populations
        for popName in self.__populations:
            self.__populations[popName].timestep(N,a,b,c)
        # self.__populations = pop
        self.__updateComSize()
        return None

    """
    Other functions
    """

    def infection(self):
        return None 
    
    def __updateComSize(self):
        """Re-total community size across all populations"""
        populations = self.__populations
        comSize = 0

        for i in populations:
            comSize += populations[i].popSize()
        self.__totalComSize = comSize
        return None

    
