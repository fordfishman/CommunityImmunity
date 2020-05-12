## Ford Fishman

class PhageReceptor():
    """
    PhageReceptor
    Protein receptors on cell membranes to which phage bind to insert DNA
    ---------------------------
    Attributes:
    active (bool): 
    function (float): relative fitness of this specific receptor
    """
    def __init__(self, name:str, fitness:float=1.0, expressed:bool=True):
        self.isExpressed = expressed
        self.name = name
        self.fitness = fitness

    """
    Attribute functions
    """

    # def isExpressed(self):
    #     return self.expressed
    
    # def name(self):
    #     return self.name

    # def fitness(self):
    #     return self.__fitness