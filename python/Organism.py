## Organism.py
## Ford Fishman

import numpy as np

##########################################################################################################

class Organism():

    """
    Organism
    -----------------
    name (str)
    pop (float)
    """

    def __init__(self, name:str, pop:float, type_:str,evoTraits:list=None):

        self.name = name
        self.pop = pop
        self.evoTraits = evoTraits
        self.type = type_
