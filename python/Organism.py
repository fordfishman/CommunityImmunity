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

    def __init__(self, name:str, pop:float):

        self.name = name
        self.pop = pop