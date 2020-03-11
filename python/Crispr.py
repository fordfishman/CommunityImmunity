## Ford Fishman

import Spacer
import numpy as np

class Crispr():
    """
    Crispr
    Represents a CRISPR locus
    -------------------------------------
    Attributes:
    spacers (set) - sequences to which the spacer corresponds; might change if order matters eventually
    """
    def __init__(self, spacers:set=None, spacerLength:int = 20):
        self.__spacers = spacers
        self.__spacerLength = spacerLength

    """
    Attribute Functions
    """

    def hasSpacer(self, genome:str):
        """Only checks for single match"""
        match = False
        seq = ''

        for i in range(0, len(genome) - self.__spacerLength):
            seq = genome[i,i + self.__spacerLength] 

            if seq in self.__spacers:
                match = True
                break

        return match

    """
    Other functions
    """

    def newSpacer(self, genome:str):
        """Make a new spacer based on the incoming DNA"""
        inds = [ i for i in range( 0, len(genome) - self.__spacerLength ) ] # all possible starting spots for a new spacer
        i = np.random.choice(inds)
        self.__spacers.add(genome[i : i + self.__spacerLength])

        return None

