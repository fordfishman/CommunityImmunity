## Ford Fishman

import Spacer
import numpy as np

##########################################################################################################

class Crispr():

    """
    Crispr
    Represents a CRISPR locus
    -------------------------------------
    Attributes:
    spacers (set) - sequences to which the spacer corresponds; might change if order matters eventually
    spacerLength (int) - all spacer of this CRISPR have this defined length
    hasCost (bool) - does this unit have a current cost (can change scenarios with active cost)
    isActive (bool) - is the crispr system on or off
    """
    def __init__(self, spacers:set=None, spacerLength:int = 20, isActive:bool=True):
        if not spacers is None:
            self.__spacers = spacers
        else:
            self.__spacers = set()
        self.__spacerLength = spacerLength
        self.__isActive = isActive
        self.__updateCost()
        
##########################################################################################################

    """
    Attribute Functions
    """

    def hasSpacer(self, genome:str):
        """Only checks for single match"""
        match = False
        seq = ''

        if len(self.__spacers) != 0: # if system has spacers
            for i in range(0, len(genome) - self.__spacerLength): # check entire genome for spacers
                seq = genome[i:i + self.__spacerLength] 

                if seq in self.__spacers:
                    match = True
                    break

        return match
    
    def hasCost(self):
        return self.__hasCost

    def isActive(self):
        return (self.__isActive)

##########################################################################################################

    """
    Other functions
    """

    def newSpacer(self, genome:str):
        """Make a new spacer based on the incoming DNA"""
        inds = [ i for i in range( 0, len(genome) - self.__spacerLength ) ] # all possible starting spots for a new spacer
        i = int(np.random.choice(inds))
        self.__spacers.add(genome[i : (i + self.__spacerLength)])

        self.__updateCost() # check to see if there is a cost to the CRISPR now

        return None

    def removeSpacer(self):
        """Remove a spacer randomly"""
        self.__spacers.pop()

        self.__updateCost() # check to see if there is a cost to the CRISPR now

        return None
    
    def activate(self):
        """Turn system on"""
        self.__isActive = True

        self.__updateCost() # check to see if there is a cost to the CRISPR now

        return None

    def deactivate(self):
        """Turn system off"""
        self.__isActive = False

        self.__updateCost() # check to see if there is a cost to the CRISPR now

        return None    

##########################################################################################################

    """
    Private
    Methods
    """

    def __updateCost(self):
        """Has cost if has spacer and is active"""
        self.__hasCost = self.__isActive and len(self.__spacers)>0  
        return None