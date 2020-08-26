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
    def __init__(self, spacers:set=None, spacerLength:int = 30, isActive:bool=True):
        if not spacers is None:
            self.spacers = spacers
        else:
            self.spacers = set()
        self.spacerLength = spacerLength
        self.isActive = isActive
        self.__updateCost()
        
##########################################################################################################

    """
    Attribute Functions
    """
    def __len__(self):
        return len(self.spacers)

    def hasSpacer(self, genome:str):
        """Only checks for single match"""
        match = False
        seq = ''

        if len(self.spacers) != 0: # if system has spacers
            for i in range(0, len(genome) - self.spacerLength): # check entire genome for spacers
                seq = genome[i:i + self.spacerLength] 

                if seq in self.spacers:
                    match = True
                    break

        return match
    

##########################################################################################################

    """
    Other functions
    """

    def makeSpacer(self, genome:str)->str:
        """Make a new spacer based on the incoming DNA"""
        inds = [ i for i in range( 0, len(genome) - self.spacerLength ) ] # all possible starting spots for a new spacer
        i = int(np.random.choice(inds))
        return genome[i : (i + self.spacerLength)]

    def addSpacer(self, spacer:str):
        
        self.spacers.add(spacer)

        self.__updateCost() # check to see if there is a cost to the CRISPR now

        return None

    def removeSpacer(self):
        """Remove a spacer randomly"""
        self.spacers.pop()

        self.__updateCost() # check to see if there is a cost to the CRISPR now

        return None
    
    def activate(self):
        """Turn system on"""
        self.isActive = True

        self.__updateCost() # check to see if there is a cost to the CRISPR now

        return None

    def deactivate(self):
        """Turn system off"""
        self.isActive = False

        self.__updateCost() # check to see if there is a cost to the CRISPR now

        return None    

##########################################################################################################

    """
    Private
    Methods
    """

    def __updateCost(self):
        """Has cost if has spacer and is active"""
        self.hasCost = self.isActive and len(self.spacers)>0  
        return None