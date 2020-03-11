## Ford Fishman

class Spacer():

    """
    Spacer
    Represents a spacer in a CRISPR locus
    -------------------------------------
    Attributes:
    seq (str) - sequences to which the spacer corresponds
    """

    def __init__(self, seq:str):
        self.__seq = seq

    """
    Attribute functions
    """

    def seq(self):
        return self.__seq
