## Ford Fishman

from enum import Enum

class Type(Enum):
    """
    Enumerated class
    Used for generating names
    """
    RECEPTOR = "r"
    STRAIN = "s"
    SPACER = "sp"
    POPULATION = "pop"
    COMMUNITY = "c"
    PHAGE = "p"
    PROTO = "proto"


class Mutation(Enum):
    """
    Different types of mutations
    """
    SNP = "snp"
    DELETION = "del"
    PROTOCHANGE = "change"
    PROTOADD = "add"
    PROTODELETE = "del"


    


