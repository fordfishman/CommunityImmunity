## Ford Fishman

"""
Modules
"""
import numpy as np
import Strain; import Population; import Community 
import Phage; import PhageReceptor; import Crispr


"""
Parameters
"""

timesteps = 10
pS = 0.001 # prob of spacer forming if infection occurs
b = 1.5 # initial max intrinsic growth
a = 0.01 # competition coefficient (inter and intraspecific)
crisprCost = 0.1 # fitness cost of having active CRISPR system


"""
Functions called by main
"""

def initialize():
    """
    Set up initial strains and phages
    Return: com (Community)
    """
    receptor1 = PhageReceptor.PhageReceptor( name = "r" + "1" ) # change numbering system

    crispr0 = Crispr.Crispr()

    strain1 = Strain.Strain(
        name = "s" + "1",
        crispr = crispr0,
        phReceptors = {
            receptor1.name():receptor1
            }
    )

    pop1 = Population.Population(
        name = "pop" + "1",
        strains = {
            strain1.name(): strain1
        } 
    )

    phage1 = Phage.Phage(
        name = "p" + "1",
        phageReceptor = receptor1
    )

    com = Community.Community(
        k=1000,
        populations = {
            pop1.name(): pop1
        },
        phages = {
            phage1.name(): phage1
        }
    )
    return com 


def timestep(community):

    return community


def generateName(type:str, num):
    """Make names for different objects of different classes"""

    return None
    
"""
Main function for running simulation
"""
def main():

    community = initialize()
    N = community.totalComSize() # community size
    for i in range(timesteps): # run for # of timesteps

        N = community.totalComSize() # community size
        # community.getPopulation("pop1").getStrain("s1").timestep(N,a,b,crisprCost)
        community.timestep(N, a, b, crisprCost)
    
    print(N)
    # print("strain nam")
    # print(community.getPopulation("pop1").getStrain("s1").name())
    return None


"""
Execute main function
"""
if __name__ == "__main__":
    main()







# import numpy as np
# import Strain; import Population; import Community 
# import Phage; import PhageReceptor; import Crispr


# """
# Parameters
# """

# timesteps = 10
# pS = 0.001 # prob of spacer forming if infection occurs
# b = 1.25 # initial max intrinsic growth
# a = 0.5 # competition coefficient (inter and intraspecific)
# crisprCost = 0.1 # fitness cost of having active CRISPR system


# """
# Functions called by main
# """

# def initialize():
#     """
#     Set up initial strains and phages
#     Return: com (Community)
#     """
#     receptor1 = PhageReceptor.PhageReceptor( name = "r" + "1" ) # change numbering system

#     crispr0 = Crispr.Crispr()

#     strain1 = Strain.Strain(
#         name = "s" + "1",
#         crispr = crispr0,
#         phReceptors = {
#             receptor1.name():receptor1
#             }
#     )

#     pop1 = Population.Population(
#         name = "pop" + "1",
#         strains = {
#             strain1.name(): strain1
#         } 
#     )

#     phage1 = Phage.Phage(
#         name = "p" + "1",
#         phageReceptor = receptor1
#     )

#     com = Community.Community(
#         k=1000,
#         populations = {
#             pop1.name(): pop1
#         },
#         phages = {
#             phage1.name(): phage1
#         }
#     )
#     return com 



# community = initialize()
# N = community.totalComSize() # community size
# print(community.)
# for i in range(timesteps): # run for # of timesteps

#     N = community.totalComSize() # community size

#     community.getPopulation("pop1").getStrain("s1").timestep(N,a,b,crisprCost)
#     # community.timestep(N, a, b, crisprCost)

# print(N)
# # print("strain nam")
# # print(commu