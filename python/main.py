## Ford Fishman

"""
Modules
"""

import numpy as np; import pandas as pd
import argparse
import Strain; import Population; import Community; import Phage; import PhageReceptor; import Crispr
import general as gen

"""
Setting up arguments
"""
parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', default = 'comim', help = "Desired name of output path")
arguments = parser.parse_args()
out = arguments.output

"""
Parameters
"""
outputMain = out + "_main.csv"
outputRichness = out + "_richness.csv"

timesteps = 1000
# host parameters
pS = 10**-6 # prob of spacer forming if infection occurs per host
bH = 1.2 # initial max intrinsic growth of hosts
aH = 10**7 # region where density dependence sets in for hosts
y = 1 # rate of density dependence setting in around a, y = 1 makes classic beverton holt 
crisprCost = 0.1 # fitness cost of having active CRISPR system

# phage parameters
bP = 10 # burst size of phage
absP = 10**-7 # absorbtion rate of phage
dP = 0.1 # natural decay rate of phage
m = 10**-5

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
            },
        pop = 100000
    )

    phage1 = Phage.Phage(
        name = "p" + "1",
        receptor = receptor1,
        pop = 100000,
        genomeLength=10**3
    )


    pop1 = Population.Population(
        name = "pop" + "1",
        strains = {
            strain1.name(): strain1
        } 
    )



    com = Community.Community(
        populations = {
            pop1.name(): pop1
        },
        phages = {
            phage1.name(): phage1
        }
    )
    return com 


# def timestep(community):

#     return community

    
"""
Main function for running simulation
"""
def main():

    """
    Implement __add__, __str__, and __repr__ at some point
    """

    community = initialize()
    N = community.totalComSize() # community size
    pRichness = [len( community.phagesPopDict() )] # phage richness over time
    sRichness = [1] # strain richness over time
    cRichness = list() # spacer richness over time (list of lists)


    for i in range(timesteps): # run for # of timesteps
        
        community.timestep(aH=aH,bH=bH, c=crisprCost, y=y, bP=bP, absP=absP, dP=dP, pS=pS, m=m)
        pRichness.append( len( community.phagesPopDict() ) )
        sRichness.append( community.richness() )
        cRichness.append( community.spacerRichness() )


    N = str(community.comSizeOverTime()[-1]) # community size
    P = str(community.phagePopOverTime()[-1])
    print("hosts:\t"+N)
    print("phages:\t"+P)
    print("Strains:")
    print(sRichness[-1])
    print("Phages:")
    print(pRichness[-1])
    # print("Spacers:")
    # print(cRichness)

    df1 = pd.DataFrame( 
        list( 
            zip(
                community.comSizeOverTime(), 
                community.phagePopOverTime(),
                community.resOverTime(),
                community.vulnOverTime(), 
                range(1,timesteps+1), 
                ) 
            ),
        columns= ['host','phage','resistant','vulnerable','time'],
        )
    # cols = [ "strain"+str(i) for i in range(0,len(cRichness)) ]
    df2 = pd.DataFrame( 
        list( 
            zip(
                sRichness,
                pRichness,
                range(1,timesteps+1)
             ) 
            ),
        columns= ['HostRichness','PhageRichness','time'],
        )

    df1.to_csv(outputMain)
    df2.to_csv(outputRichness)
    # for phage in community.phages().values():
    #     print(phage.genome())
    return None

"""
Execute main function
"""
if __name__ == "__main__":
    main()

##########################################################################################################

"""
Print tests
"""


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