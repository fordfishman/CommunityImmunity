## Ford Fishman

"""
Modules
"""
# need to see if spacers are actually getting added or are just replacing old spacer set

import numpy as np; import pandas as pd; import sys
import argparse
import Strain; import Population; import Community; import Phage; import PhageReceptor; import Crispr
import general as gen; from network import createNetwork, plotBipartite, adjacencyMatrix
import timeit; import multiprocessing as mp; from functools import partial
from tqdm import tqdm

"""
Setting up arguments
"""
parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', default='comim', type=str, help="Desired name of output path")
parser.add_argument('-t','--timesteps', default=2000, type=int, help="Number of timesteps in each simulation")
parser.add_argument('-s','--single', default=True, type=bool, help='Whether or not to run the program in single simulation mode')
parser.add_argument('-S','--sims', default=100, type=int, help="Number of simulations to run (--single must be False)")
parser.add_argument('-pS', default=None, type=float, help="probability of spacer formation per infection")
parser.add_argument('-m', default=None, type=float, help="phage mutation rate per nt")
parser.add_argument('-b', default=None, type=float, help="host per capita growth rate")
parser.add_argument('-a', default=None, type=float, help="region where density dependence sets in for hosts")
parser.add_argument('-c', default=None, type=float, help="cost of CRISPR")
parser.add_argument('-f', default=None, type=float, help="CRISPR failure rate")
parser.add_argument('--adsp', default=None, type=float, help="adsorption rate of phage")
parser.add_argument('--beta', default=None, type=int, help="burst size of phage")
parser.add_argument('-d', default=None, type=float, help="phage decay rate per timestep")
parser.add_argument('-l', default=None, type=float, help="proportion of infections that lead to bursting each timestep")
parser.add_argument('--popinit',default=None, type=float, help="initial host population")
parser.add_argument('--phageinit',default=None, type=float, help="initial phage population")

arguments = parser.parse_args()
out = arguments.output
timesteps = arguments.timesteps
single_run = arguments.single
# single_run = True
sims = arguments.sims
pS = arguments.pS
m = arguments.m
b = arguments.b
a = arguments.a
c = arguments.c
f = arguments.f
adsp = arguments.adsp
beta = arguments.beta
d = arguments.d
l = arguments.l
popinit = arguments.popinit
phageinit = arguments.phageinit


"""
Parameters
"""


# timesteps = 1000
# host parameters
# pS = 10**-6 # prob of spacer forming if infection occurs per host
# b = 1.2 # initial max intrinsic growth of hosts
# a = 10**6 # region where density dependence sets in for hosts
y = 1 # rate of density dependence setting in around a, y = 1 makes classic beverton holt 
# crisprCost = 0.1 # fitness cost of having active CRISPR system

# phage parameters
# beta = 10 # burst size of phage
# adsp = 10**-7 # adsorption rate of phage
# d = 0.1 # natural decay rate of phage
# m = 10**-6 # phage mutation rate per nt
# l = 0.5 # proportion of infections that lead to bursting each timestep

timesteps=5000
a=1e6
c=0.01
f=0
l=0.9
m=1e-6
b=1.2
pS=1e-6
d=0.1
adsp=1e-8
beta=100
popinit=1e5
phageinit=1e7

"""
Functions called by main
"""

def initialize():
    """
    Set up initial strains and phages
    Return: com (Community)
    """
    from Enums import Type

    from numpy.random import uniform as uni, randint
    params = {"pS", "b", "a", "c", "f","beta", "adsp", "d", "m", "l", "popinit", "phageinit"}
    param_dict = dict() # will contain parameters already specified at the command line
    
    for arg,val in globals().items():
        if arg in params and not val is None:
            param_dict[arg] = val


    # if parameters aren't specified, draw them from distributions
    pS = param_dict.get( "pS", 10**-( uni(6,9) ) ) # default: random float 10^-5 - 10^-9
    b = param_dict.get( "b", uni(0.9, 2) )
    a = param_dict.get( "a", 10**( uni(5,8) ) )
    c = param_dict.get( "c", 10**-( uni(1,3) ) )
    f = param_dict.get( "f", 10**-( uni(5,7) ) )
    
    beta = param_dict.get( "beta", randint(50,200) ) # default: random int from 1-200
    adsp = param_dict.get( "adsp", 10**-( uni(7,9) ) )
    d = param_dict.get( "d", uni(0.0, 0.3) )
    m = param_dict.get( "m", 10**-( uni(6,9) ) )
    l = param_dict.get( "l", uni(0.0, 1.0) )

    popinit = param_dict.get( "popinit", 10**( uni(5,7) ) )
    phageinit = param_dict.get( "phageinit", 10**( uni(5,7) ))

    receptor1 = PhageReceptor.PhageReceptor( name = "r1" ) # change numbering system

    crispr0 = Crispr.Crispr()


    strain1 = Strain.Strain(
        name = "s001",
        a=a,b=b,c=c,y=y,f=f,
        crispr = crispr0,
        phReceptors = {
            receptor1.name:receptor1
            },
        pop = popinit
    )

    # strain2 = Strain.Strain(
    #     name = "s2",
    #     a=a,b=b,c=c,y=y,f=f,
    #     crispr = crispr0,
    #     phReceptors = {
    #         receptor1.name:receptor1
    #         },
    #     pop = popinit/100
    # )

    nameGenerator = gen.NameGenerator()

    protospacers = set()

    for i in range(4):

        protospacers.add( nameGenerator.generateName(Type.PROTO))

    phage1 = Phage.Phage(
        name = "p001",
        adsp = adsp,beta = beta, d = d,
        receptor = receptor1,
        pop = phageinit,
        protospacers = protospacers
        # numProto=4

    )

    # spacer = strain2.crispr.makeSpacer("AGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGT")
    # spacer = strain2.crispr.makeSpacer(phage1.genome)

    # strain2.crispr.addSpacer(spacer)

    # pop1 = Population.Population(
    #     name = "pop1",
    #     strains = {
    #         strain1.name: strain1
    #     } 
    # )

    com = Community.Community(
        c=c,pS=pS,m=m,l=l,
        strains = {
            strain1.name: strain1,
            # strain2.name: strain2
        },
        phages = {
            phage1.name: phage1
        },
        nameGenerator=nameGenerator
    )

    com.summary = pd.Series(
        data = {
            "pop":None,
            "phage":None,
            "immune":None,
            "susceptible":None,
            "richness":None,
            "phageRichness":None,
            "pS":pS,
            "b":b,
            "a":a, 
            "c":c,
            "f":f, 
            "beta":beta, 
            "adsp":adsp, 
            "d":d, 
            "m":m, 
            "l":l, 
            "popinit":popinit, 
            "phageinit":phageinit
            }        
    )
    return com 


# def timestep(community):

#     return community

    
"""
Main function for running simulation
"""
def main():
    if single_run:
        one_sim()
    else:
        multi_sim(sims)
    # one_sim()
    # one_sim(m=0, pS=0)
    # multi_sim(1000, m=0)
    # multi_sim(1000, m = 10**-7, pS=10**-5, beta=100,adsp=10**-8, d = 0.1, c=0)
    return None

def sim_proc(community, timesteps):
    pRichness = [len( community.phagesPopDict() )] # phage richness over time
    maxPRich = pRichness[0]
    sRichness = [len( community.strains )] # strain richness over time
    maxSRich = sRichness[0]

    t0 = timeit.default_timer()

    timestep = community.timestep
    
    nameGenerator = gen.NameGenerator()

    for i in range(timesteps):
        timestep(i, nameGenerator=nameGenerator)
        pRichness.append( len( community.phagesPopDict() ) )
        maxPRich = max(pRichness)
        sRichness.append( community.richness() )
        maxSRich = max(sRichness)
        t1 = timeit.default_timer()

        if (t1 - t0 ) > 30*60: # if the simulation takes longer than 5 min, end it
            break

    community.summary["richness"] = maxSRich
    community.summary["phageRichness"] = maxPRich 
    community.summary["pop"] = community.N_tot
    community.summary["phage"] = community.PList[-1]
    community.summary["immune"] = community.IList[-1]
    community.summary["susceptible"] = community.SList[-1]

    return community


def one_sim():

    """
    Runs a single simulation with given parameters
    """
    outputMain = out + "/main.csv"
    outputRichness = out + "/richness.csv"
    outputNetwork = out + "/network.png"
    outputAdjacency = out + "/adjacency.csv"

    community = initialize()
    
    print("pS:\t%s"%community.pS)
    print("b:\t%s"%community.strains["s001"].b)
    print("a:\t%s"%community.strains["s001"].a)
    print("c:\t%s"%community.strains["s001"].c)
    print("f:\t%s"%community.strains["s001"].f)
    print("beta:\t%s"%beta)
    print("adsp:\t%s"%adsp)
    print("d:\t%s"%d)
    print("m:\t%s"%community.m)
    print("l:\t%s"%community.l)
    print()
    print("popint:\t%s"%community.N_tot)
    print("phageinit:\t%s"%community.PList[0])

    N = community.N_tot # community size
    pRichness = [len( community.phagesPopDict() )] # phage richness over time
    maxPRich = pRichness[0]
    sRichness = [len( community.strains )] # strain richness over time
    maxSRich = sRichness[0]
    cRichness = list() # spacer richness over time (list of lists)

    timestep = community.timestep2

    print('Running Simulation...')

    for i in tqdm(range(timesteps)): # run for # of timesteps
        
        timestep(i)
        pRichness.append( len( community.phagesPopDict() ) )
        maxPRich = max(pRichness)
        sRichness.append( community.richness() )
        maxSRich = max(sRichness)
        cRichness.append( community.spacerRichness() )

        community.summary["richness"] = maxSRich
        community.summary["phageRichness"] = maxPRich 

        # if i == 50 and community.N != 0:
        #     s1 = community.strains["s1"]
        #     crispr0 = Crispr.Crispr()

        #     receptor = s1.phReceptors["r1"]
        #     s2 = Strain.Strain(
        #         name = "s2",
        #         a=s1.a,b=s1.b,c=s1.c,y=y,f=f,
        #         crispr = crispr0,
        #         phReceptors = {
        #             receptor.name:receptor
        #             },
        #         pop = 1e0
        #     )
        #     spacer = s2.crispr.makeSpacer(community.phages["p1"].genome)
        #     s2.crispr.addSpacer(spacer)

        #     community.strains["s2"] = s2

    community.summary["pop"] = community.N_tot
    community.summary["phage"] = community.PList[-1]
    community.summary["immune"] = community.IList[-1]
    community.summary["susceptible"] = community.SList[-1]
    
    print('Initial conditions:')
    print(community.summary)

    # edges = community.globalInfectionEdges()
    # B = createNetwork(edges, community.trimmedStrains.keys(), community.trimmedPhages.keys())
    # plotBipartite(B, outputNetwork)
    # A = adjacencyMatrix(B)
    # A.to_csv(outputAdjacency)

    N = str(community.NList[-1]) # community size
    P = str(community.PList[-1])
    print("hosts:\t%s"%(N))
    print("phages:\t%s"%(P))
    print("Strains:")
    print(sRichness[-1])
    print("Phages:")
    print(pRichness[-1])
    # print("Spacers:")
    # print(cRichness)
    ## times
    # import statistics as stat
    # print("Strain times:")
    # print("max: %s\tmean: %s" % (max(community.strainTimes), stat.mean(community.strainTimes)))
    # print("Phage times:")
    # print("max: %s\tmean: %s" % (max(community.phageTimes), stat.mean(community.phageTimes)))
    # print("df times:")
    # print("max: %s\tmean: %s" % (max(community.dfTimes), stat.mean(community.dfTimes)))
    # print("Other times:")
    # print("max: %s\tmean: %s" % (max(community.otherTimes), stat.mean(community.otherTimes)))
    # print()
    # print("max immune: %s"%(max(community.IList)))

    # df1 = community.fullRecord()
    df1 = pd.DataFrame(
        list(
            zip(
                community.NList,
                community.PList,
                community.SList,
                community.IList,
                range(1,timesteps+1)
            )
        ),
        columns = ['N','P','S','I','t'],
    )
    # # cols = [ "strain"+str(i) for i in range(0,len(cRichness)) ]
    df2 = pd.DataFrame( 
        list( 
            zip(
                sRichness,
                pRichness,
                range(1,timesteps+1)
             ) 
            ),
        columns = ['HostRichness','PhageRichness','time'],
        )

    df1.to_csv(outputMain)
    df2.to_csv(outputRichness)

    print('Output files:')
    print(outputMain)
    # print(outputNetwork)
    print(outputRichness)
    # print(outputAdjacency)

    return None

def multi_sim(sims):

    # params = {"pS", "b", "a", "c", "f","beta", "adsp", "d", "m", "l", "popinit", "phageinit"}
    # constant_params = sorted(["_%s%s" % (param,globals()[param]) for param in params if not globals()[param] is None]) 

    communities = [ initialize() for i in range(sims) ]

    # output name for run uses provided output name, number of sims, set parameters

    # output = "%s%ssims%s.csv" % (out, sims, "".join(constant_params)) 
    output = "%s/summary.csv" % (out) 

    # timestep = [ community.timestep for community in communities ] # list of all timestep functions for all sims
    # print("initial sizes")
    # for community in communities:
    #     print(community.N)
    df = pd.DataFrame(None, 
        columns=["pop","phage","immune","susceptible","richness","phageRichness","pS", "b", "a", "c", "f","beta", "adsp", "d", "m", "l", "popinit", "phageinit"],
    )
    
    timesteps = 1000
    # times_com = list()
    # times_timestep = list()
    # times_remaining = list()
    # times_ratio = list()

    map_proc = partial(sim_proc, timesteps=timesteps)

    pool = mp.Pool(mp.cpu_count())

    print("%s Processors\n" % mp.cpu_count() )

    final_communities = pool.map(map_proc, communities, chunksize=1)

    pool.close()

    for i,community in enumerate(final_communities):
        df.loc[i] = df.columns.map( community.summary )

    df.to_csv(output)
    


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


