## Ford Fishman

"""
Modules
"""
# need to see if spacers are actually getting added or are just replacing old spacer set

import numpy as np; import pandas as pd; import sys
import argparse
import Strain; import Population; import Community; import Phage; import PhageReceptor; import Crispr
import general as gen
import timeit; import multiprocessing as mp; from functools import partial
from progressbar import progressbar

"""
Setting up arguments
"""
parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', default='comim', type=str, help="Desired name of output path")
parser.add_argument('-t','--timesteps', default=2000, type=int, help="Number of timesteps in each simulation")
parser.add_argument('-s','--single', default=False, type=bool, help='Whether or not to run the program in single simulation mode')
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

"""
Functions called by main
"""

def initialize():
    """
    Set up initial strains and phages
    Return: com (Community)
    """
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
        name = "s1",
        a=a,b=b,c=c,y=y,f=f,
        crispr = crispr0,
        phReceptors = {
            receptor1.name:receptor1
            },
        pop = popinit
    )

    phage1 = Phage.Phage(
        name = "p1",
        adsp = adsp,beta = beta, d = d,
        receptor = receptor1,
        pop = phageinit,
        genomeLength=1000
    )

    # spacer = strain1.crispr.makeSpacer("AGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGT")
    # strain1.crispr.addSpacer(spacer)

    # pop1 = Population.Population(
    #     name = "pop1",
    #     strains = {
    #         strain1.name: strain1
    #     } 
    # )

    com = Community.Community(
        pS=pS,m=m,l=l,
        strains = {
            strain1.name: strain1
        },
        phages = {
            phage1.name: phage1
        }
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
    
    for i in range(timesteps):
        timestep(i)
        pRichness.append( len( community.phagesPopDict() ) )
        maxPRich = max(pRichness)
        sRichness.append( community.richness() )
        maxSRich = max(sRichness)
        t1 = timeit.default_timer()

        if (t1 - t0 ) > 30*60: # if the simulation takes longer than 5 min, end it
            break

    community.summary["richness"] = maxSRich
    community.summary["phageRichness"] = maxPRich 
    community.summary["pop"] = community.totalComSize
    community.summary["phage"] = community.phagePopOverTime[-1]
    community.summary["immune"] = community.imOverTime[-1]
    community.summary["susceptible"] = community.susOverTime[-1]

    return community


def one_sim():

    """
    Runs a single simulation with given parameters
    """
    outputMain = out + "/main.csv"
    outputRichness = out + "/richness.csv"

    community = initialize()
    
    # params = {"pS", "b", "a", "c", "beta", "adsp", "d", "m", "l", "popinit", "phageinit"}
    print("pS:\t%s"%community.pS)
    print("b:\t%s"%community.strains["s1"].b)
    print("a:\t%s"%community.strains["s1"].a)
    print("c:\t%s"%community.strains["s1"].c)
    print("f:\t%s"%community.strains["s1"].f)
    print("beta:\t%s"%community.phages["p1"].beta)
    print("adsp:\t%s"%community.phages["p1"].adsp)
    print("d:\t%s"%community.phages["p1"].d)
    print("m:\t%s"%community.m)
    print("l:\t%s"%community.l)
    print()
    print("popint:\t%s"%community.totalComSize)
    print("phageinit:\t%s"%community.phagePopOverTime[0])




    N = community.totalComSize # community size
    pRichness = [len( community.phagesPopDict() )] # phage richness over time
    maxPRich = pRichness[0]
    sRichness = [len( community.strains )] # strain richness over time
    maxSRich = sRichness[0]
    cRichness = list() # spacer richness over time (list of lists)

    timestep = community.timestep

    for i in progressbar(range(timesteps)): # run for # of timesteps
        
        timestep(i)
        pRichness.append( len( community.phagesPopDict() ) )
        maxPRich = max(pRichness)
        sRichness.append( community.richness() )
        maxSRich = max(sRichness)
        cRichness.append( community.spacerRichness() )

        community.summary["richness"] = maxSRich
        community.summary["phageRichness"] = maxPRich 

        # if i == 600 and community.totalComSize != 0:
        #     s1 = community.strains["s1"]
        #     spacer = s1.crispr.makeSpacer(community.phages["p1"].genome)
        #     receptor = s1.phReceptors["r1"]
        #     s2 = Strain.Strain(
        #         name = "s2",
        #         a=s1.a,b=s1.b,c=s1.c,y=y,
        #         crispr = Crispr.Crispr(),
        #         phReceptors = {
        #             receptor.name:receptor
        #             },
        #         pop = 1
        #     )

        #     s2.crispr.addSpacer(spacer)

        #     community.strains["s2"] = s2

    community.summary["pop"] = community.totalComSize
    community.summary["phage"] = community.phagePopOverTime[-1]
    community.summary["immune"] = community.imOverTime[-1]
    community.summary["susceptible"] = community.susOverTime[-1]

    print(community.summary)

    N = str(community.comSizeOverTime[-1]) # community size
    P = str(community.phagePopOverTime[-1])
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
    # print("max immune: %s"%(max(community.imOverTime)))
    # df1 = pd.DataFrame( 
    #     list( 
    #         zip(
    #             community.comSizeOverTime, 
    #             community.phagePopOverTime,
    #             community.imOverTime,
    #             community.susOverTime, 
    #             range(1,timesteps+1), 
    #             ) 
    #         ),
    #     columns= ['host','phage','immune','susceptible','time'],
    #     )
    df1 = community.fullRecord()
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
    # for phage in community.phages.values:
    #     print(phage.genome)
    return None

def multi_sim(sims):

    params = {"pS", "b", "a", "c", "f","beta", "adsp", "d", "m", "l", "popinit", "phageinit"}
    # constant_params = sorted(["_%s%s" % (param,globals()[param]) for param in params if not globals()[param] is None]) 

    communities = [ initialize() for i in range(sims) ]

    # output name for run uses provided output name, number of sims, set parameters

    # output = "%s%ssims%s.csv" % (out, sims, "".join(constant_params)) 
    output = "%s/summary.csv" % (out) 

    # timestep = [ community.timestep for community in communities ] # list of all timestep functions for all sims
    # print("initial sizes")
    # for community in communities:
    #     print(community.totalComSize)
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


