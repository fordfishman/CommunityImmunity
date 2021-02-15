## Ford Fishman

"""
Modules
"""
import os; import uuid; import warnings
from datetime import datetime
from subprocess import check_output
import numpy as np; import pandas as pd; import sys
import argparse
import Strain; import Community; import Phage; import PhageReceptor; import Crispr
import general as gen; from network import createNetwork, plotBipartite, adjacencyMatrix, plotNQ
import timeit; import multiprocessing as mp; from functools import partial
from tqdm import tqdm

"""
Setting up arguments
"""
parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', default='comim', type=str, help="Desired name of output path")
parser.add_argument('-t','--timesteps', default=2000, type=int, help="Number of timesteps in each simulation")
parser.add_argument('-s','--single', dest='single', action='store_true', help='Single simulation mode')
parser.add_argument('-M','--multi', dest='single', action='store_false', help='Multi simulation mode')
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
parser.set_defaults(single=True)

arguments = parser.parse_args()
out = arguments.output
timesteps = arguments.timesteps
single_run = arguments.single
# single_run = False
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

params_list =  ["pS", "b", "a", "c", "f","beta", "adsp", "d", "m", "l", "popinit", "phageinit"]
params = set(params_list)


path = os.getcwd()

"""
Parameters
"""

y = 1 # rate of density dependence setting in around a, y = 1 makes classic beverton holt 

# timesteps=5000
# a=1e6
# c=0.01
# f=0
# l=0.9
# m=1e-6
# b=1.2
# pS=1e-6
# d=0.1
# adsp=1e-8
# beta=100
# popinit=1e5
# phageinit=1e7

"""
Functions called by main
"""

def initialize(param_dict:dict=None):
    """
    Set up initial strains and phages
    Return: com (Community)
    """
    from Enums import Type

    from numpy.random import uniform as uni, randint
    
    if not param_dict is None:
        param_dict = dict() # will contain parameters already specified at the command line
        
        for arg in params:
            if arg in globals() and not globals()[arg] is None:
                param_dict[arg] = globals()[arg]

    # for arg,val in globals().items():
    #     if arg in params and not val is None:
    #         param_dict[arg] = val


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
        name = "s0001",
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

    for i in range(8):

        protospacers.add( nameGenerator.generateName(Type.PROTO))

    phage1 = Phage.Phage(
        name = "p0001",
        adsp = adsp,beta = beta, d = d,
        receptor = receptor1,
        pop = phageinit,
        protospacers = protospacers

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
            "id": uuid.uuid1(),
            "pop":np.nan,
            "phage":np.nan,
            "immune":np.nan,
            "susceptible":np.nan,
            "richness":np.nan,
            "phageRichness":np.nan,
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
            "phageinit":phageinit,
            "nodf":np.nan,
            "Q":np.nan,
            }        
    )
    return com 


# def timestep(community):

#     return community

    
"""
Main function for running simulation
"""
def main():

    now = datetime.now()
    print("Sim Start Time:",now.strftime("%Y-%d-%m %H:%M:%S"))
    print()

    unset_params = "Unset Parameters: "
    set_params = dict()
    for arg in params:
        if arg in globals():
            if globals()[arg] is None:
                unset_params += arg + " "
            else:
                set_params[arg] = globals()[arg]
    
    print(unset_params)
    print()
    
    if single_run:
        one_sim(set_params)
    else:
        
        print("Set parameters:")
        for param in params_list:
            if param in set_params:
                print("%s:\t%s" % (param, set_params[param]))

        multi_sim(sims, set_params, unset_params)

    now = datetime.now()
    print("Sim End Time:",now.strftime("%Y-%d-%m %H:%M:%S"),'\n')

    return None

def sim_proc(community, timesteps):
    np.random.seed() # ensures parallel runs have different seeds

    N = community.N_tot # community size
    pRichness = [ community.phageRichness() ] # phage richness over time
    maxPRich = pRichness[0]
    sRichness = [ community.richness() ] # strain richness over time
    maxSRich = sRichness[0]
    cRichness = list() # spacer richness over time (list of lists)


    t0 = timeit.default_timer()

    timestep = community.timestep

    for i in range(timesteps):
        timestep(i)
        pRichness.append( community.phageRichness() )
        maxPRich = max(pRichness)
        sRichness.append( community.richness() )
        maxSRich = max(sRichness)
        cRichness.append( community.spacerRichness() )
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


def one_sim(params):

    """
    Runs a single simulation with given parameters
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    outputMain = out + "/main.csv"
    outputRichness = out + "/richness.csv"
    outputFull = out + "/full.csv"
    outputNetwork = out + "/network.png"
    outputAdjacency = out + "/adjacency.csv"
    networkR = '%s/r/multinetwork.R' % (path)

    community = initialize(params)
    
    print("pS:\t%s"%community.pS)
    print("b:\t%s"%community.strains["s0001"].b)
    print("a:\t%s"%community.strains["s0001"].a)
    print("c:\t%s"%community.strains["s0001"].c)
    print("f:\t%s"%community.strains["s0001"].f)
    print("beta:\t%s"%beta)
    print("adsp:\t%s"%adsp)
    print("d:\t%s"%d)
    print("m:\t%s"%community.m)
    print("l:\t%s"%community.l)
    print()
    print("popint:\t%s"%community.N_tot)
    print("phageinit:\t%s"%community.PList[0])

    N = community.N_tot # community size
    pRichness = [ community.phageRichness() ] # phage richness over time
    maxPRich = pRichness[0]
    sRichness = [ community.richness() ] # strain richness over time
    maxSRich = sRichness[0]
    cRichness = list() # spacer richness over time (list of lists)

    timestep = community.timestep

    print('Running Simulation...')

    for i in tqdm(range(timesteps)): # run for # of timesteps
        
        timestep(i)
        pRichness.append( community.phageRichness() )
        maxPRich = max(pRichness)
        sRichness.append( community.richness() )
        maxSRich = max(sRichness)
        cRichness.append( community.spacerRichness() )

        community.summary["richness"] = maxSRich
        community.summary["phageRichness"] = maxPRich 

    community.summary["pop"] = community.N_tot
    community.summary["phage"] = community.PList[-1]
    community.summary["immune"] = community.IList[-1]
    community.summary["susceptible"] = community.SList[-1]


    net = community.net
    plotBipartite(net, outputNetwork)
    strainIDS, phageIDS = community.strainIDS, community.phageIDS
    strainIDS.sort(); phageIDS.sort()

    A = pd.DataFrame(community.A, columns=phageIDS, index=strainIDS)
    A.to_csv(outputAdjacency)

    if community.A.shape[0] > 1 and community.A.shape[1] > 1:

        output = check_output([networkR, '-f', outputAdjacency]).decode('utf-8')
        items = output.split(" ")
        
        nodf = items[0]
        community.summary['nodf'] = nodf

        if len(items) > 1:
            Q = items[1]
            community.summary['Q'] = Q

    print('Initial conditions:')
    print(community.summary)
    print()

    N = str(community.NList[-1]) # community size
    P = str(community.PList[-1])
    print("hosts:\t%s"%(N))
    print("phages:\t%s"%(P))
    print("Strains:")
    print(sRichness[-1])
    print("Phages:")
    print(pRichness[-1])
    print()

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

    df3 = community.fullDF()

    df3.drop(df3[df3['pop'] <= 0].index, inplace=True)

    df1.to_csv(outputMain)
    df2.to_csv(outputRichness)
    df3.to_csv(outputFull)

    print('Output files:')
    print(outputMain)
    print(outputFull)
    print(outputNetwork)
    print(outputRichness)
    print(outputAdjacency)
    print()

    return None

def map_store(community):
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    temp = '%s/temp/adjacency_%s.csv' % (out, community.summary['id'])
    networkR = '%s/r/multinetwork.R' % (path)
    strainIDS, phageIDS = community.strainIDS, community.phageIDS
    strainIDS.sort(); phageIDS.sort()

    A = pd.DataFrame(community.A, columns=phageIDS, index=strainIDS)
    A.to_csv(temp)
    
    if community.A.shape[0] > 1 and community.A.shape[1] > 1:
        output = check_output([networkR, '-f', temp]).decode('utf-8')
        items = output.split(" ")

        nodf = items[0]
        community.summary['nodf'] = nodf

        if len(items) > 1:
            Q = items[1]
            community.summary['Q'] = Q

    return community.summary

def multi_sim(sims, set_params, unset_params):

    communities = [ initialize(set_params) for i in range(sims) ]

    # output name for run uses provided output name, number of sims, set parameters
    # temp = '%s/temp/adjacency.csv' % (out)
    # networkR = '%s/r/multinetwork.R' % (path)
    output = "%s/summary.csv" % (out) 
    nq_path = "%s/nq.png" % (out) 

    df = pd.DataFrame(None, 
        columns=['id',"pop","phage","immune","susceptible","richness","phageRichness","pS", "b", "a", "c", "f","beta", "adsp", "d", "m", "l", "popinit", "phageinit","nodf","Q"],
    )
    
    timesteps = 5000

    map_proc = partial(sim_proc, timesteps=timesteps)

    pool = mp.Pool(mp.cpu_count())

    tasks = len(communities)

    print("%s Processors\n" % mp.cpu_count() )

    final_communities = pool.map(map_proc, communities, chunksize=1)

    df_list = pool.map(map_store, final_communities)
    pool.close()
    
    df = pd.DataFrame(df_list)

    df.to_csv(output)
    plotNQ(df, nq_path, unset=unset_params)

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


