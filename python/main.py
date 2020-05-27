## Ford Fishman

"""
Modules
"""

import numpy as np; import pandas as pd; import sys
import argparse
import Strain; import Population; import Community; import Phage; import PhageReceptor; import Crispr
import general as gen

"""
Setting up arguments
"""
parser = argparse.ArgumentParser()
# parser.add_argument('-pS', default=None)
parser.add_argument('-o','--output', default = 'comim', help = "Desired name of output path")
arguments = parser.parse_args()
out = arguments.output

"""
Parameters
"""


timesteps = 1000
# host parameters
# pS = 10**-6 # prob of spacer forming if infection occurs per host
# b = 1.2 # initial max intrinsic growth of hosts
# a = 10**6 # region where density dependence sets in for hosts
y = 1 # rate of density dependence setting in around a, y = 1 makes classic beverton holt 
# crisprCost = 0.1 # fitness cost of having active CRISPR system

# phage parameters
# beta = 10 # burst size of phage
# absp = 10**-7 # absorbtion rate of phage
# d = 0.1 # natural decay rate of phage
# m = 10**-6 # phage mutation rate per nt
# l = 0.5 # proportion of infections that lead to bursting each timestep

"""
Functions called by main
"""

def initialize(**kw):
    """
    Set up initial strains and phages
    Return: com (Community)
    """
    from numpy.random import uniform as uni, randint
    params = {"pS", "b", "a", "c", "beta", "absp", "d", "m", "l", "popinit", "phageinit"}
    
    for arg in kw:
        if arg not in params:
            print("Unknown parameter passed")
            sys.exit()

    # if parameters aren't specified, draw them from distributions
    pS = kw.get( "pS", 10**-( uni(6,9) ) ) # default: random float 10^-5 - 10^-9
    b = kw.get( "b", uni(0.9, 2) )
    a = kw.get("a", 10**( uni(5,8) ) )
    c = kw.get( "c", 10**-( uni(0,3) ) )
    
    beta = kw.get( "beta", randint(1,200) ) # default: random int from 1-200
    absp = kw.get( "absp", 10**-( uni(7,9) ) )
    d = kw.get( "d", uni(0.0, 0.3) )
    m = kw.get( "m", 10**-( uni(6,9) ) )
    l = kw.get( "l", uni(0.0, 1.0) )

    popinit = kw.get( "popinit", 10**( uni(5,7) ) )
    phageinit = kw.get( "phageinit", 10**( uni(5,7) ))

    receptor1 = PhageReceptor.PhageReceptor( name = "r1" ) # change numbering system

    crispr0 = Crispr.Crispr()

    strain1 = Strain.Strain(
        name = "s1",
        a=a,b=b,c=c,y=y,
        crispr = crispr0,
        phReceptors = {
            receptor1.name:receptor1
            },
        pop = popinit
    )

    phage1 = Phage.Phage(
        name = "p1",
        absp = absp,beta = beta, d = d,
        receptor = receptor1,
        pop = phageinit,
        genomeLength=10**3
    )


    pop1 = Population.Population(
        name = "pop1",
        strains = {
            strain1.name: strain1
        } 
    )



    com = Community.Community(
        pS=pS,m=m,l=l,
        populations = {
            pop1.name: pop1
        },
        phages = {
            phage1.name: phage1
        }
    )

    com.summary = pd.Series(
        data = {
            "pop":None,
            "phage":None,
            "resistant":None,
            "vulnerable":None,
            "richness":None,
            "phageRichness":None,
            "pS":pS,
            "b":b,
            "a":a, 
            "c":c, 
            "beta":beta, 
            "absp":absp, 
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
    # one_sim()
    # one_sim(m=0, pS=0)
    multi_sim(1)
    return None

def one_sim(**kw):

    """
    Runs a single simulation with given parameters
    """
    outputMain = out + "_main.csv"
    outputRichness = out + "_richness.csv"

    community = initialize(**kw)
    
    # params = {"pS", "b", "a", "c", "beta", "absp", "d", "m", "l", "popinit", "phageinit"}
    print("pS:\t%s"%community.pS)
    print("b:\t%s"%community.strains["s1"].b)
    print("a:\t%s"%community.strains["s1"].a)
    print("c:\t%s"%community.strains["s1"].c)
    print("beta:\t%s"%community.phages["p1"].beta)
    print("absp:\t%s"%community.phages["p1"].absp)
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

    for i in range(timesteps): # run for # of timesteps
        
        timestep(i)
        pRichness.append( len( community.phagesPopDict() ) )
        maxPRich = max(pRichness)
        sRichness.append( community.richness() )
        maxSRich = max(sRichness)
        cRichness.append( community.spacerRichness() )

        community.summary["richness"] = maxSRich
        community.summary["phageRichness"] = maxPRich 

    community.summary["pop"] = community.totalComSize
    community.summary["phage"] = community.phagePopOverTime[-1]
    community.summary["resistant"] = community.resOverTime[-1]
    community.summary["vulnerable"] = community.vulnOverTime[-1]

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
    # print("max resistant: %s"%(max(community.resOverTime)))
    df1 = pd.DataFrame( 
        list( 
            zip(
                community.comSizeOverTime, 
                community.phagePopOverTime,
                community.resOverTime,
                community.vulnOverTime, 
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
        columns = ['HostRichness','PhageRichness','time'],
        )

    df1.to_csv(outputMain)
    df2.to_csv(outputRichness)
    # for phage in community.phages.values:
    #     print(phage.genome)
    return None

def multi_sim(sims, **kw):


    communities = [ initialize(**kw) for i in range(sims) ]

    # output name for run uses provided output name, number of sims, set parameters
    constant_params = sorted(["_%s%s" % (param,val) for param, val in kw.items()]) 
    output = "%s_%ssims%s.csv" % (out, sims, "".join(constant_params)) 

    timestep = [ community.timestep for community in communities ] # list of all timestep functions for all sims
    # print("initial sizes")
    # for community in communities:
    #     print(community.totalComSize)
    df = pd.DataFrame(None, 
        columns=["pop","phage","resistant","vulnerable","richness","phageRichness","pS", "b", "a", "c", "beta", "absp", "d", "m", "l", "popinit", "phageinit"],
    )
    
    timesteps = 1000
    
    for ind in range( len(communities) ):
        community = communities[ind]
        pRichness = [len( community.phagesPopDict() )] # phage richness over time
        maxPRich = pRichness[0]
        sRichness = [len( community.strains )] # strain richness over time
        maxSRich = sRichness[0]

        for i in range(timesteps):

            timestep[ind](i)
            pRichness.append( len( community.phagesPopDict() ) )
            maxPRich = max(pRichness)
            sRichness.append( community.richness() )
            maxSRich = max(sRichness)

            community.summary["richness"] = maxSRich
            community.summary["phageRichness"] = maxPRich 

        community.summary["pop"] = community.totalComSize
        community.summary["phage"] = community.phagePopOverTime[-1]
        community.summary["resistant"] = community.resOverTime[-1]
        community.summary["vulnerable"] = community.vulnOverTime[-1]

        df.loc[ind] = df.columns.map( community.summary )



    # print("after sizes")
    # for community in communities:
    #     print(community.totalComSize)

    # print(df)
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


