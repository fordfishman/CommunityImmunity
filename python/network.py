## network.py
## Ford Fishman

import networkx as nx
from networkx.algorithms import bipartite
import numpy as np; import pandas as pd
from scipy.sparse import csr_matrix
import matplotlib as mpl 
import matplotlib.pyplot as plt

mpl.use('Agg') # no UI backend


def createNetwork(edges:list, strainNames:list, phageNames:list):
    """
    Creates a bipartite infection matrix
    -----------------------------------
    phageNames: phage names
    strainNames (list): strain names
    edges (list of tuples): tuples should contain names of connected strain to phage

    """
    # Create bipartite graph
    B = nx.Graph()

    B.add_nodes_from( strainNames, bipartite=0 ) 
    B.add_nodes_from( phageNames, bipartite=1 )

    B.add_edges_from(edges)

    return B

def addToNetwork(B, edges:list, strainNames:list=None, phageNames:list=None):

    if not strainNames is None:
        B.add_nodes_from( strainNames, bipartite=0 ) 
    
    if not phageNames is None:
        B.add_nodes_from( phageNames, bipartite=1 )

    B.add_edges_from(edges)

    return B

def bipartNodes(graph):
    x = {n for n,d in graph.nodes(data=True) if d["bipartite"]==0}
    y = set(graph) - x
    return (x, y)

def adjacencyMatrix(graph):
    x, y = bipartNodes(graph)
    A = nx.bipartite.biadjacency_matrix(graph, row_order=x, column_order=y)
    df = pd.DataFrame.sparse.from_spmatrix(A, index=x, columns = y)
    return df

def plotBipartite(B:nx.Graph, path)->None:
    """
    Plots a bipartite graph using matplotlib
    B: bipartite infection graph
    """
    # 
    # x = {n for n,d in B.nodes(data=True) if d["bipartite"]==0}
    # y = set(B) - x
    x, y = bipartNodes(B)
    pos = dict()
    pos.update( (n, (1, i)) for i,n in enumerate(x))
    pos.update( (n, (2, i)) for i,n in enumerate(y))

    # plt.subplot(121)
    nx.draw(B, pos=pos, with_labels=True, font_weight='bold')
    plt.savefig(path)
    return None

def nestedness(G):

    """
    Calculates NODF
    adapted from R vegan NODF sourcecode: https://rdrr.io/cran/vegan/src/R/nestednodf.R
    """

    # Create adjacency matrix
    x, y = bipartNodes(G)
    A = nx.bipartite.biadjacency_matrix(G, row_order=x, column_order=y)
    m, n = A.shape

    # Sum rows and columns and reorder by sums
    rowsum = A.sum(1)
    colsum = A.sum(0).reshape(m, 1) 

    rowsum_np = np.array([ i[0,0] for i in rowsum ])
    colsum_np = np.array([ i[0,0] for i in colsum ])

    rorder = np.argsort(-1*rowsum_np)[:m]
    corder = np.argsort(-1*colsum_np)[:n]

    B = A.toarray()[rorder,...]
    B = B[...,corder]
    
    rowsum = np.sort(rowsum_np)[::-1]
    colsum = np.sort(colsum_np)[::-1] 

    A = csr_matrix(B)
    print(A.toarray())
    nrowpairs = list()
    ncolpairs = list()
    denom = (n*(n-1)/2) + (m*(m-1)/2) 

    # calculate Nrow

    for i in range(m-1):

        row1 = A[i,:].toarray()
        rowsum1 = rowsum[i]
        # print(rowsum1)

        for j in range(i+1, m): # don't want to repeat over pairs already done

            row2 = A[j,:].toarray()
            rowsum2 = rowsum[j]

            if rowsum1 <= rowsum2 or rowsum1==0 or rowsum2==0:
                nrowpairs.append(0)

            else:
                # print(np.sum((row1 + row2) == 2))
                # print(row2)

                N = np.sum((row1 + row2) == 2)/rowsum2
                nrowpairs.append(N)

    # calculate Ncol

    for i in range(n):

        col1 = A[:,i].toarray()
        colsum1 = colsum[i]

        for j in range(i+1, n):

            col2 = A[:,j].toarray()
            colsum2 = colsum[j]

            if colsum1 <= colsum2 or colsum1==0 or colsum2==0:
                ncolpairs.append(0) 

            else:
                N = sum((col1 + col2) == 2)/colsum2
                ncolpairs.append(N[0])

    Npair = np.array([*ncolpairs, *nrowpairs]) 

    nodf = np.sum(Npair)*100/denom

    return nodf

def plotNQ(summary, path,unset:list=None):

    # summary = summary.dropna(how='all')
    nodf = summary.nodf
    Q = summary.Q

    num_unset = 0
    
    if isinstance(unset, list):
        num_unset = len(unset)
    
    panels = 1 + num_unset*2

    if panels==1:
        fig, axs = plt.subplots(figsize=(8, 4), dpi=80)
        axs.scatter(Q, nodf, c='black');axs.set_xlabel('Modularity (Q)'); axs.set_ylabel('Nestedness (NODF)');
        return None
    
    fig, axs = plt.subplots(panels, 1, figsize=(8, 4*panels), dpi=80)
    axs[0].scatter(Q, nodf, c='black');axs[0].set_xlabel('Modularity (Q)'); axs[0].set_ylabel('Nestedness (NODF)');

    ind = 0
    for i in range(num_unset):
        var = unset[i]
        data = summary[var]

        if np.mean(data) > 1e4 or np.mean(data) < 1e-4:
            data = np.log10(data)

        ind += 1
        axs[ind].scatter(data, Q, c='black');axs[ind].set_xlabel(var); axs[ind].set_ylabel('Modularity (Q)');
        ind +=1
        axs[ind].scatter(data, nodf, c='black');axs[ind].set_xlabel(var); axs[ind].set_ylabel('Nestedness (NODF)');

    plt.savefig(path)

    return None

######### TESTING ##########

a = [1,2,3,4,5]
b = ['a','b','c','d','e']
# e = [('a',1), ('d', 4), ('c',2), ('b',5), ('e',3)]
e = [
    # ('a',1), 
    # ('b',1), ('b',2), 
    # ('c',1), ('c',2), ('c',3), 
    # ('d',1), ('d',2), ('d',3), ('d', 4),
    ('e',1), ('e',2),('e',3), ('e', 4), ('e', 5)
]
B = createNetwork(e, a, b)
addToNetwork(B, [('a',6)],[6])
# print(adjacencyMatrix(B))
# x = {n for n,d in B.nodes(data=True) if d["bipartite"]==0}
# y = set(B) - x
# print(x)
# print(y)
# nestedness(B)
# print(nestedness(B))
# modularity(B)
# plotBipartite(B, '/mnt/c/Users/fordf/Downloads/network.png')



