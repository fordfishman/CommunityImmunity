## network.py
## Ford Fishman

import networkx as nx
from networkx.algorithms import bipartite
import numpy as np
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

def bipartNodes(graph):
    x = {n for n,d in graph.nodes(data=True) if d["bipartite"]==0}
    y = set(graph) - x
    return (x, y)

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

    x, y = bipartNodes(G)
    A = nx.bipartite.biadjacency_matrix(G, row_order=x, column_order=y)
    m, n = A.shape

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

    print(rowsum)

    A = csr_matrix(B)
    print(A.toarray())
    nrowpairs = list()
    ncolpairs = list()
    denom = (n*(n-1)/2) + (m*(m-1)/2) 

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
    # print(ncolpairs)
    # print(nrowpairs)

    return nodf
    # return None


def modularity(G):

    x, y = bipartNodes(G)
    A = nx.bipartite.biadjacency_matrix(G, row_order=x, column_order=y)

######### TESTING ##########

a = [1,2,3,4,5]
b = ['a','b','c','d','e']
# e = [('a',1), ('d', 4), ('c',2), ('b',5), ('e',3)]
e = [
    ('a',1), 
    ('b',1), ('b',2), 
    ('c',1), ('c',2), ('c',3), 
    ('d',1), ('d',2), ('d',3), ('d', 4),
    ('e',1), ('e',2),('e',3), ('e', 4), ('e', 5)
]
B = createNetwork(e, a, b)
# x = {n for n,d in B.nodes(data=True) if d["bipartite"]==0}
# y = set(B) - x
# print(x)
# print(y)
# nestedness(B)
print(nestedness(B))
plotBipartite(B, '/mnt/c/Users/fordf/Downloads/network.png')

# print(rowsum.reshape(5,1))


# G = nx.Graph()

# G.add_node(1)

# G.add_nodes_from( [2,3] )

# G.add_nodes_from([
#     (4, {'color': 'red'}),
#     (5, {'color': 'green'}),
# ])

# G.add_edge(1,2)

# e = (2,3)

# G.add_edge(*e)  

# G.add_edges_from( [ (1,2), (1,3) ] )

# G.add_node('spam')

# # print(G.nodes)

# H = nx.DiGraph(G)

# H = nx.petersen_graph()

