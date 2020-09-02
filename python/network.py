## network.py
## Ford Fishman

import networkx as nx
from networkx.algorithms import bipartite

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

def plotBipartite(B:nx.Graph)->None:
    """
    Plots a bipartite graph using matplotlib
    B: bipartite infection graph
    """
    # 
    x = {n for n,d in B.nodes(data=True) if d["bipartite"]==0}
    y = set(B) - x
    pos = dict()
    pos.update( (n, (1, i)) for i,n in enumerate(x))
    pos.update( (n, (2, i)) for i,n in enumerate(y))

    # plt.subplot(121)
    nx.draw(B, pos=pos, with_labels=True, font_weight='bold')
    plt.savefig('/mnt/c/Users/fordf/OneDrive/Documents/GitHub/CommunityImmunity/output/test.png')
    return None


######### TESTING ##########

# a = [1,2,3,4,5]
# b = ['a','b','c','d','e']
# e = [('a',1), ('d', 4), ('c',2), ('b',5), ('e',3)]
# B = createNetwork(e, a, b)



# # print(nx.is_connected(B))

# # x, y = bipartite.sets(B)
# x = {n for n,d in B.nodes(data=True) if d["bipartite"]==0}
# y = set(B) - x
# pos = dict()
# pos.update( (n, (1, i)) for i,n in enumerate(x))
# pos.update( (n, (2, i)) for i,n in enumerate(y))
# # G = nx.Graph()

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

# H = nx.DiGraph(G)

# H = nx.petersen_graph()

