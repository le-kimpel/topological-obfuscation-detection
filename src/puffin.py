import angr
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from simplicial_complex import SimplicialComplex
from itertools import chain, combinations, permutations
from operator import itemgetter

'''
Main code for supporting topological binary analysis.
TODO get more complex data/binaries. 
'''
def ordered_powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(permutations(s, r) for r in range(len(s)+1))

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def filter_cfg(cfg, k, metric="distance", undirected=True):
    '''
    Filtration of binary data using a distance metric over the graph representation of the CFG, using seed data chosen based on vertex importance.
    We can interpret this data as either the abstract simplicial complex over an undirected graph, 
    or we can build a path complex from the CFG digraph. 
    '''
    # compute G
    if (undirected == True):
        G = cfg.graph.to_undirected()
    else:
        G = cfg.graph
    l = []
    central = []
    V = nx.degree_centrality(G)
    for item in V:
        n = V[item]
        central.append((n, item))
    central.sort(key=itemgetter(0), reverse=True)
    seed = central[0][1]

    # starting with the seed node, compute all nodes within distance k of seed.
    ans = list(nx.bfs_edges(G, source=seed, depth_limit=k))
    for path in ans:
        l.append(path)
    
    return l

def check_faces(cfg, d1, simplices, dimension, undirected=True):
    '''
    Ensure that an edge exists in the o.g. graph if it's in the set of simplices
    '''
    if (undirected==True):
        G = cfg.graph.to_undirected()
    else:
        G = cfg.graph
    M = []
    if (dimension < 2):
        return simplices
    if (dimension == 2):
        for s in simplices:
            if G.has_edge(s[0],s[1]) and [s[0]] in d1 and [s[1]] in d1:
                M.append(s)
        return M
    elif (dimension == 3):
        for s in simplices:
            if G.has_edge(s[0], s[1]) and [s[0]] in d1 and [s[1]] in d1:
                if G.has_edge(s[1], s[2]) and [s[2]] in d1:
                    if (undirected == True):
                        if (G.has_edge(s[2], s[0])):
                            M.append(s)
                    else:
                        M.append(s)
        return M
def build_simplex(l):
    '''
    Returns 1, 2, and 3-simplices from the filtered CFG data.
    '''
    # 1-simplices
    E = []
    for i in l:
        for v in i:
            E.append(v)

    E = list(set(E))
    d1 = []
    for vertex in E:
        d1.append([vertex])
    
    # 2-simplices
    d2 = []
    for s in powerset(E):
        if len(s) == 2:
            d2.append(s)
    # 3-simplices
    d3 = []
    for s in powerset(E):
        if len(s) == 3:
            d3.append(s)
    # ensure we get rid of duplicates
    return [d1,list(set(d2)),list(set(d3))]

# build a networkx graph from the 2d edges in a simplex.
def graph_from_simplex(l):
    G = nx.Graph()
    for edge in l:
        G.add_edge(edge[0], edge[1])
    return G

if __name__ == "__main__":

    # first, build the angr CFG 
    filename = '../binaries/bin/obfuscated/helloobf'
    blob = angr.Project(filename, load_options={'auto_load_libs':False})
    cfg = blob.analyses.CFGEmulated(keep_state=True)

    nx.draw(cfg.graph.to_undirected())
    plt.show()

    # now get all nodes within distance k

    distances = []
    H0_hlist = []
    H1_hlist = []
    H2_hlist = []
    for i in range(1,19):
        print("DISTANCE = " + str(i))
        l = filter_cfg(cfg,i, undirected=False)    
        A = build_simplex(l)
        distances.append(i)
        Cp = [check_faces(cfg, A[0], A[indx], indx+1, undirected=False) for indx in range(0,len(A))]
        #N = graph_from_simplex(Cp[1])
        #nx.draw(N)
        #plt.show()
    
        # build the complex
        SC = SimplicialComplex(Cp)
        H0 = SC.compute_homologies(1)
        H1 = SC.compute_homologies(2)
        H2 = SC.compute_homologies(3)
        
        H0_hlist.append(SC.compute_homology_rank(1))
        H1_hlist.append(SC.compute_homology_rank(2))
        H2_hlist.append(SC.compute_homology_rank(3))
    
        print("--------------- S T A T S ---------------") 
        print("Dimension of SC: " + str(SC.dimension))
        print("Rank H0: " + str(SC.compute_homology_rank(1)))
        print("Rank H1: " + str(SC.compute_homology_rank(2)))
        print("Rank H2: " + str(SC.compute_homology_rank(3)))
        print("Cyclomatic Complexity: " + str(SC.compute_cyclomatic_complexity()))
        print("--------------- S T A T S ---------------")

    # try graphing the persistent homologies!
    plt.scatter(distances, H0_hlist, color='r')
    plt.scatter(distances, H1_hlist, color='g')
    plt.scatter(distances, H2_hlist,  color='k')
    plt.xlabel('Path Distance')
    plt.ylabel('Ranks of H_0, H_1, and H_2')
    labels = ['H_0', 'H_1', 'H_2']
    plt.legend(labels)
    plt.show()

