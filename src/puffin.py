import angr
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from simplicial_complex import SimplicialComplex
from itertools import chain, combinations

'''
Main code for supporting topological binary analysis.
'''
def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
def get_union(A,B):
    '''
    Gets the union of tuples
    '''
    return tuple(set(A).union(set(B)))
def get_intersection(A, B):
    '''
    Gets the intersection of tuples
    '''
    return tuple(set(A) & set(B))
def filter_cfg(cfg, k, metric="centrality"):
    '''
    Filtration of binary data using a distance metric over the adjacency matrix of a CFG.
    '''
    # compute G
    G = cfg.graph
    l = []
    
    # for each node n in G, compute the top k central nodes, and return them in a list.
    V = nx.degree_centrality(G)
    for item in V:
        n = V[item]
        if n <= k:
            l.append(item)
    return l

def check_faces(cfg, simplices, dimension):
    '''
    Ensure that an edge exists in the o.g. graph if it's in the set of simplices
    '''
    G = cfg.graph
    if (dimension < 2):
        return simplices
    return 

def build_simplex(l):
    '''
    Returns 1, 2, and 3-simplices from the filtered CFG data.
    '''
    # 1-simplices
    d1 = []
    for i in l:
        d1.append([i])
    # 2-simplices
    d2 = []
    for s in powerset(l):
        if len(s) == 2:
            d2.append(s)
    
    # 3-simplices
    d3 = []
    for s in powerset(l):
        if len(s) == 3:
            d3.append(s)
    return [d1,d2,d3]

if __name__ == "__main__":

    # first, build the angr CFG 
    filename = '../binaries/bin/hello'
    blob = angr.Project(filename, load_options={'auto_load_libs':False})
    cfg = blob.analyses.CFGEmulated(keep_state=True)
    l = filter_cfg(cfg,0.2)
    A = build_simplex(l)
    #check_faces(cfg, A, 2)
    SC = SimplicialComplex(A)
    H0 = SC.compute_homologies(1)
    H1 = SC.compute_homologies(2)
    H2 = SC.compute_homologies(3)

    print("--------------- S T A T S ---------------") 
    print("Dimension of SC: " + str(SC.dimension))
    print("Rank H0: " + str(SC.compute_homology_rank(1)))
    print("Rank H1: " + str(SC.compute_homology_rank(2)))
    print("Rank H2: " + str(SC.compute_homology_rank(3)))
    print("Euler Characteristic: " + str(SC.compute_euler_characterisic()))
    print("--------------- S T A T S ---------------")
