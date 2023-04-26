import angr
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

'''
Main code for supporting topological binary analysis.
'''
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

def build_simplex(G):
    '''
    Returns 1, 2, and 3-simplices from the filtered CFG data.
    '''
    return 

if __name__ == "__main__":

    # first, build the angr CFG 
    filename = '../binaries/bin/hello'
    blob = angr.Project(filename, load_options={'auto_load_libs':False})
    cfg = blob.analyses.CFGEmulated(keep_state=True)
    l = filter_cfg(cfg,0.2)
    
