import angr
import numpy as np
import matplotlib.pyplot as plt
from sklearn import ensemble
from sklearn.datasets import make_classification
from sklearn import model_selection
import networkx as nx
from simplicial_complex import SimplicialComplex
from itertools import chain, combinations, permutations
import pandas as pd
from operator import itemgetter
from scipy.stats import pearsonr
'''
Main code for supporting topological binary analysis.
TODO: 

1. Replace test data with actual data
2. Polish off the computation of the path complexes
'''
def is_subset(t1, t2):
    return set(t1).issubset(t2)
def get_intersection(A, B):
    '''
    Gets the intersection of tuples
    '''
    return tuple(set(A) & set(B))

def get_pval(df):
    col = pd.DataFrame(columns=df.columns)
    pval = col.transpose().join(col, how='outer')
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            pvalues[r][c] = round(pearsonr(tmp[r], tmp[c])[1], 4)
    return pvalues
def ordered_powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(permutations(s, r) for r in range(len(s)+1))

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def filter_cfg_new(cfg, k):
    '''
    Use the directed distance metric as defined by Chartand to identify paths within a distance 
    of k. The lengths of these paths correspond to n-simplices within a path complex.
    '''
    G = cfg.graph
    l = []

    # use the most central node as our starting point
    central = []
    V = nx.degree_centrality(G)
    for item in V:
        n = V[item]
        central.append((n, item))
    central.sort(key=itemgetter(0), reverse=True)
    seed = central[0][1]

    for v in G.nodes:
        A = nx.all_simple_paths(G, seed, v, cutoff=k)
        B = nx.all_simple_paths(G, v, seed, cutoff=k)
        
        # lengths of all the possible paths between nodes
        d1 = [(path, len(path)) for path in A]
        d2 = [(path, len(path)) for path in B]
        
        # sort the distances and paths
        d1.sort(key = itemgetter(1), reverse=True)
        d2.sort(key = itemgetter(1), reverse=True)
        
        if (d1 == [] and d2 == []):
            l.append([])
        elif (d1 == [] and d2 != []):
            l.append(d2[0][0])
        elif (d2 == [] and d1 != []):
            l.append(d1[0][0])
        else:
            # get the longest path and stuff that into a simplex.
            if (d1[0][1] > d2[0][1]):
                l.append(d1[0][0])
            elif (d2[0][1] > d1[0][1]):
                l.append(d2[0][0])
            else:
                l.append(d1[0][0])
    return l

def build_simplex(paths, cfg):
    '''
    Returns 1, 2, and 3-simplices from the filtered CFG data.
    '''

    G = cfg.graph
    
    # 0-simplices
    E = []
    for path in paths:
        for v in path:
            E.append(v)

    E = list(set(E))
    d0 = []
    for vertex in E:
        d0.append([vertex])
    
    # 1-simplices
    d1 = []
    for path in paths:
        for i in range(0, len(path)-1):
            if (len(path) >= 2):
                j = i+1
                if (G.has_edge(path[i], path[j])):
                    t = (path[i], path[j])
                    if t in powerset(path):
                        d1.append(t)
                        
    d1 = [t for t in (set(tuple(i) for i in d1))]
 
    # 2-simplices
    d2 = []
    for path in paths:
        if (len(path) >= 3):
            for i in range(0, len(path)-2):
                j = i+1
                k = j+1
                if (G.has_edge(path[i], path[j])) and (G.has_edge(path[j], path[k])):
                    t = (path[i], path[j], path[k])
                    if t in powerset(path):
                        d2.append(t)
                    
    d2 = [t for t in (set(tuple(i) for i in d2))]

    # 3-simplices
    d3 = []
    for path in paths:
        if (len(path) >= 4):
            for i in range(0, len(path)-3):
                j = i + 1
                k = j + 1
                l = k + 1
                if (G.has_edge(path[i], path[j])) and (G.has_edge(path[j], path[k])) and (G.has_edge(path[k], path[l])):
                    t = (path[i], path[j], path[k], path[l])
                    if t in powerset(path):
                        d3.append(t)
    d3 = [t for t in (set(tuple(i) for i in d3))]
    return [d0, d1, d2, d3]

# build a networkx graph from the 2d edges in a simplex.
def graph_from_simplex(l):
    G = nx.Graph()
    for edge in l:
        G.add_edge(edge[0], edge[1])
    return G

if __name__ == "__main__":

    '''
    Ci = [[0],[1],[2],[3],[4],[5],[6],[7],[8], [(0,1), (0,2), (0,3), (0,4), (0,5), (0,6), (0,7), (0,8), (1,2), (3,4), (3,5), (4,5), (6,7), (6,8), (7,8)], [(0,1,2), (0,3,4), (0,3,5), (0,4,5), (0,6,7), (0,6,8), (3,4,5),(6,7,8)], [(0,6,7,8), (0,3,4,5)]]
    A = SimplicialComplex(Ci)

    rankH0 = A.compute_homology_rank(1)
    rankH1 = A.compute_homology_rank(2)
    rankH2 = A.compute_homology_rank(3)
    print("Rank H0: " + str(rankH0))
    print("Rank H1: " + str(rankH1))
    print("Rank H2: " + str(rankH2))
    print(A.compute_euler_characteristic())
    
    # A debugger simplicial complex
    Ci = [[[1],[2],[3]], [(1,2),(1,3), (2,3)], [(1,2,3)]]
    A = SimplicialComplex(Ci)
        
    
    H0 = A.compute_homologies(1)
    H1 = A.compute_homologies(2)
    H2 = A.compute_homologies(3)
    
    rankH0 = A.compute_homology_rank(1)
    rankH1 = A.compute_homology_rank(2)
    rankH2 = A.compute_homology_rank(3)
    print("Rank H0: " + str(rankH0))
    print("Rank H1: " + str(rankH1))
    print("Rank H2: " + str(rankH2))
    print(A.compute_euler_characteristic())
    '''

    filename_list = ['../binaries/bin/obfuscated/helloobf', '../binaries/bin/orig/hello', '../binaries/bin/obfuscated/t1obf', '../binaries/bin/orig/t1', '../binaries/bin/obfuscated/t2obf', '../binaries/bin/orig/t2', '../binaries/bin/obfuscated/t3obf', '../binaries/bin/orig/t3', '../binaries/bin/obfuscated/t4obf', '../binaries/bin/orig/t4', '../binaries/bin/obfuscated/t5obf', '../binaries/bin/orig/t5']
    # first, build the angr CFG 
    H0_hlist = []
    H1_hlist = []
    H2_hlist = []
    H3_hlist = []
    
    is_obf = []
    iota = []
    for filename in filename_list:
        
        blob = angr.Project(filename, load_options={'auto_load_libs':False})
        cfg = blob.analyses.CFGEmulated(keep_state=True)
        print(filename[-3:])
        if (filename[-3:] == "obf"):
            obf = 1
        else:
            obf = 0
        # now get all nodes within distance k
        distances = []
        H0 = []
        H1 = []
        H2 = []
        H3 = []
        
        for i in range(0, 40):
            print("DISTANCE = " + str(i))
            is_obf.append(obf)
            paths = filter_cfg_new(cfg, i)
            Cp = build_simplex(paths, cfg)
            distances.append(i)
            
            # build the complex
            SC = SimplicialComplex(Cp)

            H0_hlist.append(SC.compute_homology_rank(1))
            H1_hlist.append(SC.compute_homology_rank(2))
            H2_hlist.append(SC.compute_homology_rank(3))
            H3_hlist.append(SC.compute_homology_rank(4))
            
            H0.append(SC.compute_homology_rank(1))
            H1.append(SC.compute_homology_rank(2))
            H2.append(SC.compute_homology_rank(3))
            H3.append(SC.compute_homology_rank(4))
            
            iota.append(SC.compute_cyclomatic_complexity())

            print("--------------- S T A T S ---------------") 
            print("Dimension of SC: " + str(SC.dimension))
            print("Rank H0: " + str(SC.compute_homology_rank(1)))
            print("Rank H1: " + str(SC.compute_homology_rank(2)))
            print("Rank H2: " + str(SC.compute_homology_rank(3)))
            print("Rank H3: " + str(SC.compute_homology_rank(4)))
            print("Cyclomatic Complexity: " + str(SC.compute_cyclomatic_complexity()))
            print("Euler Characteristic: " + str(SC.compute_euler_characteristic()))
            print("--------------- S T A T S ---------------")

        # try graphing the persistent homologies!
        plt.plot(distances, distances, color='r')

        labels = []
        labels.append('distances')
        # if the homologies are 0-vectors...get rid of them.
        if (sum(np.array(H0)) != 0):
            for n in range(0, len(H0)):
                if (H0[n] == 0):
                    H0[n] = np.nan
            plt.scatter(H0, distances, color='r')
            labels.append('H_0')
        if (sum(np.array(H1)) != 0):
            for n in range(0, len(H1)):
                if (H1[n] == 0):
                    H1[n] = np.nan
            plt.scatter(H1, distances, color='g')
            labels.append('H_1')
        if (sum(np.array(H2)) != 0):
            for n in range(0,len(H2)):
                if (H2[n] == 0):
                    H2[n] = np.nan
            plt.scatter(H2, distances, color='k')
            labels.append('H_2')
        if (sum(np.array(H3)) != 0):
            for n in range(0, len(H3)):
                if (H3[n] == 0):
                    H3[n] = np.nan
            plt.scatter(H3, distances, color='b')
            labels.append('H_3')
        plt.xlabel('Birth')
        plt.ylabel('Death')
        plt.legend(labels)
        plt.show()
    
       
    # write these out to a dataframe
    df = pd.DataFrame({'H0': H0_hlist, 'H1': H1_hlist, 'H2': H2_hlist, 'H3': H3_hlist, 'iota': iota, 'obf' : is_obf })
    print(df)

    
    X = df.drop("obf", axis=1).values
    y = df["obf"]
    model = ensemble.RandomForestClassifier(n_estimators=100, criterion="entropy", random_state=0)
    model.fit(X,y)

    blob = angr.Project('../binaries/bin/orig/test', load_options={'auto_load_libs':False})
    cfg = blob.analyses.CFGEmulated(keep_state=True)

    H0_hlist = []
    H1_hlist = []
    H2_hlist = []
    is_obf = []
    iota = []
    
    # now get all nodes within distance k
    distances = []

    H0 = []
    H1 = []
    H2 = []
    H3 = []
    for i in range(0,40):
        print("DISTANCE = " + str(i))
        
        paths = filter_cfg_new(cfg,i)    
        Cp = build_simplex(paths, cfg)
        distances.append(i)
        
        # build the complex
        SC = SimplicialComplex(Cp)
    
        H0.append(SC.compute_homology_rank(1))
        H1.append(SC.compute_homology_rank(2))
        H2.append(SC.compute_homology_rank(3))
        H3.append(SC.compute_homology_rank(4))

        iota.append(SC.compute_cyclomatic_complexity())
        print("--------------- S T A T S ---------------") 
        print("Dimension of SC: " + str(SC.dimension))
        print("Rank H0: " + str(SC.compute_homology_rank(1)))
        print("Rank H1: " + str(SC.compute_homology_rank(2)))
        print("Rank H2: " + str(SC.compute_homology_rank(3)))
        print("Rank H3: " + str(SC.compute_homology_rank(4)))
        print("Cyclomatic Complexity: " + str(SC.compute_cyclomatic_complexity()))
        print("--------------- S T A T S ---------------")
            
       
    # write these out to a dataframe
    df = pd.DataFrame({'H0': H0, 'H1': H1, 'H2': H2, 'H3': H3, 'iota': iota})
    pred = model.predict(df)
    print(pred)

