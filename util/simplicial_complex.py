import numpy as np
from sympy import Matrix, linsolve, symbols

'''
simplicial_complex.py - a script for computing representations of chain groups in chain complexes
specificially. 

DISCLAIMER: 
It is typical in mathematical literature to refer to the dimensions of p-simplices by 0,..p.
Here, I label them by 1,...p instead. So a p-chain dimension of 1 is a 0-chain, 2 is a 1-chain,
etc. I'll try to avoid doing that in the future. It's pretty confusing. And also a little wrong.

'''
class pchain:
    '''
    Python representation of a pchain object
    '''
    def __init__(self, data):
        self.mdata = data
        self.dimension = 0
        return
    def compute_boundary(self):
        # recall the equation for computing the boundary of a pchain!
        total = []
        
        eqn = ''
        if (self.dimension > 1):
            for i in range(0,self.mdata.size):
                arr = np.delete(self.mdata, i)
               
                if (i is not self.mdata.size-1):
                    if (i%2 == 0):
                        eqn += "(p)" + str(arr) + " + "
                        total.append((1, arr))
                    else:
                        eqn += "(n)" + str(arr) + " + "
                        total.append((-1, arr))
                else:
                    if (i%2 == 0):
                        eqn += "(p)" + str(arr)
                        total.append((1, arr))
                    else:
                        eqn += "(n)" + str(arr)
                        total.append((-1, arr))
            return total, "Boundary of " + str(self.mdata) + ": " + eqn
        else:
            return total, "Boundary of " + str(self.mdata) + ": 0"
class SimplicialComplex:
    '''
    Python representation of a generic simplicial complex.
    '''
    def __init__(self, Cp):
        self.dimension = len(Cp)+1
        self.Cp = Cp
        for i in Cp:
            if i == []:
                self.dimension -= 1
        self.pchains = self.init_pchains()
        return
    def init_pchains(self):
        '''
        Iteratively build  out pchain objects for each of the p-dimensional chains in deltas
        '''
        plist = []
        max_dimension = 0
        for chain in self.Cp:
            for data in chain:
                # get the dimension of the pchain
                if not isinstance(data, int) or not isinstance(data, list):
                    dim = len(data)
                else:
                    dim = 1
                p = pchain(np.array(data))
                p.dimension = dim
                if (max_dimension < p.dimension):
                    max_dimension = p.dimension
                plist.append(p)
        if max_dimension != self.dimension:
            self.dimension = max_dimension
        return plist 
    def compute_boundary_matrix(self, dimension):
        '''
        Returns the boundary matrix representation of the chains that span Cp
        and Cp-1
        '''
        if (dimension == 1 or dimension >= self.dimension):
            return []
        Cp = self.get_pchains(dimension)
        C_ = self.get_pchains(dimension - 1)
    
        # build an m x n numpy matrix
        D = np.zeros((len(C_), len(Cp)))
        
        # set an index equal to 1 if Cp-1 belongs to the boundary of Cp, 0 if not
        for i in range(0, len(C_)):
            for j in range(0, len(Cp)): 
                b,p = Cp[j].compute_boundary()
                res = p.split(":")[1:]
                index = res[0].find(str(C_[i].mdata))
                n = 0
                for k in b:
                    if C_[i].mdata in k[1]:
                        break
                    else:
                        n+=1
                if n == len(b):
                    D[i][j] = 0
                elif b[n][0] == 1:
                    D[i][j] = 1
                elif b[n][0] == -1:
                    D[i][j] = -1
        return D

    def compute_cycles(self, dimension):
        '''
        Row-reduce to find the kernel of the boundary map as the kernel of a linear map.
        '''

        # C0 = spanZ0!
        if (dimension == 1):
            return self.get_pchains(dimension)
        
        # first, compute the boundary matrix
        M = self.compute_boundary_matrix(dimension)
        M = Matrix(M)
        kernel = M.nullspace()
        for basis in kernel:
            basis = np.array(basis)
        # now write the actual basis in terms of the simplices
        Cp = self.get_pchains(dimension)
        KERNEL = ''
        indx = 0
        for basis in kernel:
            for i in range (0, len(basis)):
                if (int(basis[i]) > 0):
                    KERNEL += str(Cp[i].mdata)
                    if (i+1 < len(basis)):
                        KERNEL += ' + '
                elif (int(basis[i]) < 0):
                    KERNEL += '-' + str(Cp[i].mdata)
                    if (i+1 < len(basis)):
                        KERNEL += ' + '
                if (i+1 == len(basis)):
                    KERNEL += ','
        # format this output so we can produce the kernel generators
        return KERNEL
        
    def get_pchains(self, p):
        '''
        Return the p-chains of dimension p
        '''
        res = []
        for pchain in self.pchains:
            if(pchain.dimension == p):
                res.append(pchain)
        return res

    def compute_homologies(self, dimension):
        '''
        Returns the homology representation of the simplicial complex.
        '''
        Bp = []
        Zp = []
        res = '<'
        res2 = '<'

        if (dimension >= 1 and dimension < self.dimension):
            boundary_pchains = self.get_pchains(dimension+1)

            for chain in boundary_pchains:
                b,p = chain.compute_boundary()
                p = p.split(":")[1:]
                Bp.append(p)
    
        kernel_pchains = self.get_pchains(dimension)
        kernel = self.compute_cycles(dimension)
        if (dimension == 1):
            for i in kernel_pchains:
                Zp.append(i.mdata)
        elif (dimension > 1 and dimension < self.dimension):
            Zp.append(kernel)
        else:
            return 0

        # stuff every possible chain into the kernel equation and then ensure that the members of Bp do not belong to the vector spanned by result
        indx = 0
        for image in Bp:
            res += str(image) + ","
            if indx+1 == len(Bp):
                res += str(image) + ">"
            indx+=1

        indx = 0
        for k in Zp:
            if indx+1 == len(Zp):
                res2 += str(k) + ">"
            else:
                indx+=1
                if (dimension == 1):
                    res2 += str(k) + ",  "
                else:
                    res2 = Zp
                    continue
        if (len(res) == 1):
            res += "0>"
        return res2 + " / " + res
    
    def compute_boundary_rank(self, dimension):
        '''
        Compute the ranks of the boundaries
        '''
        if dimension > self.dimension:
            return 0
        M = self.compute_boundary_matrix(dimension)
        rank = np.linalg.matrix_rank(M)
        return rank

    def compute_cycle_rank(self, dimension):
        '''
        Compute the ranks of the cycles: from the Rank-Nullity Theorem,
        we have that 

        rank Zp = col(M) - rank(M)
        '''
        if (dimension == 1):
            return len(self.get_pchains(1))
       
        boundary_rank = self.compute_boundary_rank(dimension)
        M = Matrix(self.compute_boundary_matrix(dimension))
        M_rref = M.rref()[0]
        return M_rref.shape[1] - M.rank()

    def compute_homology_rank(self, dimension):
        '''
        Compute the rank of the homologies:
        first, take the rank of the cycles, and then the rank of the boundaries
        of the next dimension, and then do the arithmetic.
        E.g.:
        rank Hp = rank Zp - rank Bp. 
        '''
        Zp = self.compute_cycle_rank(dimension)
        Bp = self.compute_boundary_rank(dimension+1)    
        return Zp - Bp

    def compute_euler_characteristic(self):
        '''
        Chi = Vertices - Edges + Faces
        '''
        C0 = self.get_pchains(1)
        C1 = self.get_pchains(2)
        C2 = self.get_pchains(3)
        print(len(C0))
        print(len(C1))
        Chi = len(C0) - len(C1) + len(C2)
        return Chi

    def compute_cyclomatic_complexity(self):
        '''
        Iota = Edges - Vertices + 1
        '''
        C0 = self.get_pchains(1)
        C1 = self.get_pchains(2)
        Iota = len(C1) - len(C0) + 1
        return Iota
        
def compute_boundary_with_matrix(M):
    '''
    Pass in a boundary matrix;  
    then use that matrix to produce the boundaries of each relevant chain.

    Right now we're just going to get the generators for this group, 
    returning them as columns of M.
    '''
   
    boundaries = []
    rows, cols = M.shape
    for i in range(0, cols):
        boundaries.append(M[:i])
    return boundaries

if __name__ == "__main__":

    # an integer representation
    cow = 1
    rabbit = 2
    horse = 3
    dog = 4
    fish = 5
    dolphin = 6
    oyster = 7
    broccoli = 8
    fern = 9
    onion = 10
    apple = 11
    '''
    # try not to neglect the vertices here either
    C0 = [(horse), (cow), (rabbit), (dog), (fish), (oyster), (dolphin), (broccoli), (fern), (onion), (apple)]
    
    # manual representation grabbed from the painstaking labor in coding hw #1
    C1 =  [ (cow,rabbit),
      (cow, horse),
      (cow, dog),
      (rabbit, horse),
      (rabbit, dog),
      (horse, dog),
      (fish, dolphin),
      (fish, oyster),
      (dolphin, oyster),
      (broccoli, fern),
      (fern, onion),
      (fern,apple),
    (onion, apple),
    (broccoli, fern),
    (broccoli, apple),
    (broccoli, onion)]

    C2 = [(cow,rabbit, horse), (cow, rabbit, dog), (cow, horse, dog), (rabbit, horse, dog), (fish, dolphin, oyster), (broccoli, fern, onion), (broccoli, fern, apple), (broccoli, onion, apple), (fern, onion, apple)]

    # this is essentially Cp 
    Cp = [C0, C1, C2]
    
    A = SimplicialComplex(Cp)

    print("---------------------------------------------------------")
    print(" C O M P L E X   A " )
    print("---------------------------------------------------------")

    
    print("Cycle rank Z0: " + str(A.compute_cycle_rank(1)))
    print("Boundary rank B0: " + str(A.compute_boundary_rank(2)))
    print("Homology rank H0: " + str(A.compute_homology_rank(1)))

    print("")

    print("Boundary matrix (dell1): ")
    print(A.compute_boundary_matrix(2))
    
    print("Cycle rank Z1: " + str(A.compute_cycle_rank(2)))
    print("Boundary rank B1: " + str(A.compute_boundary_rank(3)))
    print("Homology rank H1: " + str(A.compute_homology_rank(2)))

    print("")
        
    print("Boundary matrix (dell2): ")
    print(A.compute_boundary_matrix(3))
    
    print("Cycle rank Z2: " + str(A.compute_cycle_rank(3)))
    print("Boundary rank B2: " + str(A.compute_boundary_rank(4)))
    print("Homology rank H2: " + str(A.compute_homology_rank(3)))

    print("H0: " + A.compute_homologies(1))
    print("H1:  " + A.compute_homologies(2))
    print("H2: " + A.compute_homologies(3))

    '''
    D0 = [(cow), (rabbit), (fish), (oyster), (broccoli), (onion), (apple), (dog), (horse), (fern), (dolphin)]
    
    D1 = [(cow, rabbit), (cow, fish), (cow, oyster), (cow, oyster), (cow, broccoli), (cow, onion), 
(cow, apple), (rabbit, fish), (rabbit, oyster), (rabbit, broccoli), (rabbit, onion),  (rabbit, apple), (fish, oyster), (fish, broccoli), (fish, onion), (fish, apple), (oyster, broccoli), (oyster, onion), (oyster, apple), (broccoli, onion), (broccoli, apple), (onion, apple), (horse, dog), (horse, dolphin), (horse, fern), (dog, dolphin), (dog, fern), (dolphin, fern)]

    D2 = [(cow, broccoli, apple), (cow, onion, apple), (rabbit, broccoli, apple),  (rabbit, onion, apple), (fish, broccoli, apple), (fish, onion, apple),  
(oyster, broccoli, apple), (oyster, onion, apple)]


    Dp = [D0, D1, D2]
    
    B = SimplicialComplex(Dp)

    print("---------------------------------------------------------")
    print(" C O M P L E X   B " )
    print("---------------------------------------------------------")

    print("Cycle rank Z0: " + str(B.compute_cycle_rank(1)))
    print("Boundary rank B0: " + str(B.compute_boundary_rank(2)))
    print("Homology rank H0: " + str(B.compute_homology_rank(1)))

    print("")

    print("Boundary matrix (dell1): ")
    print(B.compute_boundary_matrix(2))

    print("Cycle rank Z1: " + str(B.compute_cycle_rank(2)))
    print("Boundary rank B1: " + str(B.compute_boundary_rank(3)))
    print("Homology rank H1: " + str(B.compute_homology_rank(2)))

    print("")

        
    print("Boundary matrix (dell2): ")
    print(B.compute_boundary_matrix(3))
       
    print("Cycle rank Z2: " + str(B.compute_cycle_rank(3)))
    print("Boundary rank B2: " + str(B.compute_boundary_rank(4)))
    print("Homology rank H2: " + str(B.compute_homology_rank(3)))

    print("H0: " + B.compute_homologies(1))
    print("H1:  " + B.compute_homologies(2))
    print("H2: " + B.compute_homologies(3))
