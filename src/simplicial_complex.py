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
                if (i !=  self.mdata.size-1):
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
                if not isinstance(data, int):
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
        if (dimension == 1 or dimension > self.dimension):
            return []
        Cp = self.get_pchains(dimension)
        C_ = self.get_pchains(dimension - 1)

        K = [chain.mdata.tolist() for chain in C_]

        if Cp == [] or C_ == []:
            return []
        
        # do some preliminary checks...
        for i in range(0, len(C_)):
            for j in range(0, len(Cp)):
                b,p = Cp[j].compute_boundary()
                for k in b:
                    # if an item in the computed boundary does not belong..get rid of it
                    if k[1].tolist() not in K:
                        if Cp[j] in self.pchains:
                            self.pchains.remove(Cp[j])
        # build an m x n numpy matrix
        D = np.zeros((len(C_), len(Cp)))
        
        # set an index equal to 1 if Cp-1 belongs to the boundary of Cp, 0 if not
        for i in range(0, len(C_)):
            for j in range(0, len(Cp)): 
                b,p = Cp[j].compute_boundary()
                n = 0
                for k in b:
                    if C_[i].mdata.size == 1:
                        if C_[i].mdata in k[1]:
                            break
                    if np.array_equal(C_[i].mdata,  k[1]):
                        break
                    else:
                        n+=1
                if n == len(b):
                    D[i][j] = 0
                elif b[n][0] == 1:
                    D[i][j] = 1
                elif b[n][0] == -1:
                    D[i][j] = -1
        if D == []:
            self.dimension -=1
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
        if M == []:
            return 0
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

        if M == []:
            return 0
        
        M_rref = np.array(M.rref()[0])
        print("SHAPE: " + str(M_rref.shape[1]))
        print("RANK: " + str(M.rank()))
        
        # get rid of zero-valued columns
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

        print("Z" + str(dimension) + ":" + " " + str(Zp))
        print("B" + str(dimension+1) + ":" + " " + str(Bp))
        
        return Zp - Bp

    def compute_euler_characteristic(self):
        '''
        Chi = Vertices - Edges + Faces
        '''
        Chi = 0
        C0 = self.get_pchains(1)
        C1 = self.get_pchains(2)
        C2 = self.get_pchains(3)
        Chi = len(C0) - len(C1) + len(C2)
        return Chi

    def compute_cyclomatic_complexity(self):
        C0 = self.get_pchains(1)
        C1 = self.get_pchains(2)
        C2 = self.get_pchains(3)
        
        iota = len(C1) - len(C0) + 1
        return iota
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
