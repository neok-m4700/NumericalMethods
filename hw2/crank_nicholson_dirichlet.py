import numpy as np
import scipy as sp
import scipy.sparse
import scipy.linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import common
from read_params import read_params
import conditions



def main():
   #Read parameters
   P = read_params("params.txt")

   # Initialize empty domain vector
   dom = np.zeros([P.N, 1])

   # Initialze with inital conditions
   conditions.initial_domain(dom, P)

   #Set C
   C = P.Alpha*P.dt/(2*(P.dx*P.dx))

   # Offsets of i, j, and k for the various terms in FTCS 
   # [0,0,0] -> [i,j,k]
   # [0,-1,0] -> [i,j-1,k] etc
   offsets_B = np.array(np.mat( """0  0  0;
           -1  0  0; 0 -1  0; 0  0 -1;
            1  0  0; 0  1  0; 0  0  1"""))
   coeffs_B = [1-6*C, C, C, C, C, C, C]

   # Initialize B matrix, such that A*dom(m+1) = B*dom(m)
   B = sp.sparse.eye(P.N).tolil() #Diagonal values are 1 by default
   B = common.gen_matrix(B, offsets_B, coeffs_B, P, wrap=False)
   B = sp.sparse.csr_matrix(B)

   offsets_A = np.array(np.mat( """0  0  0;
           -1  0  0; 0 -1  0; 0  0 -1;
            1  0  0; 0  1  0; 0  0  1"""))
   coeffs_A = [1 + 6*C, -C, -C, -C, -C, -C, -C]

   # Initialize B matrix, such that A*dom(m+1) = B*dom(m)
   A = sp.sparse.eye(P.N).tolil() #Diagonal values are 1 by default
   A = common.gen_matrix(A, offsets_A, coeffs_A, P, wrap=False)

   # Decompose A
   PLU = sp.linalg.lu_factor(A.toarray())

   # Initialize dirichlet mask
   M = conditions.boundary_mask(P)  # zero if boundary, 1 otherwise
   Minv = 1 - M  # zero if non-boundary, 1 otherwise
   Boundary = conditions.boundary_vector(P, 0)  #Assuming constant boundary conditions
   #These can be made time dependant if so desired


   #Initialize plot
   fig = plt.figure()
   ax = fig.gca(projection="3d")

   # Plot initial domain
   common.plot_slice(ax, 'g', dom, P, i_slice=P.Nx/2)

   #evolve the system
   for i in range(P.nSteps):
       #solve for x, where Ax = B*dom
       dom = sp.linalg.lu_solve(PLU, B*dom)
       dom = dom + conditions.source(P, i) # Add heat source
       dom = dom * M + Boundary * Minv # Apply boundary conditions

   common.plot_slice(ax, 'r', dom, P, i_slice=P.Nx/2)
   plt.show()


if __name__ == "__main__":
    main()
