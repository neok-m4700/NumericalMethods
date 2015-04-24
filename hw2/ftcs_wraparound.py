import numpy as np
import scipy as sp
import scipy.sparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import common
from read_params import read_params
from initial_conditions import initial_domain



def main():
   #Read parameters
   P = read_params("params.txt")

   # Initialize empty domain vector
   dom = np.zeros([P.N, 1])

   # Initialze with inital conditions
   initial_domain(dom, P)

   #Set C
   C = P.Alpha*P.dt/(P.dx*P.dx)

   # Offsets of i, j, and k for the various terms in FTCS 
   # [0,0,0] -> [i,j,k]
   # [0,-1,0] -> [i,j-1,k] etc
   offsets = np.array(np.mat( """0  0  0;
           -1  0  0; 0 -1  0; 0  0 -1;
            1  0  0; 0  1  0; 0  0  1"""))
   coeffs = [1-6*C, C, C, C, C, C, C]

   # Initialize B matrix, such that A*dom(m+1) = B*dom(m)
   B = sp.sparse.lil_matrix((P.N,P.N))
   B = common.gen_matrix(B, offsets, coeffs, P, wrap=True)
   B = sp.sparse.csr_matrix(B)

   #Initialize plot
   fig = plt.figure()
   ax = fig.gca(projection="3d")

   # Plot initial domain
   common.plot_slice(ax, 'g', dom, P, i_slice=P.Nx/2)

   #evolve the system
   for i in range(P.nSteps/40):
       dom = B*dom

   common.plot_slice(ax, 'r', dom, P, i_slice=P.Nx/2)
   plt.show()


if __name__ == "__main__":
    main()
