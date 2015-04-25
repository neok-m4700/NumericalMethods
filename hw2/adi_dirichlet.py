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
   C = P.Alpha*P.dt/(3*(P.dx*P.dx))

   data = np.ones([3, P.N])*C

   Ai = sp.sparse.spdiags(data , [0 , -P.Ny*P.Nz , P.Ny*P.Nz ] , P.N , P.N)
   Aj = sp.sparse.spdiags(data , [0 , -P.Nz      , P.Nz      ] , P.N , P.N)
   Ak = sp.sparse.spdiags(data , [0 , -1         , 1         ] , P.N , P.N)
   I  = sp.sparse.eye(P.N)

   # Decomposing left side matrices
   LUi = sp.linalg.lu_factor((I-Ai).toarray())
   LUj = sp.linalg.lu_factor((I-Aj).toarray())
   LUk = sp.linalg.lu_factor((I-Ak).toarray())

   # Right side matrics
   Bi = I + Aj + Ak
   Bj = I + Ak + Ai
   Bk = I + Ai + Aj

   # Initialize dirichlet mask
   M = conditions.boundary_mask(P)  # zero if boundary, 1 otherwise
   Minv = 1 - M  # zero if non-boundary, 1 otherwise
   Boundary = conditions.boundary_vector(P, 0)  #Assuming constant boundary conditions
   #These can be made time dependant if so desired


   #Initialize plot
   fig = plt.figure()
   ax = fig.gca(projection="3d")

   # Plot initial domain
   common.plot_slice(ax, 'b', dom, P, i_slice=P.Nx/2)

   #evolve the system
   for i in range(P.nSteps):
      dom = sp.linalg.lu_solve(LUi, Bi*dom)
      dom = sp.linalg.lu_solve(LUj, Bj*dom)
      dom = sp.linalg.lu_solve(LUk, Bk*dom)
      #Apply boundary conditions
      dom = dom * M + Boundary * Minv # Apply boundary conditions


   common.plot_slice(ax, 'r', dom, P, i_slice=P.Nx/2)
   plt.show()

if __name__ == "__main__":
    main()




