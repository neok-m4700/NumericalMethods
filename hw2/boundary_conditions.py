import numpy as np
import common

def initial_domain(dom, Params):
   """Returns a guassian centered around the centre of the hypercube"""
   P = Params
   for i in range(P.Nx):
      x = 1.0*i/P.Nx;
      for j in range(P.Ny):
         y = 1.0*j/P.Nx
         for k in range(P.Nz):
            z = 1.0*k/P.Nx
            dom[common.ijk2idx([i,j,k],[P.Nx,P.Ny,P.Nz])] = P.A*np.exp(
                  -1.0/(2*P.sig*P.sig)*(x*x + y*y + z*z))
