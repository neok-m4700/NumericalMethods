## Initial, boundary and source equations
import numpy as np
import common

def initial_domain(P):
   """Returns a guassian centered around the centre of the hypercube"""
   dom = np.zeros([P.N, 1])
   for i in range(P.Nx):
      x = 1.0*i/P.Nx - 0.5;
      for j in range(P.Ny):
         y = 1.0*j/P.Nx - 0.5;
         for k in range(P.Nz):
            z = 1.0*k/P.Nx - 0.5;
            dom[common.ijk2idx([i,j,k],P)] = (
               P.A*np.exp( -1.0/(2*P.sig*P.sig)*(x*x + y*y + z*z)) +
               np.random.rand()*P.A*P.Randomness)
   return dom

def boundary_vector(P, t):
   """Returns a vector representing boundary conditions"""
   bv = np.zeros([P.N, 1])
   return bv

def source(P, tstep):
    """Heat source"""
    return P.S*P.dt

