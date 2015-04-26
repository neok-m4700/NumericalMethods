## Initial, boundary and source equations
import numpy as np
import common
from params import Params

def initial_domain():
   """Returns a guassian centered around the centre of the hypercube"""
   P = Params.instance
   dom = np.zeros([P.N, 1])
   for i in range(P.Nx):
      x = 1.0*i/P.Nx - 0.5;
      for j in range(P.Ny):
         y = 1.0*j/P.Nx - 0.5;
         for k in range(P.Nz):
            z = 1.0*k/P.Nx - 0.5;
            dom[common.ijk2idx([i,j,k])] = (
               P.A*np.exp( -1.0/(2*P.sig*P.sig)*(x*x + y*y + z*z)) +
               np.random.rand()*P.A*P.Randomness)
   return dom

def boundary_vector(t):
   """Returns a vector representing boundary conditions"""
   P = Params.instance
   bv = np.zeros([P.N, 1])
   return bv

def source(tstep):
    """Heat source"""
    P = Params.instance
    return P.S*P.dt

