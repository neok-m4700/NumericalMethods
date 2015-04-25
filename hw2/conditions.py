## Initial, boundary and source equations
import numpy as np
import common

def initial_domain(dom, Params):
   """Returns a guassian centered around the centre of the hypercube"""
   P = Params
   for i in range(P.Nx):
      x = 1.0*i/P.Nx - 0.5;
      for j in range(P.Ny):
         y = 1.0*j/P.Nx - 0.5;
         for k in range(P.Nz):
            z = 1.0*k/P.Nx - 0.5;
            dom[common.ijk2idx([i,j,k],[P.Nx,P.Ny,P.Nz])] = (
               P.A*np.exp( -1.0/(2*P.sig*P.sig)*(x*x + y*y + z*z)) +
               np.random.rand()*P.A*P.Randomness)



def boundary_mask(P):
   """Mask the boundary values on the domain vector
   value is zero if the point is on the boudary"""
   dims_size = np.array([P.Nx, P.Ny, P.Nz])
   mask = np.zeros([P.N, 1])
   for idx in range(P.N):
      ijk = common.idx2ijk(idx, dims_size)
      if common.is_inner(ijk, dims_size):
         mask[idx] = 1
   return mask

def boundary_vector(P, t):
   """Returns a boundary vector"""
   bv = np.zeros([P.N, 1])
   return bv

def source(P, tstep):
    """Heat source"""
    return P.S*P.dt

