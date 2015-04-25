import numpy as np
import scipy as sp
import scipy.sparse
import common
import conditions

class Dirichlet(object):
   """Implements the Wraparound boundary conditions
   Actual conditions specified as a function in module 'conditions' 
   This module enforces those conditions"""

   def __init__(self, P):
      # Initialize dirichlet mask
      self.P        = P
      self.mask     = self.gen_mask()  # zero if boundary, 1 otherwise
      self.mask_inv = 1 - self.mask  # zero if non-boundary, 1 otherwise
      #Assuming constant boundary conditions 
      #These can be made time dependant if so desired
      self.boundary = conditions.boundary_vector(P, 0) 

   def gen_mask(self):
      """Mask the boundary values on the domain vector
      value is zero if the point is on the boudary"""
      P    = self.P
      mask = np.zeros([P.N, 1])
      for idx in range(P.N):
         ijk = common.idx2ijk(idx, P)
         if common.is_inner(ijk, P):
            mask[idx] = 1
      return mask

   def update(self, dom):
      return dom*self.mask + self.boundary*self.mask_inv
   

class Wraparound(object):
   """Implements the Wraparound boundary conditions"""

   def __init__(self, P):
      # Initialize dirichlet mask
      self.P = P
      wrapper = sp.sparse.lil_matrix((P.N,P.N))
      for i in range(P.Nx):
         iw   = (i-1) % (P.Nx - 2) + 1
         for j in range(P.Ny):
            jw    = (j-1) % (P.Ny - 2) + 1
            for k in range(P.Nz):
               kw   = (k-1) % (P.Nz - 2) + 1
               idx  = common.ijk2idx([i  , j  , k]  , P)
               idxw = common.ijk2idx([iw , jw , kw] , P)
               wrapper[idx,idxw] = 1
      self.wrapper = wrapper
               
   def update(self, dom):
      return self.wrapper*dom
