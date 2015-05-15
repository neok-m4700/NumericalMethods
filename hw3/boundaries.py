# Defines the various boundary condition implementations
import numpy as np
import scipy as sp
import scipy.sparse
import common
import conditions
from params import Params

class Dirichlet(object):

   def __init__(self):
      # Initialize dirichlet mask
      self.mask = self.gen_mask()  # zero if boundary, 1 otherwise

   def gen_mask(self):
      """Mask the boundary values on the domain vector
      value is zero if the point is on the boudary"""
      P = Params.instance
      mask = np.zeros([P.N, 1])
      for idx in range(P.N):
         ijk = common.idx2ijk(idx)
         if common.is_inner(ijk):
            mask[idx] = 1
      return mask

   def update(self, dom):
      return np.multiply(dom, self.mask)


