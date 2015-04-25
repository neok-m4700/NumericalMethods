import numpy as np
import common

class Wraparound(object):
   """Implements the Wraparound boundary conditions"""

   def __init__(self, P):
      # Initialize dirichlet mask
      self.P = P
      wrapper = sp.sparse.eye(P.N).tolil()
      for idx in range(P.N):
         if not common.is_inner(common.idx2ijk(idx, P)):
            wrapper[idx,:] = 0
            wrapper[idx, self.wrap(idx)] = 1
      self.wrapper = sp.sparse.csr(wrapper)

   def wrap(self, idx):
      ijk = common.idx2ijk(idx,self.P)
      ijk_wrap = [(l-1) % (Nl - 2) + 1 for l,Nl in zip(ijk, P.dims_size)]
      return common.ijk2idx(ijk_wrap, self.P)

   def update(self, dom):
      return self.wrapper*dom
