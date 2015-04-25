import numpy as np
import scipy as sp
import scipy.sparse
import common

class FTCS(object):
   def __init__(self,P):
      C = P.Alpha * P.dt / (P.dx * P.dx)
      offsets = np.array(np.mat("""0  0  0;
           -1  0  0; 0 -1  0; 0  0 -1;
            1  0  0; 0  1  0; 0  0  1"""))
      coeffs = [1 - 6 * C, C, C, C, C, C, C]
      data = np.zeros([7,P.N])
      for i in range(7):
         data[i,:] = coeffs[i]
      diags = [common.ijk2idx(offset, P) for offset in offsets]
      self.B = sp.sparse.spdiags(data, diags, P.N, P.N)

   def step(self, dom):
      dom = self.B*dom
      return dom
