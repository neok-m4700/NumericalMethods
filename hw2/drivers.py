# Implement the 3 drivers: FTCS, CN, and ADI

import numpy as np
import scipy as sp
import scipy.sparse
import scipy.linalg
import scipy.sparse.linalg
from params import Params

class FTCS(object):
   """Forward time Center Space Method"""
   def __init__(self):
      P = Params.instance
      C = P.Alpha * P.dt / (P.dx * P.dx)
      coeffs = [1 - 6 * C, C, C, C, C, C, C]
      data = np.tile(np.reshape(coeffs,[7,1]),[1,P.N])
      diags = [0, P.Ny*P.Nz, P.Nz, 1, -P.Ny*P.Nz, -P.Nz, -1]
      self.B = sp.sparse.spdiags(data, diags, P.N, P.N)

   def step(self, dom):
      return self.B * dom


class CN(object):
   """Crank Nicholson Method"""
   def __init__(self):
      P = Params.instance
      C = P.Alpha * P.dt / (2*P.dx * P.dx)
      coeffs = [1 - 6 * C, C, C, C, C, C, C]
      data = np.tile(np.reshape(coeffs,[7,1]),[1,P.N])
      diags = [0, P.Ny*P.Nz, P.Nz, 1, -P.Ny*P.Nz, -P.Nz, -1]
      self.B = sp.sparse.spdiags(data, diags, P.N, P.N)

      coeffs = [1 + 6 * C, -C, -C, -C, -C, -C, -C]
      data = np.tile(np.reshape(coeffs,[7,1]),[1,P.N])
      diags = [0, P.Ny*P.Nz, P.Nz, 1, -P.Ny*P.Nz, -P.Nz, -1]
      A = sp.sparse.spdiags(data, diags, P.N, P.N).toarray()
      self.LU = sp.linalg.lu_factor(A)

   def step(self, dom):
      return sp.linalg.lu_solve(self.LU, self.B*dom)


class ADI(object):
   """Aternate Directions Implicit method"""
   def __init__(self):
      P = Params.instance
      C = P.Alpha*P.dt/(3*(P.dx*P.dx))
      data = np.tile(np.array([-2*C, C, C]).reshape(3,1), [1, P.N])
   
      Ai = sp.sparse.spdiags(data , [0 , -P.Ny*P.Nz , P.Ny*P.Nz ] , P.N , P.N)
      Aj = sp.sparse.spdiags(data , [0 , -P.Nz      , P.Nz      ] , P.N , P.N)
      Ak = sp.sparse.spdiags(data , [0 , -1         , 1         ] , P.N , P.N)
      I  = sp.sparse.eye(P.N)
   
      # Decomposing left side matrices
      self.LU1 = sp.linalg.lu_factor((I-Ai).toarray())
      self.LU2 = sp.linalg.lu_factor((I-Aj).toarray())
      self.LU3 = sp.linalg.lu_factor((I-Ak).toarray())
   
      # Right side matrics
      self.B1 = I + Aj + Ak
      self.B2 = I + Ak + Ai
      self.B3 = I + Ai + Aj
   
   def step(self, dom):
      dom = sp.linalg.lu_solve(self.LU1, self.B1*dom)
      dom = sp.linalg.lu_solve(self.LU2, self.B2*dom)
      dom = sp.linalg.lu_solve(self.LU3, self.B3*dom)
      return dom

