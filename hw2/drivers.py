import numpy as np
import scipy as sp
import scipy.sparse
import scipy.linalg
import scipy.sparse.linalg
import common

class FTCS(object):
   """Forward time Center Space Method"""
   def __init__(self,P):
      C = P.Alpha * P.dt / (P.dx * P.dx)
      coeffs = [1 - 6 * C, C, C, C, C, C, C]
      data = np.tile(np.reshape(coeffs,[7,1]),[1,P.N])
      diags = [0, P.Ny*P.Nz, P.Nz, 1, -P.Ny*P.Nz, -P.Nz, -1]
      self.B = sp.sparse.spdiags(data, diags, P.N, P.N)

   def step(self, dom):
      dom = self.B*dom
      return dom


class CN(object):
   """Crank Nicholson Method"""
   def __init__(self,P):
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
   def __init__(self, P):
      C = P.Alpha*P.dt/(3*(P.dx*P.dx))
      data = np.ones([3, P.N])*C
   
      Ai = sp.sparse.spdiags(data , [0 , -P.Ny*P.Nz , P.Ny*P.Nz ] , P.N , P.N)
      Aj = sp.sparse.spdiags(data , [0 , -P.Nz      , P.Nz      ] , P.N , P.N)
      Ak = sp.sparse.spdiags(data , [0 , -1         , 1         ] , P.N , P.N)
      I  = sp.sparse.eye(P.N)
   
      # Decomposing left side matrices
      self.LUi = sp.linalg.lu_factor((I-Ai).toarray())
      self.LUj = sp.linalg.lu_factor((I-Aj).toarray())
      self.LUk = sp.linalg.lu_factor((I-Ak).toarray())
   
      # Right side matrics
      self.Bi = I + Aj + Ak
      self.Bj = I + Ak + Ai
      self.Bk = I + Ai + Aj
   
   def step(self, dom):
      dom = sp.linalg.lu_solve(self.LUi, self.Bi*dom)
      dom = sp.linalg.lu_solve(self.LUj, self.Bj*dom)
      dom = sp.linalg.lu_solve(self.LUk, self.Bk*dom)
      return dom

