# Implement the 3 drivers: FTCS, CN, and ADI

import numpy as np
import scipy as sp
import scipy.sparse
import scipy.linalg
import scipy.sparse.linalg
from params import Params
import solvers


class CN(object):
   """Crank Nicholson Method"""

   def __init__(self, solver):
      P = Params.instance
      C = P.Alpha * P.dt / (2 * P.dx * P.dx)
      coeffs = [1 - 6 * C, C, C, C, C, C, C]
      data = np.tile(np.reshape(coeffs, [7, 1]), [1, P.N])
      diags = [0, P.Ny * P.Nz, P.Nz, 1, -P.Ny * P.Nz, -P.Nz, -1]
      self.B = sp.sparse.spdiags(data, diags, P.N, P.N)

      coeffs = [1 + 6 * C, -C, -C, -C, -C, -C, -C]
      data = np.tile(np.reshape(coeffs, [7, 1]), [1, P.N])
      diags = [0, P.Ny * P.Nz, P.Nz, 1, -P.Ny * P.Nz, -P.Nz, -1]
      A = sp.sparse.spdiags(data, diags, P.N, P.N).toarray()

      # Assign solver
      self.solver = solver(A)

   def step(self, dom):
      return self.solver.solve(self.B * dom)

