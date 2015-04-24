#GENMATRIX Generate matrix for Heat equation time step 
#  Given a relationship like X(i,j,k) = SumOverL[Y(i+d_li, j+d_lj, k+d_lk)*C_l]
# L is number of terms
# This will return a matrix A, such that AY = X
# Y and X are vectors, that range from 000 to IJK
# IJK = dims_size(1) * dims_size(2) * dims_size(3)
# off_coeff stores the d's and the C's as [ di dj dk C]
# off_coeff has L rows

import numpy as np
import scipy as sp
import scipy.sparse

def genMatrix(off_coeffs, dims_size):
   nterms = np.shape(off_coeffs)[0]
   N = np.prod(dims_size) # n_x * n_y * n_z
   B = np.zeros([nterms, N]) # columns of B are diagonals of A
   d = np.zeros(nterms) # specifies distance to central diaganoal 

   for t in range(nterms):
      off_coeff = off_coeffs[t,:]
      coeff = off_coeff[-1]
      B[t,:] = coeff
      offset_nd = off_coeff[:-1]
      offset_1d = ijk2idx(offset_nd, dims_size)
      d[t] = offset_1d

   return sp.sparse.spdiags(B, d, N, N)

def ijk2idx(indices, dims_size):
   idx = 0
   for d in range(np.size(dims_size)):
      idx = idx*dims_size[d] + indices[d]
   return idx

if __name__ == "__main__":
   off = np.vstack([np.zeros(3), np.eye(3)*(-1), np.eye(3)])
   coeff = np.array([[1, 2, 3, 4, -2, -3, -4]]).reshape((7,1))
   off_coeff = np.hstack([off, coeff])
   dims_size = np.array([2,2,3])
   print genMatrix(off_coeff, dims_size).toarray()
