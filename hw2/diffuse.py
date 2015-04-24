from read_params import read_params
from genMatrix import genMatrix
import scipy as sp
import scipy.linalg as linalg

import crank_nicholson

def evolve():
   P = read_params("params.txt")
   A_off_coeffs = crank_nicholson.A_off_coeffs
   B_off_coeffs = crank_nicholson.B_off_coeffs
   dims_size = [P["Nx"], P["Ny"], P["Nz"]]
   A = genMatrix(A_off_coeffs, np.array(dims_size))
   B = genMatrix(B_off_coeffs, np.array(dims_size))
   LU = linalg.lu_factor(A)






