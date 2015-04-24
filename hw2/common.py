import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def ijk2idx(coords, dims_size):
    """Given the spatial coordinates of a point (i,j,k) 
    and given the size of the domain (Nx, Ny, Nz), 
    return a single scalar value representative of the
    point."""
    idx = 0
    for d in range(len(dims_size)):
       idx = idx*dims_size[d] + coords[d]
    return idx

def idx2ijk(idx, dims_size):
    """Inverse of ijk2idx"""
    l = len(dims_size)
    coords = np.zeros(l)
    for d in range(l):
        coords[l-d-1] = idx % dims_size[l-d-1]
        idx = idx / dims_size[l-d-1]
    return coords

def wrap_around(coords, dims_size):
    """Apply wrap around boundary conditions on a tuple of
    coordinates"""
    for d in range(len(dims_size)):
        coords[d] = coords[d] % dims_size[d]
    return coords

def in_limits(coords, dims_size):
   """Return true if all coordinates respect spatial boundaries"""
   lower = all([coord >= 0 for coord in coords]) 
   upper = all([coord < dim_size for coord,dim_size in zip(coords, dims_size)])
   return lower and upper

def gen_matrix(B, offsets, coeffs, P, wrap = False):
   """Construct sparse coefficient matrix using offsets
   and wrap around boundary conditions"""
   dims_size = np.array([P.Nx, P.Ny, P.Nz])
   for row_idx in range(P.N):
      ijk = idx2ijk(row_idx, dims_size)
      for term in range(len(offsets)):
         ijk_offset = ijk + offsets[term] 
         if wrap: #apply wrap_around bondary conditions
            ijk_offset = wrap_around(ijk_offset, dims_size)
         else: #ijk + offset should not exeed spatial boundaries
            assert in_limits(ijk_offset, dims_size)
         col_idx = ijk2idx(ijk_offset, dims_size)
         B[row_idx, col_idx] = coeffs[term]
   return B

def plot_slice(ax, color, dom, P, i_slice = 0):
   dims_size = np.array([P.Nx, P.Ny, P.Nz])
   ijkvec = [idx2ijk(idx, dims_size) for idx in range(P.N)]
   slice_i = [[j,k,d] for [i,j,k],d in zip(ijkvec, dom) if i == i_slice]
   slice_i = np.array(slice_i).T
   ax.scatter(*slice_i, c=color)
