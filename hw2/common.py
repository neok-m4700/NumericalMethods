# Common functions
import numpy as np

def ijk2idx(coords, P):
    """Given the spatial coordinates of a point (i,j,k) 
    and given the size of the domain (Nx, Ny, Nz), 
    return a single scalar value representative of the
    point."""
    idx = 0
    for d in range(len(P.dims_size)):
       idx = idx*P.dims_size[d] + coords[d]
    return idx

def idx2ijk(idx, P):
    """Inverse of ijk2idx"""
    l = len(P.dims_size)
    coords = np.zeros(l)
    for d in range(l):
        coords[l-d-1] = idx % P.dims_size[l-d-1]
        idx = idx / P.dims_size[l-d-1]
    return coords

def is_inner(coords, P):
   """Return true if point is not on boundary"""
   dims_size = P.dims_size
   lower = all([coord > 0 for coord in coords]) 
   upper = all([coord < dim_size-1 for coord,dim_size in zip(coords, dims_size)])
   return lower and upper
