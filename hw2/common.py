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

def wrap_around(coords, P):
    """Apply wrap around boundary conditions on a tuple of
    coordinates"""
    for d in range(len(P.dims_size)):
        coords[d] = coords[d] % P.dims_size[d]
    return coords

def is_inner(coords, P):
   """Return true if point is not on boundary"""
   dims_size = P.dims_size
   lower = all([coord > 0 for coord in coords]) 
   upper = all([coord < dim_size-1 for coord,dim_size in zip(coords, dims_size)])
   return lower and upper

def plot_slice(ax, color, dom, P, i_slice = None):
   if i_slice is None:
      i_slice = P.Nx/2
   dims_size = P.dims_size
   ijkvec = [idx2ijk(idx, P) for idx in range(P.N)]
   slice_i = [[j,k,d] for [i,j,k],d in zip(ijkvec, dom) if i == i_slice]
   slice_i = np.array(slice_i).T
   ax.scatter(*slice_i, c=color)
