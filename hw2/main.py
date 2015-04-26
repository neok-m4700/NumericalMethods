import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import common
from read_params import read_params
import conditions
import boundaries
import drivers
import sys

def main():
   P = read_params("params.txt")
   #Set number of steps in each direction plus boundaries
   P.Nx = int(P.Lx/P.dx) + 2
   P.Ny = int(P.Ly/P.dx) + 2
   P.Nz = int(P.Lz/P.dx) + 2
   #Set total number of steps
   P.N =  P.Nx * P.Ny * P.Nz
   #Set total number of time steps
   P.nSteps = int(P.T/P.dt)
   P.dims_size = [P.Nx, P.Ny, P.Nz]
   #Initialize Domain
   dom_initial = conditions.initial_domain(P)
   if len(sys.argv) > 1:       
       #Override params.txt if a driver is provided on the command line
       P.driver = sys.argv[1]
   driver = drivers.__dict__[P.driver](P)
   boundary = boundaries.__dict__[P.boundary](P)

   #Run it!
   dom_final = evolve(dom_initial, driver, boundary, P)

   #Plot it!
   if(P.plot):
       plot(dom_initial, dom_final, P)

def evolve(dom, driver, boundary, P):
   #MAIN LOOP
   for i in range(P.nSteps):
      dom = driver.step(dom) #step
      dom = boundary.update(dom) #boundary
      dom = dom + conditions.source(P, i) #source
   return dom 

def plot(dom1, dom2, P):
   fig = plt.figure()
   i_slice = P.Nx/2
   ijkvec = [common.idx2ijk(idx, P) for idx in range(P.N)]
   def plot_slice(dom, plotnum):
       ax  = fig.add_subplot(2,1,plotnum,projection="3d")
       ax.set_zlim3d(0,2)
       X = np.arange(1,P.Nx-1)
       Y = np.arange(1,P.Ny-1)
       X,Y = np.meshgrid(X,Y)
       indices = np.array(X*P.Ny*P.Nz + Y*P.Nz + P.Nz/2)
       Z = dom.T[0][indices]
       surf = ax.plot_surface(X*P.dx,Y*P.dx, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
               linewidth=0, antialiased=False)
       fig.colorbar(surf)

   plot_slice(dom1,1)
   plot_slice(dom2,2)
   plt.show()

if __name__ == "__main__":
    main()

