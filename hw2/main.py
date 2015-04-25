import numpy as np
import scipy as sp
import scipy.sparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import common
from read_params import read_params
import conditions
import boundaries
import drivers
import sys


def main():
   # Read parameters
   P = read_params("params.txt")
   # Initialze domain with inital conditions
   dom = conditions.initial_domain(P)
   dom_orig = np.array(dom)
   # Initialize driver
   if len(sys.argv) > 1:
      P.driver = sys.argv[1]
   driver = drivers.__dict__[P.driver](P)
   # Initialize Boundary conditions
   boundary = boundaries.__dict__[P.boundary](P)

   #evolve the system
   for i in range(P.nSteps):
      #step
      dom = driver.step(dom) 
      #apply boundary conditions
      dom = boundary.update(dom)
      #Add heat source
      #dom = dom + conditions.source(P, i)

   if(P.plot):
       #Initialize plot
       fig = plt.figure()
       ax  = fig.gca(projection = "3d")
       # Plot initial domain
       common.plot_slice(ax , 'g' , dom_orig , P)
       common.plot_slice(ax , 'r' , dom      , P)
       plt.show()


if __name__ == "__main__":
    main()
