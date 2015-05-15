import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import argparse

import common
from params import Params
import conditions
import drivers
import boundaries
import solvers

import sys
import time


def main():
    # Read Params
    parser = argparse.ArgumentParser()
    parser.add_argument("-s" , "--solver", default="Jacobi", help="one of Jacobi, Gauss, SOR")
    parser.add_argument("-w" , "--omega" , default="1.0", help="w for SOR")
    parser.add_argument("-x" , "--length" , default="1.0", help="Length of one side of domain")
    args = parser.parse_args()
    P = Params("params.txt")
    P.solver = args.solver
    P.Lx = P.Ly = P.Lz = float(args.length)
    P.omega = args.omega
    P.set_dependent()

    # Initialize Domain
    dom_initial = conditions.initial_domain()

    # Assign Solver
    solver = {'Jacobi': solvers.Jacobi, 'Gauss': solvers.Gauss, 'SOR': solvers.SOR}[P.solver]

    # Assign Driver
    driver = drivers.CN(solver)

    # Assign Boundary
    boundary = boundaries.Dirichlet()

    # Run it!
    tic = time.clock()
    dom_final, meaniters = evolve(dom_initial, driver, boundary)
    print (time.clock() - tic) / P.nSteps, meaniters

    # Plot it!
    if (P.plot):
        plot(dom_initial, dom_final)


def evolve(dom, driver, boundary):
    # MAIN LOOP
    P = Params.instance
    sumiters = 0
    for i in range(P.nSteps):
        dom, iters = driver.step(dom)  # step
        sumiters += iters
        dom = boundary.update(dom)
    return dom, sumiters/P.nSteps


def plot(dom1, dom2):
    P = Params.instance
    fig = plt.figure()
    i_slice = P.Nx / 2
    ijkvec = [common.idx2ijk(idx) for idx in range(P.N)]

    def plot_slice(dom, plotnum):
        ax = fig.add_subplot(2, 1, plotnum, projection="3d")
        ax.set_zlim3d(0, 2)
        X = np.arange(1, P.Nx - 1)
        Y = np.arange(1, P.Ny - 1)
        X, Y = np.meshgrid(X, Y)
        indices = np.array(X * P.Ny * P.Nz + Y * P.Nz + P.Nz / 2)
        Z = dom.flat[indices]
        surf = ax.plot_surface(X * P.dx, Y * P.dx, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        fig.colorbar(surf)

    plot_slice(dom1, 1)
    plot_slice(dom2, 2)
    plt.show()


if __name__ == "__main__":
    main()

