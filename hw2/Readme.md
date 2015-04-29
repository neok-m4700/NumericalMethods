# 3D Diffusion Solver

## Dependencies
1. Python 2.7
2. Numpy
3. Scipy
4. Matplotlib (for plotting)

## Configuration
Edit params.txt and set values as desired.

## Execution
    $python main.py [method] [domain length]
If the 'method' and 'domain length' parameters are provided on the command line then the values of the same set in parmas.txt will be overridden.

## Results
Data can be collected by running run.sh. 
Turn plotting off in 'params.txt' if using this script to run many simulations
Data is already stored for a run in data/
You can find plotter performance in data/performance.png
You can generate this file given the gnpuplot script in data/

## Python Files
* main.py - Main executable. Puts all the pieces together
* params.py - Read the params file and export the Params object
* drivers.py - The code for the three drivers - FTCS, CN and ADI
* boundaries.py - The code for the two boundary conditions - Wraparound and Dirichlet
* common.py - Common functions used by more than one module
* conditions.py - Initial and Boundary condiditions specified here

