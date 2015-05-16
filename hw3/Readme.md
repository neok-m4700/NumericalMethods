# 3D Diffusion Solver

## Dependencies
1. Python 2.7
2. Numpy
3. Scipy
4. Matplotlib (for plotting)

## Configuration
Edit params.txt and set values as desired.

## Execution
$ python main.py
optional arguments:
  -h, --help            show this help message and exit
  -s SOLVER, --solver SOLVER
                        one of Jacobi, Gauss, SOR
  -w OMEGA, --omega OMEGA
                        w for SOR
  -x LENGTH, --length LENGTH
                        Length of one side of domain

## Results
Data can be collected by running run.sh. 
Turn plotting off in 'params.txt' if using this script to run many simulations
Data is already stored for a run in data/
You can find plot of performance in data/performance.png
You can find a plot of the number of iterations for each method in data/iterations.png
You can generate this file given the gnpuplot script in data/

The "length per side" is indicated in metres. dx is set to 0.1 metres by default. So a lenght per side of "2" indicates 20 slices per side. 

Performance of all 3 methos is comparable, with Jacobi and Gauss Seidel both convergnig in around 3 iterations. With SOR the number of iterations is very variable, ranging between 4 and 10. I tried with various values of w but there is no universal values for which it gives a best performance, it changes with size of domain. I settled on 1.0 as w for a uniform comparison. 

The criterion for for convergence is the absolute mean error being less than 1E-6. 

## Cubic Spline

The derivation of the cubic spline is based on the provided reference with the "Natural Splines" method. You can find the plot in Cubic/. The number of points plottet can be changed by setting the argument to the cubic(n) function. 
