#Compilation

Just run make

#Execution

To change parameters, you can edit main.c, OR you can pass them on the command line:
$ main.out -p v
where p is the parameter symbol and v is it's value
The following symbols are accepted:

 'n': size of matrix is n x n 
 'N': Number of time steps
 'd': dx = dy
 't': dt
 'a': diffusivity
 'A': amplitude of Gaussian
 's': sigma of Gaussian
 'r': randomization
 'S': source term
 'c': ncycles for MG
 'm': method. 0 = MG, 1 = CG, 2 = MGCG, 3 = SOR

#Plotting

The results of a simulation are saved in out.data. You can plot this with running "gnuplot plot.gnuplot". Then the plot will be saved in out.png

#Results

//dx = 1.0; 
//dt = 7.0; alpha = 0.1; A = 100; sig = 5; r = 10; source = 0.001; ncycles = 2;
Method   n     time
MG       33    0.029101
MG       65    0.126472
MG       129   0.487210
MG       257   2.199205
MG       513   11.978596
CG       33    0.029336
CG       65    0.123792
CG       129   0.529666
CG       257   2.186132
CG       513   9.673456
MG+CG    33    0.146938
MG+CG    65    0.572467
MG+CG    129   2.548330
MG+CG    257   11.479111
MG+CG    513   62.704048
SOR      33    0.074579
SOR      65    0.289886
SOR      129   1.147321
SOR      257   4.317905
SOR      513   16.954153

#Comments:
Of MG, CG and SOR, CG ranks the best, followed by MG and SOR. MG+CG took the most amount of time, probably because there is many cycles of MG running (one for each iteratin of CG). 

We can study of effects of ncycles on MG by running just MG with different ncycles:

//dx = 1.0; dt = 7.0; alpha = 0.1; A = 100; sig = 5; r = 10; source = 0.001;
Method   n     time       ncycles
MG       33    0.029958   2
MG       65    0.124927   2
MG       129   0.520377   2
MG       33    0.059304   4
MG       65    0.229842   4
MG       129   0.990300   4
MG       33    0.110204   8
MG       65    0.452684   8
MG       129   1.899379   8
MG       33    0.208820   16
MG       65    0.876607   16
MG       129   3.672616   16
