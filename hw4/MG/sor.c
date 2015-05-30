#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "mg.h"

void sor(double **xold, int n, double w)
{
    double **x = dmatrix(1,n,1,n);
    copy(x, xold, n);
    double a = 1.0/(1 + 4*C);
    double rs,xlast;
    int k;
    for (k = 0; k < 1e3; k++) { 
        rs = 0;
        for(int i = 2; i < n ; ++i) 
            for(int j = 2; j < n; ++j) {
                xlast = x[i][j];
                x[i][j] = w*a*((1 - 4*C)*xold[i][j]
                        + C*(xold[i+1][j] + xold[i-1][j] + xold[i][j+1] + xold[i][j-1]))
                        + (1-w)*x[i][j] + w*a*C*(x[i+1][j] + x[i-1][j] + x[i][j+1] + x[i][j-1]) ;
                rs += fabs(xlast - x[i][j]);
            }
        if (rs/(n*n) < 1e-10) break;
    }
    copy(xold, x, n);
    free_dmatrix(x, 1, n, 1, n);
}
