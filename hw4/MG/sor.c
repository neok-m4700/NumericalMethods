#include "nrutil.h"
#include "mg.h"

void sor(double **xold, int n, double w)
{
    double **x = dmatrix(1,n,1,n);
    copy(x, xold, n);
    double a = 1/(1 + 4*C);
    double rs,r,b;

    for (int k = 0; k < 1e6; k++) { 
        rs = 0;
        for(int i = 2; i < n ; ++i) 
            for(int j = 2; j < n; ++j) {
                x[i][j] = w*a*((1 - 4*C)*xold[i][j]
                        + C*(xold[i+1][j] + xold[i-1][j] + xold[i][j+1] + xold[i][j-1] + xold[i][j] + xold[i][j]))
                        + (1-w)*x[i][j] + w*a*C*(x[i+1][j] + x[i-1][j] + x[i][j+1] + x[i][j-1] + x[i][j] + x[i][j]);
                //resid = b - Ax;
                b  = ((1 - 4*C)*xold[i][j] + C*(xold[i+1][j] + xold[i-1][j] + xold[i][j+1] + xold[i][j-1] + xold[i][j] + xold[i][j]));
                r = b - ((1 + 4*C)*x[i][j]    - C*(x[i+1][j] + x[i-1][j] + x[i][j+1] + x[i][j-1] + x[i][j] + x[i][j]));
            }
        if (r*r < 1e-10) break;
    }
    copy(xold, x, n);
}
