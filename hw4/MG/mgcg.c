#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"
#include "mg.h"

void mgcg(double **xold, double **x, int n) {
    double alpha, beta, pr, pAp, rAp;
    double **r = dmatrix(1,n,1,n);
    double **p = dmatrix(1,n,1,n);
    double **Ap = dmatrix(1,n,1,n);

    //r = b - Ax
    for(int i = 2; i < n; ++i) 
        for(int j = 2; j < n; ++j) {
            double b = (1-4*C)*xold[i][j] + 
                C*(xold[i-1][j] + xold[i+1][j] + xold[i][j-1] + xold[i][j+1]);
            r[i][j] = b - (1+4*C)*x[i][j] + 
                C*(x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1]);
        }

    //boundaries of r
    for(int i = 1; i <= n; ++i)
        r[i][1] = r[i][n] = r[1][i] = r[n][i] = 0;

    //precondition with mg
    mglin(r, n, 1); 
    //p = r
    copy(p, r, n);

    //main iterative loop of cg
    for (int k = 0; k < n; ++k) {

        //Ap = A*p
        for(int i = 2; i < n; ++i) 
            for(int j = 2; j < n; ++j) 
                Ap[i][j] = (1+4*C)*p[i][j] 
                    - C*(p[i-1][j] + p[i+1][j] + p[i][j-1] + p[i][j+1]);
        for(int i = 1; i <= n; ++i)
            Ap[i][1] = Ap[i][n] = Ap[1][i] = Ap[n][i] = 0;

        pr = dotprod(n, p, r);
        pAp = dotprod(n, p, Ap);
        alpha = pr / pAp;

        for(int i = 2; i < n; ++i) 
            for(int j = 2; j < n; ++j) {
                x[i][j] += alpha*p[i][j];
                r[i][j] -= alpha*Ap[i][j];
            }

        double rsnew = dotprod(n, r, r);
        if (rsnew < 1e-10) break;

        //precondition with mg
        mglin(r, n, 1);
        rAp = dotprod(n, r, Ap);
        beta = rAp/pAp;

        //p = r + Bp
        for(int i = 2; i < n; ++i) 
            for(int j = 2; j < n; ++j) 
                p[i][j] = r[i][j] + beta*p[i][j];
    }

    copy(xold, x, n); //copy solution into xold (implicit return)

    free_dmatrix(r, 1, n, 1, n);
    free_dmatrix(p, 1, n, 1, n);
    free_dmatrix(Ap, 1, n, 1, n);
}

