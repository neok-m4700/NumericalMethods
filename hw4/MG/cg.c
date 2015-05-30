#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"
#include "mg.h"

double dotprod(int n, double **r1, double **r2) {
    double rdot = 0;
    for (int i = 1; i <= n; ++i)
       for(int j = 1; j <= n; ++j) 
            rdot += r1[i][j] * r2[i][j];
    return rdot;
}

void linear_comb(int n, double **r1, double a, double **r2, double b, double **r3){
    //r1 = a*r2 + b*r3
    for (int i = 1; i <= n; ++i) 
       for(int j = 1; j <= n; ++j) 
          r1[i][j] = a*r2[i][j] + b*r3[i][j];
}

void cg(double **b, double **x, int n) {
    double alpha;
    double beta;
    double **r = dmatrix(1,n,1,n);
    double **p = dmatrix(1,n,1,n);
    double **Ap = dmatrix(1,n,1,n);
    double **Ax = dmatrix(1,n,1,n);
    double **aux = dmatrix(1,n,1,n);

    //r = b - Ax
    copy(Ax, x, n);
    bi_diagonal_prod(-C, n, aux, Ax); //multiplies argument by A
    linear_comb(n, r, 1, b, -1, Ax);

    //p = r
    copy(p, r, n);

    double rsold, rsnew;
    rsold = dotprod(n, r, r);

    //main iterative loop of cg
    int k;
    for (k = 0; k < n; ++k) {
        //Ap = A*p
        copy(Ap, p, n);
        bi_diagonal_prod(-C, n, aux, Ap);

        double pAp = dotprod(n, p, Ap);
        alpha = rsold / pAp;

        linear_comb(n, x, 1.0, x, alpha, p);
        linear_comb(n, r, 1.0, r, -alpha, Ap);

        rsnew = dotprod(n, r, r);
        if (rsnew < 1e-10) break;

        beta = rsnew / rsold;
        linear_comb(n, p, 1.0, r, beta, p);
        rsold = rsnew;
    }

    printf("%d\n", k);
    copy(b, x, n); //copy solution into b (implicit return)

    free_dmatrix(r, 1, n, 1, n);
    free_dmatrix(p, 1, n, 1, n);
    free_dmatrix(Ap, 1, n, 1, n);
    free_dmatrix(Ax, 1, n, 1, n);
    free_dmatrix(aux, 1, n, 1, n);
}
