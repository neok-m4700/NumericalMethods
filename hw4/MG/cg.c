#include <stdlib.h>
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

void cg(double **x, int n) {
    double alpha;
    double beta;
    double r[n][n];
    double p[n][n];
    double Ap[n][n];
    double Ax[n][n];
    double aux[n][n];

    //r = b - Ax
    bi_diagonal_prod(-C, n, aux, x);
    linear_comb(N, r, 1, b, -1, Ax);

    //p = r
    copy_vector(N, p, r);

    double rsold, rsnew;
    rsold = dotprod(N, r, r);

    //main iterative loop of cg
    for (int k = 0; k < N; ++k) {
        //Ap = A*p
        bi_diagonal_prod(-C, n, Ap, p);

        double pAp = dotprod(N, p, Ap);
        alpha = rsold / pAp;

        linear_comb(N, x, 1.0, x, alpha, p);
        linear_comb(N, r, 1.0, r, -alpha, Ap);

        rsnew = dotprod(N, r, r);
        if (rsnew < 1e-10) return;

        beta = rsnew / rsold;
        linear_comb(N, p, 1.0, r, beta, p);
        rsold = rsnew;
    }
}
