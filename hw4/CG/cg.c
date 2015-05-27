#include <stdlib.h>
#include "common.h"


void solve(int n, double *x, double C, double *b) {
    int N = n*n;
    double alpha;
    double beta;
    double r[N];
    double p[N];
    double Ap[N];
    double Ax[N];

    //r = b - Ax
    bi_diagonal_prod(-C, n, Ax, x);
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
