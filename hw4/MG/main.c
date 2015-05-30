#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "mg.h"
#include "nrutil.h"


typedef struct _Opt_ { //Options
    int n; //number of steps in a direction
    int N; //number of time steps
    double dx; //size of step in a direction
    double dt; //time step
    double alpha; //diffusivity
    double A; // Gaussian A
    double sig; // Gaussian sigma
    double r; //noise scaling
    double source; //source
} Opt;

void gaussian(Opt *o, double **b);
void save_output(Opt *o, double **b);
void noise(Opt *o, double **b);
void bi_diagonal_prod(double C, int n, double **aux, double **dom);

int main(int argc, char **argv) {
    int c;
    Opt *o = malloc(sizeof(struct _Opt_));
    int ncycles; //for mutligrid

    //specify default values
    o->n = 129;
    o->N = 100;
    o->dx = 1.0;
    o->dt = 5.0;
    o->alpha = 0.1;
    o->A = 100;
    o->sig = 5;
    o->r = 0.0;
    o->source = 0.000;
    ncycles = 1;

    //get command line values
    while ((c = getopt (argc, argv, "n:d:N:t:a:A:s:r:S:c:")) != -1)
        switch(c)
        {
            case 'n': o->n = atoi(optarg); break;
            case 'N': o->N = atoi(optarg); break;
            case 'd': o->dx = atof(optarg); break;
            case 't': o->dt = atof(optarg); break;
            case 'a': o->alpha = atof(optarg); break;
            case 'A': o->A = atof(optarg); break;
            case 's': o->sig = atof(optarg); break;
            case 'r': o->r = atof(optarg); break;
            case 'S': o->source = atof(optarg); break;
            case 'c': ncycles = atof(optarg); break;
            default: abort();
        }

    int n = o->n;

    double **dom = dmatrix(1, n, 1, n);
    double **x = dmatrix(1,n,1,n);
    double **b = dmatrix(1,n,1,n);
    double **aux = dmatrix(1,n,1,n);

    //initialize gaussian
    gaussian(o, dom);

    //add noise
    noise(o, dom);

    //caculate C
    C = o->alpha * o->dt / (2.0 * o->dx * o->dx);

    for (int k = 0; k < o->N; ++k) { //time stepping

        //calculate B*dom
        bi_diagonal_prod(C, n,aux, dom);

        //solve Ax = b
        mglin(b, n, ncycles);

        //add source and reset boundaries
        for (int i = 2; i < n; ++i)
            for (int j = 2; j < n; ++j)
                dom[i][j] += o->source;

    }

    //save the domain
    save_output(o, dom);

    //free stuff
    free_dmatrix(dom, 1, n, 1, n);
    free_dmatrix(x, 1, n, 1, n);
    free_dmatrix(b, 1, n, 1, n);
    free_dmatrix(aux, 1, n, 1, n);
    free(o);
    //return
    return 0;
}



void noise(Opt *o, double **b) {
   int n = o->n;
    for (int i = 2; i <= n-1; ++i) {
        for (int j = 2; j <= n-1; ++j) {
            b[i][j] += random()*o->r/RAND_MAX;
        }
    }

}


void gaussian(Opt *o, double **b){
    int n = o->n;
    double x0 = o->dx * n / 2.0;
    double y0 = o->dx * n / 2.0;
    double sig2 = 2*o->sig*o->sig;
    double x,y,x2,y2;


    for(int i = 1; i <= n; ++i) {
        b[i][1] = 0;
        b[i][n] = 0;
        b[1][i] = 0;
        b[n][i] = 0;
    }

    for (int i = 2; i <= n-1; ++i) {
        x = i * o->dx - x0;
        x2 = x*x/sig2;
        for (int j = 2; j <= n-1; ++j) {
            y = j * o->dx - y0;
            y2 = y*y/sig2;
            b[i][j] = o->A * exp(-(x2 + y2));
        }
    }
    return;
}


void save_output(Opt *o, double **b) {
   int n = o->n;
    FILE *f = fopen("out.data", "w");
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            fprintf(f, "%lf %lf %lf\n", i*o->dx, j*o->dx, b[i][j]);
        }
    }
    fclose(f);
}

void bi_diagonal_prod(double C, int n, double **aux, double **dom) {

   copy(aux, dom, n);

    //boundary condition = 0
    for(int i = 1; i <= n; ++i) {
      dom[i][1] = 0;
      dom[i][n] = 0;
      dom[1][i] = 0;
      dom[n][i] = 0;
    }
    for (int i = 2; i <= n-1; ++i) { //iterating over domain
        for (int j = 2; j <= n-1; ++j) {
            dom[i][j] = (1.0 - 4*C) * aux[i][j]
                         +  C * aux[i-1][j]
                         +  C * aux[i+1][j]
                         +  C * aux[i][j-1]
                         +  C * aux[i][j+1];
        }
    }

}

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
