#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

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
    double boundary; //boundary value
} Opt;


void gaussian(Opt *o, double *b);
void save_output(Opt *o, double *b);
void noise(Opt *o, double *b);
void cg(Opt *o, double *x, double C, double *b); //solves Ax = B
double dotprod(int N, const double *r1, const double *r2);
void copy_vector(int N, double *r1, double *r2) ;
void linear_comb(int N, double *r1, double a, double *r2, double b, double *r3);
void bi_diagonal_prod(double C, int l, int h, int n, double *product, const double *multiplicant);

int main(int argc, char **argv) {
    int c;
    Opt *o = malloc(sizeof(struct _Opt_));

    //specify default values
    o->n = 20;
    o->N = 5000;
    o->dx = 0.1;
    o->dt = 0.001;
    o->alpha = 1;
    o->A = 2;
    o->sig = 1;
    o->r = 0.1;
    o->source = 0.001;
    o->boundary = 0;

    //get command line values
    while ((c = getopt (argc, argv, "n:d:N:t:a:A:s:r:S:b:")) != -1)
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
            case 'b': o->boundary = atof(optarg); break;
            default: abort();
        }

    int n = o->n;

    double *dom = malloc(sizeof(double) * n*n); //domain vector
    double *x   = malloc(sizeof(double) * n*n); //solution vector
    double *b   = malloc(sizeof(double) * n*n); //rhs vector

    //initialize
    for (int i = 0; i < n*n ; ++i) {
        dom[i] = o->boundary;
        b[i]   = o->boundary;
        x[i]   = o->boundary;
    }

    //initialize gaussian
    gaussian(o, dom);

    //add noise
    noise(o, dom);

    //caculate C
    double C = o->alpha * o->dt / (2.0 * o->dx * o->dx);

    for (int k = 0; k < o->N; ++k) { //time stepping

        //calculate B*dom
        bi_diagonal_prod(-C, 1, n-1, n, b, dom);

        //solve Ax = b
        cg(o, x, C, b);

        //add source and reset boundaries
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j){
                x[n*i + j] += o->source;
                if (i==0 || i==n-1 || j==0 || j==n-1)
                    x[n*i + j] = o->boundary;
            }

        //swap pointers
        double *tmp = dom;
        dom = x;
        x = tmp;
    }

    //save the domain
    save_output(o, dom);

    //free stuff
    free(dom); free(x); free(b); free(o);
    //return
    return 0;
}

void cg(Opt *o, double *x, double C, double *b) {
    int n = o->n;
    int N = n*n;
    double alpha;
    double beta;
    double r[N];
    double p[N];
    double Ap[N];
    double Ax[N];

    //r = b - Ax
    bi_diagonal_prod(C, 0, n, n, Ax, x);
    linear_comb(N, r, 1, b, -1, Ax);

    //p = r
    copy_vector(N, p, r);

    double rsold, rsnew;
    rsold = dotprod(N, r, r);

    //main iterative loop of cg
    for (int k = 0; k < N; ++k) {
        //Ap = A*p
        bi_diagonal_prod(C, 0, n, n, Ap, p);

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


void bi_diagonal_prod(double C, int l, int h, int n, double *product, const double *multiplicant) {
    for (int i = l; i < h; ++i) { //iterating over domain
        for (int j = l; j < h; ++j) {
            int idx = i*n + j;
            product[idx] =         (1.0 + 4*C) * multiplicant[idx]
                         - ( i-1 == -1 ? 0 : C * multiplicant[idx - n])
                         - ( i+1 ==  n ? 0 : C * multiplicant[idx + n])
                         - ( j-1 == -1 ? 0 : C * multiplicant[idx - 1])
                         - ( j+1 ==  n ? 0 : C * multiplicant[idx + 1]);
        }
    }
}

double dotprod(int N, const double *r1, const double *r2) {
    double rdot = 0;
    for (int i = 0; i < N; ++i)
            rdot += r1[i] * r2[i];
    return rdot;
}

void copy_vector(int N, double *r1, double *r2) {
    for (int i = 0; i < N; ++i) {
        r1[i] = r2[i];
    }
}

void linear_comb(int N, double *r1, double a, double *r2, double b, double *r3){
    //r1 = a*r2 + b*r3
    for (int i = 0; i < N; ++i) {
       r1[i] = a*r2[i] + b*r3[i];
    }
}

void noise(Opt *o, double *b) {
    for (int i = 1; i < o->n -1; ++i) {
        for (int j = 1; j < o->n - 1; ++j) {
            b[i*o->n + j] += random()*o->r/RAND_MAX;
        }
    }

}


void gaussian(Opt *o, double *b){
    double x0 = o->dx * (o->n) / 2.0;
    double y0 = o->dx * (o->n) / 2.0;
    double sig2 = 2*o->sig*o->sig;
    double x,y,x2,y2;

    for (int i = 1; i < o->n-1; ++i) {
        x = i * o->dx - x0;
        x2 = x*x/sig2;
        for (int j = 1; j < o->n-1; ++j) {
            y = j * o->dx - y0;
            y2 = y*y/sig2;
            b[i*o->n + j] = o->A * exp(-(x2 + y2));
        }
    }
    return;
}


void save_output(Opt *o, double *b) {
    FILE *f = fopen("out.data", "w");
    for (int i = 0; i < o->n; ++i) {
        for (int j = 0; j < o->n; ++j) {
            fprintf(f, "%lf %lf %lf\n", i*o->dx, j*o->dx, b[i*o->n + j]);
        }
    }
    fclose(f);
}


