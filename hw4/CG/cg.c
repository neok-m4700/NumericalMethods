#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

typedef struct _Opt_ {
    int n; //number of steps in a direction
    int N; //number of time steps
    double dx; //size of step in a direction
    double dt; //time step
    double alpha; //diffusivity
    double A; // Gaussian A
    double sig; // Gaussian sigma
    double r; //noise scaling
} Opt;


void gaussian(Opt *o, double *b);
void save_output(Opt *o, double *b);
void noise(Opt *o, double *b);
void cg(Opt *o, double *x, double C, double *b); //solves Ax = B
double dotprod(int N, const double *r1, const double *r2);
void copy_vector(int N, double *r1, double *r2) ;
void linear_comb(int N, double *r1, double a, double *r2, double b, double *r3);

int main(int argc, char **argv) {
    int c;
    Opt *o = malloc(sizeof(struct _Opt_));

    //specify default values
    o->n = 20;
    o->N = 1000;
    o->dx = 0.1;
    o->dt = 0.001;
    o->alpha = 1;
    o->A = 2;
    o->sig = 1;
    o->r = 0.001;

    //get command line values
    while ((c = getopt (argc, argv, "n:d:N:t:a:A:s:r:")) != -1)
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
            default: abort();
        }

    int n = o->n;

    //make a domain (b)
    double *dom = malloc(sizeof(double) * o->n * o->n);
    //make a solution domain
    double *x = malloc(sizeof(double) * o->n * o->n);
    //make a rhs vector
    double *b = malloc(sizeof(double) * o->n * o->n);
    for (int i = 0; i < o->n*o->n; ++i) {
        dom[i] = 0;
        b[i] = 0;
    }
    //initialize gaussian
    gaussian(o, dom);

    //add noise
    noise(o, dom);

    //caculate C
    double C = o->alpha * o->dt / (2.0 * o->dx * o->dx);

    for (int k = 0; k < o->N; ++k) {
        //calculate B*dom
        for (int i = 1; i < n-1; ++i) {
            for (int j = 1; j < n-1; ++j) {
                int idx = i*n + j;
                b[idx] = (1.0 - 4.0 * C) * dom[idx]
                       + C * dom[idx+n]
                       + C * dom[idx-n]
                       + C * dom[idx+1]
                       + C * dom[idx-1];
            }
        }

        cg(o, x, C, b);
        //copy_vector(n*n, x, b);
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
    int n = o->n; int N = n*n;
    double *r = malloc(sizeof(double)*N);
    double *p = malloc(sizeof(double)*N);
    double *Ap = malloc(sizeof(double)*N);
    double alpha;
    double beta;
    for (int i = 0; i < N; ++i) {
        x[i] = 0;
    }
    copy_vector(N, r, b);
    copy_vector(N, p, b);
    double rsold, rsnew;
    rsold = dotprod(N, r, r);

    for (int k = 0; k < N; ++k) { //main iterative loop of cg
        for (int i = 0; i < n; ++i) { //iterating over domain
            for (int j = 0; j < n; ++j) {
                int idx = i*n + j;
                Ap[idx] = (1.0 + 4*C) * p[idx]
                        - ( i-1 == -1 ? 0 : C * p[idx - n])
                        - ( i+1 ==  n ? 0 : C * p[idx + n])
                        - ( j-1 == -1 ? 0 : C * p[idx - 1])
                        - ( j+1 ==  n ? 0 : C * p[idx + 1]);
            }
        }
        double pAp = dotprod(N, p, Ap);
        alpha = rsold / pAp;

        linear_comb(N, x, 1.0, x, alpha, p);
        linear_comb(N, r, 1.0, r, -alpha, Ap);

        rsnew = dotprod(N, r, r);
        if (sqrt(rsnew) < 1e-10) return;

        beta = rsnew / rsold;
        linear_comb(N, p, 1.0, r, beta, p);
        rsold = rsnew;
    }
    free(r); free(p); free(Ap);
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


