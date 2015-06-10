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

enum { MG, CG, MG_CG, SOR } methods;

void gaussian(Opt *o, double **b);
void save_output(Opt *o, double **b);
void print_output(Opt *o, double **b);
void noise(Opt *o, double **b);

int main(int argc, char **argv) {
    int c;
    Opt *o = malloc(sizeof(struct _Opt_));
    int ncycles; //for mutligrid
    int method; //which method of the 4

    //specify default values
    o->n = 33;
    o->N = 100;
    o->dx = 1.0;
    o->dt = 7.0;
    o->alpha = 0.1;
    o->A = 100;
    o->sig = 5;
    o->r = 2;
    o->source = 0.001;
    ncycles = 2;
    method = 0;

    //get command line values
    while ((c = getopt (argc, argv, "n:d:N:t:a:A:s:r:S:c:m:")) != -1)
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
            case 'm': method = atoi(optarg); break;
            default: abort();
        }

    int n = o->n;

    double **dom = dmatrix(1, n, 1, n);
    double **aux = dmatrix(1,n,1,n);
    double **x = dmatrix(1,n,1,n);

    //initialize gaussian
    gaussian(o, dom);

    //add noise
    noise(o, dom);
    //print_output(o, dom);

    //caculate C
    C = o->alpha * o->dt / (2.0 * o->dx * o->dx);

    switch(method){
        case(MG): printf("Running MG\n"); break;
        case(CG): printf("Running CG\n"); break;
        case(MG_CG): printf("Running MG+CG\n"); break;
        case(SOR): printf("Running SOR\n"); break;
    }


    for (int k = 0; k < o->N; ++k) { //time stepping

        switch(method) //run solver
        {
            case (MG): 
                mglin(dom, n, ncycles); 
                break;
            case (CG):
                copy(x, dom, n); //guess goes in x
                cg(dom, x, n); //run cg
                break;
            case (MG_CG):
                copy(x, dom, n); //guess goes in x
                mglin(x, n, 1); //run one cycle of mg
                cg(dom, x, n); //run cg
                break;
            case(SOR):
                sor(dom, n, 0.9);
        }

        //add source
        for (int i = 2; i < n; ++i)
            for (int j = 2; j < n; ++j)
                dom[i][j] += o->source;

    }

    //save the domain
    save_output(o, dom);

    //free stuff
    free_dmatrix(dom, 1, n, 1, n);
    free_dmatrix(aux, 1, n, 1, n);
    free_dmatrix(x, 1, n, 1, n);
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
        b[i][1] = b[i][n] = b[1][i] = b[n][i] = 0;
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

void print_output(Opt *o, double **b) {
   int n = o->n;
    FILE *f = fopen("out.data", "w");
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            printf("%0.5lf ", b[i][j]);
        }
        printf("\n");
    }
    fclose(f);
}

