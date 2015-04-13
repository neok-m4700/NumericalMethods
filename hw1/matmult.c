#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

double **matrix_new(int n){
   int i;
   double **m = malloc(sizeof(double*)*n);
   m[0] = malloc(sizeof(double)*n*n);
   for(i = 0; i < n; ++i){
      m[i] = &m[0][i*n];
   }
   return m;
}

void matrix_free(double **m){
   free(m[0]);
   free(m);
}

double **generate_matrix(double max, int n){
   int i,j;
   double *a = malloc(sizeof(double)*n);
   double **m = matrix_new(n);

   for(i = 0; i < n; ++i) 
      a[i] = max/n*i;
   for(i = 0; i < n; ++i){
      for(j = 0; j < n; ++j) 
         m[i][j] = a[i]*a[j];
   }
   free(a);
   return m;
}

double **matmult(double **m1, double **m2, int n){
   int i,j,k;
   double **m = matrix_new(n);
   for(i = 0; i < n; ++i){
      for(j = 0; j < n; ++j){
         m[i][j] = 0;
         for(k = 0; k < n; ++k){
            m[i][j] += m1[i][k]*m2[k][j];
         }
      }
   }
}

float time_it(int n){
   clock_t time;
   double **m1 = generate_matrix(1,n);
   double **m2 = generate_matrix(2,n);
   time = clock();
   double **m3 = matmult(m1,m2, n);
   time = clock() - time;
   matrix_free(m1);
   matrix_free(m2);
   return (float)time/CLOCKS_PER_SEC;
}

int main(int argc, char **argv){
   int i;
   float runt;

   if (argc > 1){ //if n is provided, just run for that n
      i = atoi(argv[1]);
      runt = time_it(i);
      printf("%4.4e\n", runt);
   } else { //run for a range of n's
      for(i = 100; i <= 1000; i += 100){
         runt = time_it(i);
         printf("%d  %4.4e\n", i, runt);
      }
   }

   return 0;
} 


