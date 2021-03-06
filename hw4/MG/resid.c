#include "mg.h"

void resid(double **res, double **u, double **xold, int n)
/*
Returns minus the residual for the model problem. Input quantities are u[1..n][1..n] and
rhs[1..n][1..n], while res[1..n][1..n] is returned.
*/
{
  int i,j;
  /* Interior points.*/
  for (j=2;j<n;j++) 
    for (i=2;i<n;i++)
      res[i][j] = 
         (1-4*B)*xold[i][j] + B*(xold[i-1][j] + xold[i+1][j] + xold[i][j-1] + xold[i][j+1])
         + C*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])
         - (1+4*C)*u[i][j];
  /* Boundary points are copied.*/
  for (i=1;i<=n;i++) 
    res[i][1]=res[i][n]=res[1][i]=res[n][i]=0.0;
}
