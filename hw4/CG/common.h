double dotprod(int N, const double *r1, const double *r2);
void copy_vector(int N, double *r1, double *r2) ;
void linear_comb(int N, double *r1, double a, double *r2, double b, double *r3);
void bi_diagonal_prod(double C, int n, double *product, const double *multiplicant);
void solve(int n, double *x, double C, double *b); //solves Ax = B
