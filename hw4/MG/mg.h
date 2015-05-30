#define NPRE  1
#define NPOST 1
#define NGMAX 15

double C;
void bi_diagonal_prod(double C, int n, double **aux, double **dom);
void addint(double **uf, double **uc, double **res, int nf);
void copy(double **aout, double **ain, int n);
void fill0(double **u, int n);
void interp(double **uf, double **uc, int nf);
void relax(double **u, double **rhs, int n);
void resid(double **res, double **u, double **rhs, int n);
void rstrct(double **uc, double **uf, int nc);
void slvsml(double **u, double **rhs);
void mglin(double **u, int n, int ncycle);
void cg(double **rhs, double **guess, int n);
