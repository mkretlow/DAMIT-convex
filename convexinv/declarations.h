
void trifac(int nrows, int **ifp);
void areanorm(double t[], double f[], int ndir, int nfac, int **ifp, 
              double at[], double af[]);
void sphfunc(int ndir, double at[], double af[]);
void ellfit(double r[], double a, double b, double c, 
            int ndir, int ncoef, double at[], double af[]);
void lubksb(double **a, int n, int indx[], double b[]);
void ludcmp(double **a, int n, int indx[], double d[]);
void mrqmin(double **x1, double **x2, double x3[], double y[], 
            double sig[], double a[], int ia[], int ma, 
	    double **covar, double **alpha, double (*funcs)());
double mrqcof(double **x1, double **x2, double x3[], double y[], 
              double sig[], double a[], int ia[], int ma, 
	      double **alpha, double beta[], double (*funcs)());
void curv(double cg[]);
void blmatrix(double bet, double lam);
double conv(int nc, double dres[], int ma);
void gauss(double **aa, int n, double b[]);
void covsrt(double **covar, int ma, int ia[], int mfit);
void phasec(double dcdp[], double alpha, double p[]);
void matrix(double omg, double t, double tmat[][4], double dtm[][4][4]);
double bright(double ee[], double ee0[], double t, double cg[], 
            double dyda[], int ncoef);

double *vector_double(int length);
int *vector_int(int length);
double **matrix_double(int rows, int columns);
int **matrix_int(int rows, int columns);
double ***matrix_3_double(int n_1, int n_2, int n_3);
void deallocate_vector(void *p_x);
void deallocate_matrix(void **p_x, int rows);
void deallocate_matrix_3(void ***p_x, int n_1, int n_2);

double dot_product(double a[], double b[]);
