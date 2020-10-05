void trifac(int nrows, int **ifp);
double conv_cg(int nc, double p[]);
double conv_cg_deriv(int nc, double dres[], double cg[]);
double dot_product(double a[], double b[]);

double *vector_double(int length);
int *vector_int(int length);
double **matrix_double(int rows, int columns);
int **matrix_int(int rows, int columns);
double ***matrix_3_double(int n_1, int n_2, int n_3);
void deallocate_vector(void *p_x);
void deallocate_matrix(void **p_x, int rows);
void deallocate_matrix_3(void ***p_x, int n_1, int n_2);

void areanorm_cg(double t[], double f[], int ndir, int nfac, int **ifp, 
              double at[], double af[], double a_ell, double b_ell, double c_ell);
double bright_all(double a[]);
void dbright_all(double a[], double dyda[]);
double bright_cg(double ee[], double ee0[], double t, double cg[]);
double bright_cg_deriv(double ee[], double ee0[], double t, double cg[], double dyda[]);
double phasec_cg(double alpha);
void matrix_cg(double t, double tmat[][4]);
void blmatrix_cg();

void frprmn(double p[], int n, int *iter, double *fret, double (*func)(), void (*dfunc)(), int itmax,
            double conw, int ndata, int verbose);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)());
double dbrent(double ax, double bx, double cx, double (*f)(), double (*df)(), double tol, double *xmin);
void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []), void (*dfunc)(double [], double []));
double f1dim(double x);
double df1dim(double x);

