void blmatrix_direct(double bet, double lam);
double phasec_direct(double alpha, double p[]);
void matrix_direct(double omg, double t, double tmat[][4]);
double bright_direct(double ee[], double ee0[], double t, double omg, double par[], double cl, double cls);

double *vector_double(int length);
int *vector_int(int length);
double **matrix_double(int rows, int columns);
int **matrix_int(int rows, int columns);
double ***matrix_3_double(int n_1, int n_2, int n_3);
void deallocate_vector(void *p_x);
void deallocate_matrix(void **p_x, int rows);
void deallocate_matrix_3(void ***p_x, int n_1, int n_2);

double dot_product(double a[], double b[]);
double cross_product(double a[], double b[], double c[]);
