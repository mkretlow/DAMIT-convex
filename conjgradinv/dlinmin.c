/* Copyright - Press, Teukolsky, Vetterling, and Flannery - Numerical Recipes in C, CUP 1992 */

#define NRANSI
#define TOL 2.0e-4

#include "declarations.h"

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);
void (*nrdfun)(double [], double []);

void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []), void (*dfunc)(double [], double []))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector_double(n);
	xicom=vector_double(n);
	nrfunc=func;
	nrdfun=dfunc;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	deallocate_vector((void *) xicom);
	deallocate_vector((void *) pcom);
}
#undef TOL
/* (C) Copr. 1986-92 Numerical Recipes Software .2#5n2. */
