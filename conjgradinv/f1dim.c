/* Copyright - Press, Teukolsky, Vetterling, and Flannery - Numerical Recipes in C, CUP 1992 */

#define NRANSI
#include "declarations.h"

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);

double f1dim(double x)
{
	int j;
	double f,*xt;

	xt=vector_double(ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=nrfunc(xt);
	deallocate_vector((void *) xt);
	return f;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software .2#5n2. */
