/* Copyright - Press, Teukolsky, Vetterling, and Flannery - Numerical Recipes in C, CUP 1992 */

#include "declarations.h"

extern int ncom;	/* defined in DLINMIN */
extern double *pcom,*xicom,(*nrfunc)();
extern void (*nrdfun)();

double df1dim(double x)
{
	int j;
	double df1=0.0;
	double *xt,*df;

	xt=vector_double(ncom);
	df=vector_double(ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	nrdfun(xt,df);
	for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
	deallocate_vector((void *) df);
	deallocate_vector((void *) xt);
	return df1;
}
