/* Copyright - Press, Teukolsky, Vetterling, and Flannery - Numerical Recipes in C, CUP 1992 */
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "declarations.h"
#include "globals.h"

void frprmn(double p[], int n, int *iter, double *fret, double (*func)(), void (*dfunc)(), int itmax, double conw, int ndata, int verbose)
{
	int j, its, i;
	double gg, gam, fp, dgg, totarea, dev, expp, dark;
	double *g, *h, *xi, res[4];
	
	g = vector_double(n);
	h = vector_double(n);
	xi = vector_double(n);

	fp = func(p);
	dfunc(p, xi);
	
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}

	/* iteration loop */
	for (its=1;its<=itmax;its++) {
		*iter=its;
		dlinmin(p,xi,n,fret,func,dfunc);
		
		fp=func(p);
		dfunc(p,xi);

                if (verbose == 1)
                {
		   totarea = res[1] = res[2] = res[3] = 0;
		   for (i = 1; i <= Numfac; i++)
		   {
		      expp = exp(p[i]);
		      totarea += expp;
		      for (j = 1; j <= 3; j++)
			 res[j] += expp * Nor[i][j];
		   }	 
		   dev = sqrt(Rchisq / (ndata - 3));
		   dark = sqrt(pow(res[1],2) + pow(res[2],2) + pow(res[3],2)) / totarea * 100;
		   printf("%d   chi2 %f   dev %f   dark facet %4.2f%%\n", its, Rchisq, dev, dark);
                }

		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) 
			return;
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	return;
}

