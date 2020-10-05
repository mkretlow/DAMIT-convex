/* Find the curv. fn. Laplace series for given ellipsoid */

/*
Copyright (C) 2006  Mikko Kaasalainen, Josef Durech

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <math.h>
#include "globals.h"
#include "declarations.h"

void ellfit(double cg[], double a, double b, double c,  
            int ndir, int ncoef, double at[], double af[])
{
   int i, m, l, n, j, k;
   int *indx;
   
   double sum, st;
   double *fitvec, *d, *er,
          **fmat, **fitmat;

   indx = vector_int(ncoef);

   fitvec = vector_double(ncoef);
   er = vector_double(ndir);
   d = vector_double(1);

   fmat = matrix_double(ndir, ncoef);
   fitmat = matrix_double(ncoef, ncoef);

   /* Compute the LOGcurv.func. of the ellipsoid */
   for (i = 1; i <= ndir; i++)
   {
      st = sin(at[i]);
      sum = pow(a * st * cos(af[i]), 2) + pow(b * st * sin(af[i]), 2) + 
            pow(c * cos(at[i]), 2);
      er[i] = 2 * (log(a * b * c) - log(sum));
   }
   /* Compute the sph. harm. values at each direction and
      construct the matrix fmat from them */
   for (i = 1; i <= ndir; i++)
   {
      n = 0;
      for (m = 0; m <= Mmax; m++)
         for (l = m; l <= Lmax; l++)
	 {
            n++;
            if (m != 0) 
	    {
               fmat[i][n] = Pleg[i][l][m] * cos(m * af[i]);
               n++;
               fmat[i][n] = Pleg[i][l][m] * sin(m * af[i]);
            }
            else
               fmat[i][n] = Pleg[i][l][m];
         }
   }

   /* Fit the coefficients r from fmat[ndir,ncoef]*r[ncoef]=er[ndir] */
   for (i = 1; i <= ncoef; i++)
   {
      for (j = 1; j <= ncoef; j++)
      {
         fitmat[i][j]=0;

         for (k = 1; k <= ndir; k++)
            fitmat[i][j] = fitmat[i][j] + fmat[k][i] * fmat[k][j];
	    
      }
      fitvec[i]=0;
      
      for (j = 1; j <= ndir; j++)
         fitvec[i] = fitvec[i] + fmat[j][i] * er[j];
   }


   ludcmp(fitmat,ncoef,indx,d);
   lubksb(fitmat,ncoef,indx,fitvec);
   
   for (i = 1; i <= ncoef; i++)
      cg[i] = fitvec[i];

   deallocate_matrix((void *) fitmat, ncoef);
   deallocate_matrix((void *) fmat, ndir);
   deallocate_vector((void *) fitvec);
   deallocate_vector((void *) d);
   deallocate_vector((void *) indx);
   deallocate_vector((void *) er);
         
}

