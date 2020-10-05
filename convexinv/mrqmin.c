/* this is a modified version of the routine 
   from Press, Teukolsky, Vetterling, and Flannery - Numerical Recipes in C, CUP 1992 */
   
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
   

#include "globals.h"
#include "declarations.h"
#include "constants.h"

void mrqmin(double **x1, double **x2, double x3[], double y[], 
            double sig[], double a[], int ia[], int ma, 
            double **covar, double **alpha, double (*funcs)())
{     

   int j, k, l;
   static int mfit; /* it is set in the firs call */

   static double *atry, *beta, *da;

   double temp;

   /* dealocates memory when usen in period_search */
   if (Deallocate == 1)
   {
      deallocate_vector((void *) atry);
      deallocate_vector((void *) beta);
      deallocate_vector((void *) da); 
      return;
   }
         
   if (Lastcall != 1)
   {
      if(Alamda < 0)
      {
         atry = vector_double(ma);
         beta = vector_double(ma);
         da = vector_double(ma);
   
         /* number of fitted parameters */
         mfit = 0;
         for (j = 1; j <= ma; j++)  
            if (ia[j]) mfit++;

         Alamda = ALAMDA_START; /* initial alambda */

         temp = mrqcof(x1,x2,x3,y,sig,a,ia,ma,alpha,beta,funcs);

         Ochisq = temp;
         for (j = 1; j <= ma; j++)
            atry[j] = a[j];
      }
      for (j = 1; j <= mfit; j++)
      {
         for (k = 1; k <= mfit; k++)
            covar[j][k] = alpha[j][k];
         covar[j][j] = alpha[j][j] * (1 + Alamda);
         da[j] = beta[j];
      }

      gauss(covar,mfit,da);

      if (Alamda == 0)
      {
         covsrt(covar,ma,ia,mfit);
         return;
      }
      
      j = 0;
      for (l = 1; l <= ma; l++)
        if(ia[l]) 
	{
           j++;
           atry[l] = a[l] + da[j];
        }
   } /* Lastcall != 1 */
   
   if (Lastcall == 1)
      for (l = 1; l <= ma; l++)
         atry[l] = a[l];
  
   temp = mrqcof(x1,x2,x3,y,sig,atry,ia,ma,covar,da,funcs);
   Chisq = temp;


   if (Lastcall == 1) 
   {
      deallocate_vector((void *) atry);
      deallocate_vector((void *) beta);
      deallocate_vector((void *) da); 
      return;
   }
   
   if (temp < Ochisq)
   {
      Alamda = Alamda / ALAMDA_COEFF;
      for (j = 1; j <= mfit; j++)
      {
         for (k = 1; k <= mfit; k++)
            alpha[j][k] = covar[j][k];
         beta[j] = da[j];
      }
      for (l = 1; l <= ma; l++)
         a[l] = atry[l];
   }
   else
   {
      Alamda = ALAMDA_COEFF * Alamda;
      Chisq = Ochisq;
   }

}

