/* This program take the input lightcurves, scans over the
   given period range starting from 6 different poles and finds the best 
   period+pole+shape+scattering 
   solution. All but the period is forgotten. The period and the rms residual 
   of the fit given to the output.
   
   syntax:
   cat input_lcs | period_scan [-v] input_parameters output_periods
*/

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
   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <memory.h>

#include "declarations.h"
#include "constants.h"

/* global parameters */
int Lmax, Mmax, Niter, Lastcall, Ncoef, Numfac, Lcurves, Nphpar,
    Lpoints[MAX_LC+1], Inrel[MAX_LC+1], Deallocate;
    
double Ochisq, Chisq, Alamda, Alamda_incr, Alamda_start, Phi_0, Scale,
       Area[MAX_N_FAC+1], Darea[MAX_N_FAC+1], Sclnw[MAX_LC+1],
       Yout[MAX_N_OBS+1],
       Fc[MAX_N_FAC+1][MAX_LM+1], Fs[MAX_N_FAC+1][MAX_LM+1], 
       Tc[MAX_N_FAC+1][MAX_LM+1], Ts[MAX_N_FAC+1][MAX_LM+1], 
       Dsph[MAX_N_FAC+1][MAX_N_PAR+1], Dg[MAX_N_FAC+1][MAX_N_PAR+1],    
       Nor[MAX_N_FAC+1][4], Blmat[4][4],
       Pleg[MAX_N_FAC+1][MAX_LM+1][MAX_LM+1],
       Dblm[3][4][4],
       Weight[MAX_N_OBS+1];
    
/*--------------------------------------------------------------*/

FILE *f_1, *f_per;
  
int main(int argc, char *argv[])
{
   int i, j, l, m, k, n, nrows, ndata, k2, ndir, i_temp, niter_best,
       n_iter_max, n_iter_min, ind_par_file, ind_out_file, verbose, arg_shift, 
       *ia,  ial0, ia_par[4], ia_cl,
       **ifp; 
  
   double per_start, per_step_coef, per_end,
          freq, freq_start, freq_step, freq_end, jd_min, jd_max,
          dev_old, dev_new, iter_diff, iter_diff_max, stop_condition,
          jd_0, fi_0=0, conw, a=1.05, b=1, c=0.95, prd, cl, al0, ave, e0len, elen, cos_alpha,
	  totarea, sum, dark, dth, dph, rfit, escl, dev_best, per_best, chisq_best, dark_best,
          *t, *f, *at, *af,   
          *brightness, e[4], e0[4], **ee,
          **ee0, *cg, *cg_first, **covar,
          **aalpha, *sig, chck[4], lambda_start[7], beta_start[7],
          *tim, *al, par[4], rchisq;

   char *str_temp;

   str_temp = (char *) malloc (MAX_LINE_LENGTH);

   ee = matrix_double(MAX_N_OBS,3);
   ee0 = matrix_double(MAX_N_OBS,3);
   covar = matrix_double(MAX_N_PAR,MAX_N_PAR);
   aalpha = matrix_double(MAX_N_PAR,MAX_N_PAR);
   ifp = matrix_int(MAX_N_FAC,4);
   
   tim = vector_double(MAX_N_OBS);
   brightness = vector_double(MAX_N_OBS);
   sig = vector_double(MAX_N_OBS);
   cg = vector_double(MAX_N_PAR);
   cg_first = vector_double(MAX_N_PAR);
   t = vector_double(MAX_N_FAC);
   f = vector_double(MAX_N_FAC);
   at = vector_double(MAX_N_FAC);
   af = vector_double(MAX_N_FAC);
      
   ia = vector_int(MAX_N_PAR);

   /* recognizes arguments */
   ind_par_file = 1;
   ind_out_file = 2;
   arg_shift = 0;
   verbose = 0;
   for (i = 1; i <= argc-1; i++)
   {
      if (strcmp(argv[i], "-v") == 0)
      {
         verbose = 1;
	 arg_shift++;
      }      
   }      
   ind_par_file += arg_shift;
   ind_out_file += arg_shift;
   
   /* parameters file */
   if ((f_1 = fopen(argv[ind_par_file], "r")) == NULL)
      {printf("cannot open 'input parameters' file \n"); fflush(stdout); exit(2);}

   /* period interval (hrs) */
   fscanf(f_1, "%lf %lf %lf", &per_start, &per_end, &per_step_coef); fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* the weight factor for conv. reg. */
   fscanf(f_1, "%lf", &conw);                                 fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* degree and order of the Laplace series */
   fscanf(f_1, "%d %d", &Lmax, &Mmax);                        fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* nr. of triangulation rows per octant */   
   fscanf(f_1, "%d", &nrows);                                 fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* Initial guesses for phase funct. params. */   
   fscanf(f_1, "%lf %d", &par[1], &ia_par[1]);                fgets(str_temp, MAX_LINE_LENGTH, f_1);
   fscanf(f_1, "%lf %d", &par[2], &ia_par[2]);                fgets(str_temp, MAX_LINE_LENGTH, f_1);
   fscanf(f_1, "%lf %d", &par[3], &ia_par[3]);                fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* Initial Lambert coeff. (L-S=1) */
   fscanf(f_1, "%lf %d", &cl, &ia_cl);                        fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* maximum number of iterations (when > 1) or
      minimum difference in dev to stop (when < 1) */
   fscanf(f_1, "%lf", &stop_condition);                       fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* minimum number of iterations when stop_condition < 1 */
   fscanf(f_1, "%d", &n_iter_min);                            fgets(str_temp, MAX_LINE_LENGTH, f_1);
   fclose(f_1);

   /*
   For low Lmax, Mmax the results strongly depends on
   the initial epocj jd_0, and rotation angle fi_0 !!!
   */
   if (verbose == 1)
   {
      printf("%g %g   period interval (hrs) \n", per_start, per_end);   
      printf("%g  the weight factor for conv. reg.\n", conw);  
      printf("%d %d  degree and order of the Laplace series\n", Lmax, Mmax);  
      printf("%d  nr. of triangulation rows per octant\n", nrows);    
      printf("%g %g %g  initial guesses for phase funct. params. (%d,%d,%d)\n", par[1], par[2], par[3], ia_par[1], ia_par[2], ia_par[3]);  
      printf("%g  initial Lambert coeff. (L-S=1) (%d)\n", cl, ia_cl);   
      printf("%g  stop condition\n", stop_condition);   
      printf("%d  minimum number of iterations\n", n_iter_min);   
   }

   /* lightcurves + geometry file */   
   fscanf(stdin, "%d", &Lcurves);
      
   if (Lcurves > MAX_LC)
   {
      fprintf(stderr, "\nError: Number of lcs  is greater than MAX_LC = %d\n", MAX_LC);
      fflush(stderr); exit(1);
   }

   al = vector_double(Lcurves);   
   
   ndata = 0; /* total number of data */
   k2 = 0;   /* index */
   al0 = PI; /* the smallest solar phase angle */
   ial0 = -1; /* initialization, index of al0 */
   jd_min = 1e20; /* initial minimum and minimum JD */
   jd_max = -1e20; 
  
   /* loop over lightcurves */   
   for (i = 1; i <= Lcurves; i++)
   {
      ave = 0; /* average */
      fscanf(stdin, "%d %d", &Lpoints[i], &i_temp); /* points in this lightcurve */
      Inrel[i] = 1 - i_temp;
     
      if (Lpoints[i] > POINTS_MAX)
      {
         fprintf(stderr, "\nError: Number of lc points is greater than POINTS_MAX = %d\n", POINTS_MAX);
         fflush(stderr); exit(1);
      }

      /* loop over one lightcurve */
      for (j = 1; j <= Lpoints[i]; j++)
      {
         ndata++;
	 
         if (ndata > MAX_N_OBS)
         {
            fprintf(stderr, "\nError: Number of data is greater than MAX_N_OBS = %d\n", MAX_N_OBS);
            fflush(stderr); exit(1);
         }

	 fscanf(stdin, "%lf %lf", &tim[ndata], &brightness[ndata]); /* JD, brightness */	 
	 fscanf(stdin, "%lf %lf %lf", &e0[1], &e0[2], &e0[3]); /* ecliptic astr_tempocentric coord. of the Sun in AU */
	 fscanf(stdin, "%lf %lf %lf", &e[1], &e[2], &e[3]); /* ecliptic astrocentric coord. of the Earth in AU */	 
	 /* selects the minimum and maximum JD */
	 if (tim[ndata] < jd_min) jd_min = tim[ndata];
	 if (tim[ndata] > jd_max) jd_max = tim[ndata];
	 
	 /* normals of distance vectors */
         e0len=sqrt(e0[1]*e0[1] + e0[2]*e0[2] + e0[3]*e0[3]);
         elen=sqrt(e[1]*e[1] + e[2]*e[2] + e[3]*e[3]);

         ave = ave + brightness[ndata];

	 /* normalization of distance vectors */
         for (k = 1; k <= 3; k++)
	 {
            ee[ndata][k] = e[k] / elen;
            ee0[ndata][k] = e0[k] / e0len;
         }

         cos_alpha = dot_product(e, e0) / (elen * e0len); 
         al[i] = acos(cos_alpha); /* solar phase angle */
         /* Find the smallest solar phase al0 (not important, just for info) */
         if (al[i] < al0) 
         {
            al0 = al[i];
            ial0 = ndata;
         }
      } /* j, one lightcurve */           
   
      ave = ave / Lpoints[i];

      /* Mean brightness of lcurve
         Use the mean brightness as 'sigma' to renormalize the
         mean of each lightcurve to unity */

      for (j = 1; j <= Lpoints[i]; j++)
      {
         k2++;
         sig[k2] = ave;
      }
   } /* i, all lightcurves */        

   /* jd_0 is set to the day before the lowest JD in the data */
   jd_0 = (int) jd_min;
   if (verbose == 1)
      printf("\nThe used epoch of zero time  %f\n\n", jd_0);
       
   /* loop over data - subtraction of jd_0 */   
   for (i = 1; i <= ndata; i++)
      tim[i] = tim[i] - jd_0;         
   
   fi_0 = fi_0 * DEG2RAD;

   /* Initial shape guess:
      N.B. the initial guess should not be exactly a sphere (for which
      l=1-coeffs. are redundant, resulting in a bad first step)
      This wild rfit guess is for L-S at opposition (usually irrelevant) */
   rfit = sqrt(2 * sig[ial0] / (0.5 * PI * (1+cos(al0))));
   escl = rfit / sqrt((a * b + b * c + a * c) / 3);
   a = a * escl;
   b = b * escl;
   c = c * escl;
 
   /* Convexity regularization: make one last 'lightcurve' that
      consists of the three comps. of the residual nonconv. vect.
      that should all be zero */
   Lcurves = Lcurves + 1;
   Lpoints[Lcurves] = 3;
   Inrel[Lcurves] = 0;
   conw /= escl * escl;
   for (j = 1; j <= 3; j++)
   {
      ndata++;
      brightness[ndata] = 0;
      sig[ndata] = 1 / conw;
   }

   /* the ordering of the coeffs. of the Laplace series */
   Ncoef = 0; /* number of coeffs. */
   for (m = 0; m <= Mmax; m++)
      for (l = m; l <= Lmax; l++)
      {
         Ncoef++;
         if (m != 0) Ncoef++;
      }

   /*  Fix the directions of the triangle vertices of the Gaussian image sphere 
       t = theta angle, f = phi angle */
   dth = PI / (2 * nrows); /* step in theta */
   k = 1;
   t[1] = 0;
   f[1] = 0;
   for (i = 1; i <= nrows; i++)
   {
      dph = PI / (2 * i); /* step in phi */
      for (j = 0; j <= 4 * i - 1; j++)
      {
         k++;
         t[k] = i * dth;
         f[k] = j * dph;
      }
   }

   /* go to south pole (in the same rot. order, not a mirror image) */
   for (i = nrows - 1; i >= 1; i--)
   {
      dph = PI / (2 * i);
      for (j = 0; j <= 4 * i - 1; j++)
      {
         k++;
         t[k] = PI - i * dth;
         f[k] = j * dph;
      }
   }      

   ndir = k + 1; /* number of vertices */

   t[ndir] = PI;
   f[ndir] = 0;
   Numfac = 8 * nrows * nrows;
   
   if (Numfac > MAX_N_FAC)
   {   
      fprintf(stderr, "\nError: Number of facets is greater than MAX_N_FAC!\n"); 
      fflush(stdout); exit(1);
   }   

   /* makes indices to triangle vertices */
   trifac(nrows, ifp);
   /* areas and normals of the triangulated Gaussian image sphere */
   areanorm(t,f,ndir,Numfac, ifp, at, af);
   /* Precompute some function values at each normal direction*/
   sphfunc(Numfac, at, af);
   
   ellfit(cg_first, a, b, c, Numfac, Ncoef, at, af);

   freq_start = 1 / per_start;
   freq_end = 1 / per_end;
   freq_step = 0.5 / (jd_max - jd_min) / 24 * per_step_coef;

   /* pole and period are free */
   ia[Ncoef+1] = 1;
   ia[Ncoef+2] = 1;
   ia[Ncoef+3] = 1;
   /* phase function parameters */
   Nphpar = 3;
   /* shape is free to be optimized */
   for (i = 2; i <= Ncoef; i++)
      ia[i] = 1;
   /* The first shape param. fixed for relative br. fit */
   ia[1] = 0;
   /* is there any absolute lc.? */
   for (i = 1; i <= Lcurves - 1; i++)
      if (Inrel[i] == 0)
         ia[1] = 1; 
	 
   ia[Ncoef+3+Nphpar+1] = ia_cl;
   /* Lommel-Seeliger part is fixed */
   ia[Ncoef+3+Nphpar+2] = 0;
 
   if ((Ncoef+3+Nphpar+1) > MAX_N_PAR)
   {
      fprintf(stderr, "\nError: Number of parameters is greater than MAX_N_PAR = %d\n", MAX_N_PAR);
      fflush(stderr); exit(1);
   }

   if ((f_per = fopen(argv[ind_out_file], "w")) == NULL)
   {
      fprintf(stderr, "\nError: cannot open 'output periods' file \n"); 
      fflush(stderr); exit(2);
   }
 
   /* six initial poles */
   lambda_start[1] = 0;   beta_start[1] = 0;
   lambda_start[2] = 180; beta_start[2] = 0;   
   lambda_start[3] = 90;  beta_start[3] = 60;
   lambda_start[4] = 240; beta_start[4] = 60;
   lambda_start[5] = 90;  beta_start[5] = -60;
   lambda_start[6] = 240; beta_start[6] = -60;

   printf("period   rms      chi2      iter. dark area %%\n");
   
   for (n = 1; n <= (int) ((freq_start - freq_end) / freq_step) + 1; n++)
   {
       dev_best = 1e10;
       per_best = 0;
       chisq_best = 0;
       niter_best = 0;
       dark_best = 0;
      /* loop over poles */
      for (k = 1; k <= 6; k++)
      {
         /* starts from the initial ellipsoid */
         for (i = 1; i <= Ncoef; i++)
            cg[i] = cg_first[i];
   
         freq = freq_start - (n - 1) * freq_step;
         prd = 1 / freq;

         cg[Ncoef+1] = beta_start[k];
         cg[Ncoef+2] = lambda_start[k];
   
         /* The formulas use beta measured from the pole */
         cg[Ncoef+1] = 90 - cg[Ncoef+1];
         /* conversion of lambda, beta to radians */
         cg[Ncoef+1] = DEG2RAD * cg[Ncoef+1];
         cg[Ncoef+2] = DEG2RAD * cg[Ncoef+2];
      
         /* Use omega instead of period */
         cg[Ncoef+3] = 24 * 2 * PI / prd;

         for (i = 1; i <= Nphpar; i++)
         {
            cg[Ncoef+3+i] = par[i];
            ia[Ncoef+3+i] = ia_par[i];
         }
         /* Lommel-Seeliger part */
         cg[Ncoef+3+Nphpar+2] = 1;
         /* Use logarithmic formulation for Lambert to keep it positive */
         cg[Ncoef+3+Nphpar+1] = log(cl);

         /* Levenberg-Marquardt loop */
         n_iter_max = 0; 
         iter_diff_max = -1;
         rchisq = -1;
         if (stop_condition > 1)
         { 
            n_iter_max = (int) stop_condition;      
            iter_diff_max = 0;
            n_iter_min = 0; /* to not overwrite the n_iter_max value */
         }      
         if (stop_condition < 1)
         { 
            n_iter_max = MAX_N_ITER; /* to avoid neverending loop */     
            iter_diff_max = stop_condition;
         }      
         Alamda = -1;
         Niter = 0;
         iter_diff = 1e40;
         dev_old = 1e30;
         dev_new = 0;
         Lastcall = 0; 
      
         while (((Niter < n_iter_max) && (iter_diff > iter_diff_max)) || (Niter < n_iter_min))
         {
            mrqmin(ee,ee0,tim,brightness,sig,cg,ia,Ncoef+5+Nphpar,covar,aalpha,bright);
            Niter++;
    
            if ((Niter == 1) || (Chisq < Ochisq))
            {
               Ochisq = Chisq;
               curv(cg);
               for (i = 1; i <= 3; i++)
               {
                  chck[i] = 0;
		  totarea = 0;
                  for (j = 1; j <= Numfac; j++)
		  {
                     chck[i] += Area[j] * Nor[j][i];
		     totarea += Area[j];
		  }
               }
               rchisq = Chisq - (pow(chck[1] / totarea, 2) + pow(chck[2] / totarea, 2) + pow(chck[3] / totarea, 2)) * pow(conw,2);
            }
            dev_new = sqrt(rchisq / (ndata - 3));
            /* only if this step is better than the previous,
               1e-10 is for numeric errors */
            if (dev_old - dev_new > 1e-10)
            {
               iter_diff = dev_old - dev_new;
               dev_old = dev_new;
            }
/*	    printf("%d %f %f %f \n", Niter, rchisq, dev_new, Alamda);*/

         }

         totarea = 0;
         for (i = 1; i <= Numfac; i++)
            totarea = totarea + Area[i];
         sum = pow(chck[1],2) + pow(chck[2],2) + pow(chck[3],2);
         dark = sqrt(sum);
       
         /* period solution */
         prd = 2 * PI / cg[Ncoef+3];

         if (dev_new < dev_best)
	 {
	    dev_best = dev_new;
	    per_best = 24 * prd;
	    chisq_best = dev_new * dev_new * (ndata - 3);
	    niter_best = Niter;
	    dark_best = dark / totarea * 100;
         }	    
	 
         /* deallocates variables used in mrqmin */
         Deallocate = 1;
         mrqmin(ee,ee0,tim,brightness,sig,cg,ia,Ncoef+5+Nphpar,covar,aalpha,bright);
         Deallocate = 0;
      } /* pole loop */

      
      if (verbose == 1)
         printf("%f %f %f  %d   %4.1f\n", per_best, dev_best, chisq_best, niter_best, dark_best);

      /* output file */
      fprintf(f_per, "%f %f %f %d %4.1f\n", per_best, dev_best, chisq_best, niter_best, dark_best);
      
   } /* period loop */      

   fclose(f_per);
   
   return(0);
}


      
   

