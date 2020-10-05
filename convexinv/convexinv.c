/*     A procedure for shape+spin+scattering modelling of lightcurves
       (calibrated, uncalibrated, or sparse). Scattering law+shape representation 
       are simple and robustly converging.
       
       Original code written in Fortran by Mikko Kaasalainen, converted to C by Josef Durech.
       
       syntax:
       
       cat lcs_file | convexinv [-v] [-s] [-o output_file] [-p output_param_file] input_parameters output_lcs
      
       options:	 -v verbose
    		 -o write results to the output_file
	         -p write parameters to output file
                 -s write areas to standard output
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
int Lmax, Mmax, Niter, Lastcall, Ncoef, Numfac, Lcurves, Nphpar, Deallocate = 0,
    Lpoints[MAX_LC+1], Inrel[MAX_LC+1];
    
double Ochisq, Chisq, Alamda, Phi_0, Scale,
       Area[MAX_N_FAC+1], Darea[MAX_N_FAC+1], Sclnw[MAX_LC+1], Yout[MAX_N_OBS+1],
       Fc[MAX_N_FAC+1][MAX_LM+1], Fs[MAX_N_FAC+1][MAX_LM+1], 
       Tc[MAX_N_FAC+1][MAX_LM+1], Ts[MAX_N_FAC+1][MAX_LM+1], 
       Dsph[MAX_N_FAC+1][MAX_N_PAR+1], Dg[MAX_N_FAC+1][MAX_N_PAR+1],    
       Nor[MAX_N_FAC+1][4], Pleg[MAX_N_FAC+1][MAX_LM+1][MAX_LM+1],
       Blmat[4][4], Dblm[3][4][4];
    
FILE *f_1, *f_lc_out, *f_areas, *f_par;

/* -------------------------------------------------------------------------------*/
  
int main(int argc, char *argv[])
{
   int i, j, l, m, k, nrows, ndata, k2, ndir,  
       *ia,  ial0, ial0_abs, ia_beta_pole, ia_lambda_pole, ia_prd, ia_par[4], ia_cl,
       ind_par_file, ind_lc_out_file, ind_out_file, arg_shift, verbose, output, onlyrel,
       i_temp, ind_parout_file, paramout, std_out, n_iter_max, **ifp;
          
   double jd_0, conw, a, b, c, prd, cl, al0, al0_abs, ave, *e0len, *elen, cos_alpha,
          dth, dph, rfit, escl, totarea, sum, dark, beta_pole, lambda_pole, par[4], rchisq,
          dev_old, dev_new, iter_diff, iter_diff_max, stop_condition, jd_min, 
	  *brightness, e[4], e0[4], **ee,
          **ee0, *cg, **covar, *t, *f, *at, *af, **aalpha, *sig, chck[4],
          *tim, *al, *alpha;

   char *str_temp;

   str_temp = (char *) malloc (MAX_LINE_LENGTH);

   Lastcall = 0; 

   ee = matrix_double(MAX_N_OBS,3);
   ee0 = matrix_double(MAX_N_OBS,3);
   covar = matrix_double(MAX_N_PAR,MAX_N_PAR);
   aalpha = matrix_double(MAX_N_PAR,MAX_N_PAR);
   ifp = matrix_int(MAX_N_FAC,4);
   
   tim = vector_double(MAX_N_OBS);
   brightness = vector_double(MAX_N_OBS);
   alpha = vector_double(MAX_N_OBS);
   sig = vector_double(MAX_N_OBS);
   e0len = vector_double(MAX_N_OBS);
   elen = vector_double(MAX_N_OBS);
   cg = vector_double(MAX_N_PAR);
   t = vector_double(MAX_N_FAC);
   f = vector_double(MAX_N_FAC);
   at = vector_double(MAX_N_FAC);
   af = vector_double(MAX_N_FAC);
      
   ia = vector_int(MAX_N_PAR);

   /* command line arguments */
   ind_par_file = 1;
   ind_lc_out_file = 2;
   ind_out_file = 0;
   ind_parout_file = 0;
   arg_shift = 0;
   verbose = 0;
   output = 0;
   paramout = 0;
   std_out = 0;
   for (i = 1; i <= argc-1; i++)
   {
      if (strcmp(argv[i], "-v") == 0)
      {
         verbose = 1;
	 arg_shift++;
      }      
      if (strcmp(argv[i], "-s") == 0)
      {
         std_out = 1;
	 arg_shift++;
      }      
      if (strcmp(argv[i], "-o") == 0)
      {
         output = 1;
	 arg_shift += 2;
	 ind_out_file = i + 1;
      }      
      if (strcmp(argv[i], "-p") == 0)
      {
         paramout = 1;
	 arg_shift += 2;
	 ind_parout_file = i + 1;
      }      
   }      
   ind_par_file += arg_shift;
   ind_lc_out_file += arg_shift;

   /* input parameters file */
   if ((f_1 = fopen(argv[ind_par_file], "r")) == NULL)
      {fprintf(stderr, "\nError: cannot open 'input parameters' file \n"); fflush(stdout); exit(2);}
   /* initial lambda fixed or free  */
   fscanf(f_1, "%lf %d", &lambda_pole, &ia_lambda_pole);      fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* initial beta fixed or free  */
   fscanf(f_1, "%lf %d", &beta_pole, &ia_beta_pole);          fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* initial period (hrs) fixed or free */
   fscanf(f_1, "%lf %d", &prd, &ia_prd);                      fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* epoch of zero time */
   fscanf(f_1, "%lf", &jd_0);                                 fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* initial rotation angle fi0 */   
   fscanf(f_1, "%lf", &Phi_0);                                 fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* the weight factor for conv. reg. */
   fscanf(f_1, "%lf", &conw);                                 fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* degree and order of the Laplace series */
   fscanf(f_1, "%d %d", &Lmax, &Mmax);                        fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* nr. of triangulation rows per octant */   
   fscanf(f_1, "%d", &nrows);                                 fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* initial values for phase funct. params., fixed or free */   
   fscanf(f_1, "%lf %d", &par[1], &ia_par[1]);                fgets(str_temp, MAX_LINE_LENGTH, f_1);
   fscanf(f_1, "%lf %d", &par[2], &ia_par[2]);                fgets(str_temp, MAX_LINE_LENGTH, f_1);
   fscanf(f_1, "%lf %d", &par[3], &ia_par[3]);                fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* Initial Lambert coeff., fixed or free (Lommel-Seeliger=1) */
   fscanf(f_1, "%lf %d", &cl, &ia_cl);                        fgets(str_temp, MAX_LINE_LENGTH, f_1);
   /* maximum number of iterations (when > 1) or minimum difference in dev to stop (when < 1) */
   fscanf(f_1, "%lf", &stop_condition);                       fgets(str_temp, MAX_LINE_LENGTH, f_1);
   fclose(f_1);

   if (verbose == 1)
   {
      printf("%g %g   initial lambda, beta (%d,%d)\n", lambda_pole, beta_pole, ia_lambda_pole, ia_beta_pole);  
      printf("%f  initial period (hrs) (%d)\n", prd, ia_prd);   
      printf("%g epoch of zero time t0\n", jd_0);  
      printf("%g  initial fixed rotation angle fi0\n", Phi_0);  
      printf("%g  the weight factor for conv. reg.\n", conw);  
      printf("%d %d  degree and order of the Laplace series\n", Lmax, Mmax);  
      printf("%d  nr. of triangulation rows per octant\n", nrows);    
      printf("%g %g %g  initial guesses for phase funct. params. (%d,%d,%d)\n", par[1], par[2], par[3], ia_par[1], ia_par[2], ia_par[3]);  
      printf("%g  initial Lambert coeff. (Lommel-Seeliger part = 1) (%d)\n", cl, ia_cl);   
      printf("%g  stop condition\n", stop_condition);   
      printf("\n");   
   }      

   /* lightcurves + geometry file */   
   /* number of lightcurves */
   fscanf(stdin, "%d", &Lcurves);

   if (Lcurves > MAX_LC)
   {
      fprintf(stderr, "\nError: Number of lcs is greater than MAX_LC = %d\n", MAX_LC);
      fflush(stderr); exit(1);
   }

   al = vector_double(Lcurves);   
   
   ndata = 0; /* total number of data */
   k2 = 0;   /* index */
   al0 = al0_abs = PI; /* the smallest solar phase angle */
   ial0 = ial0_abs = -1; /* initialization, index of al0 */
   jd_min = 1e20; /* initial minimum JD */
   onlyrel = 1; 
       
   /* loop over lightcurves */   
   for (i = 1; i <= Lcurves; i++)
   {
      ave = 0; /* average */
      fscanf(stdin, "%d %d", &Lpoints[i], &i_temp); /* number of points in this lc., is it absolute or relative */
      Inrel[i] = 1 - i_temp;
      /* are there some calibrated lcs.? */
      if (Inrel[i] == 0)
         onlyrel = 0;

      if (Lpoints[i] > POINTS_MAX)
      {
         fprintf(stderr, "\nError: Number of lc. points is greater than POINTS_MAX = %d\n", POINTS_MAX);
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

	 /* select the minimum JD */
	 if (tim[ndata] < jd_min)
	    jd_min = tim[ndata];
	 
	 /* normals of distance vectors */
         e0len[ndata] = sqrt(e0[1]*e0[1] + e0[2]*e0[2] + e0[3]*e0[3]);
         elen[ndata] = sqrt(e[1]*e[1] + e[2]*e[2] + e[3]*e[3]);

         ave = ave + brightness[ndata];

	 /* normalization of distance vectors */
         for (k = 1; k <= 3; k++)
	 {
            ee[ndata][k] = e[k] / elen[ndata];
            ee0[ndata][k] = e0[k] / e0len[ndata];
         }
	 cos_alpha = dot_product(e, e0) / (elen[ndata] * e0len[ndata]); 
	 alpha[ndata] =	acos(cos_alpha); /* solar phase angle */
         if (j == 1)
	 {
            al[i] = alpha[ndata];
            /* Find the smallest solar phase al0 */
            if (al[i] < al0)  
	    {
               al0 = al[i];
               ial0 = ndata;
            }
            /* the same for abs. lcs. */
            if ((al[i] < al0_abs) && (Inrel[i] == 0))
	    {
	       al0_abs = al[i];
	       ial0_abs = ndata;
	     } 
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

   /* If input jd_0 <= 0 then the jd_0 is set to the epoch of the lowest JD in the data */
   if (jd_0 <= 0)
   {
      jd_0 = (int) jd_min;
      if (verbose == 1)
         printf("the used epoch of zero time  %f\n", jd_0);
   }      
      
   /* loop over data - subtraction of jd_0 */   
   for (i = 1; i <= ndata; i++)
      tim[i] = tim[i] - jd_0;         

   Phi_0 = Phi_0 * DEG2RAD;

   /* use calibrated data if possible */
   if (onlyrel == 0)
   {
      al0 = al0_abs;
      ial0 = ial0_abs;
   }      
   
   /* Initial shape guess:
      N.B. the initial guess should not be exactly a sphere (for which
      l=1-coeffs. are redundant, resulting in a bad first step)
      This wild rfit guess is for L-S at opposition (usually irrelevant) */
   a = A_ELL_INIT;
   b = B_ELL_INIT;
   c = C_ELL_INIT;
   rfit = sqrt(2 * sig[ial0] / (0.5 * PI * (1+cos(al0))));
   escl = rfit / sqrt((a * b + b * c + a * c) / 3);     
   if (onlyrel == 0)  /* for some reason it helps */
      escl *= 0.8;
   a = a * escl;
   b = b * escl;
   c = c * escl;
   
   /* Convexity regularization: make one last 'lightcurve' that
      consists of the three comps. of the residual nonconv. vect.
      that should all be zero */
   Lcurves = Lcurves + 1;
   Lpoints[Lcurves] = 3;
   Inrel[Lcurves] = 0;
   conw /= escl * escl; /* to make 'conw' parameter size independent */
   for (j = 1; j <= 3; j++)
   {
      ndata++;
      brightness[ndata] = 0;
      sig[ndata] = 1 / conw; /* the weight of conv. reg. */
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

   /* make indices to triangle vertices */
   trifac(nrows, ifp);
   /* areas and normals of the triangulated Gaussian image sphere */
   areanorm(t,f,ndir,Numfac, ifp, at, af);
   /* precompute some function values at each normal direction*/
   sphfunc(Numfac, at, af);
   /* sph. harm. coeffs. of the initial ell. */
   ellfit(cg,a,b,c,Numfac,Ncoef, at, af);

   cg[Ncoef+1] = beta_pole;
   cg[Ncoef+2] = lambda_pole;
   
   /* the formulas use beta measured from the pole */
   cg[Ncoef+1] = 90 - cg[Ncoef+1];
   /* use omega instead of period */
   cg[Ncoef+3] = 24 * 2 * PI / prd;
   /* conversion of lambda, beta to radians */
   cg[Ncoef+1] = DEG2RAD * cg[Ncoef+1];
   cg[Ncoef+2] = DEG2RAD * cg[Ncoef+2];
   /* give ia the value 0/1 if it's fixed/free */
   ia[Ncoef+1] = ia_beta_pole;
   ia[Ncoef+2] = ia_lambda_pole;
   ia[Ncoef+3] = ia_prd;
   /* phase function parameters */
   Nphpar = 3;
   for (i = 1; i <= Nphpar; i++)
   {
      cg[Ncoef+3+i] = par[i];
      ia[Ncoef+3+i] = ia_par[i];
   }
   /* Lommel-Seeliger part */
   cg[Ncoef+3+Nphpar+2] = 1;
   /* Use logarithmic formulation for Lambert to keep it positive */
   cg[Ncoef+3+Nphpar+1] = log(cl);
   /* shape is free to be optimized */
   for (i = 2; i <= Ncoef; i++)
      ia[i] = 1;
   /* The first shape param. fixed for relative br. fit */
   ia[1] = 0;
   /* is there any absolute lc.?, if yes, then the first shape param is free */
   if (onlyrel == 0)
      ia[1] = 1;
   ia[Ncoef+3+Nphpar+1] = ia_cl;
   /* Lommel-Seeliger part is fixed */
   ia[Ncoef+3+Nphpar+2] = 0;

   if ((Ncoef+3+Nphpar+1) > MAX_N_PAR)
   {
      fprintf(stderr, "\nError: Number of parameters is greater than MAX_N_PAR = %d\n", MAX_N_PAR);
      fflush(stderr); exit(1);
   }

   /* when to stop iterations? */
   n_iter_max = 0; 
   iter_diff_max = -1;
   rchisq = -1;
   if (stop_condition > 1)
   { 
      n_iter_max = (int) stop_condition;      
      iter_diff_max = 0;
   }      
   if (stop_condition < 1)
   { 
      n_iter_max = MAX_N_ITER; /* to avoid neverending loop */     
      iter_diff_max = stop_condition;
   }      

   if (verbose == 1)  printf("\n");

   /* initialization */   
   Alamda = -1;
   Niter = 0;
   iter_diff = 1e40;
   dev_old = 1e30;

   /* Levenberg-Marquardt loop */
   while ((Niter < n_iter_max) && (iter_diff > iter_diff_max))
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
          rchisq = Chisq - (pow(chck[1], 2) + pow(chck[2], 2) + pow(chck[3], 2)) * pow(conw,2);
      }
      dev_new = sqrt(rchisq / (ndata - 3));
      /* only if this step is better than the previous, 1e-10 is for numeric errors */
      if (dev_old - dev_new > 1e-10)
      {
	 iter_diff = dev_old - dev_new;
         dev_old = dev_new;
      }
      if (verbose == 1)
         printf("%d  chi2 %f  dev %f  alambda %f\n", Niter, rchisq, dev_new, Alamda);
   }

   /* period solution */
   prd = 2 * PI / cg[Ncoef+3];
   if (verbose == 1)
   {
      printf("\nlambda, beta and period (hrs): %f %+f %f\n", cg[Ncoef+2] * RAD2DEG , 90 - (cg[Ncoef+1] * RAD2DEG), 24 * prd);
      printf("phase function parameters: ");
      for (i = 1; i <= Nphpar; i++)
         printf("%f ", cg[Ncoef+3+i]);
   }	 
   if (verbose == 1)      
      printf("\nLambert coefficient: %g\n", exp(cg[Ncoef+Nphpar+4]));
   Alamda = 0;
   mrqmin(ee,ee0,tim,brightness,sig,cg,ia,Ncoef+5+Nphpar,covar,aalpha,bright);

   /* Write the model brightnesses to output file */
   Lastcall = 1;

   mrqmin(ee,ee0,tim,brightness,sig,cg,ia,Ncoef+5+Nphpar,covar,aalpha,bright);

   if (output == 1)
      if ((f_areas = fopen(argv[ind_out_file], "w")) == NULL)
         {fprintf(stderr, "\nError: cannot open 'output areas' file\n"); fflush(stdout); exit(2);}

   if (output == 1)
      fprintf(f_areas, "%d\n",  Numfac+1);

   if (std_out == 1)
      fprintf(stdout, "%d\n",  Numfac+1);
      
   totarea = 0;
   for (i = 1; i <= Numfac; i++)
   {
      if (output == 1)
      {
         fprintf(f_areas, "%14.12e\n", Area[i]);
         fprintf(f_areas, "%14.12e %14.12e %14.12e\n", Nor[i][1], Nor[i][2], Nor[i][3]);
      }
      if (std_out == 1)
      {
         fprintf(stdout, "%14.12e\n", Area[i]);
         fprintf(stdout, "%14.12e %14.12e %14.12e\n", Nor[i][1], Nor[i][2], Nor[i][3]);
      }	 
      totarea = totarea + Area[i];
   }
   sum = pow(chck[1],2) + pow(chck[2],2) + pow(chck[3],2);
   
   /* An additional facet to make the collection convex */
   dark = sqrt(sum);
   if (output == 1)
   {
      fprintf(f_areas, "%14.12e\n", dark);
      fprintf(f_areas, "%14.12e %14.12e %14.12e\n", -chck[1]/dark,-chck[2]/dark,-chck[3]/dark);
      fclose(f_areas);
   }

   if (std_out == 1)
   {
      fprintf(stdout, "%14.12e\n", dark);
      fprintf(stdout, "%14.12e %14.12e %14.12e\n", -chck[1]/dark,-chck[2]/dark,-chck[3]/dark);
   }
   
   if (verbose == 1)
      printf("plus a dark facet with area %3.2f%%\n\n", dark / totarea * 100);
   
   if ((f_lc_out = fopen(argv[ind_lc_out_file], "w")) == NULL)
      {fprintf(stderr, "\nError: cannot open 'output lc' file \n"); fflush(stdout); exit(4);}

   /* modelled brightness - output lcs. */
   k = 0;
   /* loop over lightcurves */   
   for (i = 1; i <= Lcurves - 1; i++)
   {
      /* loop over one lightcurve */
      for (j = 1; j <= Lpoints[i]; j++)
      {
         k++;
	 if (Inrel[i] == 1)    
            fprintf(f_lc_out, "%g\n", Yout[k] * Sclnw[i] / (cg[Ncoef+4] * exp(-alpha[k] / cg[Ncoef+5]) + cg[Ncoef+6] * alpha[k] + 1));
         else
            fprintf(f_lc_out, "%g\n", Yout[k]);
       }
    }
	 
   fclose(f_lc_out);
   
   /* output parameters */
   if (paramout == 1)
   {
      if ((f_par = fopen(argv[ind_parout_file], "w")) == NULL)
         {fprintf(stderr, "\nError: cannot open 'output par' file \n"); fflush(stdout); exit(2);}

      fprintf(f_par, "%5.1f %+5.1f %12.8f\n", cg[Ncoef+2] * RAD2DEG , 90 - (cg[Ncoef+1] * RAD2DEG), 24 * prd);
      fprintf(f_par, "%f %g\n", jd_0, Phi_0 * RAD2DEG);
      for (i = 1; i <= Nphpar; i++)
         fprintf(f_par, "%g ", cg[Ncoef+3+i]);
      fprintf(f_par, "\n%g\n", exp(cg[Ncoef+Nphpar+4]));
      
      fclose(f_par);
   }

   return(0);
}

   
