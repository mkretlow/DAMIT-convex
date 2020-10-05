/*     A procedure for shape modelling of asteroids.
       
       syntax:
       
       cat lcs_file | conjgrad [-v] [-s] [-o output_file] input_parameters input_rot_param out_lcs
       
       -v verbose
       -s write results to the standard output
       -o write results to the output_file
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
int Niter, Numfac, Lcurves, 
    Lpoints[MAX_LC+1], Inrel[MAX_LC+1];
       
double Phi_0, Scale, Rchisq, Conw, 
       Area[MAX_N_FAC+1], Darea[MAX_N_FAC+1],
       Yout[MAX_N_OBS+1],
       Nor[MAX_N_FAC+1][4], Blmat[4][4],
       Sc_par[4], Cl, Cls,
       *Brightness, **Ee, **Ee0, *Sig, *Tim, Beta, Lambda, Omg;
    
/*--------------------------------------------------------------*/

FILE *f_in, *f_par, *f_out, *f_lc_out;
  
int main(int argc, char *argv[])
{
   int i, j, k, nrows, ndata, k2, ndir, ial0, ial0_abs, 
       ind_param_file, arg_shift, verbose, output, std_out,
       i_temp, ind_input_file, ind_out_file, ind_lc_out_file, onlyrel,  
       n_iter_max, **ifp;
          
   double *t, *f, *at, *af, *weight_lc, darea[MAX_N_FAC+1],
          jd_0, conw, a = 1.05, b = 1.00, c = 0.95, prd, al0, al0_abs, ave, *e0len, *elen, cos_alpha,
          dth, dph, rfit, escl, totarea, sum, dark,
          e[4], e0[4], *cg, chck[4], *al, *alpha,
	  beta_pole, lambda_pole, rchisq, fret;

   char *str_temp;

   str_temp = (char *) malloc (MAX_LINE_LENGTH);

   Ee = matrix_double(MAX_N_OBS,3);
   Ee0 = matrix_double(MAX_N_OBS,3);
   ifp = matrix_int(MAX_N_FAC,4);
   
   Tim = vector_double(MAX_N_OBS);
   Brightness = vector_double(MAX_N_OBS);
   alpha = vector_double(MAX_N_OBS);
   Sig = vector_double(MAX_N_OBS);
   e0len = vector_double(MAX_N_OBS);
   elen = vector_double(MAX_N_OBS);
   cg = vector_double(MAX_N_FAC);
   t = vector_double(MAX_N_FAC);
   f = vector_double(MAX_N_FAC);
   at = vector_double(MAX_N_FAC);
   af = vector_double(MAX_N_FAC);
      
   /* recognizes arguments */
   ind_input_file = 1;
   ind_param_file = 2;
   ind_lc_out_file = 3;
   ind_out_file = 0;
   arg_shift = 0;
   verbose = 0;
   output = 0;
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
   }      
   ind_input_file += arg_shift;
   ind_param_file += arg_shift;
   ind_lc_out_file += arg_shift;

   /* input file */
   if ((f_in = fopen(argv[ind_input_file], "r")) == NULL)
      {fprintf(stderr, "\nError: cannot open 'input' file \n"); fflush(stderr); exit(1);}
      
   /* the weight factor for conv. reg. */
   fscanf(f_in, "%lf", &conw);                                 fgets(str_temp, MAX_LINE_LENGTH, f_in);
   /* nr. of triangulation rows per octant */   
   fscanf(f_in, "%d", &nrows);                                 fgets(str_temp, MAX_LINE_LENGTH, f_in);
   /* number of iterations */
   fscanf(f_in, "%d", &n_iter_max);                            fgets(str_temp, MAX_LINE_LENGTH, f_in);
   fclose(f_in);

   /* reads spin parameters from another input file */
   if ((f_par = fopen(argv[ind_param_file], "r")) == NULL)
      {fprintf(stderr, "\nError: cannot open 'input parameters' file \n"); fflush(stderr); exit(1);}

   /* lambda, beta, period  */
   fscanf(f_par, "%lf %lf %lf", &lambda_pole, &beta_pole, &prd);
   /* epoch of zero time t0, phi0 */
   fscanf(f_par, "%lf %lf", &jd_0, &Phi_0);
   /* phase funct. params. */   
   fscanf(f_par, "%lf %lf %lf", &Sc_par[1], &Sc_par[2], &Sc_par[3]);
   /* Lambert coeff. (L-S=1) */
   fscanf(f_par, "%lf", &Cl);
   fclose(f_par);
   
   if (verbose == 1)
   {
      printf("\n%g %g  lambda, beta\n", lambda_pole, beta_pole);  
      printf("%f  period (hrs)\n", prd);   
      printf("%f  epoch of zero time t0\n", jd_0);  
      printf("%g  initial rotation angle fi0\n", Phi_0);  
      printf("%g  the weight factor for conv. reg.\n", conw);  
      printf("%d  nr. of triangulation rows per octant\n", nrows);    
      printf("%g %g %g  phase function parameters\n", Sc_par[1], Sc_par[2], Sc_par[3]);  
      printf("%g  Lambert coeff. (Lommel-Seeliger part = 1) \n", Cl);   
      printf("%d  number of iterations\n", n_iter_max);   
      printf("\n");   
   }      

   /* lightcurves + geometry file */   
   fscanf(stdin, "%d", &Lcurves);

   if (Lcurves > MAX_LC)
   {  fprintf(stderr, "\nError: Number of lcs  is greater than MAX_LC = %d\n", MAX_LC); fflush(stderr); exit(2); }

   al = vector_double(Lcurves);   
   weight_lc = vector_double(Lcurves);   
   
   ndata = 0; /* total number of data */
   k2 = 0;   /* index */
   al0 = al0_abs = PI; /* the smallest solar phase angle */
   ial0 = ial0_abs = -1; /* initialization, index of al0 */
   onlyrel = 1;
         
   /* loop over lightcurves */   
   for (i = 1; i <= Lcurves; i++)
   {
      ave = 0; /* average */
      fscanf(stdin, "%d %d", &Lpoints[i], &i_temp); /* points in this lightcurve, is it absolute or relative? */
      Inrel[i] = 1 - i_temp;
      /* are there some calibrated lcs? */
      if (Inrel[i] == 0)
         onlyrel = 0;

      if (Lpoints[i] > POINTS_MAX)
      {  fprintf(stderr, "\nError: Number of lc points is greater than POINTS_MAX = %d\n", POINTS_MAX); fflush(stderr); exit(2); }

      /* loop over one lightcurve */
      for (j = 1; j <= Lpoints[i]; j++)
      {
         ndata++;
	 
         if (ndata > MAX_N_OBS)
         {  fprintf(stderr, "\nError: Number of data is greater than MAX_N_OBS = %d\n", MAX_N_OBS); fflush(stderr); exit(2); }

	 fscanf(stdin, "%lf %lf", &Tim[ndata], &Brightness[ndata]);
	 fscanf(stdin, "%lf %lf %lf", &e0[1], &e0[2], &e0[3]); /* ecliptic astrocentric coord. of the Sun in AU */
	 fscanf(stdin, "%lf %lf %lf", &e[1], &e[2], &e[3]); /* ecliptic astrocentric coord. of the Earth in AU */	 
         
	 /* normals of distance vectors */
         e0len[ndata] = sqrt(e0[1]*e0[1] + e0[2]*e0[2] + e0[3]*e0[3]);
         elen[ndata] = sqrt(e[1]*e[1] + e[2]*e[2] + e[3]*e[3]);

         ave += Brightness[ndata];

	 /* normalization of distance vectors */
         for (k = 1; k <= 3; k++)
	 {
            Ee[ndata][k] = e[k] / elen[ndata];
            Ee0[ndata][k] = e0[k] / e0len[ndata];
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
         Sig[k2] = ave ;
      }
   } /* i, all lightcurves */        

   /* loop over data - subtraction of jd_0 */   
   for (i = 1; i <= ndata; i++)
      Tim[i] = Tim[i] - jd_0;         

   Phi_0 = Phi_0 * DEG2RAD;

   /* use calibrated data if possible */
   if (onlyrel == 0)
   {
      al0 = al0_abs;
      ial0 = ial0_abs;
   }
         
   /* Initial shape guess */
   rfit = sqrt(2 * Sig[ial0] / (0.5 * PI * (1+cos(al0))));
   escl = rfit / sqrt((a * b + b * c + a * c) / 3);
   if (onlyrel == 0)
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
   conw *= 100; 
   Conw = conw;
   for (j = 1; j <= Lpoints[Lcurves]; j++)
   {
      ndata++;
      Brightness[ndata] = 0;
      Sig[ndata] = 1 / conw;
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
   {  fprintf(stderr, "\nError: Number of facets is greater than MAX_N_FAC!\n"); fflush(stderr); exit(2); }   

   /* makes indices to triangle vertices */
   trifac(nrows, ifp);
   /* areas and normals of the triangulated Gaussian image sphere */
   areanorm_cg(t,f,ndir,Numfac, ifp, at, af, a, b, c);

   Beta = beta_pole;
   Lambda = lambda_pole;
   
   /* The formulas use beta measured from the pole */
   Beta = 90 - Beta;
   /* Use omega instead of period */
   Omg = 24 * 2 * PI / prd;
   /* conversion of lambda, beta to radians */
   Beta *= DEG2RAD;
   Lambda *= DEG2RAD;
   /* phase function parameters */

   /* Lommel-Seeliger part */
   Cls = 1;

   Niter = 0;

   for (i = 1; i <= Numfac; i++)
      darea[i] = log(Darea[i]);

   blmatrix_cg();
   
   frprmn(darea, Numfac, &Niter, &fret, bright_all, dbright_all, n_iter_max, conw, ndata, verbose);

   totarea = 0;    
   for (j = 1; j <= Numfac; j++)
      totarea += exp(darea[j]);

   for (i = 1; i <= 3; i++)
   {
      chck[i] = 0;
      for (j = 1; j <= Numfac; j++)
         chck[i] = chck[i] + exp(darea[j]) * Nor[j][i];
   }

   if (verbose == 1)
   {
      rchisq = fret - (pow(chck[1],2) + pow(chck[2],2) + pow(chck[3],2)) * pow(conw,2) / totarea / totarea;
      printf("\nchisq %f,  rchisq %f,  dev %f\n", fret, rchisq, sqrt(rchisq / (ndata - 3))); 
   }

   /* output areas and normals */
   if (output == 1)
      if ((f_out = fopen(argv[ind_out_file], "w")) == NULL)
         {fprintf(stderr, "\nError: cannot open 'output areas' file \n"); fflush(stderr); exit(1);}

   if (output == 1)
      fprintf(f_out, "%d\n",  Numfac+1);

   if (std_out == 1)
      fprintf(stdout, "%d\n",  Numfac+1);

   totarea = 0;
   for (i = 1; i <= Numfac; i++)
   {
      if (output == 1)
      {
         fprintf(f_out, "%14.12e\n", exp(darea[i]));
         fprintf(f_out, "%14.12e %14.12e %14.12e\n", Nor[i][1], Nor[i][2], Nor[i][3]);
      }
      if (std_out == 1)
      {
         fprintf(stdout, "%14.12e\n", exp(darea[i]));
         fprintf(stdout, "%14.12e %14.12e %14.12e\n", Nor[i][1], Nor[i][2], Nor[i][3]);
      }

      totarea = totarea + exp(darea[i]);
   }
   sum = pow(chck[1],2) + pow(chck[2],2) + pow(chck[3],2);
   
   /* An additional facet to make the collection convex */
   dark = sqrt(sum);
   if (output == 1)
   {
      fprintf(f_out, "%14.12e\n", dark);
      fprintf(f_out, "%14.12e %14.12e %14.12e\n", -chck[1]/dark,-chck[2]/dark,-chck[3]/dark);
      fclose(f_out);
   }
   if (std_out == 1)
   {
      fprintf(stdout, "%14.12e\n", dark);
      fprintf(stdout, "%14.12e %14.12e %14.12e\n", -chck[1]/dark,-chck[2]/dark,-chck[3]/dark);
   }

   if (verbose == 1)
      printf("\nPlus a dark facet with area %3.2f%%\n\n", dark / totarea * 100);
      
   /* output lightcurves */
   if ((f_lc_out = fopen(argv[ind_lc_out_file], "w")) == NULL)
      {fprintf(stderr, "\nError: cannot open 'output lc' file \n"); fflush(stderr); exit(1);}

   k = 0;
   for (i = 1; i <= Lcurves - 1; i++)
      for (j = 1; j <= Lpoints[i]; j++)
      {
         k++;
         fprintf(f_lc_out, "%g\n", Yout[k]);
      }
   	 
   fclose(f_lc_out);

   return(0);
}

