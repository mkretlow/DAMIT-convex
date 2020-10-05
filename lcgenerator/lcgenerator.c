/*     A procedure that computes lightcurves for a given shape+spin+scattering model.
       It does not compute shadowing, the shape has to be convex!!!
       
       syntax:
       
       cat lcs_file | lcgenerator [-v] input_parameters shape output_lcs
      
       options:	 -v verbose
*/       

/*
Copyright (C) 2006  Josef Durech, Mikko Kaasalainen

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
int Numfac, Lcurves, Lpoints[MAX_LC+1], Inrel[MAX_LC+1];
    
double Phi_0, Area[MAX_N_FAC+1], Nor[MAX_N_FAC+1][4], Blmat[4][4];
    
FILE *f_lc_out, *f_par, *f_shape;

/* -------------------------------------------------------------------------------*/
  
int main(int argc, char *argv[])
{
   int i, j, k, ndata, k2,   
       ind_par_file, ind_lc_out_file, ind_shape_file, arg_shift, verbose, 
       i_temp, n_ver, **fac;
          
   double jd_0, prd, cl, ave, *e0len, *elen, 
          beta_pole, lambda_pole, par[4], 
          *brightness, e[4], e0[4], **ee,
          **ee0, *tim, 
	  *x, *y, *z, g[4], h[4],
	  **d_vert, **e_vert, **f_vert,  **normal, *br_comp;

   ee = matrix_double(MAX_N_OBS,3);
   ee0 = matrix_double(MAX_N_OBS,3);
   
   tim = vector_double(MAX_N_OBS);
   brightness = vector_double(MAX_N_OBS);
   br_comp = vector_double(MAX_N_OBS);
   e0len = vector_double(MAX_N_OBS);
   elen = vector_double(MAX_N_OBS);

   /* command line arguments */
   ind_par_file = 1;
   ind_shape_file = 2;
   ind_lc_out_file = 3;
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
   ind_shape_file += arg_shift;
   ind_lc_out_file += arg_shift;

   /* input parameters file */
   if ((f_par = fopen(argv[ind_par_file], "r")) == NULL)
      {fprintf(stderr, "\nError: cannot open 'input parameters' file \n"); fflush(stdout); exit(2);}
   /* lambda, beta, P */
   fscanf(f_par, "%lf %lf %lf", &lambda_pole, &beta_pole, &prd);
   /* epoch of zero time and initial rotation angle */
   fscanf(f_par, "%lf %lf", &jd_0, &Phi_0);                     
   /* values for phase funct. params. */   
   fscanf(f_par, "%lf %lf %lf", &par[1], &par[2], &par[3]);     
   /* Lambert coeff. (Lommel-Seeliger=1) */
   fscanf(f_par, "%lf" , &cl);
   fclose(f_par);

   if (verbose == 1)
   {
      printf("%g %g lambda, beta\n", lambda_pole, beta_pole);  
      printf("%f  period (hrs)\n", prd);   
      printf("%f epoch of zero time t0\n", jd_0);  
      printf("%g  initial fixed rotation angle fi0\n", Phi_0);  
      printf("%g %g %g  phase funct. params.\n", par[1], par[2], par[3]);  
      printf("%g  Lambert coeff. (Lommel-Seeliger part = 1)\n", cl);   
      printf("\n");   
   }      

   /* lightcurves + geometry file */   
   /* number of lightcurves */
   fscanf(stdin, "%d", &Lcurves);

   if (Lcurves > MAX_LC)
   {  fprintf(stderr, "\nError: Number of lcs is greater than MAX_LC = %d\n", MAX_LC); fflush(stderr); exit(1); }

   ndata = 0; /* total number of data */
       
   /* loop over lightcurves */   
   for (i = 1; i <= Lcurves; i++)
   {
      fscanf(stdin, "%d %d", &Lpoints[i], &i_temp); /* number of points in this lc., is it absolute or relative */
      Inrel[i] = 1 - i_temp;

      if (Lpoints[i] > POINTS_MAX)
      {  fprintf(stderr, "\nError: Number of lc. points is greater than POINTS_MAX = %d\n", POINTS_MAX); fflush(stderr); exit(1); }

      /* loop over one lightcurve */
      for (j = 1; j <= Lpoints[i]; j++)
      {
         ndata++;
	 
         if (ndata > MAX_N_OBS)
         {  fprintf(stderr, "\nError: Number of data is greater than MAX_N_OBS = %d\n", MAX_N_OBS); fflush(stderr); exit(1); }

	 fscanf(stdin, "%lf %lf", &tim[ndata], &brightness[ndata]); /* JD, brightness */	 
	 fscanf(stdin, "%lf %lf %lf", &e0[1], &e0[2], &e0[3]); /* ecliptic astr_tempocentric coord. of the Sun in AU */
	 fscanf(stdin, "%lf %lf %lf", &e[1], &e[2], &e[3]); /* ecliptic astrocentric coord. of the Earth in AU */	 

	 /* normals of distance vectors */
         e0len[ndata] = sqrt(e0[1]*e0[1] + e0[2]*e0[2] + e0[3]*e0[3]);
         elen[ndata] = sqrt(e[1]*e[1] + e[2]*e[2] + e[3]*e[3]);

	 /* normalization of distance vectors */
         for (k = 1; k <= 3; k++)
	 {
            ee[ndata][k] = e[k] / elen[ndata];
            ee0[ndata][k] = e0[k] / e0len[ndata];
         }
       } /* j, one lightcurve */           
   } /* i, all lightcurves */        

   /* loop over data - subtraction of jd_0 */   
   for (i = 1; i <= ndata; i++)
      tim[i] = tim[i] - jd_0;         

   Phi_0 = Phi_0 * DEG2RAD;

   /* input shape */  
   if ((f_shape = fopen(argv[ind_shape_file], "r")) == NULL)
      {fprintf(stderr, "\nError: cannot open 'input shape' file \n"); fflush(stdout); exit(2);}

   /* number of vertices and facets */
   fscanf(f_shape, "%d %d\n", &n_ver, &Numfac);

   if (Numfac > MAX_N_FAC)
   {  fprintf(stderr, "\nError: Number of facets is greater than MAX_N_FAC = %d\n", MAX_N_FAC); fflush(stderr); exit(1); }

   fac = matrix_int(Numfac, 3);
   x = vector_double(n_ver);
   y = vector_double(n_ver);
   z = vector_double(n_ver);
   d_vert = matrix_double(Numfac,3);
   e_vert = matrix_double(Numfac,3);
   f_vert = matrix_double(Numfac,3);
   normal = matrix_double(Numfac,3);

   /* vertices */       
   for (i = 1; i <= n_ver; i++)
      fscanf(f_shape,"%lf %lf %lf\n", &x[i], &y[i], &z[i]);

   /* facets */
   for (i = 1; i <= Numfac; i++)
   {
      fscanf(f_shape,"%d %d %d\n", &fac[i][1], &fac[i][2], &fac[i][3]);
      
      /* triangles */  
      d_vert[i][1] = x[fac[i][1]];
      d_vert[i][2] = y[fac[i][1]];
      d_vert[i][3] = z[fac[i][1]];
      e_vert[i][1] = x[fac[i][2]];
      e_vert[i][2] = y[fac[i][2]];
      e_vert[i][3] = z[fac[i][2]];
      f_vert[i][1] = x[fac[i][3]];
      f_vert[i][2] = y[fac[i][3]];
      f_vert[i][3] = z[fac[i][3]];

      /* triangle edges */
      g[1] = e_vert[i][1] - d_vert[i][1];
      g[2] = e_vert[i][2] - d_vert[i][2];
      g[3] = e_vert[i][3] - d_vert[i][3];
      h[1] = f_vert[i][1] - d_vert[i][1];
      h[2] = f_vert[i][2] - d_vert[i][2];
      h[3] = f_vert[i][3] - d_vert[i][3];

      /* normals - right hand rule !!!*/
      cross_product(g, h, normal[i]);
      /* triangle area */
      Area[i] = sqrt(normal[i][1] * normal[i][1] + normal[i][2] * normal[i][2] + normal[i][3] * normal[i][3]) / 2;
      /* normalization of normals */
      Nor[i][1] = normal[i][1] / (2 * Area[i]);
      Nor[i][2] = normal[i][2] / (2 * Area[i]);
      Nor[i][3] = normal[i][3] / (2 * Area[i]);
   }      
     
   fclose(f_shape);

   /* beta lambda rotation matrix */
   blmatrix_direct(DEG2RAD * (90 - beta_pole), DEG2RAD * lambda_pole);

   /* computed lightcurves file */
   if ((f_lc_out = fopen(argv[ind_lc_out_file], "w")) == NULL)
      {fprintf(stderr, "\nError: cannot open 'output lightcurves' file \n"); fflush(stdout); exit(2);}

   /* compute lightcurves */   
   ndata = 0;
   k2 = 0;
   for (i = 1; i <= Lcurves; i++)
   {
      /* loop over one lightcurve */
      ave = 0;
      for (j = 1; j <= Lpoints[i]; j++)
      {
         ndata++;
         br_comp[ndata] = bright_direct(ee[ndata], ee0[ndata], tim[ndata], prd, par, cl, 1);
	 ave += br_comp[ndata];
      }
      ave /= Lpoints[i];
      for (j = 1; j <= Lpoints[i]; j++)
      {
         k2++;
	 /* normalization of relative lcs. to unit mean */
	 if (Inrel[i] == 1)
	    br_comp[k2] /= ave;
	 fprintf(f_lc_out, "%f\n", br_comp[k2]);    
      }
   } /* i, all lightcurves */        

   fclose(f_lc_out);
   
   return(0);
}

