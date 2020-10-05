/* this computes brightness of a convex(!) shape  */

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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "globals.h"
#include "declarations.h"
#include "constants.h"

double bright_direct(double ee[], double ee0[], double t, double per, double par[], double cl, double cls)
{
   int i, j;
 
   double cos_alpha, br, alpha, dnom, scale, omg,
          e[4], e0[4],
	  mu[MAX_N_FAC+1], mu0[MAX_N_FAC+1+1], s[MAX_N_FAC+1], 
          tmat[4][4];
   
   cos_alpha = dot_product(ee, ee0);
   alpha = acos(cos_alpha);

   scale = phasec_direct(alpha,par); 

   omg = 2 * PI / (per / 24);
   
   matrix_direct(omg,t,tmat);

   /* Directions in the rotating system */
   for (i = 1; i <= 3; i++)
   {
      e[i] = 0;
      e0[i] = 0;
      for (j = 1; j <= 3; j++)
      {
         e[i] = e[i] + tmat[i][j] * ee[j];
         e0[i] = e0[i] + tmat[i][j] * ee0[j];
      }
   } 

   br = 0;
   /* Integrated brightness (phase coeff. used later) */
   for (i = 1; i <= Numfac; i++)
   {
      mu[i] = e[1] * Nor[i][1] + e[2] * Nor[i][2] + e[3] * Nor[i][3];
      mu0[i] = e0[1] * Nor[i][1] + e0[2] * Nor[i][2] + e0[3] * Nor[i][3];
      if((mu[i] > 0) && (mu0[i] > 0)) 
      {
         dnom = mu[i] + mu0[i];
         s[i] = mu[i] * mu0[i] * (cl + cls / dnom);
         br = br + Area[i] * s[i];
      }
    }

   /* phase angle function */
   br = br * scale;
   
   return(br);
}
