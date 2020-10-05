/* Curvature function (and hence facet area) from Laplace series */

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
#include "constants.h"

void curv(double cg[])
{
   int i, m, n, l, k;
   
   double fsum,
          g[MAX_N_FAC];
   
   for (i = 1; i <= Numfac; i++)
   {
      g[i] = 0;
      n = 0;
      for (m = 0; m <= Mmax; m++)
         for (l = m; l <= Lmax; l++)
	 {
            n++;
            fsum = cg[n] * Fc[i][m];
            if (m != 0) 
            {
	       n++;
               fsum = fsum + cg[n] * Fs[i][m];
            }
            g[i] = g[i] + Pleg[i][l][m] * fsum;
          }
      g[i] = exp(g[i]);
      Area[i] = Darea[i] * g[i];
      for (k = 1; k <= n; k++)
         Dg[i][k] = g[i] * Dsph[i][k];
   }

}
