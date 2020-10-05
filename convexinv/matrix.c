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


void matrix(double omg, double t, double tmat[][4], double dtm[][4][4])
{
   double f, cf, sf, dfm[4][4], fmat[4][4]; 
   
   int i, j, k;
   
   /* phase of rotation */
   f = omg * t + Phi_0;
   f = fmod(f, 2 * PI); 
   cf = cos(f);
   sf = sin(f);
   /* rotation matrix, Z axis, angle f */ 
   fmat[1][1] = cf;
   fmat[1][2] = sf;
   fmat[1][3] = 0;
   fmat[2][1] = -sf;
   fmat[2][2] = cf;
   fmat[2][3] = 0;
   fmat[3][1] = 0;
   fmat[3][2] = 0;
   fmat[3][3] = 1;
   /* Ders. w.r.t omg */
   dfm[1][1] = -t * sf;
   dfm[1][2] = t * cf;
   dfm[1][3] = 0;
   dfm[2][1] = -t * cf;
   dfm[2][2] = -t * sf;
   dfm[2][3] = 0;
   dfm[3][1] = 0;
   dfm[3][2] = 0;
   dfm[3][3] = 0;
   /* Construct tmat (complete rotation matrix) and its derivatives */
   for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++)
      {
         tmat[i][j] = 0;
         dtm[1][i][j] = 0;
         dtm[2][i][j] = 0;
         dtm[3][i][j] = 0;
	 for (k = 1; k <= 3; k++)
	 {
            tmat[i][j] = tmat[i][j] + fmat[i][k] * Blmat[k][j];
            dtm[1][i][j] = dtm[1][i][j] + fmat[i][k] * Dblm[1][k][j];
            dtm[2][i][j] = dtm[2][i][j] + fmat[i][k] * Dblm[2][k][j];
            dtm[3][i][j] = dtm[3][i][j] + dfm[i][k] * Blmat[k][j];
          }
      }
}
