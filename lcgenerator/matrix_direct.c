/* rotation matrix */

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
#include "globals.h"
#include "constants.h"

void matrix_direct(double omg, double t, double tmat[][4])
{
   double f, cf, sf, fmat[4][4]; 
   
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
   /* Construct tmat (complete rotation matrix) */
   for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++)
      {
         tmat[i][j] = 0;
	 for (k = 1; k <= 3; k++)
	    tmat[i][j] = tmat[i][j] + fmat[i][k] * Blmat[k][j];
      }
}
