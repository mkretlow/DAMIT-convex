/* beta and lambda rotation matrix */

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

void blmatrix(double bet, double lam)
{
   double cb, sb, cl, sl;

   cb = cos(bet);
   sb = sin(bet);
   cl = cos(lam);
   sl = sin(lam);
   Blmat[1][1] = cb * cl;
   Blmat[1][2] = cb * sl;
   Blmat[1][3] = -sb;
   Blmat[2][1] = -sl;
   Blmat[2][2] = cl;
   Blmat[2][3] = 0;
   Blmat[3][1] = sb * cl;
   Blmat[3][2] = sb * sl;
   Blmat[3][3] = cb;
   /* Ders. of Blmat w.r.t. bet */
   Dblm[1][1][1] = -sb * cl;
   Dblm[1][1][2] = -sb * sl;
   Dblm[1][1][3] = -cb;
   Dblm[1][2][1] = 0;
   Dblm[1][2][2] = 0;
   Dblm[1][2][3] = 0;
   Dblm[1][3][1] = cb * cl;
   Dblm[1][3][2] = cb * sl;
   Dblm[1][3][3] = -sb;
   /* Ders. w.r.t. lam */
   Dblm[2][1][1] = -cb * sl;
   Dblm[2][1][2] = cb * cl;
   Dblm[2][1][3] = 0;
   Dblm[2][2][1] = -cl;
   Dblm[2][2][2] = -sl;
   Dblm[2][2][3] = 0;
   Dblm[2][3][1] = -sb * sl;
   Dblm[2][3][2] = sb * cl;
   Dblm[2][3][3] = 0;
}
