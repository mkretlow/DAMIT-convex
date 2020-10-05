/* beta and lambda rotation matrix */

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

void blmatrix_direct(double bet, double lam)
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
}
