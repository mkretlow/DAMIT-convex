/* phase function */

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
#include <stdio.h>
#include "globals.h"

void phasec(double dcdp[], double alpha, double p[])
{
   double e, c;

   /* Exp-lin model (const.term=1.) */
   e = exp(-alpha / p[2]);
   c = 1 + p[1] * e + p[3] * alpha;
   /* derivatives */
   dcdp[1] = e;
   dcdp[2] = p[1] * e * alpha / (p[2] * p[2]);
   dcdp[3] = alpha;
   if (c < 0)
   {
      fprintf(stderr, "f<0 at alpha %f, ", alpha * RAD2DEG);
      fprintf(stderr, "phase param. %f %f %f \n", p[1], p[2], p[3]);
   }

   Scale = c;
}
