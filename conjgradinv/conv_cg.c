/*       Convexity regularization function */

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
#include <stdlib.h>
#include <stdio.h>
#include "globals.h"
#include "declarations.h"

double conv_cg(int nc, double p[])
{
   int i;
   
   double res, totarea, area;
   
   res = 0;
   totarea = 0;

   for (i = 1; i <= Numfac; i++)
   {
      area = exp(p[i]);
      res += area * Nor[i][nc];
      totarea += area;
   }

   /* normalization to the total surface area */
   res /= totarea;
      
   return(res);
}   
