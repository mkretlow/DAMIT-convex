/* phase function */

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
#include <stdio.h>
#include "globals.h"

double phasec_direct(double alpha, double p[])
{
   double e, c;

   /* Exp-lin model (const.term=1.) */
   e = exp(-alpha / p[2]);
   c = 1 + p[1] * e + p[3] * alpha;

   return(c);
}
