/* Areas and normals of the triangulated ellipsoid */

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
#include "declarations.h"

void areanorm_cg(double t[], double f[], int ndir, int nfac, int **ifp, 
              double at[], double af[], double a_ell, double b_ell, double c_ell)
{
   int i, j;	  

   double  st, clen2, clen;

   double c[4], vx[4], vy[4], vz[4],
          *x, *y, *z;

   x = vector_double(ndir);
   y = vector_double(ndir);
   z = vector_double(ndir);
   
   for (i = 1; i <= ndir; i++)
   {
      st = sin(t[i]);
      x[i] = st * cos(f[i]) * a_ell;
      y[i] = st * sin(f[i]) * b_ell;
      z[i] = cos(t[i]) * c_ell;
   }
   for (i = 1; i <= nfac; i++)
   {
      /* vectors of triangle edges */
      for (j = 2; j <= 3; j++)
      {
         vx[j]=x[ifp[i][j]] - x[ifp[i][1]];
         vy[j]=y[ifp[i][j]] - y[ifp[i][1]];
         vz[j]=z[ifp[i][j]] - z[ifp[i][1]];
      }
      /* The cross product for each triangle */
      c[1]=vy[2]*vz[3]-vy[3]*vz[2];
      c[2]=vz[2]*vx[3]-vz[3]*vx[2];
      c[3]=vx[2]*vy[3]-vx[3]*vy[2];
      /* Areas (on the unit sphere) and normals */
      clen2=c[1]*c[1] + c[2]*c[2] + c[3]*c[3];
      clen=sqrt(clen2);
      /* normal */
      Nor[i][1]=c[1]/clen;
      Nor[i][2]=c[2]/clen;
      Nor[i][3]=c[3]/clen;
      /* direction angles of normal */
      at[i]=acos(Nor[i][3]);
      af[i] = atan2(Nor[i][2], Nor[i][1]);
      /* triangle area */
      Darea[i]= 0.5 * clen;
   }
   
   deallocate_vector((void *) x);
   deallocate_vector((void *) y);
   deallocate_vector((void *) z);
}

