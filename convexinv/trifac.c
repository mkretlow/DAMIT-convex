/* Form the vertex triplets of standard triangulation facets */

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

#include "globals.h"
#include "declarations.h"

void trifac(int nrows, int **ifp)
{
   int nnod, i, j, j0, j1, j2, j3, ntri;
   
   int **nod;

   nod = matrix_int(2*nrows, 4*nrows);

   nnod = 1; /* index from 1 to number of vertices */
   nod[0][0] = nnod;
   for (i = 1; i <= nrows; i++)
      for (j = 0; j <= 4 * i - 1; j++)
      {
         nnod++;
         nod[i][j] = nnod;
         if (j == 0) nod[i][4*i] = nnod;
      }
   for (i = nrows - 1; i >= 1; i--)
      for (j = 0; j <= 4 * i - 1; j++)
      {
         nnod++;
         nod[2*nrows-i][j] = nnod;
         if (j == 0) nod[2*nrows-i][4*i] = nnod;
       }

   nod[2*nrows][0] = nnod + 1;
   ntri = 0;

   for (j1 = 1; j1 <= nrows; j1++)
      for (j3 = 1; j3 <= 4; j3++)
      {
         j0 = (j3-1) * j1;
         ntri++;
         ifp[ntri][1] = nod[j1-1][j0-(j3-1)];
         ifp[ntri][2] = nod[j1][j0];
         ifp[ntri][3] = nod[j1][j0+1];                            
         for (j2 = j0 + 1; j2 <= j0 + j1 - 1; j2++)
	 {
            ntri++;
            ifp[ntri][1] = nod[j1][j2];
            ifp[ntri][2] = nod[j1-1][j2-(j3-1)];
            ifp[ntri][3] = nod[j1-1][j2-1-(j3-1)];
            ntri++;
            ifp[ntri][1] = nod[j1-1][j2-(j3-1)];
            ifp[ntri][2] = nod[j1][j2];
            ifp[ntri][3] = nod[j1][j2+1];
          }
       }

   /* Do the lower hemisphere */
   for (j1 = nrows + 1; j1 <= 2 * nrows; j1++)
      for (j3 = 1; j3 <= 4; j3++)
      {
         j0 = (j3 - 1) * (2 * nrows - j1);
         ntri++;
         ifp[ntri][1] = nod[j1][j0];
         ifp[ntri][2] = nod[j1-1][j0+1+(j3-1)];
         ifp[ntri][3] = nod[j1-1][j0+(j3-1)];               
         for (j2 = j0 + 1; j2 <= j0 + (2 * nrows - j1); j2++)
	 {
            ntri++;
            ifp[ntri][1] = nod[j1][j2];
            ifp[ntri][2] = nod[j1-1][j2+(j3-1)];
            ifp[ntri][3] = nod[j1][j2-1];
            ntri++;
            ifp[ntri][1] = nod[j1][j2];
            ifp[ntri][2] = nod[j1-1][j2+1+(j3-1)];
            ifp[ntri][3] = nod[j1-1][j2+(j3-1)];
         }
      }
      
   deallocate_matrix((void *) nod, 2*nrows);      
}
