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

#include <stdio.h>
#include <stdlib.h>
#include "globals.h"
#include "declarations.h"
#include "constants.h"

double bright_all(double a[])
{
   int i,j,l, np, np1, np2, jp, ic;

   double xx1[4], xx2[4],dy,sig2i, ymod, 
          ytemp[POINTS_MAX+1], dave[MAX_N_FAC+1], 
	  coef, ave = 0, trial_chisq, rchisq;

   trial_chisq = rchisq = 0;
   np = 0;
   np1 = 0;
   np2 = 0;

   for (i = 1; i <= Lcurves; i++)
   {
      if (Inrel[i] == 1) /* is the LC relative? */
      {
         ave = 0;
         for (l = 1; l <= Numfac; l++)
            dave[l]=0;
      }
      for (jp = 1; jp <= Lpoints[i]; jp++)
      {
         np++;
         for (ic = 1; ic <= 3; ic++) /* position vectors */
         {
            xx1[ic] = Ee[np][ic];
            xx2[ic] = Ee0[np][ic];
         }
	    
         if (i < Lcurves) 
            ymod = bright_cg(xx1,xx2,Tim[np],a);
         else
            ymod = conv_cg(jp,a);

         ytemp[jp] = ymod;
	    
         if (Inrel[i] == 1)
            ave = ave + ymod;
      } /* jp, lpoints */

      for (jp = 1; jp <= Lpoints[i]; jp++)
      {
         np1++;
         if (Inrel[i] == 1) 
         {
            coef = Sig[np1] * Lpoints[i] / ave;
            ytemp[jp] = coef * ytemp[jp];
         }
      }
      for (jp = 1; jp <= Lpoints[i]; jp++)
      {
         ymod = ytemp[jp];
         np2++;
         sig2i = 1 / (Sig[np2] * Sig[np2]);
         dy = Brightness[np2] - ymod;
         Yout[np2] = ymod;
         j = 0;
         trial_chisq = trial_chisq + dy * dy * sig2i;
         if (i < Lcurves)
            rchisq += dy * dy * sig2i;
        } /* jp */
        
   } /* i,  lcurves */

   Rchisq = rchisq;
   
   return trial_chisq;
}

