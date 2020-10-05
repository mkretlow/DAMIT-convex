#include <stdio.h>

#include "constants.h"

extern int Niter, Numfac, Lcurves,
           Lpoints[MAX_LC+1], Inrel[MAX_LC+1];
	   
    
extern double Phi_0, Scale, Rchisq, Conw,
              Area[MAX_N_FAC+1], Darea[MAX_N_FAC+1], Sclnw[MAX_LC+1], 
	      Nor[MAX_N_FAC+1][4], Blmat[4][4], Yout[MAX_N_OBS+1],
	      Sc_par[4], Cl, Cls, Beta, Lambda, Omg, **Ee, **Ee0,
	      *Brightness, *Tim, *Sig;
    
	   

