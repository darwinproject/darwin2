C $Header$
C $Name$

#include "DARWIN_OPTIONS.h"

c MONOD_FIELDS.h
c --------------------------------------------------------------------
c  Description:  fields for MONOD pkg
c --------------------------------------------------------------------      

C Chl_phy :: diagnostic Chlorophyll to be used in next time step
C            for attenuation, ....

#ifndef DYNAMIC_CHL
#ifdef GEIDER
      COMMON /DARWIN_CHL/
     &      Chl_phy          
         _RL Chl_phy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy,npmax)
#endif
#endif

