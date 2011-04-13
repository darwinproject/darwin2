C $Header$
C $Name$

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef GEIDER
#ifndef DYNCHL
c
c ====================================================================
c SUBROUTINE GEIDER97
c ====================================================================
c
       SUBROUTINE GEIDER97(
     I             PARlocal,
     I             nutrient,
     I             phytotempfunction,
     O             pcarbon, ! photosynthetic carbon uptake, or...
     O             Chl2C,
     I             myThid)
       IMPLICIT NONE

#ifdef ALLOW_QUOTA
#include "QUOTA_SIZE.h"
#include "QUOTA.h"
#else
#include "MONOD_SIZE.h"
#include "MONOD.h"
#endif
c
       INTEGER myThid
       INTEGER ii,jp
c
       _RL PARlocal
       _RL phytotempfunction(npmax)
       _RL nutrient(iimax)
c
       _RL pcarbon(npmax)
       _RL Chl2C(npmax)
c
       _RL E0
c
       _RL PCmax
       _RL PCPhot(npmax)
c
       _RL pday
       pday = 86400. _d 0
c
c Geider, MacIntyre and Kana - photosynthesis models
c GMK98 for Dynamic Quota

c convert units to match GMK98
c
c Variables...
c W/m^2 --> µmol photons/m^2/s
       E0 = PARlocal * 5. _d 0
c
       print*,"PARlocal = ",PARlocal
       print*,"E0 = ",E0
c
       do jp=1,npmax
c
         if (pcref(jp).gt.0. _d 0) then
           PCmax = pcref(jp)
     &           * nutrient(iNO3)
     &           / (nutrient(iNO3)+kn(iNO3,jp))
     &           * phytotempfunction(jp)
         else
           PCmax = 0. _d 0
         endif
c
         if (PCmax.gt.0. _d 0) then
c
           Chl2C(jp) = Chl2Cmax(jp)
     &               / (1 + Chl2Cmax(jp)*alphachl(jp)*E0/(2*PCmax))
c
           Chl2C(jp)=min(Chl2C(jp),chl2cmax(jp))
           Chl2C(jp)=max(Chl2C(jp),chl2cmin(jp))
c
           if (E0.gt.5. _d -1) then
             PCPhot(jp) = PCmax
     &                  * (1-exp(-alphachl(jp)*Chl2C(jp)*E0/PCmax))
c
             if (pcarbon(jp).lt. 0. _d 0)
     &          print*,'QQ ERROR pc=',jp,pcarbon(jp)
           else
             PCPhot=0. _d 0
           endif
c
           print*,"PCmax = ",PCmax
           print*,"PCref = ",PCref(jp)
           print*,"nitrogen = ",nutrient(iNO3)
           print*,"kn = ",kn(iNO3,jp)
           print*,"phyto T func = ",phytotempfunction(jp)
           print*,"Chl2C = ",Chl2C(jp)
           print*,"PCPhot = ",PCPhot(jp)
           print*," "
         else
           PCPhot(jp) = 0. _d 0
           Chl2C(jp)  = chl2cmin(jp)
         endif
c --------------------------------------------------------------------
c convert units for output back to ecomodel
c d-1 --> s^-1
         pcarbon(jp) = PCPhot(jp) / pday
c gChl/gC --> mgChl/mmolC
         Chl2C(jp) = Chl2C(jp) * 12. _d 0
       enddo

       RETURN
       END
#endif /*not DYNCHL*/
#endif /*GEIDER*/
#endif /*ALLOW_DARWIN*/
#endif /*ALLOW_PTRACERS*/
c ==================================================================
