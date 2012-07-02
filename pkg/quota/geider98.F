C $Header$
C $Name$

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"
c
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef GEIDER
#ifdef DYNCHL
c
c ====================================================================
c SUBROUTINE GEIDER98
c ====================================================================
       SUBROUTINE GEIDER98(
     I             PARlocal,
     I             biomass,
     I             qlimit,
#ifdef FQUOTA
     I             felimit,
#endif
#ifdef QUOTA_DIAG_LIMIT
     O             Ilim,
#endif
     I             up_inorg,
     O             PP,
     I             photo_Tempfunction,
     O             dchldt, ! chlorophyll synthesis rate
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
       _RL biomass(iomax,npmax)
       _RL qlimit(npmax)
       _RL felimit(npmax)
       _RL alpha_fe
       _RL up_inorg(iimax,npmax)
       _RL N_uptake(npmax)
       _RL photo_Tempfunction
c
       _RL dchldt(npmax)
c
       _RL E0
       _RL C_biomass
       _RL chl
       _RL VCref
c
       _RL Chl2C
       _RL PCmax
       _RL PCPhot(npmax)
       _RL rhochl
       _RL chlsynth(npmax)
       _RL VCmax
       _RL VCN
       _RL PP
#ifdef QUOTA_DIAG_LIMIT
       _RL Ilim(npmax)
#endif
c
c Geider, MacIntyre and Kana (1998) photosynthesis model
c adapted for multiple nutrients following Moore et al (2002)
c
       E0 = PARlocal ! muEin m^-2 s^-1
       PP = 0.0 _d 0
       do jp=1,npmax
#ifdef FQUOTA
         alpha_fe = alphachl(jp) * felimit(jp)
#else
         alpha_fe = alphachl(jp)
#endif
         if (autotrophy(jp).gt.0. _d 0) then
           C_biomass = biomass(iCarb,jp) ! mmol C m^-3
           chl = biomass(iChlo,jp)       ! mg Chl m^-3
c
c N uptake is sum of NO3, NO2 and NH4 uptake
           VCN = 0.0 _d 0
           do ii=2,iimax
             if (ii.eq.iNO3.or.ii.eq.iNO2.or.ii.eq.iNH4) then
               VCN = VCN + up_inorg(ii,jp) ! mmol N (mmol C)^-1 s^-1
             endif
           enddo
c
c-----------------------------------------------------------------
           if (E0.gt.1. _d -1 .and.
     &       vmaxi(iDIC,jp).gt.0. _d 0 .and.
     &       C_biomass.gt.0. _d 0 .and.
     &       chl.gt.0. _d 0) then
             Chl2C = chl / C_biomass                   ! mg chl (mmol C)^-1
c
             PCmax = vmaxi(iDIC,jp)                    ! s^-1
     &             * qlimit(jp)
     &             * photo_Tempfunction
c
             if (PCmax.gt.0. _d 0.and.alpha_fe.gt.0. _d 0) then
               PCPhot(jp) = PCmax                      ! s^-1
     &          *(1. _d 0-exp(-alpha_fe*Chl2C*E0/PCmax))
c
#ifdef QUOTA_DIAG_LIMIT
               Ilim(jp)=1. _d 0
     &                 -exp(-alphachl(jp)*Chl2C*E0/vmaxi(iDIC,jp))
#endif
c
               rhochl = Chl2Nmax * PCPhot(jp)          ! mg chl (mmol N)^-1
     &                /(alphachl(jp)*Chl2C*E0)
               chlsynth(jp) = rhochl * VCN * C_biomass ! mg chl m^-3 s^-1
             else
               PCPhot(jp)   = 0. _d 0                  ! s^-1
               chlsynth(jp) = 0. _d 0                  ! mg chl m^-3 s^-1
             endif
c
           else ! else if insufficient light, max phot=0 or C biomass=0...
             PCPhot(jp)   = 0. _d 0                    ! s^-1
             chlsynth(jp) = 0. _d 0                    ! mg chl m^-3 s^-1
           endif
c-----------------------------------------------------------------
           ! for passing back to quota_plankton
           PP = PP + PCPhot(jp) * C_biomass
           up_inorg(iDIC,jp) = PCPhot(jp) - biosynth*VCN  ! s^-1
           dchldt(jp)        = chlsynth(jp)  ! mg chl m^-3 s^-1
         else ! else if jp is a heterotroph
           up_inorg(iDIC,jp) = 0. _d 0       ! s^-1
           dchldt(jp)        = 0. _d 0       ! mg chl m^-3 s^-1
c-----------------------------------------------------------------
         endif
       enddo

       RETURN
       END
c
#endif /*DYNCHL*/
#endif /*GEIDER*/
#endif /*ALLOW_DARWIN*/
#endif /*ALLOW_PTRACERS*/
C ==================================================================
