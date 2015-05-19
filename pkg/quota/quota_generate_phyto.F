c $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/quota/quota_generate_phyto.F,v 1.4 2015/05/19 14:32:43 benw Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA

c ==========================================================
c SUBROUTINE QUOTA_GENERATE_PHYTO
c generate parameters for "Operational Taxonomic Units" of plankton (index jp)
c using an allometric approach
c
c         Ben Ward 2009/10
c ==========================================================
        SUBROUTINE QUOTA_GENERATE_PHYTO(myThid)

        implicit none
#include "EEPARAMS.h"
#include "DARWIN_PARAMS.h"
#include "QUOTA_SIZE.h"
#include "QUOTA.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
        INTEGER myThid

C     === Functions ===
      _RL DARWIN_RANDOM
      EXTERNAL DARWIN_RANDOM
      _RL DARWIN_RANDOM_NORMAL
      EXTERNAL DARWIN_RANDOM_NORMAL


C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf    - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf

        _RL RandNo
        _RL mortdays
        _RL year
        _RL rtime
        _RL standin
        _RL tmpsrt
        _RL tmpend
        _RL tmprng
        _RL iimaxm1
        _RL npmaxm1
        _RL komaxm1
        _RL prd_pry
        _RL factor
#ifdef ALLOWPFT
        _RL taxon_mu(npmax)
#endif
        _RL a,b,p,error
        _RL heterotrophy(npmax)
        _RL tau1,tau2
        _RL ESD1,pi
        _RL logvol(npmax)
        INTEGER ii,io,jp,ko,ni
        INTEGER jp2,icount,ntroph
        INTEGER signvar
CEOP
c
        standin= 0. _d 0
        pi     = 4. _d 0 * datan(1. _d 0)
c each time generate another functional group add one to ngroups
        ngroups = ngroups + 1

        iimaxm1 = float(iimax-1)
        npmaxm1 = float(npmax-1)
        komaxm1 = float(komax-1)
c
c..........................................................
c Generate plankton volumes and stochastic parameters
c..........................................................
#ifdef ALLOWPFT
        factor = 2. _d 0
        tau1=1.0 _d 0
        tau2=1.0 _d 0
c       Allocate Phytoplankton Taxa
c Prochloro
        do jp=1,2
          biovol(jp)    = 0.125 _d 0 * factor**(jp-1)
          autotrophy(jp)= 1.00 _d 0
          use_NO3(jp)   = 1
          use_Si(jp)    = 0
          taxon_mu(jp)  = 1.00 _d 0
          pft(jp)       = 1
        enddo
c Synnecho
        do jp=3,5
          biovol(jp)    = 0.50 _d 0 * factor**(jp-3)
          autotrophy(jp)= 1.00 _d 0
          use_NO3(jp)   = 1
          use_Si(jp)    = 0
          taxon_mu(jp)  = 1.40 _d 0
          pft(jp)       = 2
        enddo
c Small Euk
        do jp=6,9
          biovol(jp)    = 4.00 _d 0 * factor**(jp-6)
          autotrophy(jp)= 1.0 _d 0
          use_NO3(jp)   = 1
          use_Si(jp)    = 0
          taxon_mu(jp)  = 2.10 _d 0
          pft(jp)       = 3
        enddo
c Diatoms
        do jp=10,15
          biovol(jp)    = 128.0 _d 0 * factor**(jp-10)
          autotrophy(jp)= 1.0 _d 0
          use_NO3(jp)   = 1
          use_Si(jp)    = 0
          taxon_mu(jp)  = 3.8 _d 0
          pft(jp)       = 4
        enddo
c Specialist grazers
        do jp=16,16
          biovol(jp)    = 8.0 _d 0 * factor**(jp-16)
          autotrophy(jp)= 0.00 _d 0
          use_NO3(jp)   = 0
          use_Si(jp)    = 0
          taxon_mu(jp)  = 0.00 _d 0
          pft(jp)       = 6
        enddo
c
        do jp=1,16
          heterotrophy(jp)=1.0 _d 0 - autotrophy(jp)
        enddo
#else
        ESD1   = 0.5 _d 0 ! minimum plankton ESD
        ni     = 2 ! ni size classes within diameter gaps of * 10
        factor = 1000. _d 0 ** (1. _d 0 / float(ni))
        tau1   = 1.25 _d 0
        tau2   = 1.0 _d 0 / tau1
        ntroph = 2
        if (ntroph.eq.1) then
          tau1   = 0.0 _d 0
          tau2   = 0.0 _d 0
        endif
c Allocate plankton traits
        icount=0
        do jp2=1,ntroph
          do jp=1,npmax/ntroph
            icount = icount + 1
            biovol(icount)      = pi*ESD1**3/6. _d 0 *factor**(jp-1)
            logvol(icount)      = log10(biovol(icount))
            use_NO3(icount)     = 1
            use_Si(icount)      = 0
            if (ntroph.gt.1) then
              autotrophy(icount)  = 1.0 _d 0 - float(jp2-1)
     &                                       / float(ntroph - 1)
              heterotrophy(icount)= 1.0 _d 0 - autotrophy(icount)
            else
              autotrophy(icount)  = 0.5 _d 0
              heterotrophy(icount)= 0.5 _d 0
            endif
          enddo
        enddo
#endif
c ----------------------------------------------------------------------
c Allometry
#ifdef UNCERTAINTY
        error = 1.0 _d 0
#else
        error = 0.0 _d 0
!       set stdev of allometric parameters to zero
#endif
c ----------------------------------------------------------------------
        do jp=1,npmax
          ! parameters independent of nutrient element
c CARBON CONTENT
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_qcarbon),log10(ae_qcarbon)*error)
          call invnormal(b,p,b_qcarbon,be_qcarbon*error)
          qcarbon(jp)    = 10. _d 0**a * biovol(jp) ** b
c INITIAL SLOPE P-I
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_alphachl),log10(ae_alphachl)*error)
          call invnormal(b,p,b_alphachl,be_alphachl*error)
          alphachl(jp)   = 10. _d 0**a * biovol(jp) ** b
c RESPIRATION RATE
          p = darwin_random(myThid)
          IF (a_respir.NE.0. _d 0) THEN
          call invnormal(a,p,
     &         log10(a_respir),log10(ae_respir)*error)
          call invnormal(b,p,b_respir,be_respir*error)
          respiration(jp) =  10. _d 0**a * biovol(jp) ** b
          ELSE
            respiration(jp) = 0.0 _d 0
          ENDIF
c GRAZING SIZE PREFERENCE RATIO
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_prdpry),log10(ae_prdpry)*error)
          call invnormal(b,p,b_prdpry,be_prdpry*error)
          pp_opt(jp)      = 10. _d 0**a * biovol(jp) ** b
c MAXIMUM GRAZING RATE + WIDTH OF GRAZING KERNEL
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_graz),log10(ae_graz)*error)
          call invnormal(b,p,b_graz,be_graz*error)
!#ifdef ONEGRAZER
!! if only one grazer, set max grazing by 1 * prey size
!! only P are grazed
!          graz(jp) = 10. _d 0**a * (biovol(jp)*pp_opt(jp)) ** b
!     &             * autotrophy(jp) ** tau2
!! temp fix to remove size dep.
!          graz(jp) = 10. _d 0**a * (biovol(npmax)) ** b
!     &             * autotrophy(jp) ** tau2
!! temp fix to remove size dep.
!          pp_sig(jp) = 1. _d 10
!#else
! set grazing rate by grazer size, non-grazers to zero
          graz(jp) = 10. _d 0**a *  biovol(jp)             ** b
     &             * heterotrophy(jp) ** tau2
          pp_sig(jp)      = 2. _d 0
!#endif
c FRACTION GRAZED TAND MORTALITY TO DOM
          do io=1,iomax-iChl
#ifdef ALLOWPFT
            if (pft(jp).lt.3)  beta_graz(io,jp)=0.8
            if (pft(jp).gt.2)  beta_graz(io,jp)=0.5
#else
            beta_graz(io,jp) =  0.9 _d 0 - 0.7 _d 0
     &                       / (1.0 _d 0 + exp(-logvol(jp)+2.0 _d 0))
#endif
            beta_mort(io,jp) = beta_graz(io,jp)
          enddo
c GRAZING HALF-SATURATION
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_kg),log10(ae_kg)*error)
          call invnormal(b,p,b_kg,be_kg*error)
          kg(jp)         = 10. _d 0**a * biovol(jp) ** b
#ifdef DIFFLIMIT
     &                   * heterotrophy(jp) ** tau1
#endif
c PHYTOPLANKTON SINKING
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_biosink),log10(ae_biosink)*error)
          call invnormal(b,p,b_biosink,be_biosink*error)
          biosink(jp)    =  (10.0 _d 0**a) * biovol(jp) ** b
#ifdef ALLOWPFT
          if (pft(jp).eq.6) biosink(jp) = 0. _d 0
#endif
c MORTALITY
          ! constant background mortality
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_mort),log10(ae_mort)*error)
          call invnormal(b,p,b_mort,be_mort*error)
          kmort(jp)    =  (10.0 _d 0**a) * biovol(jp) ** b
! parameters relating to inorganic nutrients
          do ii=1,iimax
c MAXIMUM NUTRIENT UPTAKE RATE
            p = darwin_random(myThid)
            call invnormal(a,p,
     &           log10(a_vmaxi(ii)),log10(ae_vmaxi(ii))*error)
            call invnormal(b,p,b_vmaxi(ii),be_vmaxi(ii)*error)
            if (ii.eq.iDIC) then
#ifdef ALLOWPFT
              vmaxi(ii,jp)= 10. _d 0**a * biovol(jp) ** b
     &                    * taxon_mu(jp)
#else
              vmaxi(ii,jp)= (3.1 _d 0+logvol(jp))
     & / (5.0 _d 0 - 3.8 _d 0*logvol(jp) + logvol(jp)**2)
     &                    / 86400. _d 0
     &                    * autotrophy(jp) ** tau1
#endif
            else
               vmaxi(ii,jp)= 10. _d 0**a * biovol(jp) ** b
     &                     * autotrophy(jp) ** tau1
c NUTRIENT HALF-SATURATION CONSTANT
               p = darwin_random(myThid)
              call invnormal(a,p,
     &             log10(a_kn(ii)),log10(ae_kn(ii))*error)
              call invnormal(b,p,b_kn(ii),be_kn(ii)*error)
              kn(ii,jp) = 10. _d 0**a * biovol(jp) ** b
#ifdef DIFFLIMIT
     &                  * autotrophy(jp) ** tau1
#endif
            endif
          enddo
#ifdef SQUOTA
! Silicate parameters to zero for non-diatoms
          vmaxi(iSi,jp) = vmaxi(iSi,jp) * float(use_Si(jp))
#endif
c
          if (use_NO3(jp).eq.0) then
            ! prochlorocococcus can't use NO3
            vmaxi(iNO3,jp) = 0.0 _d 0
            ! but have higher NH4 affinity
            vmaxi(iNH4,jp) = vmaxi(iNH4,jp) * 2.0 _d 0
          endif
! parameters relating to quota nutrients
          do io=1,iomax-iChl
c EXCRETION
            if ((io.eq.iCarb.or.io.eq.iNitr.or.io.eq.iPhos)
     &          .and.a_kexc(io).NE.0. _d 0
     &          .and.ae_kexc(io).NE.0. _d 0) then
              p = darwin_random(myThid)
              call invnormal(a,p,
     &             log10(a_kexc(io)),log10(ae_kexc(io))*error)
              call invnormal(b,p,b_kexc(io),be_kexc(io)*error)
              kexc(io,jp) = 10. _d 0**a * biovol(jp) ** b
            else
              kexc(io,jp) = 0. _d 0
            endif
            if (io.ne.iCarb) then
c MINIMUM QUOTA
              p = darwin_random(myThid)
              call invnormal(a,p,
     &             log10(a_qmin(io)),log10(ae_qmin(io))*error)
              call invnormal(b,p,b_qmin(io),be_qmin(io)*error)
              qmin(io,jp) = 10. _d 0**a * biovol(jp) ** b
!     &                    * (autotrophy(jp)   ** tau2
!     &                     + heterotrophy(jp) ** tau2)
c MAXIMUM QUOTA
              p = darwin_random(myThid)
              call invnormal(a,p,
     &             log10(a_qmax(io)),log10(ae_qmax(io))*error)
              call invnormal(b,p,b_qmax(io),be_qmax(io)*error)
              qmax(io,jp) = 10. _d 0**a * biovol(jp) ** b
            endif
          enddo
#ifdef SQUOTA
          ! Silicate parameters to zero for non-diatoms
          qmin(iSili,jp) = qmin(iSili,jp) * float(use_Si(jp))
          qmax(iSili,jp) = qmax(iSili,jp) * float(use_Si(jp))
#endif
#ifdef ALLOWPFT
c Zooplankton have approximately Redfieldian N:C ratio
          if (pft(jp).eq.6) then
            qmin(iNitr,jp) = 0.0755 _d 0
            qmax(iNitr,jp) = 0.1510 _d 0
          endif
#endif
c PREFERENCE FUNCTION
          ! assign grazing preference according to predator/prey radius ratio
          do jp2=1,npmax ! jp2 denotes prey
#ifdef ALLOWPFT
            if (heterotrophy(jp).gt.0. _d 0.and.pft(jp2).ne.6) then
#else
            if (heterotrophy(jp).gt.0. _d 0) then
#endif
              prd_pry = biovol(jp) / biovol(jp2)
              graz_pref(jp,jp2) =
#ifdef SWITCHING
     &        1.0 _d 0
#else
     &        exp(-log(prd_pry/pp_opt(jp))**2 / (2*pp_sig(jp)**2))
#endif
              if (graz_pref(jp,jp2).lt.1. _d -4) then
                graz_pref(jp,jp2)=0. _d 0
              endif
              assim_graz(jp,jp2) = ass_eff
            else
              graz_pref(jp,jp2) = 0. _d 0
            endif
          enddo
c
c..........................................................
c generate phyto Temperature Function parameters
c.......................................................
          phytoTempCoeff(jp) = tempcoeff1
          phytoTempExp1(jp) = tempcoeff3
          phytoTempExp2(jp) = tempcoeff2_small
     &                      + (tempcoeff2_big-tempcoeff2_small)
     &                      * float(jp-1)/npmaxm1
          phytoTempOptimum(jp) = 2. _d 0
          phytoDecayPower(jp) = tempdecay

c..........................................................
        enddo


        RETURN
        END
#endif  /*ALLOW_QUOTA*/
#endif  /*ALLOW_DARWIN*/
#endif  /*ALLOW_PTRACERS*/

c ===========================================================
