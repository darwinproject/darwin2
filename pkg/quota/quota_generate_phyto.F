Clphachl
c $Header$
C $Name$

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
        _RL month
        _RL fiveday
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
        _RL taxon_mu(npmax)
        _RL a,b,p,error
        _RL heterotrophy(npmax)
        INTEGER ii,io,jp,ko
        INTEGER jp2,icount
        INTEGER signvar
CEOP
c
        standin=0. _d 0

c each time generate another functional group add one to ngroups
        ngroups = ngroups + 1

        iimaxm1 = float(iimax-1)
        npmaxm1 = float(npmax-1)
        komaxm1 = float(komax-1)
c
c..........................................................
c Generate plankton volumes and stochastic parameters
c..........................................................
        factor = 2. _d 0
c       Allocate Phytoplankton Taxa
c Prochloro
        do jp=1,2
          biovol(jp)    = 1.25 _d -1 * factor**(jp-1)
          autotrophy(jp)= 1.00 _d 0
          use_NO3(jp)   = 0
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
          taxon_mu(jp)  = 1.10 _d 0
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
          biovol(jp)    = 64.0 _d 0 * factor**(jp-10)
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
c
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
          respiration(jp) = 10. _d 0**a
     &                    * (12. _d 9 * qcarbon(jp)) ** b
     &                    / qcarbon(jp)
          if (pft(jp).eq.6) then
            respiration(jp) = respiration(jp) * 0.50 _d 0
          endif
          ELSE
            respiration(jp) = 0.0 _d 0
          ENDIF
c MAXIMUM GRAZING RATE
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_graz),log10(ae_graz)*error)
          call invnormal(b,p,b_graz,be_graz*error)
          graz(jp)       = 10. _d 0**a * biovol(jp) ** b
     &                   * heterotrophy(jp)
c GRAZING SIZE PREFERENCE RATIO
          if (pft(jp).eq.5) then ! dinoflagellates prey upon similar sized plankton
            pp_opt(jp)      = 1.0 _d 0
          else ! other types follow standard relationship
            p = darwin_random(myThid)
            call invnormal(a,p,
     &           log10(a_prdpry),log10(ae_prdpry)*error)
            call invnormal(b,p,b_prdpry,be_prdpry*error)
            pp_opt(jp)      = 10. _d 0**a * biovol(jp) ** b
          endif
c STANDARD DEVIATION OF SIZE PREFERENCE
            pp_sig(jp) = 0.50 _d 0
c FRACTION GRAZED TO DOM
          do io=1,iomax-iChl
            p = darwin_random(myThid)
            call invnormal(a,p,
     &           log10(a_beta_graz(io)),log10(ae_beta_graz(io))*error)
            call invnormal(b,p,b_beta_graz(io),be_beta_graz(io)*error)
            beta_graz(io,jp) = 10. _d 0**a * biovol(jp) ** b
            beta_graz(io,jp)=min(max(
     &                       beta_graz(io,jp),0.1 _d 0),0.9 _d 0)
c FRACTION MORTALITY TO DOM
            p = darwin_random(myThid)
            call invnormal(a,p,
     &           log10(a_beta_mort(io)),log10(ae_beta_mort(io))*error)
            call invnormal(b,p,b_beta_mort(io),be_beta_mort(io)*error)
            beta_mort(io,jp)= 10. _d 0**a * biovol(jp) ** b
            beta_mort(io,jp)=min(max(
     &                       beta_mort(io,jp),0.1 _d 0),0.9 _d 0)
          enddo
c GRAZING HALF-SATURATION
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_kg),log10(ae_kg)*error)
          call invnormal(b,p,b_kg,be_kg*error)
          kg(jp)         = 10. _d 0**a * biovol(jp) ** b
c PHYTOPLANKTON SINKING
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_biosink),log10(ae_biosink)*error)
          call invnormal(b,p,b_biosink,be_biosink*error)
          if (pft(jp).eq.6) then
            biosink(jp)    =  0.0 _d 0 ! grazers don't sink
          else
            biosink(jp)    =  (10.0 _d 0**a) * biovol(jp) ** b
          endif
c SWIMMING
          p = darwin_random(myThid)
          IF (a_bioswim.NE.0. _d 0) THEN
          call invnormal(a,p,
     &         log10(a_bioswim),log10(ae_bioswim)*error)
          call invnormal(b,p,b_bioswim,be_bioswim*error)
          if (autotrophy(jp).eq.1) then
            bioswim(jp)    =  0.00 _d 0 ! only grazers can swim
          else
            bioswim(jp)    =-(10.0 _d 0**a) * biovol(jp) ** b
          endif
          ELSE
            bioswim(jp)    =  0.00 _d 0
          ENDIF
c MORTALITY
          p = darwin_random(myThid)
          call invnormal(a,p,
     &         log10(a_mort),log10(ae_mort)*error)
          call invnormal(b,p,b_mort,be_mort*error)
          kmort(jp)      = 10. _d 0**a * biovol(jp) ** b
!          if (pft(jp).eq.6) then
!            kmort(jp) = kmort(jp) / 8.0 _d 0
!          endif
! parameters relating to inorganic nutrients
          do ii=1,iimax
c MAXIMUM NUTRIENT UPTAKE RATE
            p = darwin_random(myThid)
            call invnormal(a,p,
     &           log10(a_vmaxi(ii)),log10(ae_vmaxi(ii))*error)
            call invnormal(b,p,b_vmaxi(ii),be_vmaxi(ii)*error)
            if (ii.eq.iDIC) then
              vmaxi(ii,jp)= 10. _d 0**a * biovol(jp) ** b
     &                    * taxon_mu(jp)
            else
               vmaxi(ii,jp)= 10. _d 0**a * biovol(jp) ** b
     &                     * autotrophy(jp)
               p = darwin_random(myThid)
              call invnormal(a,p,
     &             log10(a_kn(ii)),log10(ae_kn(ii))*error)
              call invnormal(b,p,b_kn(ii),be_kn(ii)*error)
              kn(ii,jp) = 10. _d 0**a * biovol(jp) ** b
!     *                  * autotrophy(jp)
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
c Zooplankton have approximately Redfieldian N:C ratio
          if (pft(jp).eq.6) then
            qmin(iNitr,jp) = 0.0755 _d 0
            qmax(iNitr,jp) = 0.1510 _d 0
          endif
c PREFERENCE FUNCTION
          ! assign grazing preference according to predator/prey radius ratio
          do jp2=1,npmax ! jp2 denotes prey
            if (heterotrophy(jp).gt.0. _d 0.and.pft(jp2).ne.6) then
              prd_pry = biovol(jp) / biovol(jp2)
              graz_pref(jp,jp2) =
#ifdef SWITCH3
              ! lower preference for larger P
!     &        1.0 _d 0
     &        biovol(jp2)**(-0.16 _d 0)
#else
     &        exp(-(log(prd_pry/pp_opt(jp))**2) / (2*pp_sig(jp)**2))
     &                                          / pp_sig(jp)/2. _d 0
#endif
c
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
