C $Header$
C $Name$
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA

c ====================================================================
c SUBROUTINE QUOTA_PLANKTON
c 0. New version for mixotrophic QUOTA model
c
c 1. Local ecological interactions for models with many phytoplankton
c    "functional groups"
c 2. Timestep plankton and nutrients locally
c 3. Same equations for DOM and POM
c 4. Remineralization of detritus also determined in routine
c 5. Sinking particles and plankton
c 6. NOT in this routine: iron chemistry
c
c       Mick Follows, Scott Grant, Fall/Winter 2005
c       Stephanie Dutkiewicz Spring/Summer 2006
c       Ben Ward, 2009/10
c
c - Dropped in initial quota version...
c   - R* diagnostics (#define DAR_DIAG_RSTAR)
c   - diazotrophy (#define ALLOW_DIAZ)
c   - mutation code (#define ALLOW_MUTANTS)
c   - new nitrogen limiting scheme (#undef OLD_NSCHEME)
c   - diversity diagnostics
c   - waveband dependence of light attenuation and absorption
c ====================================================================

         SUBROUTINE QUOTA_PLANKTON(
     I                       biomass, orgmat, nutrient,
     O                       PP,
     I                       bioabove,
     I                       orgabove,
#ifdef FQUOTA
     I                       freefelocal, inputFelocal,
#endif
#ifdef QUOTA_DIAG_LIMIT
     O                       AP, HP,
     O                       Rlim, Ilim, photo_Tempfunction,
#endif
     I                       PARlocal, Tlocal, Slocal,
     I                       bottom, surface, dzlocal,
     O                       dbiomassdt,dorgmatdt, dnutrientdt,
     I                       debug,
     I                       runtim,
     I                       MyThid)


         implicit none
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DARWIN_PARAMS.h"
#include "QUOTA_SIZE.h"
#include "QUOTA.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP
c === GLOBAL VARIABLES =====================
c iimax     = no of nutrient species (1=Carbon,2=Nitrate, ...)
c npmax     = no of plankton species
c komax     = no of organic matter classes
c biomass  = plankton biomass (auto/mixo/hetero)
c x_num    = cell density
c quota    = (average) cell quota
c nutrient = ambient inorganic nutrient concntration
c orgmat   = organic matter biomass (dissolved and particulate)
         _RL biomass(iomax,npmax)
         _RL nutrient(iimax)
         _RL orgmat(iomax-iChl,komax)
         _RL quota(iomax-iChl,npmax)
         _RL PP
         _RL bioabove(iomax,npmax)
         _RL orgabove(iomax-iChl,komax)
#ifdef FQUOTA
         _RL freefelocal
         _RL inputFelocal
#endif
#ifdef QUOTA_DIAG_LIMIT
         _RL Rlim(iomax-iChl-1,npmax),Ilim(npmax),Tlim
         _RL AP(iomax,npmax),HP(iomax,npmax)
#endif
         _RL m_ref(npmax)
         _RL PARlocal
         _RL Tlocal
         _RL Slocal
         INTEGER bottom
         INTEGER surface
         _RL dzlocal
         INTEGER debug
         _RL dbiomassdt(iomax,npmax)
         _RL dchloro(npmax)
         _RL dnutrientdt(iimax)
         _RL dorgmatdt(iomax-iChl,komax)
#ifdef ALLOW_PAR_DAY
         _RL PARdaylocal
#endif
         _RL runtim

c LOCAL VARIABLES
c -------------------------------------------------------------

c WORKING VARIABLES
c ii = nutrient element index
         integer ii,io
c jp = plankton index
         integer jp,jp2
c ko = organic matter index
         integer ko
c 'spare' indices
         integer ii2, ko2
         integer jpred, jprey
         integer iin
c variables for plankton growth rate limitation
         _RL limit
         _RL qlimit(npmax)
#ifdef FQUOTA
         _RL felimit(npmax)
#endif
c plankton specific nutrient limitation terms
         _RL qreg(iomax,npmax)
         _RL q_reg(iomax,npmax)
c photosynthesis light limitation term
         _RL ilimit(npmax)
c temperature limitation terms
         _RL photo_Tempfunction
         _RL activ_Tempfunction
c uptake of inorganic nutrients
         _RL up_inorg(iimax,npmax)
c plankton grazing rates
         _RL grazing(iomax,npmax,npmax)
         _RL food1,food2,refuge(npmax)
         _RL biomass2(npmax)
c plankton respiration rates (carbon only)
         _RL C_resp(npmax)
c varible for mimumum phyto
         _RL biomassmin(npmax)
c variables for remineralization of organic matter
c
c variables for sinking
        _RL bvert(iomax,npmax)
        _RL pomsink(iomax-iChl)
        _RL fSurf,fBott
c variables for sums of plankton and organic matter
         _RL totplankton(iomax)
         _RL totorgmat(iomax-iChl)
c ?
         _RL facpz
         _RL kpar, kinh
         _RL tmp, tmpr,tmpz, tmpgrow, tmp1, tmp2
         integer ITEST
c
c *****************************************************************
c ******************** Evaluate Quota Status **********************
c *****************************************************************
c qlimit --> Function of most limiting nutrient
c            1 = replete, 0 = no C assimilation
c            N,Si   - linear from quota
c            P & Fe - Droop from quota
c
c qreg   --> individual nutrient status for uptake regualtion
c            0 = quota full, 1 = quota empty
c            linear for all elements (not needed for Si)
c            qreg(C,jp) = inverse of qreg for most limiting element
c
      do jp = 1,npmax
c       set diagnostic to zero
        qlimit(jp)     = 1. _d 0
        qreg(iCarb,jp) = 1. _d 0
c
        do io = 2,iomax-iChl ! skip carbon index; quota = X:C biomass ratio
c quota = nutrient biomass to carbon biomass ratio
          if (biomass(iCarb,jp).gt.0. _d 0) then
            quota(io,jp) = biomass(io,jp) / biomass(iCarb,jp)
          else
            quota(io,jp) = qmin(io,jp)
          endif
c limit ranges from 1 to 0 between qmin and qmax
          if (quota(io,jp).le.qmin(io,jp)) then ! if quota empty...
            limit       = 0. _d 0
            qreg(io,jp) = 1. _d 0
          elseif (quota(io,jp).ge.qmax(io,jp)) then ! if quota full...
            limit       = 1. _d 0
            qreg(io,jp) = 0. _d 0
          else ! if quota somewhere in between...
            if (io.eq.iNitr.or.io.eq.iSili) then ! linear
              limit = (quota(io,jp) - qmin(io,jp))
     &              / ( qmax(io,jp) - qmin(io,jp))
            else ! normalised Droop
              limit = (1. _d 0 - qmin(io,jp)/quota(io,jp))
     &              / (1. _d 0 - qmin(io,jp)/qmax(io,jp) )
            endif
            ! regulation term is always linear
            qreg(io,jp) = (qmax(io,jp) - quota(io,jp))
     &                  / (qmax(io,jp) - qmin(io,jp) )
          endif
          
#ifdef QUOTA_DIAG_LIMIT
          if (io.eq.iNitr) Rlim(iNitr-1,jp) = limit
#ifdef PQUOTA
          if (io.eq.iPhos) Rlim(iPhos-1,jp) = limit
#endif
#ifdef FQUOTA
          if (io.eq.iIron) Rlim(iIron-1,jp) = limit
#endif
#endif
#ifdef FQUOTA
          if (io.eq.iIron) then
            felimit(jp) = limit
            limit      = 1. _d 0
          endif
#endif
#ifdef SQUOTA
          ! non-diatoms are not Si limited
          if (io.eq.iSili.and.use_Si(jp).eq.0) then
            limit             = 1. _d 0
            qreg(iSili,jp)    = 0. _d 0
            biomass(iSili,jp) = 0. _d 0
            bioabove(iSili,jp)= 0. _d 0
            quota(iSili,jp)   = 0. _d 0
          endif
#endif
          qlimit(jp)     = min(qlimit(jp),limit)
          qreg(iCarb,jp) = min(qreg(iCarb,jp),1. _d 0 - qreg(io,jp))
c
          q_reg(io,jp) = qreg(io,jp) ** hill
        enddo ! io
        q_reg(iCarb,jp) = qreg(iCarb,jp) ** hill
#ifdef QUOTA_DIAG_LIMIT
        ! use monod limitation diagnostic where types are extinct
        if (biomass(iCarb,jp).le.0. _d 0) then
          Rlim(iNitr-1,jp) = nutrient(iNO3)/(nutrient(iNO3)+kn(iNO3,jp))
#ifdef PQUOTA
          Rlim(iPhos-1,jp) = nutrient(iPO4)/(nutrient(iPO4)+kn(iPO4,jp))
#endif
#ifdef FQUOTA
          Rlim(iIron-1,jp) = nutrient(iFeT)/(nutrient(iFeT)+kn(iFeT,jp))
#endif
        endif
#endif
        if (autotrophy(jp).eq. 0. _d 0) then
          biomass(iChlo,jp)  = 0. _d 0 ! pure heterotroph, so chl is zero
          bioabove(iChlo,jp) = 0. _d 0
        endif
      enddo ! jp
c
c ****************************************************************
c * Determine temperature Dependent Growth function for Plankton *
c ****************************************************************
       call quota_tempfunc(
     I             Tlocal,
     O             photo_Tempfunction,
     O             activ_Tempfunction,
     I             myThid)
c
c *****************************************************************
c ******************** Resource Acquisition ***********************
c *****************************************************************
      do jp=1,npmax
        if (autotrophy(jp).gt.0.0 _d 0) then
          do ii=2,iimax ! not carbon...
            if (ii.eq.iNO3.or.ii.eq.iNO2.or.ii.eq.iNH4) io=iNitr
                                        if (ii.eq.iPO4) io=iPhos
                                        if (ii.eq.iFeT) io=iIron
                                        if (ii.eq.iSi)  io=iSili
c C-specific nutrient uptake, modulated by quota fullness and temperature
            if (nutrient(ii).gt.0. _d 0) then
              up_inorg(ii,jp) = vmaxi(ii,jp)                          ! maximum uptake rate
     &                        * nutrient(ii)/(nutrient(ii)+kn(ii,jp)) ! ambient nutrients
     &                        * q_reg(io,jp)                          ! quota satiation
     &                        * activ_Tempfunction                    ! temperature effects
#ifdef AMMON
c             apply ammonium inhibition to NO3 and NO2
              if (ii.eq.iNO3.or.ii.eq.iNO2) then
                up_inorg(ii,jp) = up_inorg(ii,jp)
     &                          * exp(-amminhib*nutrient(iNH4))
              endif
#endif
            else
              up_inorg(ii,jp) = 0. _d 0
            endif
          enddo ! ii
        else ! if autotrophy(jp).eq.0
          do ii=1,iimax
            up_inorg(ii,jp) = 0. _d 0
          enddo
        endif ! autotrophy
      enddo ! jp
c
c ****************************************************************
c ************* Photosynthetic Carbon Assimilation ***************
c ****************************************************************
      PP = 0.0 _d 0
      call GEIDER98(
     I              PARlocal,
     I              biomass,
     I              qlimit,
#ifdef FQUOTA
     I              felimit,
#endif
#ifdef QUOTA_DIAG_LIMIT
     O              Ilim,
#endif
     U              up_inorg,
     O              PP,
     I              photo_Tempfunction,
     O              dchloro, ! chlorophyll synthesis rate
     I              myThid)
c
c ****************************************************************
c ********************* Heterotrophic Grazing ********************
c ****************************************************************
c PRE-ASSIMILATION grazing of type jpredator by type jprey
      do jpred=1,npmax ! loop predators
        if (autotrophy(jpred).lt.1. _d 0) then ! not for pure autotrophs
          food1 = 0.0 _d 0
          food2 = 0.0 _d 0
          do jprey=1,npmax ! sum all the prey carbon of predator, weighted by availability (preference)
            if (graz_pref(jpred,jprey).gt.0.0 _d 0) then
              food1 = food1
     &              + graz_pref(jpred,jprey)*biomass(iCarb,jprey)
#ifdef SWITCHING
              food2 = food2
     &              + (graz_pref(jpred,jprey)*biomass(iCarb,jprey))**ns
#endif
            endif
          enddo
          ! calculate grazing effort
          if (food1.gt.0. _d 0) then
            refuge(jpred) = (1.0 _d 0 - exp(Lambda * food1))
            tmp1  = activ_Tempfunction ! saturated grazing
     &            * food1 / (food1 + kg(jpred))      ! grazing effort
     &            * refuge(jpred)                    ! grazing refuge
          else
            tmp1  = 0. _d 0
          endif
          do jprey=1,npmax ! loop prey carbon consumption
            if (food1.gt.0. _d 0) then
              grazing(iCarb,jpred,jprey)                 ! d^-1
     &          = tmp1 ! grazing effort
#ifdef ONEGRAZER
     &          * graz(jprey) ! prey dependent maximum rate
#else
     &          * graz(jpred) ! predator dependent maximum rate
#endif
#ifdef SWITCHING
     &          *(graz_pref(jpred,jprey)*biomass(iCarb,jprey))**ns/food2
#else
     &          * graz_pref(jpred,jprey)*biomass(iCarb,jprey)     /food1
#endif
            else
              grazing(iCarb,jpred,jprey) = 0. _d 0
            endif
! other organic elements (+ chlorophyll) are grazed in stoichiometric relation to carbon
            if (grazing(iCarb,jpred,jprey).gt.0. _d 0
     &           .and.biomass(iCarb,jprey).gt.0. _d 0) then
              do io=2,iomax
                grazing(io,jpred,jprey) = grazing(iCarb,jpred,jprey) ! uptake of prey carbon
     &                                  * biomass(io,jprey)          ! *
     &                                  / biomass(iCarb,jprey)       ! biomass ratio of prey
              enddo
            else
              do io=1,iomax
                grazing(io,jpred,jprey) = 0. _d 0
              enddo
            endif
          enddo ! jprey
        else ! if pure autotrophs (i.e. autotrophy(jpred).eq.1)
          do io=1,iomax
            do jprey=1,npmax
              grazing(io,jpred,jprey) = 0. _d 0
            enddo
          enddo
        endif
      enddo ! jpred
c
c ************************************************************
c end evaluate biological process terms
c -----------------------------------------------------------------
c
c -----------------------------------------------------------------
c evaluate vertical sink terms
c ************************************************************
c     biosink is +ve downwards
c     (upstream - downstream) * vertical velocity
c     bvert is a gain term

!     factors to avoid sinking into surface layer, or out of bottom
      fSurf=float(1-surface)
      fBott=float(1-bottom)
!
      do io=1,iomax
!       plankton sinking
        do jp=1,npmax
            bvert(io,jp) = (fSurf*bioabove(io,jp)-fBott*biomass(io,jp))
     &                   *  biosink(jp) / dzlocal           ! + sinking in - sinking out
        enddo
!       organic matter sinking
        if (io.ne.iChlo) then
          pomsink(io) = (fSurf*orgabove(io,2)-fBott*orgmat(io,2))
     &                *  orgsink(2) / dzlocal
        endif
      enddo
c ************************************************************
c end evaluate vertical sink terms
c -----------------------------------------------------------------
c
c -------------------------------------------------------------------
c calculate tendency terms (and some diagnostics)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BIOMASS
c inorganic uptake
      do jp=1,npmax
        dbiomassdt(iCarb,jp) = biomass(iCarb,jp)*up_inorg(iDIC,jp)
        dbiomassdt(iNitr,jp) = biomass(iCarb,jp)*up_inorg(iNO3,jp)
#ifdef NITRITE
     &                       + biomass(iCarb,jp)*up_inorg(iNO2,jp)
#endif
#ifdef AMMON
     &                       + biomass(iCarb,jp)*up_inorg(iNH4,jp)
#endif
#ifdef PQUOTA
        dbiomassdt(iPhos,jp) = biomass(iCarb,jp)*up_inorg(iPO4,jp)
#endif
#ifdef SQUOTA
        dbiomassdt(iSili,jp) = biomass(iCarb,jp)*up_inorg(iSi,jp)
#endif
#ifdef FQUOTA
        dbiomassdt(iIron,jp) = biomass(iCarb,jp)*up_inorg(iFeT,jp)
#endif
#ifdef QUOTA_DIAG_LIMIT
        do io=1,iomax
          AP(io,jp) = dbiomassdt(io,jp)
        enddo
#endif
c
        dbiomassdt(iChlo,jp) = dchloro(jp)
c
! respiration
        dbiomassdt(iCarb,jp) = dbiomassdt(iCarb,jp)
     &                       - respiration(jp) *  biomass(iCarb,jp)
     &                       * activ_Tempfunction
c
c grazing uptake
        do io=1,iomax
#ifdef QUOTA_DIAG_LIMIT
          HP(io,jp) = 0.0 _d 0
#endif
          if (io.ne.iSili.and.io.ne.iChlo) then ! don't take up silicate or chlorophyll
c           Grazing uptake of everything but silicate and chlorophyll
            do jprey=1,npmax
              dbiomassdt(io,jp) = dbiomassdt(io,jp)
     &                          + biomass(iCarb,jp)   ! carbon biomass of predator
     &                          * grazing(io,jp,jprey)! * carbon specific rate
     &                          * assim_graz(jp,jprey) * q_reg(io,jp)
#ifdef QUOTA_DIAG_LIMIT
              HP(io,jp)         = HP(io,jp)
     &                          + biomass(iCarb,jp)   ! carbon biomass of predator
     &                          * grazing(io,jp,jprey)! * carbon specific rate
     &                          * assim_graz(jp,jprey) * q_reg(io,jp)
#endif
            enddo ! jprey
c           Exudation of elemental reservoirs
            if (io.ne.iChl.or.io.ne.iSili) then
              dbiomassdt(io,jp) = dbiomassdt(io,jp)
     &                          - kexc(io,jp) * biomass(io,jp)
            endif
          endif
!
! calculate (temperature adjusted?) mortality rates
          m_ref(jp) = kmort(jp) !* activ_Tempfunction
!
! Loss and sinking terms - include silicate and chlorophyll
            dbiomassdt(io,jp) = dbiomassdt(io,jp)
     &                        - biomass(io,jp)
     &                        * m_ref(jp)
     &                        + bvert(io,jp)
          do jpred=1,npmax
            ! - losses to other predators
            dbiomassdt(io,jp) = dbiomassdt(io,jp)
     &                        - biomass(iCarb,jpred)      ! carbon biomass of predator
     &                        * grazing(io,jpred,jp)      ! * carbon specific rate
          enddo ! jpred
        enddo ! io
      enddo ! jp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c NUTRIENTS
! uptake by phytoplankton
      do ii=1,iimax
        dnutrientdt(ii) = 0. _d 0
        do jp=1,npmax
          dnutrientdt(ii) = dnutrientdt(ii)
     &                    - biomass(iCarb,jp)*up_inorg(ii,jp) ! - uptake of inorganic nutrients
        enddo ! jp
      enddo ! ii
c
! remineralisation of organic matter
      do ko=1,komax
        dnutrientdt(iDIC) = dnutrientdt(iDIC)
     &                    + orgmat(iCarb,ko) * remin(iCarb,ko)
     &                                       * activ_Tempfunction
#ifndef AMMON
        dnutrientdt(iNO3) = dnutrientdt(iNO3)
     &                    + orgmat(iNitr,ko) * remin(iNitr,ko) ! straight to NO3, if no NH4 ...
     &                                       * activ_Tempfunction
#else
        dnutrientdt(iNH4) = dnutrientdt(iNH4)
     &                    + orgmat(iNitr,ko) * remin(iNitr,ko) ! ... or to NH4
     &                                       * activ_Tempfunction
#endif
#ifdef PQUOTA
        dnutrientdt(iPO4) = dnutrientdt(iPO4)
     &                    + orgmat(iPhos,ko) * remin(iPhos,ko)
     &                                       * activ_Tempfunction
#endif
#ifdef SQUOTA
        dnutrientdt(iSi)  = dnutrientdt(iSi)
     &                    + orgmat(iSili,ko) * remin(iSili,ko)
     &                                       * activ_Tempfunction
#endif
#ifdef FQUOTA
        dnutrientdt(iFeT) = dnutrientdt(iFeT)
     &                    + orgmat(iIron,ko) * remin(iIron,ko)
     &                                       * activ_Tempfunction
#endif
      enddo !ko
#ifdef FQUOTA
      dnutrientdt(iFeT) = dnutrientdt(iFeT)
     &                  - scav  * freefelocal          ! scavenging of free iron
     &                  + alpfe * inputFelocal/dzlocal ! atmospheric input
#endif

! respiration
      do jp=1,npmax
        dnutrientdt(iDIC) = dnutrientdt(iDIC)
     &                    + respiration(jp) *  biomass(iCarb,jp)
     &                    * activ_Tempfunction
      enddo

! oxidation of NH4 and NO2 compounds
#ifdef AMMON
      dnutrientdt(iNH4) = dnutrientdt(iNH4)
     &                  - amm2nrite * nutrient(iNH4)     ! ammonium to (nitrite or nitrate)
#ifdef NITRITE
      dnutrientdt(iNO2) = dnutrientdt(iNO2)
     &                  + amm2nrite * nutrient(iNH4)     ! nitrite from ammonium
     &                  - nrite2nrate * nutrient(iNO2)   ! nitrite to nitrate
      dnutrientdt(iNO3) = dnutrientdt(iNO3)
     &                  + nrite2nrate * nutrient(iNO2)   ! nitrate from nitrite
#else
      dnutrientdt(iNO3) = dnutrientdt(iNO3)              ! or
     &                  + amm2nrite * nutrient(iNH4)     ! nitrate from ammonium
#endif
#endif
c
c********************************************************************************
c organic matter
      do io=1,iomax-iChl
        dorgmatdt(io,1) = 0. _d 0    ! dissolved
        dorgmatdt(io,2) = pomsink(io) ! particulate
        do jp=1,npmax
!         mortality & excretion
          dorgmatdt(io,1) = dorgmatdt(io,1)
     &                    + biomass(io,jp)
     &                    * m_ref(jp)
     &                    * beta_mort(io,jp)
          if (io.eq.iCarb.or.io.eq.iNitr.or.io.eq.iPhos) then
            dorgmatdt(io,1) = dorgmatdt(io,1)
     &                      + kexc(io,jp) * biomass(io,jp)
          endif
          dorgmatdt(io,2) = dorgmatdt(io,2)
     &                    + biomass(io,jp)
     &                    * m_ref(jp)
     &                    *(1. _d 0-beta_mort(io,jp))
          do jprey=1,npmax
!           unassimilated grazing
            dorgmatdt(io,1) = dorgmatdt(io,1)
     &                      + biomass(iCarb,jp)
     &                      * grazing(io,jp,jprey)
     &                      *(1. _d 0-assim_graz(jp,jprey)*q_reg(io,jp))
     &                      * beta_graz(io,jprey)
            dorgmatdt(io,2) = dorgmatdt(io,2)
     &                      + biomass(iCarb,jp)
     &                      * grazing(io,jp,jprey)
     &                      *(1. _d 0-assim_graz(jp,jprey)*q_reg(io,jp))
     &                      *(1. _d 0-beta_graz(io,jprey))
          enddo ! jprey
        enddo ! jp
!       remineralisation of organic matter
        do ko=1,komax
           dorgmatdt(io,ko) = dorgmatdt(io,ko)
     &                      - orgmat(io,ko) * remin(io,ko)
     &                                       * activ_Tempfunction
        enddo ! ko
      enddo ! io
c********************************************************************************
c -------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c --------------------------------------------------------------------------
         RETURN
         END
#endif  /*ALLOW_QUOTA*/
#endif  /*ALLOW_DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================
