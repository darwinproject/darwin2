C $Header$
C $Name$

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA

c===============================================================================
C===============================================================================
CStartofinterface
      SUBROUTINE QUOTA_INIT_FIXED(myThid)
C     =============== Global data ==========================================
C     === Global variables ===
         implicit none
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "GCHEM.h"
#include "DARWIN_PARAMS.h"
#include "QUOTA_SIZE.h"
#include "QUOTA.h"

        INTEGER myThid
C============== Local variables ============================================
C     msgBuf    - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL pday
      INTEGER i,j,k,bi,bj,ii,jp,ko
      INTEGER tmp
      INTEGER prec
      CHARACTER*(MAX_LEN_MBUF) fn
C     /--------------------------------------------------------------\
C     |  initialise common block biochemical parameters               |
C     \--------------------------------------------------------------/

      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') '// Quota loading parameters'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )


c c define 1 day in seconds
      pday = 86400.0 _d 0
c par parameters
      parfrac=  0.4 _d 0                    !fraction Qsw that is PAR
      parconv= 1. _d 0/0.2174 _d 0          !conversion from W/m2 to uEin/m2/s

      istar=90. _d 0                        ! w/m2
c
c
c set temperature function
       tempcoeff1 = 1. _d 0/3. _d 0
       tempcoeff2_small = 0.001 _d 0
       tempcoeff2_big   = 0.0003 _d 0
       tempcoeff3 = 1.04 _d 0
       tempmax = 20. _d 0      ! 32. _d 0
       temprange = 22. _d 0    ! 30. _d 0
       tempnorm = 0.3 _d 0  ! 1. _d 0
       tempdecay = 4. _d 0

c set up diagnostics
#ifdef ALLOW_MNC
      IF ( useMNC ) THEN
        CALL DARWIN_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

COJ set up diagnostics
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDIAGNOSTICS ) THEN
        CALL QUOTA_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
COJ
c-----------------------------------------------------
c "Standard" parameters
      ! maximum assimilation efficiency
      ass_eff          =  0.70 _d 0
      ! organic matter sinking rate
      orgsink(1)       =  0.00 _d 0  /pday
      orgsink(2)       =  10.0 _d 0  /pday
      ! cost of biosynthesis
      biosynth         =  2.33 _d 0
      ! nitrogen uptake
      amminhib         =  4.60 _d 0
      ! nitrification
      amm2nrite        =  2.00 _d 0  /pday
      nrite2nrate      =  0.10 _d 0  /pday
      PAR0             =  10.0 _d 0
      ! iron
!      alpfe            =  1.00 _d 0 ! 3D model
      alpfe            =  0.05 _d 0 ! 1D model
!      scav             =  1.10 _d -3 /pday
      scav             =  4.40 _d -3 /pday

      ligand_stab      =  2.00 _d +5
      ligand_tot       =  1.00 _d -3
      freefemax        =  0.10 _d -3
#ifdef IRON_SED_SOURCE
c iron sediment source
      depthfesed=1000.0 _d 0       !depth above which to add sediment source
      fesedflux =1.0 _d 0 * 1.0 _d -3 / (86400.0 _d 0) !iron flux (mmol/m2/s)
      fesedflux_pcm =0.68 _d 0 * 1.0 _d -3 !iron flux (mmol/m3/s) per
c                                            mmol POC/m3/s
#endif
      ! basic light
      k_w              =  0.04 _d 0
      k_chl            =  0.03 _d 0
      chl2nmax         =  3.00 _d 0 ! Moore 2002 (labelled as GD98)
!      chl2nmax         =  5.6 _d 0  ! Geider 1998
c-----------------------------------------------------
c ALLOMETRIC PARAMETERS
c parameter = aV^b
c if errors are relative (*/) then state as "log(error)"
c if errors are absolute (+-) then state as "error"
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Initial slope of PI curve
      a_alphachl       =  3.83 _d -7        ! mmol C (mg chl a)^-1 (mu Ein m^-2)^-1
      ae_alphachl      =  1.00 _d 0
      b_alphachl       =  0.00 _d 0
      be_alphachl      =  0.00 _d 0
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! maximum specific grazing rate
      a_graz           =  21.9 _d 0 / pday  ! d^-1
      ae_graz          =  1.00 _d 0
      b_graz           = -0.16 _d 0
      be_graz          =  0.02 _d 0
! half saturation grazing prey carbon concentration
!      a_kg             =  2.24 _d +1       ! mmol C m^-3
      a_kg             =  1.00 _d 0       ! mmol C m^-3
      ae_kg            =  3.80 _d 0
      b_kg             =  0.00 _d 0
      be_kg            =  0.03 _d 0
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! sinking (enter as positive downwards)
      a_biosink        =  0.28 _d -1 / pday ! m d^-1
      ae_biosink       =  1.00 _d 0
      b_biosink        =  0.39 _d 0
      be_biosink       =  0.00 _d 0
! swimming velocity (enter as positive upwards - converted in quota_generate_phyto.F)
      a_bioswim        =  0.00 _d 0 / pday ! m d^-1
      ae_bioswim       =  1.00 _d 0
      b_bioswim        =  0.18 _d 0
      be_bioswim       =  0.00 _d 0
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! mortality
      a_mort           =  0.02 _d 0 / pday  ! d^-1
      ae_mort          =  1.00 _d 0
      b_mort           =  0.00 _d 0
      be_mort          =  0.00 _d 0
! predator prey preference distribution parameters
      a_prdpry         =  1024. _d 0  ! dimensionless
      ae_prdpry        =  1.00 _d 0
      b_prdpry         =  0.00 _d 0
      be_prdpry        =  0.00 _d 0
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! carbon
      ! max photosynthetic rate (modified in quota_generate_plankton.F)
      a_vmaxi(iDIC)    =  1.00 _d 0 / pday   ! d^-1
      ae_vmaxi(iDIC)   =  1.00 _d 0
      b_vmaxi(iDIC)    = -0.15 _d 0
      be_vmaxi(iDIC)   =  0.05 _d 0
      ! cellular carbon content
      a_qcarbon        =  1.80 _d -11        ! mmol C cell^-1
      ae_qcarbon       =  1.36 _d 0
      b_qcarbon        =  0.94 _d 0
      be_qcarbon       =  0.04 _d 0
      ! respiration (Note function of cellular C --> aC^b!)
!      a_respir         =  3.21 _d -11 / pday ! mmol C cell^-1 d^-1
      a_respir         =  0.00 _d 0
      ae_respir        =  1.00 _d 0
      b_respir         =  0.93 _d 0
      be_respir        =  0.00 _d 0
      ! carbon excretion
!      a_kexc(iCarb)    =  0.32 _d -1 / pday  ! d^-1
      a_kexc(iCarb)    =  0.00 _d 0
      ae_kexc(iCarb)   =  1.36 _d 0
      b_kexc(iCarb)    = -0.33 _d 0
      be_kexc(iCarb)   =  0.09 _d 0
      ! fraction grazing to DOC
      a_beta_graz(iCarb)  =  1.00 _d 0
      ae_beta_graz(iCarb) =  1.00 _d 0
      b_beta_graz(iCarb)  = -0.40 _d 0
      be_beta_graz(iCarb) =  0.00 _d 0
      ! fraction mortality to DOC
      a_beta_mort(iCarb)  =  1.00 _d 0
      ae_beta_mort(iCarb) =  1.00 _d 0
      b_beta_mort(iCarb)  = -0.40 _d 0
      be_beta_mort(iCarb) =  0.00 _d 0
      ! carbon remineralisation rate
      remin(iCarb,1)   =  0.02 _d 0 /pday
      remin(iCarb,2)   =  0.04 _d 0 /pday
! nitrogen & nitrate
      ! maximum NO3 uptake rate
      a_vmaxi(iNO3)    =  0.51 _d 0 / pday   ! mmol N (mmol C)^-1 d^-1
      ae_vmaxi(iNO3)   =  1.00 _d 0
      b_vmaxi(iNO3)    = -0.27 _d 0
      be_vmaxi(iNO3)   =  0.00 _d 0
      ! NO3 half-saturation
      a_kn(iNO3)       =  0.17 _d 0         ! (mmol N m^-3)
      ae_kn(iNO3)      =  1.36 _d 0
      b_kn(iNO3)       =  0.27 _d 0
      be_kn(iNO3)      =  0.07 _d 0
      ! N minimum quota
      a_qmin(iNitr)    =  0.07 _d 0       ! mmol N (mmol C)^-1
      ae_qmin(iNitr)   =  1.00 _d 0
      b_qmin(iNitr)    = -0.17 _d 0
      be_qmin(iNitr)   =  0.00 _d 0
      ! N maximum quota
      a_qmax(iNitr)    =  0.25 _d 0        ! mmol N (mmol C)^-1
      ae_qmax(iNitr)   =  1.00 _d 0
      b_qmax(iNitr)    = -0.13 _d 0
      be_qmax(iNitr)   =  0.00 _d 0
      ! nitrogen excretion
      a_kexc(iNitr)    =  0.24 _d -1 / pday ! d^-1
      a_kexc(iNitr)    =  0.00 _d 0
      ae_kexc(iNitr)   =  1.36 _d 0
      b_kexc(iNitr)    = -0.33 _d 0
      be_kexc(iNitr)   =  0.09 _d 0
      ! fraction grazing to DON
      a_beta_graz(iNitr)  =  1.00 _d 0
      ae_beta_graz(iNitr) =  1.00 _d 0
      b_beta_graz(iNitr)  = -0.40 _d 0
      be_beta_graz(iNitr) =  0.00 _d 0
      ! fraction mortality to DON
      a_beta_mort(iNitr)  =  1.00 _d 0
      ae_beta_mort(iNitr) =  1.00 _d 0
      b_beta_mort(iNitr)  = -0.40 _d 0
      be_beta_mort(iNitr) =  0.00 _d 0
      ! N remineralisation rate
      remin(iNitr,1)   =  0.02 _d 0  /pday
      remin(iNitr,2)   =  0.04 _d 0  /pday
#ifdef NITRITE
! nitrite
      ! maximum NO2 uptake rate
      a_vmaxi(iNO2)    =  0.51 _d 0 / pday  ! mmol N (mmol C)^-1 d^-1
      ae_vmaxi(iNO2)   =  1.36 _d 0
      b_vmaxi(iNO2)    = -0.27 _d 0
      be_vmaxi(iNO2)   =  0.07 _d 0
      ! NO2 half-saturation
      a_kn(iNO2)       =  0.17 _d 0         ! (mmol N m^-3)
      ae_kn(iNO2)      =  1.36 _d 0
      b_kn(iNO2)       =  0.27 _d 0
      be_kn(iNO2)      =  0.08 _d 0
#endif
#ifdef AMMON
! ammonium
      ! maximum NH4 uptake rate
      a_vmaxi(iNH4)    =  0.26 _d 0 / pday  ! mmol N (mmol C)^-1 d^-1
      ae_vmaxi(iNH4)   =  1.36 _d 0
      b_vmaxi(iNH4)    = -0.27 _d 0
      be_vmaxi(iNH4)   =  0.07 _d 0
      ! NH4 half-saturation
      a_kn(iNH4)       =  0.85 _d -1        ! (mmol N m^-3)
      ae_kn(iNH4)      =  1.36 _d 0
      b_kn(iNH4)       =  0.27 _d 0
      be_kn(iNH4)      =  0.08 _d 0
#endif
#ifdef PQUOTA
! phosphate
      ! maximum PO4 uptake rate
      a_vmaxi(iPO4)    =  0.77 _d -1 / pday ! mmol P (mmol C)^-1 d^-1
      ae_vmaxi(iPO4)   =  1.36 _d 0
      b_vmaxi(iPO4)    = -0.27 _d 0
      be_vmaxi(iPO4)   =  0.07 _d 0
      ! PO4 half-saturation
      a_kn(iPO4)       =  0.26 _d -1        ! (mmol N m^-3)
      ae_kn(iPO4)      =  1.36 _d 0
      b_kn(iPO4)       =  0.27 _d 0
      be_kn(iPO4)      =  0.08 _d 0
      ! minimum P quota
      a_qmin(iPhos)    =  2.00 _d -3       ! mmol N (mmol C)^-1
      ae_qmin(iPhos)   =  1.00 _d 0
      b_qmin(iPhos)    =  0.00 _d 0
      be_qmin(iPhos)   =  0.00 _d 0
      ! maximum P quota
      a_qmax(iPhos)    =  0.01 _d 0        ! mmol N (mmol C)^-1
      ae_qmax(iPhos)   =  1.20 _d 0
      b_qmax(iPhos)    =  0.00 _d 0
      be_qmax(iPhos)   =  0.00 _d 0
      ! P excretion
      a_kexc(iPhos)    =  0.24 _d -1 / pday ! d^-1
      ae_kexc(iPhos)   =  1.36 _d 0
      b_kexc(iPhos)    = -0.33 _d 0
      be_kexc(iPhos)   =  0.09 _d 0
      ! fraction grazing to DOP
      a_beta_graz(iPhos)  =  1.00 _d 0
      ae_beta_graz(iPhos) =  1.00 _d 0
      b_beta_graz(iPhos)  = -0.40 _d 0
      be_beta_graz(iPhos) =  0.00 _d 0
      ! fraction mortality to DOP
      a_beta_mort(iPhos)  =  1.00 _d 0
      ae_beta_mort(iPhos) =  1.00 _d 0
      b_beta_mort(iPhos)  = -0.40 _d 0
      be_beta_mort(iPhos) =  0.00 _d 0
      ! P remineralisation rate
      remin(iPhos,1)   =  0.02 _d 0  /pday
      remin(iPhos,2)   =  0.04 _d 0  /pday
#endif
#ifdef SQUOTA
! silicate
      ! maximum Si uptake rate
      a_vmaxi(iSi)     =  0.77 _d -1 / pday ! mmol Fe (mmol C)^-1 d^-1
      ae_vmaxi(iSi)    =  1.36 _d 0
      b_vmaxi(iSi)     = -0.27 _d 0
      be_vmaxi(iSi)    =  0.07 _d 0
      ! Si half-saturation
      a_kn(iSi)        =  0.24 _d -1        ! (mmol N m^-3)
      ae_kn(iSi)       =  1.36 _d 0
      b_kn(iSi)        =  0.27 _d 0
      be_kn(iSi)       =  0.08 _d 0
      ! minimum Si quota
!      a_qmin(iSili)    =  0.84 _d -1        ! mmol N (mmol C)^-1
!      ae_qmin(iSili)   =  1.00 _d 0
!      b_qmin(iSili)    = -0.17 _d 0
!      be_qmin(iSili)   =  0.00 _d 0
      a_qmin(iSili)    =  2.00 _d -3       ! mmol N (mmol C)^-1
      ae_qmin(iSili)   =  1.00 _d 0
      b_qmin(iSili)    =  0.00 _d 0
      be_qmin(iSili)   =  0.00 _d 0
      ! maximum Si quota
!      a_qmax(iSili)    =  0.30 _d 0         ! mmol N (mmol C)^-1
!      ae_qmax(iSili)   =  1.00 _d 0
!      b_qmax(iSili)    = -0.13 _d 0
!      be_qmax(iSili)   =  0.00 _d 0
      a_qmax(iSili)    =  4.00 _d -3       ! mmol N (mmol C)^-1
      ae_qmax(iSili)   =  1.20 _d 0
      b_qmax(iSili)    =  0.00 _d 0
      be_qmax(iSili)   =  0.00 _d 0
      ! Si excretion
      a_kexc(iSili)    =  0.00 _d 0  / pday ! d^-1
      ae_kexc(iSili)   =  1.00 _d 0
      b_kexc(iSili)    =  0.00 _d 0
      be_kexc(iSili)   =  0.00 _d 0
      ! fraction grazing to DOSi
      a_beta_graz(iSili)  =  0.00 _d 0
      ae_beta_graz(iSili) =  1.00 _d 0
      b_beta_graz(iSili)  =  0.00 _d 0
      be_beta_graz(iSili) =  0.00 _d 0
      ! fraction mortality to DOSi
      a_beta_mort(iSili)  =  0.00 _d 0
      ae_beta_mort(iSili) =  1.00 _d 0
      b_beta_mort(iSili)  =  0.00 _d 0
      be_beta_mort(iSili) =  0.00 _d 0
      ! POSi remineralisation rate
      remin(iSili,1)   =  0.00 _d 0 /pday
      remin(iSili,2)   =  0.33 _d -2 /pday
#endif
#ifdef FQUOTA
! iron
      ! maximum Fe uptake rate
!      a_vmaxi(iFeT)    =  96.2 _d -6 / pday ! mmol Fe (mmol C)^-1 d^-1
      a_vmaxi(iFeT)    =  14.0 _d -6 / pday
      ae_vmaxi(iFeT)   =  1.36 _d 0
      b_vmaxi(iFeT)    = -0.27 _d 0
      be_vmaxi(iFeT)   =  0.07 _d 0
      ! Fe half-saturation
!      a_kn(ifeT)       =  32.1 _d -6        ! (mmol N m^-3)
      a_kn(ifeT)       =  80.0 _d -6
      ae_kn(iFeT)      =  1.36 _d 0
      b_kn(iFeT)       =  0.27 _d 0
      be_kn(iFeT)      =  0.08 _d 0
      ! minimum Fe quota
      a_qmin(iIron)    =  1.50 _d -6       ! mmol N (mmol C)^-1 - Mongin (2006)
!      a_qmin(iIron)    =  5.00 _d -6       ! mmol N (mmol C)^-1
      ae_qmin(iIron)   =  1.00 _d 0
      b_qmin(iIron)    =  0.00 _d 0
      be_qmin(iIron)   =  0.00 _d 0
      ! maximum Fe quota
      a_qmax(iIron)    =  80.0 _d -6       ! mmol N (mmol C)^-1 - Mongin (2006)
!      a_qmax(iIron)    =  15.0 _d -6       ! mmol N (mmol C)^-1
      ae_qmax(iIron)   =  1.20 _d 0
      b_qmax(iIron)    =  0.00 _d 0
      be_qmax(iIron)   =  0.00 _d 0
      ! Fe excretion
      a_kexc(iIron)    =  0.00 _d 0  / pday ! d^-1
      ae_kexc(iIron)   =  1.00 _d 0
      b_kexc(iIron)    =  0.00 _d 0
      be_kexc(iIron)   =  0.00 _d 0
      ! fraction grazing to DOFe
      a_beta_graz(iIron)  =  1.00 _d 0
      ae_beta_graz(iIron) =  1.00 _d 0
      b_beta_graz(iIron)  = -0.40 _d 0
      be_beta_graz(iIron) =  0.00 _d 0
      ! fraction mortality to DOFe
      a_beta_mort(iIron)  =  1.00 _d 0
      ae_beta_mort(iIron) =  1.00 _d 0
      b_beta_mort(iIron)  = -0.40 _d 0
      be_beta_mort(iIron) =  0.00 _d 0
      ! Fe remineralisation rate
      remin(iIron,1)   =  0.02 _d 0  /pday
      remin(iIron,2)   =  0.04 _d 0  /pday
#endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c end allometric scaling
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

        RETURN
        END
C============================================================================
#endif /*ALLOW_QUOTA*/
#endif /*ALLOW_DARWIN*/
#endif /*ALLOW_PTRACERS*/
