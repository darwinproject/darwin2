C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/quota/quota_cons.F,v 1.2 2013/06/20 21:47:20 jahn Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA

c=============================================================
c subroutine DARWIN_cons
c check conservation in the model
C==============================================================
         SUBROUTINE QUOTA_Cons(
     &                  myTime,myIter,myThid)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM.h"

#include "QUOTA_SIZE.h"
#include "QUOTA.h"
#include "DARWIN_IO.h"

#include "DYNVARS.h"

C     === Global variables ===
       _RL myTime
       INTEGER myIter
       INTEGER myThid

#ifdef CHECK_CONS
C============== Local variables ============================================
      _RL tmptotC
      _RL tmptotN
      _RL tmptotP
      _RL tmptotF
      _RL tmptotS
      _RL vol, voltot

      INTEGER i,j,k,bi,bj
      INTEGER ii, jp, ko, npnum

      _BEGIN_MASTER(myThid)
c
      tmptotC=0. _d 0
      tmptotN=0. _d 0
#ifdef PQUOTA
      tmptotP=0. _d 0
#endif
#ifdef FQUOTA
      tmptotF=0. _d 0
#endif
#ifdef SQUOTA
      tmptotS=0. _d 0
#endif

      voltot =0. _d 0

      DO bj=myByLo(myThid),myByHi(myThid)
      DO bi=myBxLo(myThid),myBxHi(myThid)

        DO i=1,sNx
        DO j=1,sNy
          do k= 1, NR
c-----------------------------------------------------------------
            vol=rA(i,j,bi,bj)*drF(k)
            voltot=voltot+vol
c.................................................................
c Carbon budget
            tmptotC=tmptotC+Ptracer(i,j,k,bi,bj,iDIC)*vol    !DIC
            do jp=1,npmax ! loop Plankton
              npnum=ibiomass-1 + (iCarb-1)*npmax + jp
              tmptotC=tmptotC+Ptracer(i,j,k,bi,bj,npnum)*vol !C biomass
            enddo
            do ko=1,komax ! loop OM
              npnum=iorgmat-1 + (iCarb-1)*komax + ko
              tmptotC=tmptotC+Ptracer(i,j,k,bi,bj,npnum)*vol !DOC/POC
            enddo
c.................................................................
c Nitrogen budget
            tmptotN=tmptotN+Ptracer(i,j,k,bi,bj,iNO3)*vol    !NO3
#ifdef AMMON
            tmptotN=tmptotN+Ptracer(i,j,k,bi,bj,iNH4)*vol    !NH4
#endif
#ifdef NITRITE
            tmptotN=tmptotN+Ptracer(i,j,k,bi,bj,iNO2)*vol    !NO2
#endif
            do jp=1,npmax ! loop Plankton
              npnum=ibiomass-1 + (iNitr-1)*npmax + jp
              tmptotN=tmptotN+Ptracer(i,j,k,bi,bj,npnum)*vol !N biomass
            enddo
            do ko=1,komax ! loop OM
              npnum=iorgmat-1 + (iNitr-1)*komax + ko
              tmptotN=tmptotN+Ptracer(i,j,k,bi,bj,npnum)*vol !DON/PON
            enddo
c.................................................................
#ifdef PQUOTA
c Phosphorous budget
            tmptotP=tmptotP+Ptracer(i,j,k,bi,bj,iPO4)*vol    !PO4
            do jp=1,npmax ! loop Plankton
              npnum=ibiomass-1 + (iPhos-1)*npmax + jp
              tmptotP=tmptotP+Ptracer(i,j,k,bi,bj,npnum)*vol !P biomass
            enddo
            do ko=1,komax ! loop OM
              npnum=iorgmat-1 + (iPhos-1)*komax + ko
              tmptotP=tmptotP+Ptracer(i,j,k,bi,bj,npnum)*vol !DOP/POP
            enddo
#endif
c.................................................................
#ifdef FQUOTA
c Iron budget
            tmptotF=tmptotF+Ptracer(i,j,k,bi,bj,iFeT)*vol    !Fe
            do jp=1,npmax ! loop Plankton
              npnum=ibiomass-1 + (iIron-1)*npmax + jp
              tmptotF=tmptotF+Ptracer(i,j,k,bi,bj,npnum)*vol !Fe biomass
            enddo
            do ko=1,komax ! loop OM
              npnum=iorgmat-1 + (iIron-1)*komax + ko
              tmptotF=tmptotF+Ptracer(i,j,k,bi,bj,npnum)*vol !DOFe/POFe
            enddo
#endif
c.................................................................
#ifdef SQUOTA
c Silica budget
            tmptotS=tmptotS+Ptracer(i,j,k,bi,bj,iSi)*vol     !Si
            do jp=1,npmax ! loop Plankton
              npnum=ibiomass-1 + (iSili-1)*npmax + jp
              tmptotF=tmptotF+Ptracer(i,j,k,bi,bj,npnum)*vol !Si biomass
            enddo
            npnum=iorgmat-1 + (iSili-1)*komax + 2
            tmptotS=tmptotS+Ptracer(i,j,k,bi,bj,npnum)*vol   !POSi
#endif
c------------------------------------------------------------
          enddo
        ENDDO
        ENDDO
      ENDDO
      ENDDO

      _GLOBAL_SUM_RL(tmptotC,myThid)
      print*,'tmptotC',tmptotC
      write(DAR_cons_unitC,*), myIter, tmptotC, tmptotC/voltot
      _GLOBAL_SUM_RL(tmptotN,myThid)
      print*,'tmptotN',tmptotN
      write(DAR_cons_unitN,*), myIter, tmptotN, tmptotN/voltot
#ifdef PQUOTA
      _GLOBAL_SUM_RL(tmptotP,myThid)
      print*,'tmptotP',tmptotP
      write(DAR_cons_unitP,*), myIter, tmptotP, tmptotP/voltot
#endif
#ifdef FQUOTA
      _GLOBAL_SUM_RL(tmptotF,myThid)
      print*,'tmptotF',tmptotF
      write(DAR_cons_unitF,*), myIter, tmptotF, tmptotF/voltot
#endif
#ifdef SQUOTA
      _GLOBAL_SUM_RL(tmptotS,myThid)
      print*,'tmptotS',tmptotS
      write(DAR_cons_unitS,*), myIter, tmptotS, tmptotS/voltot
#endif
      print*,'-----------------------------------'
      _END_MASTER(myThid)

C--   Everyone else must wait
      _BARRIER
#endif
c
c -----------------------------------------------------
      RETURN
      END
#endif  /*ALLOW_QUOTA*/
#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/

C============================================================================
