C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/quota/quota_forcing.F,v 1.3 2013/12/27 17:29:00 jahn Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA

c=============================================================
c subroutine quota_forcing
c step forward bio-chemical tracers in time
C==============================================================
         SUBROUTINE QUOTA_FORCING(
     U                  Ptr,
     I                  bi,bj,imin,imax,jmin,jmax,
     I                  myTime,myIter,myThid)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#include "QUOTA_SIZE.h"
#include "QUOTA.h"
#include "DARWIN_IO.h"
#include "DYNVARS.h"
#ifdef USE_QSW
#include "FFIELDS.h"
#endif

C     === Global variables ===
c tracers
       _RL Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nDarwin)
       INTEGER bi,bj,imin,imax,jmin,jmax
       _RL myTime
       INTEGER myIter
       INTEGER myThid

C============== Local variables ============================================
c biomodel tracer arrays
      _RL nutrient(iimax)
      _RL biomass(iomax,npmax)
      _RL orgmat(iomax-iChl,komax)
#ifdef FQUOTA
c iron partitioning
      _RL freefe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL freefu
      _RL inputFel
#endif
c upstream arrays for sinking/swimming
      _RL bioabove(iomax,npmax)
      _RL biobelow(iomax,npmax)
      _RL orgabove(iomax-iChl,komax)
c some working variables
      _RL sumpy
      _RL sumpyup
c light variables
      _RL PAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sfac(1-OLy:sNy+OLy)
      _RL atten,lite
      _RL newtime     ! for sub-timestepping
      _RL runtim      ! time from tracer initialization
c
#ifdef ALLOW_DIAGNOSTICS
COJ for diagnostics
      _RL  PParr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
#ifdef ALLOW_TIMEAVE
#ifdef QUOTA_DIAG_LIMIT
      _RL Nlim(npmax)
      _RL Flim(npmax)
      _RL Ilim(npmax)
      _RL Tlim
#endif
#endif
c

c some local variables
       _RL Tlocal
       _RL Slocal
       _RL PARlocal
       _RL dzlocal
       _RL dtplankton
       _RL PP
c local tendencies
       _RL dbiomass(iomax,npmax)
       _RL dorgmat(iomax-iChl,komax)
       _RL dnutrient(iimax)
       _RL tmp

      INTEGER bottom
      INTEGER surface
      INTEGER i,j,k,it,itmp,ktmp
      INTEGER ii,io,jp,ko, jp2, jpsave
      INTEGER place
      INTEGER debug
      CHARACTER*8 diagname

c
c--------------------------------------------------
c initialise variables
      DO j=1-OLy,sNy+OLy
      DO i=1-OLx,sNx+OLx
       do k=1,Nr
#ifdef FQUOTA
           freefe(i,j,k) = 0.0 _d 0
# endif
           PAR(i,j,k)    = 0.0 _d 0
#ifdef ALLOW_DIAGNOSTICS
COJ for diagnostics
           PParr(i,j,k) = 0. _d 0
#endif
        enddo !k
       ENDDO !i
       ENDDO !j
c
c bio-chemical time loop
c--------------------------------------------------
      DO it=1,nsubtime
c -------------------------------------------------
COJ cannot use dfloat because of adjoint
COJ division will be double precision anyway because of dTtracerLev
         newtime=myTime-dTtracerLev(1)+
     &           float(it)*dTtracerLev(1)/float(nsubtime)
c        print*,'it  ',it,newtime,nsubtime,myTime
         runtim=myTime-float(PTRACERS_Iter0)*dTtracerLev(1)

#ifdef FQUOTA
c determine iron partitioning  - solve for free iron
         call darwin_fe_chem(bi,bj,iMin,iMax,jMin,jMax,
     &                       Ptr(1-OLx,1-OLy,1,bi,bj,iFeT), freefe,
     &                       myIter, mythid)
#endif

c find light in each grid cell
c ---------------------------
c determine incident light
#ifndef  READ_PAR
#ifdef USE_QSW
         DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           sur_par(i,j,bi,bj)=-parfrac*Qsw(i,j,bi,bj)*
     &                         parconv*maskC(i,j,1,bi,bj)
         ENDDO
         ENDDO
#else
         DO j=1-OLy,sNy+OLy
          sfac(j)=0. _d 0
         ENDDO
         call darwin_insol(newTime,sfac,bj)
         DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           sur_par(i,j,bi,bj)=sfac(j)*maskC(i,j,1,bi,bj)/86400. _d 6
c          if (i.eq.1.and.j.ge.1.and.j.le.sNy)
c    &                write(24,*) sur_par(i,j,bi,bj)
         ENDDO
         ENDDO
#endif
#endif

C.................................................................
C.................................................................


       DO j=1,sNy
        DO i=1,sNx
c surface PAR
c take ice coverage into account
#if (defined (ALLOW_SEAICE) && defined (USE_QSW))
COJ ice coverage already taken into account by seaice package
           lite=sur_par(i,j,bi,bj)
#else
#if (defined (ALLOW_SEAICE) && defined (USE_QSW))
c if using Qsw and seaice, then ice fraction is already
c taken into account
           lite=sur_par(i,j,bi,bj)
#else
           lite=sur_par(i,j,bi,bj)*(1. _d 0-fice(i,j,bi,bj))
#endif
#endif
           atten = 0. _d 0
           sumpy = 0. _d 0
c
c FOR EACH LAYER ...
        do k= 1, NR
          if (HFacC(i,j,k,bi,bj).gt.0. _d 0) then
c ---------------------------------------------------------------------
c benw
c
c Fetch biomodel variables from ptr (ptracers)
c (making sure they are .ge. 0  - brute force)
c
c (set biomodel tendencies to zero, at the same time)
c
c *********************************************************************
            place = 0
c Inorganic Nutrients
            do ii=1,iimax
              place = place + 1
c             ambient nutrients for each element (1 to iimax)
              nutrient(ii)  = max(Ptr(i,j,k,bi,bj,place),0. _d 0)
              dnutrient(ii) = 0. _d 0
            enddo ! ii
c *********************************************************************
c Unicellular biomass (including chlorophyll biomass - for non-grazers)
            do io=1,iomax
              do jp=1,npmax
                if (io.ne.iChlo.or.pft(jp).ne.6) then ! no grazer chlorophyll
                  place = place + 1
                  biomass(io,jp)  = max(Ptr(i,j,k,bi,bj,place),0. _d 0)
! biomasses above current layer for sinking
                  if (k.eq.1) then
                    bioabove(io,jp)=0. _d 0
                  endif
! biomasses below current layer for swimming
                  if (k.eq.Nr) then
                   biobelow(io,jp)=0. _d 0
                  elseif (hFacC(i,j,k+1,bi,bj).eq.0. _d 0) then
                   biobelow(io,jp)=0. _d 0
                  else
                   biobelow(io,jp)=max(Ptr(i,j,k+1,bi,bj,place),0. _d 0)
                  endif
! initialise biomass rate of change
                  dbiomass(io,jp) = 0. _d 0
                else ! if grazer, fill chl biomass with zeros
                  biomass(io,jp)  = 0. _d 0
                endif
              enddo ! jp
            enddo
c *********************************************************************
c Organic matter
            do io=1,iomax-iChl
              do ko=1,komax
c               mass of element x for all OM classes
                place = place + 1
                orgmat(io,ko) = max(Ptr(i,j,k,bi,bj,place),0. _d 0)
! biomasses above current layer for sinking
                if (k.eq.1) then
                 orgabove(io,ko) = 0. _d 0
                endif
#ifdef SQUOTA
                if (ko.and.1.and.io.eq.iSili) then
                  place = place - 1
                  orgmat(iSili,1)   = 0. _d 0
                  orgabove(iSili,1) = 0. _d 0
                endif
#endif
                dorgmat(io,ko) = 0. _d 0
              enddo ! ko
            enddo ! io
c *********************************************************************
c
c ---------------------------------------------------------------------


c find local light for level k
            sumpyup = sumpy
            sumpy = 0. _d 0
            do jp=1,npmax
#ifndef GEIDER
              ! sum nitrogen biomass
              sumpy = sumpy + biomass(iNitr,jp)
#else
              ! sum chlorophyll
              sumpy = sumpy + biomass(iChlo,jp)
#endif
            enddo

            atten= atten + (k_w + k_chl*sumpy)*5. _d -1*drF(k)
            if (k.gt.1)then
              atten = atten + (k_w+k_chl*sumpyup)*5. _d -1*drF(k-1)
            endif
            PAR(i,j,k) = lite*exp(-atten)
c
c Physical variables
            PARlocal = PAR(i,j,k)
            Tlocal = theta(i,j,k,bi,bj)
            Slocal = salt(i,j,k,bi,bj)
c Free Iron
#ifdef FQUOTA
            freefu = max(freefe(i,j,k),0. _d 0)
            if (k.eq.1) then
              inputFel = inputFe(i,j,bi,bj)
            else
              inputFel = 0. _d 0
            endif
#endif
c Layer thickness
            dzlocal = drF(k)*HFacC(i,j,k,bi,bj)
c
c set bottom=1.0 if the layer below is not ocean
            ktmp=min(nR,k+1)
            if(hFacC(i,j,ktmp,bi,bj).eq.0. _d 0.or.k.eq.Nr) then
              bottom = 1
            else
              bottom = 0
            endif
            if (k.eq.1) then
              surface = 1
            else
              surface = 0
            endif

c set other arguments to zero
            debug=0

            if (debug.eq.7) print*,'Inorganic nutrients',nutrient
            if (debug.eq.7) print*,'Plankton biomass', biomass
            if (debug.eq.7) print*,'Organic nutrients',orgmat
            if (debug.eq.8) print*,'k, PARlocal, dzlocal',
     &                                      k,PARlocal,dzlocal
c ---------------------------------------------------------------------
            CALL QUOTA_PLANKTON(
     I                       biomass, orgmat, nutrient,
     O                       PP,
     I                       bioabove, biobelow,
     I                       orgabove,
#ifdef FQUOTA
     I                       freefu, inputFel,
#endif
#ifdef ALLOW_TIMEAVE
#ifdef QUOTA_DIAG_LIMIT
     O                       Nlim, Flim, Ilim, Tlim,
#endif
#endif
     I                       PARlocal, Tlocal, Slocal,
     I                       bottom, surface, dzlocal,
     O                       dbiomass, dorgmat, dnutrient,
     I                       debug,
     I                       runtim,
     I                       MyThid)
c ---------------------------------------------------------------------
#ifdef FQUOTA
#ifdef IRON_SED_SOURCE
c only above minimum depth (continental shelf)
             if (rF(k).lt.depthfesed) then
c only if bottom layer
               if (HFacC(i,j,k+1,bi,bj).eq.0. _d 0) then
#ifdef IRON_SED_SOURCE_VARIABLE
c calculate sink of POC into bottom layer
                tmp=orgsink(2)*orgabove(iCarb,2)/dzlocal
c convert to dPOCl
                dnutrient(iFeT) = dnutrient(iFeT)
     &                          + fesedflux_pcm*tmp
#else
                dnutrient(iFeT) = dnutrient(iFeT)
     &                          + fesedflux/(drF(k)*hFacC(i,j,k,bi,bj))
#endif
               endif
             endif
#endif
#endif
c ---------------------------------------------------------------------
c save un-updated biomass as layer above
            do io=1,iomax
              do jp=1,npmax
                bioabove(io,jp)=biomass(io,jp)
              enddo
              if (io.ne.iChlo) then
                do ko=1,komax
                  orgabove(io,ko)=orgmat(io,ko)
                enddo
              endif
            enddo
c ---------------------------------------------------------------------
c now update main tracer arrays
c for timestep dtplankton
            dtplankton = dTtracerLev(k)/float(nsubtime)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
            place = 0
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c Inorganic nutrients
            do ii=1,iimax
              place = place + 1
              Ptr(i,j,k,bi,bj,place) = Ptr(i,j,k,bi,bj,place)
     &                               + dtplankton*dnutrient(ii)
            enddo ! ii
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c Biomass
            do io=1,iomax
              do jp=1,npmax
                if (io.ne.iChlo.or.pft(jp).ne.6) then ! if not a grazer
                  place = place + 1
                  Ptr(i,j,k,bi,bj,place) = Ptr(i,j,k,bi,bj,place)
     &                                   + dtplankton*dbiomass(io,jp)
                  if (pft(jp).eq.6.and.io.eq.iChlo) then
                    Ptr(i,j,k,bi,bj,place) = 0. _d 0
                  endif
                endif
              enddo ! jp
            enddo ! io
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Organic matter
            do io=1,iomax-iChl
              do ko=1,komax
                if (ko.ne.1.or.io.ne.iSili) then
                  place = place + 1
                  Ptr(i,j,k,bi,bj,place) = Ptr(i,j,k,bi,bj,place)
     &                                   + dtplankton*dorgmat(io,ko)
                endif
              enddo ! ko
            enddo ! io
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
#ifdef ALLOW_DIAGNOSTICS
COJ for diagnostics
            PParr(i,j,k) = PP
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_TIMEAVE
             PPave(i,j,k,bi,bj) = PPave(i,j,k,bi,bj)
     &                          + PP * dtplankton
            PARave(i,j,k,bi,bj) = PARave(i,j,k,bi,bj)
     &                          + PARlocal * dtplankton
c
#ifdef QUOTA_DIAG_LIMIT
            do jp=1,npmax
              Nlimave(i,j,k,bi,bj,jp) = Nlimave(i,j,k,bi,bj,jp)
     &                                + Nlim(jp) * dtplankton
              Flimave(i,j,k,bi,bj,jp) = Flimave(i,j,k,bi,bj,jp)
     &                                + Flim(jp) * dtplankton
              Ilimave(i,j,k,bi,bj,jp) = Ilimave(i,j,k,bi,bj,jp)
     &                                + Ilim(jp) * dtplankton
            enddo
            Tlimave(i,j,k,bi,bj) = Tlimave(i,j,k,bi,bj)
     &                           + Tlim * dtplankton
#endif
#endif
          endif
c end if hFac>0
         enddo ! k
c end layer loop
c
       ENDDO ! i
       ENDDO ! j
c
c
COJ fill diagnostics
#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics ) THEN
        diagname = 'PP      '
        CALL DIAGNOSTICS_FILL( PParr(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
       ENDIF
#endif
COJ

#ifdef FQUOTA
c determine iron partitioning  - solve for free iron
         call darwin_fe_chem(bi,bj,iMin,iMax,jMin,jMax,
     &                       Ptr(1-OLx,1-OLy,1,bi,bj,iFeT), freefe,
     &                       myIter, mythid)
#endif

c
#ifdef ALLOW_TIMEAVE
c save averages
       dar_timeave(bi,bj) = dar_timeave(bi,bj) + dtplankton
#endif
c
c -----------------------------------------------------
       ENDDO ! it
c -----------------------------------------------------
c end of bio-chemical time loop
c
        RETURN
        END

#endif  /*ALLOW_QUOTA*/
#endif  /*ALLOW_DARWIN*/
#endif  /*ALLOW_PTRACERS*/

C============================================================================
