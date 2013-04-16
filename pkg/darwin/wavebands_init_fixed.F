C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/wavebands_init_fixed.F,v 1.5 2013/04/16 20:21:57 jahn Exp $
C $Name:  $

c ANNA wavebands_init_fixed.F reads-in and assigns input paramters for WAVEBANDS. 

#include "DARWIN_OPTIONS.h"

CBOP
C     !ROUTINE: WAVEBANDS_INIT_FIXED
C     !INTERFACE:
       subroutine wavebands_init_fixed(myThid)

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE WAVEBANDS_INIT_FIXED
C     | o reads-in and assigns input paramters for WAVEBANDS.
C     *==========================================================*
C     \ev

C     !USES:
       implicit none
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DARWIN_SIZE.h"
#include "SPECTRAL_SIZE.h"
#include "SPECTRAL.h"
#include "WAVEBANDS_PARAMS.h"
#include "DARWIN_IO.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid     :: my Thread Id number
      integer myThid      
CEOP

#ifdef WAVEBANDS

C     !LOCAL VARIABLES:
C     == Local variables ==
c local variables      
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      character*80 title
      integer iUnit
      integer swlambda,splambda,ssflambda
      _RL sap,sap_ps,sbp,sbbp
      _RL saw,sbw
      _RL ssf
c      _RL planck, c, hc, oavo, hcoavo, rlamm
#ifdef DAR_CALC_ACDOM
      _RL rlam
#else
      _RL sacdom
#endif
c local indeces
      integer nabp,i,ilam

C fill in missing waveband information:
C "representative values" darwin_waves need not be centered within
C waveband boundaries darwin_wavebands, so both may be given
C if representative values are not given, compute from waveband
C boundaries
      do i = 1,tlam
        if (darwin_waves(i) .gt. 0) then
          pwaves(i) = darwin_waves(i)
        elseif (darwin_wavebands(i).ge.0 .and.
     &          darwin_wavebands(i+1).ge.0 ) then
          pwaves(i) = .5*(darwin_wavebands(i)+darwin_wavebands(i+1))
        else
          WRITE(msgBuf,'(3A)') 'WAVEBANDS_INIT_FIXED: ',
     &    'please provide wavelengths in darwin_waves or ',
     &    'waveband boundaries in darwin_wavebands.'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
        endif
      enddo

C if waveband boundaries not given, compute from representative values
C these will be used to compute waveband widths
      do i=1,tlam+1
        if (darwin_wavebands(i).LT.0) then
C         put boundaries half-way between central values
C         but first and last boundary are at first and last "central" value
          if (i.eq.1) then
            darwin_wavebands(i) = pwaves(1)
          elseif (i.le.tlam) then
            darwin_wavebands(i) = .5*(pwaves(i-1)+pwaves(i))
          else
            darwin_wavebands(i) = pwaves(tlam)
          endif
        endif
      enddo

C waveband widths used to compute total PAR and alpha_mean
      wb_totalWidth = 0.0
      do i=1,tlam
        wb_width(i) = darwin_wavebands(i+1) - darwin_wavebands(i)
        wb_totalWidth = wb_totalWidth + wb_width(i)
C       allow for zero-width wavebands...
        if (wb_width(i).LT.0) then
          WRITE(msgBuf,'(2A,I3)') 'WAVEBANDS_INIT_FIXED: ',
     &    'negative waveband width encountered, waveband: ',i
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: ',
     &    'Please check darwin_waves and darwin_wavebands.'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED: wavebands:'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: ',
     &    ' idx       low   rep      high    width'
          CALL PRINT_ERROR( msgBuf, myThid )
          do ilam=1,tlam
           WRITE(msgBuf,'(A,I4,F10.3,I6,F10.3,F9.3)')
     &     'WAVEBANDS_INIT_FIXED: ', ilam,darwin_wavebands(ilam),
     &     pwaves(ilam),darwin_wavebands(ilam+1),wb_width(ilam)
           CALL PRINT_ERROR( msgBuf, myThid )
          enddo
          STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
        endif
      enddo
C     ...but require at least one non-zero-width band
      if (wb_totalWidth.LE.0) then
#ifdef OASIM
       if (tlam.eq.1) then
C        in this case it doesn't matter, just make non-zero for
C        alpha_mean computation
         wb_width(1) = 1.
         wb_totalWidth = 1.
       else
#else
       if (.TRUE.) then
#endif
        WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: ',
     &  'need to provide waveband boundaries in darwin_wavebands.'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
       endif
      endif


c  Water data files
      if (darwin_waterabsorbFile .NE. ' '  ) THEN
        CALL MDSFINDUNIT( iUnit, myThid )
        open(iUnit,file=darwin_waterabsorbFile,
     &                           status='old',form='formatted')  
        do i = 1,6                 ! six lines of text for the header
         read(iUnit,'(a50)')title     ! trucates or pads (with spaces) to 50 characters length
        enddo
        do ilam = 1,tlam           
         read(iUnit,20)swlambda,saw,sbw
         if (swlambda.NE.pwaves(ilam)) then
           WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: ',
     &     "wavelength for water spectrum doesn't match darwin_waves:"
           CALL PRINT_ERROR( msgBuf, myThid )
           WRITE(msgBuf,'(2A,I3,A,I4,A,I4)') 'WAVEBANDS_INIT_FIXED: ',
     &     'ilam', ilam, ': ', swlambda, ' versus ', pwaves(ilam)
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
         endif
         aw(ilam) = saw
         bw(ilam) = sbw
        enddo
        close(iUnit)
20      format(i5,f15.4,f10.4)
      else
        WRITE(msgBuf,'(A)')
     &      'WAVEBANDS_INIT_FIXED: need to specify water absorption'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
      endif
      

c phyto data files
c ANNA phyto input data files must have a column for absorption by PS pigs
c ANNA easiest way to 'turn off' PS for growth is to put same values in both abs columns
      if (darwin_phytoabsorbFile.NE. ' '  ) THEN
        CALL MDSFINDUNIT( iUnit, myThid )
        open(iUnit,file=darwin_phytoabsorbFile,
     &                                   status='old',form='formatted')
        do i = 1,6               ! six lines of text for the header
         read(iUnit,'(a50)')title
        enddo        
        sbbp = 0. _d 0
        do nabp = 1,tnabp        
         read(iUnit,'(a50)')title   ! reads one line of text for the phytoplankton type header
         do ilam  = 1,tlam         
#ifdef DAR_NONSPECTRAL_BACKSCATTERING_RATIO
          read(iUnit,30)splambda,sap,sap_ps,sbp
#else
          read(iUnit,'(i4,3f10.0,f20.0)')splambda,sap,sap_ps,sbp,sbbp
#endif
          if (splambda.NE.pwaves(ilam)) then
           WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: ',
     &     "wavelength for phyto spectrum doesn't match darwin_waves:"
           CALL PRINT_ERROR( msgBuf, myThid )
           WRITE(msgBuf,'(2A,I3,A,I4,A,I4)') 'WAVEBANDS_INIT_FIXED: ',
     &     'ilam', ilam, ': ', splambda, ' versus ', pwaves(ilam)
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
          endif
          ap(nabp,ilam) = sap
          ap_ps(nabp,ilam) = sap_ps
          bp(nabp,ilam) = sbp
          bbp(nabp,ilam) = sbbp
         enddo
        enddo
        close(iUnit)
30      format(i4,3f10.4)
      else
        WRITE(msgBuf,'(A)')
     &      'WAVEBANDS_INIT_FIXED: need to specify phyto absorption'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
      endif


#ifndef OASIM
c QQ NEED IN HERE ifndef OASIM
c surface spectrum for initial use
      if (darwin_surfacespecFile .NE. ' '  ) THEN
       CALL MDSFINDUNIT( iUnit, myThid )
       open(iUnit,file=darwin_surfacespecFile,
     &                             status='old',form='formatted')
       do i = 1,3             ! three lines of text for the header
        read(iUnit,'(a50)')title
       enddo
       do ilam = 1,tlam
        read(iUnit,40)ssflambda,ssf
        if (ssflambda.NE.pwaves(ilam)) then
           WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: ',
     &     "wavelength for surface spectrum doesn't match darwin_waves:"
           CALL PRINT_ERROR( msgBuf, myThid )
           WRITE(msgBuf,'(2A,I3,A,I4,A,I4)') 'WAVEBANDS_INIT_FIXED: ',
     &     'ilam', ilam, ': ', ssflambda, ' versus ', pwaves(ilam)
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
        endif
        sf(ilam) = ssf
       enddo
       close(iUnit)
40     format(i5,f15.6)  
      else
        WRITE(msgBuf,'(A)')
     &      'WAVEBANDS_INIT_FIXED: need surface spectrum'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
      endif
#endif /* not OASIM */


c absorption by cdom
#ifndef DAR_CALC_ACDOM
c if no file given then CDOM is zero
      if (darwin_acdomFile.NE. ' '  ) THEN
        CALL MDSFINDUNIT( iUnit, myThid )
        open(iUnit,file=darwin_acdomFile,
     &                        status='old',form='formatted')
        do i = 1,6            ! six lines of text for the header
         read(iUnit,'(a50)')title
        enddo
        do i = 1,tlam
         read(iUnit,50)sacdom
         acdom(i) = sacdom
        enddo
        close(iUnit)
50      format(f10.4)
      else
        WRITE(msgBuf,'(A)')
     &      'WAVEBANDS_INIT_FIXED: no aCDOM'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, 1 )

        do i = 1,tlam
         acdom(i) = 0. _d 0
        enddo
      endif
#else /* DAR_CALC_ACDOM */
c for 3-D or for direct comparison to RADTRANS would need the same formulation for CDOM as in radtrans.
c   CDOM absorption exponent
      nlaCDOM = 0
      do ilam = 1,tlam
       if (pwaves(ilam) .eq. darwin_lambda_aCDOM) nlaCDOM = ilam
       rlam = float(pwaves(ilam))
       excdom(ilam) = exp(-darwin_Sdom*(rlam-darwin_lambda_aCDOM))
      enddo
      if (nlaCDOM.eq.0) then
        WRITE(msgBuf,'(A,I3,A)')
     &      'WAVEBANDS_INIT_FIXED: no waveband found at ',
     &      darwin_lambda_aCDOM, ' nm (needed for DAR_CALC_ACDOM).'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
      endif

      WRITE(msgBuf,'(A,1P1E20.12)')
     &      'WAVEBANDS_INIT_FIXED: darwin_aCDOM_fac = ',darwin_aCDOM_fac
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      WRITE(msgBuf,'(A,1P1E20.12)')
     &      'WAVEBANDS_INIT_FIXED: darwin_Sdom = ', darwin_Sdom
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      WRITE(msgBuf,'(A,I3,A,I4)')
     &      'WAVEBANDS_INIT_FIXED: nlaCDOM = ', nlaCDOM, ', lambda = ',
     &      pwaves(nlaCDOM)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
#endif /* DAR_CALC_ACDOM */

#ifdef DAR_DIAG_ACDOM
c     find waveband index for acdom diagnostic
      if (darwin_diag_acdom_ilam.GE.100) then
       do ilam = 1,tlam
        if (pwaves(ilam) .eq. darwin_diag_acdom_ilam) then
          darwin_diag_acdom_ilam = ilam
          goto 60
        endif
       enddo
       WRITE(msgBuf,'(2A,I3,A)') 'WAVEBANDS_INIT_FIXED: ',
     & 'darwin_diag_acdom_ilam =',darwin_diag_acdom_ilam,
     & ' not found in darwin_waves'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
60     continue
      endif

      WRITE(msgBuf,'(A,I3,A,I4)')
     &      'WAVEBANDS_INIT_FIXED: aCDOM diag ilam = ',
     &      darwin_diag_acdom_ilam, ', lambda = ',
     &      pwaves(darwin_diag_acdom_ilam)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
#endif

      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED: wavebands:'
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: ',
     &   ' idx       low   rep      high    width'
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      do i=1,tlam
        WRITE(msgBuf,'(A,I4,F10.3,I6,F10.3,F9.3)')
     &  'WAVEBANDS_INIT_FIXED: ', i,
     &  darwin_wavebands(i),pwaves(i),darwin_wavebands(i+1),wb_width(i)
        CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      enddo
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )

#ifdef DAR_RADTRANS
c  absorption and scattering by particles
      if (darwin_particleabsorbFile .NE. ' '  ) THEN
        CALL MDSFINDUNIT( iUnit, myThid )
        open(iUnit,file=darwin_particleabsorbFile,
     &                           status='old',form='formatted')  
        do i = 1,6                 ! six lines of text for the header
         read(iUnit,'(a50)')title     ! trucates or pads (with spaces) to 50 characters length  
        enddo
        do ilam = 1,tlam           
         read(iUnit,'(I4,3F15.0)')splambda,sap,sbp,sbbp  
         if (splambda.NE.pwaves(ilam)) then
           WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: ',
     &    "wavelength for particle spectrum doesn't match darwin_waves:"
           CALL PRINT_ERROR( msgBuf, myThid )
           WRITE(msgBuf,'(2A,I3,A,I4,A,I4)') 'WAVEBANDS_INIT_FIXED: ',
     &     'ilam', ilam, ': ', splambda, ' versus ', pwaves(ilam)
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R WAVEBANDS_INIT_FIXED'
         endif
         apart(ilam) = sap
         bpart(ilam) = sbp
         bbpart(ilam) = sbbp
         apart_P(ilam) = sap/darwin_part_size_P
         bpart_P(ilam) = sbp/darwin_part_size_P
         bbpart_P(ilam) = sbbp/darwin_part_size_P
        enddo
        close(iUnit)
      else
        do ilam = 1,tlam           
         apart(ilam) = 0. _d 0
         bpart(ilam) = 0. _d 0
         bbpart(ilam) = 0. _d 0
         apart_P(ilam) = 0. _d 0
         bpart_P(ilam) = 0. _d 0
         bbpart_P(ilam) = 0. _d 0
        enddo
      endif

c     print a summary
#ifndef OASIM
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED: surface spectrum:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      WRITE(msgBuf,'(A,A)') 'WAVEBANDS_INIT_FIXED: ',
     &     ' lam sf'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      do ilam = 1,tlam           
        WRITE(msgBuf,'(A,I4,F15.6)') 'WAVEBANDS_INIT_FIXED: ',
     &     pwaves(ilam), sf(ilam)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, 1 )
      enddo
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
#endif
c
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED: water spectra:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      WRITE(msgBuf,'(A,A)') 'WAVEBANDS_INIT_FIXED: ',
     &     ' lam aw             bw'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      do ilam = 1,tlam           
        WRITE(msgBuf,'(A,I4,F15.4,F10.4)') 'WAVEBANDS_INIT_FIXED: ',
     &     pwaves(ilam), aw(ilam), bw(ilam)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, 1 )
      enddo
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
c
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED: phyto spectra:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      do nabp = 1,tnabp        
        WRITE(msgBuf,'(A,I4)') 'WAVEBANDS_INIT_FIXED: type ',nabp
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
        WRITE(msgBuf,'(A,A)') 'WAVEBANDS_INIT_FIXED: ',
     &     ' lam ap        ap_ps     bp        bbp'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
        do ilam = 1,tlam           
          WRITE(msgBuf,'(A,I4,3F10.4,F20.9)') 'WAVEBANDS_INIT_FIXED: ',
     &     pwaves(ilam), ap(nabp,ilam), ap_ps(nabp,ilam),
     &     bp(nabp,ilam), bbp(nabp,ilam)
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, 1 )
        enddo
        WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED:'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      enddo
c
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED: particulate spectra:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      WRITE(msgBuf,'(A,A)') 'WAVEBANDS_INIT_FIXED: ',
     &     ' lam apart          bpart          bbpart'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      do ilam = 1,tlam           
        WRITE(msgBuf,'(A,I4,1P3G15.6)')'WAVEBANDS_INIT_FIXED: ',
     &     pwaves(ilam), apart(ilam), bpart(ilam), bbpart(ilam)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, 1 )
      enddo
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
c
      WRITE(msgBuf,'(2A)') 'WAVEBANDS_INIT_FIXED: particulate spectra ',
     &                    'in phosphorus units:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      WRITE(msgBuf,'(A,A)') 'WAVEBANDS_INIT_FIXED: ',
     &     ' lam apart_P        bpart_P        bbpart_P'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      do ilam = 1,tlam           
        WRITE(msgBuf,'(A,I4,2F15.9,F15.12)') 'WAVEBANDS_INIT_FIXED: ',
     &     pwaves(ilam), apart_P(ilam), bpart_P(ilam), bbpart_P(ilam)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, 1 )
      enddo
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
c
#ifndef DAR_CALC_ACDOM
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED: CDOM spectrum:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      WRITE(msgBuf,'(A,A)') 'WAVEBANDS_INIT_FIXED: ',
     &     ' lam aCDOM'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
      do ilam = 1,tlam           
        WRITE(msgBuf,'(A,I4,F10.4)') 'WAVEBANDS_INIT_FIXED: ',
     &     pwaves(ilam), acdom(ilam)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, 1 )
      enddo
      WRITE(msgBuf,'(A)') 'WAVEBANDS_INIT_FIXED:'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, 1 )
#endif

c     constants
      pid = DACOS(-1.0D0)
      rad = 180.0D0/pid
#endif

#endif /* WAVEBANDS */

      return
      end

