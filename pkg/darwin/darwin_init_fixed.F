C $Header$
C $Name$

#include "DARWIN_OPTIONS.h"
#ifdef ALLOW_EXF
#include "EXF_OPTIONS.h"
#endif

C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_INIT_FIXED(myThid )

C !DESCRIPTION:
C calls subroutines that initialized fixed variables for any tracer
C experiment

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_EXF
#include "EXF_PARAM.h"
#endif
#include "DARWIN_SIZE.h"
#include "DARWIN_IO.h"
#ifdef ALLOW_SPECTRAL
#include "SPECTRAL_SIZE.h"
#include "SPECTRAL_PARAMS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_DARWIN

C     !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef ALLOW_SPECTRAL
      INTEGER ilam
      _RL planck, c, hc, oavo, hcoavo, rlamm
#endif

      _BEGIN_MASTER(myThid)

#ifdef ALLOW_SPECTRAL
c   Quanta conversion
      planck = 6.6256 _d -34   !Plancks constant J sec
      c = 2.998 _d 8                 !speed of light m/sec
      hc = 1.0/(planck*c)
      oavo = 1.0/6.023 _d 23   ! 1/Avogadros number
      hcoavo = hc*oavo
      do ilam = 1,tlam
       rlamm = darwin_waves(ilam)*1 _d -9  !lambda in m
       WtouEins(ilam) = 1 _d 6*rlamm*hcoavo      !Watts to uEin/s conversion
      enddo
#endif

C which ice area fraction to use (by preference):
C 1. read from darwin_iceFile if given
C 2. read from areamaskfile through exf if given
C 3. use iceMask from THSICE if compiled and used
C 4. use AREA from seaice if compiled and used
      darwin_useiceFile  = .FALSE.
      darwin_useareamask = .FALSE.
      darwin_useiceMask  = .FALSE.
      darwin_useAREA     = .FALSE.
      darwin_haveIce     = .FALSE.
      IF ( darwin_iceFile .NE. ' '  ) THEN
        darwin_useiceFile  = .TRUE.
        darwin_haveIce = .TRUE.
      ENDIF
#ifdef ALLOW_EXF
      IF (.NOT.darwin_haveIce.AND.useEXF.AND.areamaskFile.NE.' ') THEN
        darwin_useareamask = .TRUE.
        darwin_haveIce = .TRUE.
#ifndef EXF_SEAICE_FRACTION
        WRITE(msgBuf,'(2A)') 'S/R DARWIN_INIT_FIXED:',
     &    ' you specified an areamaskfile but EXF_SEAICE_FRACTION'
        CALL PRINT_ERROR( msgBuf , 1)
        WRITE(msgBuf,'(2A)') 'S/R DARWIN_INIT_FIXED:',
     &    ' is not defined.  You need to define it in EXF_OPTIONS.h'
        CALL PRINT_ERROR( msgBuf , 1)
        WRITE(msgBuf,'(2A)') 'S/R DARWIN_INIT_FIXED:',
     &    ' in order to use areamaskfile.'
        CALL PRINT_ERROR( msgBuf , 1)
        STOP 'DARWIN_INIT_FIXED'
#endif
      ENDIF
#endif
#ifdef ALLOW_THSICE
      IF ( .NOT.darwin_haveIce .AND. useTHSICE ) THEN
        darwin_useiceMask = .TRUE.
        darwin_haveIce = .TRUE.
      ENDIF
#endif
#ifdef ALLOW_SEAICE
      IF ( .NOT.darwin_haveIce .AND. useSEAICE ) THEN
        darwin_useAREA = .TRUE.
        darwin_haveIce = .TRUE.
      ENDIF
#endif

      _END_MASTER(myThid)

#ifdef ALLOW_QUOTA
      call QUOTA_INIT_FIXED(myThid)
#endif
#ifdef ALLOW_MONOD
      call MONOD_INIT_FIXED(myThid)
#endif

      _BARRIER

#endif /* ALLOW_DARWIN */

      RETURN
      END
