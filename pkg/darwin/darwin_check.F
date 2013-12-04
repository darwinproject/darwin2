C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_check.F,v 1.2 2013/12/04 21:27:54 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

      SUBROUTINE DARWIN_CHECK( myThid )
C     /==========================================================\
C     | SUBROUTINE DARWIN_CHECK                                  |
C     | o Check dependences with other packages:                 |
C     |   Print darwin indices (iPO4,...) and ptracer names      |
C     |   to allow for manual check                              |
C     \==========================================================/
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_IO.h"

C     === Routine arguments ===
C     myThid -  Number of this instances
      INTEGER myThid

#ifdef ALLOW_DARWIN

C     === Local variables ===
C !LOCAL VARIABLES: ====================================================
C  iUnit   :: unit number for I/O
C  msgBuf  :: message buffer
      INTEGER iUnit,errIO
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      _BEGIN_MASTER(myThid)

C--   Print a summary of dependent parameter values:
      iUnit = standardMessageUnit
      WRITE(msgBuf,'(A)') '// ==================================='
      CALL PRINT_MESSAGE( msgBuf, iUnit, SQUEEZE_RIGHT , myThid )
      WRITE(msgBuf,'(A)') '// darwin dependent parameters'
      CALL PRINT_MESSAGE( msgBuf, iUnit, SQUEEZE_RIGHT , myThid )
      WRITE(msgBuf,'(A)') '// ==================================='
      CALL PRINT_MESSAGE( msgBuf, iUnit, SQUEEZE_RIGHT , myThid )
      IF ( darwin_useiceFile ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_CHECK:',
     &    ' using ice area from darwin_iceFile'
        CALL PRINT_MESSAGE( msgBuf, iUnit, SQUEEZE_RIGHT , myThid )
      ENDIF
      IF ( darwin_useareamask ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_CHECK:',
     &    ' using ice area from exf'
        CALL PRINT_MESSAGE( msgBuf, iUnit, SQUEEZE_RIGHT , myThid )
      ENDIF
      IF ( darwin_useiceMask ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_CHECK:',
     &    ' using ice area from THSICE pkg'
        CALL PRINT_MESSAGE( msgBuf, iUnit, SQUEEZE_RIGHT , myThid )
      ENDIF
      IF ( darwin_useAREA ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_CHECK:',
     &    ' using ice area from SEAICE pkg'
        CALL PRINT_MESSAGE( msgBuf, iUnit, SQUEEZE_RIGHT , myThid )
      ENDIF
      WRITE(msgBuf,'(A)') '// ==================================='
      CALL PRINT_MESSAGE( msgBuf, iUnit, SQUEEZE_RIGHT , myThid )

      _END_MASTER(myThid)

#ifdef ALLOW_QUOTA
      CALL QUOTA_CHECK( myThid )
#endif

#ifdef ALLOW_MONOD
      CALL MONOD_CHECK( myThid )
#endif

#endif /* ALLOW_DARWIN */
       
      RETURN
      END

