C $Header$
C $Name$

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA

      SUBROUTINE QUOTA_CHECK( myThid )
C     /==========================================================\
C     | SUBROUTINE QOUTA_CHECK                                   |
C     | o Check dependences with other packages:                 |
C     |   Print darwin indices (iPO4,...) and ptracer names      |
C     |   to allow for manual check                              |
C     \==========================================================/
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#include "QUOTA_SIZE.h"
#include "QUOTA.h"

C     === Routine arguments ===
C     myThid -  Number of this instances
      INTEGER myThid

C     === Local variables ===
C     msgBuf      - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER i,j,bi,bj,ii,io,jp,ko
      INTEGER place

      WRITE(msgBuf,'(A)') 'QUOTA_CHECK: #define ALLOW_DARWIN'
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)

      if ( nDarwin .GT. PTRACERS_num ) THEN
        WRITE(msgBuf,'(A,I5,A,I5)')
     &    'QUOTA_CHECK: ERROR: nDarwin > PTRACERS_num: ',
     &    nDarwin,' vs ',PTRACERS_num
        CALL PRINT_ERROR( msgBuf , 1)
        STOP 'ABNORMAL END: S/R QUOTA_CHECK'
      ENDIF
      WRITE(msgBuf,'(A)') 'QUOTA_CHECK: iPtr darwin PTRACERS_names'
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                         SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)')
     &  '// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
c quota tracers
      place = 0
c     inorganic nutrient tracers
      DO ii=1,iimax
        place = place + 1
        WRITE(msgBuf,'(A,I4,X,A,I4,A,A32)') 'QUOTA_CHECK: ',
     &  place,'inoganic nutrient (',ii,')  - ',PTRACERS_names(place)
          CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                                         SQUEEZE_RIGHT,myThid)
      ENDDO
      WRITE(msgBuf,'(A)')
     &  '// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
c     biomass tracers
      DO io=1,iomax
        DO jp=1,min(npmax,99)
          if (io.ne.iSili) then
            place = place + 1
            WRITE(msgBuf,'(A,I4,X,A,I1,A,I2,A,A32)') 'QUOTA_CHECK: ',
     &      place,'biomass (',io,',',jp,')    - ', PTRACERS_names(place)
              CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                                           SQUEEZE_RIGHT,myThid)
          endif
        ENDDO
      ENDDO
      WRITE(msgBuf,'(A)')
     &  '// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
c     organic matter tracers
      DO io=1,iomax-iChl
        DO ko=1,min(komax,99)
          if (io.ne.iSili.or.ko.ne.1) then
            place = place + 1
            WRITE(msgBuf,'(A,I4,X,A,I1,A,I2,A,A32)') 'QUOTA_CHECK: ',
     &      place,'orgmat(',io,',',ko,')      - ',PTRACERS_names(place)
            CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                                           SQUEEZE_RIGHT,myThid)
          endif
        ENDDO
      ENDDO


      RETURN
      END

#endif /* ALLOW_QUOTA */
#endif /* ALLOW_DARWIN */
#endif /* ALLOW_PTRACERS */

