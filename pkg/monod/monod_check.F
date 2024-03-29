C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/monod/monod_check.F,v 1.2 2013/12/04 21:17:46 jahn Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_MONOD

      SUBROUTINE MONOD_CHECK( myThid )
C     /==========================================================\
C     | SUBROUTINE MONOD_CHECK                                  |
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
#include "MONOD_SIZE.h"

C     === Routine arguments ===
C     myThid -  Number of this instances
      INTEGER myThid

C     === Local variables ===
C     msgBuf      - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER i,j,bi,bj,nz,iPtr

      WRITE(msgBuf,'(A)') 'MONOD_CHECK: #define ALLOW_MONOD'
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)

      if ( nDarwin .GT. PTRACERS_numInUse ) THEN
        WRITE(msgBuf,'(A,I5,A,I5)')
     &    'MONOD_CHECK: ERROR: nDarwin > PTRACERS_numInUse: ',
     &    nDarwin,' vs ',PTRACERS_numInUse
        CALL PRINT_ERROR( msgBuf , 1)
        STOP 'ABNORMAL END: S/R MONOD_CHECK'
      ENDIF

      WRITE(msgBuf,'(A)') 'MONOD_CHECK: iPtr darwin   PTRACERS_names'
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iPO4   ,'PO4      ', PTRACERS_names(iPO4   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)      
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iNO3   ,'NO3      ', PTRACERS_names(iNO3   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,           
     &                   SQUEEZE_RIGHT,myThid)                 
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iFeT   ,'FeT      ', PTRACERS_names(iFeT   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,           
     &                   SQUEEZE_RIGHT,myThid)                 
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iSi    ,'Si       ', PTRACERS_names(iSi    )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,           
     &                   SQUEEZE_RIGHT,myThid)                 
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iDOP   ,'DOP      ', PTRACERS_names(iDOP   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,           
     &                   SQUEEZE_RIGHT,myThid)                 
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iDON   ,'DON      ', PTRACERS_names(iDON   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,           
     &                   SQUEEZE_RIGHT,myThid)                 
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iDOFe  ,'DOFe     ', PTRACERS_names(iDOFe  )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)      

      DO nz = 1,nzmax
        iPtr = iZooP (nz)
        WRITE(msgBuf,'(A,I4,X,A,I1,A,A32)') 'MONOD_CHECK: ',
     &    iPtr, 'ZooP (',nz,') ', PTRACERS_names(iPtr)
        CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                     SQUEEZE_RIGHT,myThid)      
        iPtr = iZooN (nz)
        WRITE(msgBuf,'(A,I4,X,A,I1,A,A32)') 'MONOD_CHECK: ',
     &    iPtr, 'ZooN (',nz,') ', PTRACERS_names(iPtr)
        CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                     SQUEEZE_RIGHT,myThid)      
        iPtr = iZooFe(nz)
        WRITE(msgBuf,'(A,I4,X,A,I1,A,A32)') 'MONOD_CHECK: ',
     &    iPtr, 'ZooFe(',nz,') ', PTRACERS_names(iPtr)
        CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                     SQUEEZE_RIGHT,myThid)      
        iPtr = iZooSi(nz)
        WRITE(msgBuf,'(A,I4,X,A,I1,A,A32)') 'MONOD_CHECK: ',
     &    iPtr, 'ZooSi(',nz,') ', PTRACERS_names(iPtr)
        CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                     SQUEEZE_RIGHT,myThid)      
      ENDDO
      
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iPOP   ,'POP      ', PTRACERS_names(iPOP   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,                    
     &                   SQUEEZE_RIGHT,myThid)                          
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iPON   ,'PON      ', PTRACERS_names(iPON   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,                    
     &                   SQUEEZE_RIGHT,myThid)                          
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iPOFe  ,'POFe     ', PTRACERS_names(iPOFe  )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,                    
     &                   SQUEEZE_RIGHT,myThid)                          
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iPOSi  ,'POSi     ', PTRACERS_names(iPOSi  )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,                    
     &                   SQUEEZE_RIGHT,myThid)                          
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iNH4   ,'NH4      ', PTRACERS_names(iNH4   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,                    
     &                   SQUEEZE_RIGHT,myThid)                          
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iNO2   ,'NO2      ', PTRACERS_names(iNO2   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,                    
     &                   SQUEEZE_RIGHT,myThid)                          
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iPhy   ,'Phy(1)   ', PTRACERS_names(iPhy   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid) 
      WRITE(msgBuf,'(A,I4,X,A,I4,A,A32)') 'MONOD_CHECK: ',
     &  iPhy+npmax-1,'Phy(',npmax,')', PTRACERS_names(iPhy+npmax-1)
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid) 

#ifdef DYNAMIC_CHL
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iChl   ,'Chl(1)   ', PTRACERS_names(iChl   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A,I4,X,A,I4,A,A32)') 'MONOD_CHECK: ',
     &  iChl+npmax-1,'Chl(',npmax,')', PTRACERS_names(iChl+npmax-1)
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
#endif

#ifdef ALLOW_CARBON
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iDIC   ,'DIC      ', PTRACERS_names(iDIC   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iDOC   ,'DOC      ', PTRACERS_names(iDOC   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iPOC   ,'POC      ', PTRACERS_names(iPOC   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iPIC   ,'PIC      ', PTRACERS_names(iPIC   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iALK   ,'ALK      ', PTRACERS_names(iALK   )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A,I4,X,A,A32)') 'MONOD_CHECK: ',
     &  iO2    ,'iO2      ', PTRACERS_names(iO2    )
      CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                   SQUEEZE_RIGHT,myThid)
      DO nz = 1,nzmax
        iPtr = iZooC (nz)
        WRITE(msgBuf,'(A,I4,X,A,I1,A,A32)') 'MONOD_CHECK: ',
     &    iPtr, 'ZooC (',nz,') ', PTRACERS_names(iPtr)
        CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &                     SQUEEZE_RIGHT,myThid)
      ENDDO

#endif
       
      RETURN
      END

#endif /* ALLOW_MONOD */
#endif /* ALLOW_PTRACERS */

