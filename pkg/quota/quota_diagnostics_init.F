C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/quota/quota_diagnostics_init.F,v 1.1 2014/01/17 17:59:30 jahn Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_QUOTA

CBOP
C     !ROUTINE: QUOTA_DIAGNOSTICS_INIT
C     !INTERFACE:
      SUBROUTINE QUOTA_DIAGNOSTICS_INIT( myThid )
C     !DESCRIPTION:
C     Routine to initialize Darwin diagnostics

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "PTRACERS_SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GCHEM.h"
#include "QUOTA_SIZE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid -  Number of this instance of QUOTA_DIAGNOSTICS_INIT
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf      - Informational/error meesage buffer
c     CHARACTER*(MAX_LEN_MBUF) msgBuf

      INTEGER       diagNum
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

      INTEGER      np

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Add diagnostics to the (long) list

      diagName = 'PP      '
      WRITE(diagTitle, '(A)') 'Primary Production'
      diagUnits = 'mmol P/m^3/s    '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
C============================================================================
#endif
#endif
