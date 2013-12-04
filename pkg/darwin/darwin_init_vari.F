C $Header$
C $Name$

#include "DARWIN_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_INIT_VARI(myThid )

C !DESCRIPTION:
C calls subroutines that initialize any time dependent variables
C for any tracer experiment

C !USES: ===============================================================
      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_IO.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_DARWIN

C !FUNCTIONS: ==========================================================
      LOGICAL MASTER_CPU_THREAD
      EXTERNAL MASTER_CPU_THREAD

C !LOCAL VARIABLES: ====================================================
      INTEGER i,j,bi,bj

c set up ecosystem coefficients
c
c seed randomization
      _BEGIN_MASTER( myThid )
      CALL DARWIN_RANDOM_INIT(darwin_seed, myThid)
      _END_MASTER( myThid )
      _BARRIER

#ifdef ALLOW_QUOTA
      call QUOTA_INIT_VARI(myThid )
#endif

#ifdef ALLOW_MONOD
      call MONOD_INIT_VARI(myThid )
#endif

#ifdef CALC_RATE_TOTALS
c find unused units for darwin_rate output
      IF (MASTER_CPU_THREAD(myThid)) THEN
       CALL MDSFINDUNIT( DAR_ratetot_unit, mythid )
       open(DAR_ratetot_unit,file='darwin_ratetot.txt',status='unknown')
#ifdef ALLOW_CARBON
       CALL MDSFINDUNIT( DAR_fluxtot_unit, mythid )
       open(DAR_fluxtot_unit,file='darwin_fluxtot.txt',status='unknown')
#endif
      ENDIF
#endif

#endif /* ALLOW_DARWIN */

      RETURN
      END
