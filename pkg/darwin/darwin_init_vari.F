C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_init_vari.F,v 1.1 2011/04/13 18:56:24 jahn Exp $
C $Name:  $

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

C     !LOCAL VARIABLES:
C     === Local variables ===
CEOP

c set up ecosystem coefficients
c
c seed randomization
      CALL DARWIN_RANDOM_INIT(darwin_seed, myThid)


#ifdef ALLOW_QUOTA
      call QUOTA_INIT_VARI(myThid )
#endif

#ifdef ALLOW_MONOD
      call MONOD_INIT_VARI(myThid )
#endif

#endif /* ALLOW_DARWIN */

      RETURN
      END
