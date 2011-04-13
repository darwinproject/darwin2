
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA
c ====================================================================
c SUBROUTINE QUOTA_TEMPFUNC
c ====================================================================

       SUBROUTINE QUOTA_TEMPFUNC(
     I                    Tlocal,
     O                    photo_TempFunction,
     O                    activ_TempFunction,
     I                    myThid)

         implicit none
#include "QUOTA_SIZE.h"
#include "QUOTA.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
       INTEGER myThid
       _RL photo_TempFunction
       _RL activ_TempFunction
       _RL Tlocal
c local
       _RL Temp_P,Temp_A,Tempref
       INTEGER jp
c

!       Temp_P = 0.0637 _d 0
!       Temp_A = 0.1 _d 0
       Temp_P = 0.05 _d 0
       Temp_A = 0.05 _d 0

       Tempref= 30. _d 0
c
       photo_TempFunction = exp(Temp_P*(Tlocal - Tempref))
       photo_TempFunction = max(photo_TempFunction,1. _d -10)
c
       activ_TempFunction = exp(Temp_A*(Tlocal - Tempref))
       activ_TempFunction = max(activ_TempFunction,1. _d -10)
c
       RETURN
       END
#endif  /*ALLOW_QUOTA*/
#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================
