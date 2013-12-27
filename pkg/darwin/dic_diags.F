#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

CBOP
      SUBROUTINE DIC_DIAGS( myTime, myIter, myThid )
C     /==========================================================\
C     | SUBROUTINE DIC_BIOTIC_DIAGS                              |
C     | o Do DIC diagnostic output                               |
C     |==========================================================|
C     \==========================================================/
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_DARWIN
#include "DARWIN_SIZE.h"
#include "DARWIN_FLUX.h"
#include "DARWIN_IO.h"
#endif
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#endif
C     == Routine arguments ==
C     myTime - Current time of simulation ( s )
C     myIter - Iteration number
C     myThid -  Number of this instance of INI_FORCING
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_CARBON
#ifdef ALLOW_TIMEAVE

C     == Local variables ==
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE
      _RL  co2fluxtot
      INTEGER bi, bj, k
      CHARACTER*(MAX_LEN_MBUF) suff
      CHARACTER*(1) pf

      IF ( writeBinaryPrec .EQ. precFloat64 ) THEN
        pf(1:1) = 'D'
      ELSE
        pf(1:1) = 'R'
      ENDIF

C     Dump files and restart average computation if needed
      IF ( DIFFERENT_MULTIPLE(taveFreq,myTime,deltaTClock).
     &                               and.myIter.ne.nIter0 ) THEN

C      Normalize by integrated time
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         CALL TIMEAVE_NORMALIZE(SURave,dic_timeave, 1 ,
     &                 bi,bj,myThid)
         CALL TIMEAVE_NORMALIZE(SURCave,dic_timeave, 1 ,
     &                 bi,bj,myThid)
         CALL TIMEAVE_NORMALIZE(SUROave,dic_timeave, 1 ,
     &                 bi,bj,myThid)
#ifdef pH_3D
         CALL TIMEAVE_NORMALIZE(pCO2ave,dic_timeave, Nr,
     &                 bi,bj,myThid)
         CALL TIMEAVE_NORMALIZE(pHave,dic_timeave, Nr,
     &                 bi,bj,myThid)
#else
         CALL TIMEAVE_NORMALIZE(pCO2ave,dic_timeave, 1 ,
     &                 bi,bj,myThid)
         CALL TIMEAVE_NORMALIZE(pHave,dic_timeave, 1 ,
     &                 bi,bj,myThid)
#endif
        ENDDO
       ENDDO

#ifdef CALC_RATE_TOTALS
       call tracer_sumarea(
     &             SURCave, co2fluxtot, myThid )
c mmol/s to GtC/y
       co2fluxtot=co2fluxtot* 12. _d 0* 1. _d -3* 31536000. _d 0
     &                      *1. _d -15
       _BEGIN_MASTER( myThid )
       write(dar_fluxtot_unit,*) myIter, co2fluxtot
       _END_MASTER( myThid )
#endif


       IF ( PTRACERS_timeave_mdsio )  THEN
         WRITE(suff,'(I10.10)') myIter
         CALL WRITE_FLD_XY_RL('DIC_Surtave.',suff,SURave,
     &        myIter,myThid)
         CALL WRITE_FLD_XY_RL('DIC_SurCtave.',suff,SURCave,
     &        myIter,myThid)
         CALL WRITE_FLD_XY_RL('DIC_SurOtave.',suff,SUROave,
     &        myIter,myThid)
#ifdef pH_3D
         CALL WRITE_FLD_XYZ_RL('DIC_pHtave.',suff,pHave,
     &        myIter,myThid)
         CALL WRITE_FLD_XYZ_RL('DIC_pCO2tave.',suff,pCO2ave,
     &        myIter,myThid)
#else
         CALL WRITE_FLD_XY_RL('DIC_pHtave.',suff,pHave,
     &        myIter,myThid)
         CALL WRITE_FLD_XY_RL('DIC_pCO2tave.',suff,pCO2ave,
     &        myIter,myThid)
#endif
       ENDIF

#ifdef ALLOW_MNC
       IF ( PTRACERS_timeave_mnc ) THEN
         CALL MNC_CW_SET_UDIM('dic_tave', -1, myThid)
         CALL MNC_CW_RL_W_S('D','dic_tave',0,0,'T',myTime,myThid)
         CALL MNC_CW_SET_UDIM('dic_tave', 0, myThid)
         CALL MNC_CW_I_W_S('I','dic_tave',0,0,'iter',myIter,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_SUR_ave',SURave,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_SURC_ave',SURCave,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_SURO_ave',SUROave,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_pH_ave',pHave,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_pCO2_ave',pCO2ave,myThid)
       ENDIF
#endif

C      Reset averages to zero
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         CALL TIMEAVE_RESET(SURave,1,bi,bj,myThid)
         CALL TIMEAVE_RESET(SURCave,1,bi,bj,myThid)
         CALL TIMEAVE_RESET(SUROave,1,bi,bj,myThid)
#ifdef pH_3D
         CALL TIMEAVE_RESET(pCO2ave,Nr,bi,bj,myThid)
         CALL TIMEAVE_RESET(pHave,Nr,bi,bj,myThid)
#else
         CALL TIMEAVE_RESET(pCO2ave,1,bi,bj,myThid)
         CALL TIMEAVE_RESET(pHave,1,bi,bj,myThid)
#endif
         DIC_TimeAve(bi,bj)=0.d0
        ENDDO
       ENDDO

      ENDIF

#endif /* ALLOW_TIMEAVE */
#endif /* ALLOW_CARBON */
#endif /* ALLOW_DARWIN */
#endif /* ALLOW_PTRACERS */

      RETURN
      END

c ==================================================================
