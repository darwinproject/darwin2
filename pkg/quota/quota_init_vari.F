C $Name:  $

#include "CPP_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA

c ==========================================================
c SUBROUTINE QUOTA_INIT_VARI()
c initialize stuff for generalized quota plankton model
c adapted from NPZD2Fe - Mick Follows, Fall 2005
c               modified - Stephanie Dutkiewicz, Spring 2006
c               modified - Ben Ward, 2009/2010
c ==========================================================
c
        SUBROUTINE QUOTA_INIT_VARI(myThid)

        IMPLICIT NONE

#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "QUOTA_SIZE.h"
#include "QUOTA.h"
c#include "GCHEM.h"
#include "DARWIN_IO.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C     === Functions ===
      _RL DARWIN_RANDOM
      EXTERNAL DARWIN_RANDOM

C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf    - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(MAX_LEN_MBUF) char_str
      INTEGER char_n
      INTEGER IniUnit1, IniUnit2, IniUnit3, IniUnit4, IniUnit5

      INTEGER bi, bj, k, i, j, iPAR
      INTEGER ii,io,jp,jp2,ko
      _RL pday
c length of day (seconds)
       pday = 86400.0 _d 0
CEOP

        WRITE(msgBuf,'(A)')
     &  '// ======================================================='
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') '// Quota init variables >>> START <<<'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)')
     &  '// ======================================================='
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )

c test....................
c       write(6,*)'testing in npzd2fe_init_vari '
c test....................


c set up ecosystem coefficients
c

c initialize total number of functional groups tried
         ngroups = 0
         CALL quota_generate_phyto(MyThid)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c write out initial plankton and organic matter characteristics
         CALL MDSFINDUNIT( IniUnit1, mythid )
         open(IniUnit1,file='plankton-ini-char.dat',status='unknown')
         CALL MDSFINDUNIT( IniUnit2, mythid )
         open(IniUnit2,file='plankton_ini_char_nohead.dat',
     &        status='unknown')
c-----------------------------------------------
         char_str='         pft     bio_vol     qcarbon'
     &          //'     biosink     bioswim   mortality'
     &          //' respiration  autotrophy      pp_opt'
         char_n=108
c loop elements
         do ii=1,iimax
           WRITE(msgBuf,'(I1)'),ii
           char_str=char_str(1:char_n)//'     vmaxi_'//msgBuf(1:1)
           char_n=char_n+12
         enddo
         do ii=2,iimax ! skip carbon
           WRITE(msgBuf,'(I1)'),ii
           char_str=char_str(1:char_n)//'        kn_'//msgBuf(1:1)
           char_n=char_n+12
         enddo
c loop quotas
         do io=2,iomax-iChl ! skip carbon
           WRITE(msgBuf,'(I1)'),io
           char_str=char_str(1:char_n)//'      qmin_'//msgBuf(1:1)
     &                                //'      qmax_'//msgBuf(1:1)
           char_n=char_n+24
         enddo
         do io=1,iomax-iChl
           WRITE(msgBuf,'(I1)'),io
           if (io.eq.iCarb.or.io.eq.iNitr.or.io.eq.iPhos) then
             char_str=char_str(1:char_n)//'   excretn_'//msgBuf(1:1)
             char_n=char_n+12
           endif
           char_str=char_str(1:char_n)//' beta_mort_'//msgBuf(1:1)
     &                                //' beta_graz_'//msgBuf(1:1)
           char_n=char_n+24
         enddo
c
         char_str=char_str(1:char_n)//'    alphachl'
     &                              //'     maxgraz      k_graz'
         char_n=char_n+36
         write(IniUnit1,'(A)'),char_str(1:char_n)
c-----------------------------------------------
         do jp = 1, npmax
           write(msgBuf,120)float(pft(jp)),biovol(jp),qcarbon(jp),
     &       biosink(jp)*pday,bioswim(jp)*pday,kmort(jp)*pday,
     &       respiration(jp)*pday,autotrophy(jp),pp_opt(jp)
           char_str=msgBuf
           char_n=108
c loop elements
           do ii=1,iimax
             write(msgBuf,111)vmaxi(ii,jp)*pday
             char_str=char_str(1:char_n)//msgBuf
             char_n=char_n+12
           enddo
           do ii=2,iimax ! skip carbon
             write(msgBuf,111)kn(ii,jp)
             char_str=char_str(1:char_n)//msgBuf
             char_n=char_n+12
           enddo
c loop quotas
           do io=2,iomax-iChl ! skip carbon
             write(msgBuf,112)qmin(io,jp),qmax(io,jp)
             char_str=char_str(1:char_n)//msgBuf
             char_n=char_n+24
           enddo
           do io=1,iomax-iChl
             if (io.eq.iCarb.or.io.eq.iNitr.or.io.eq.iPhos) then
               write(msgBuf,111)kexc(io,jp)*pday
               char_str=char_str(1:char_n)//msgBuf
               char_n=char_n+12
             endif
             write(msgBuf,112)beta_mort(io,jp),beta_graz(io,jp)
             char_str=char_str(1:char_n)//msgBuf
             char_n=char_n+24
           enddo
           write(msgBuf,113)alphachl(jp),
     &                    graz(jp)*pday,kg(jp)
           char_str=char_str(1:char_n)//msgBuf
           char_n=char_n+36
c-----------------------------------------------
           write(IniUnit1,'(A)')char_str(1:char_n)
           write(IniUnit2,'(A)')char_str(1:char_n)
         enddo
c<><><><><><><><><><><><><><><><><><><><><><><><><><>
         close(IniUnit2)
         close(IniUnit1)
c-----------------------------------------------
c        write out grazing max rate and half sat matrices
         CALL MDSFINDUNIT( IniUnit3, mythid )
         open(IniUnit3,file='plankton-grazing.dat',status='unknown')
! max ingestion rates
         do jp=1,npmax
           char_n=0
           do jp2=1,npmax
             write(msgBuf,'(e9.3)')graz(jp)*pday
             char_str=char_str(1:char_n)//msgBuf(1:10)
             char_n=char_n+10
           enddo
           write(IniUnit3,'(A)')char_str(1:char_n)
         enddo
         char_n=0
         do jp2=1,npmax
           char_str=char_str(1:char_n)//'----------'
           char_n=char_n+10
         enddo
         write(IniUnit3,'(A)')char_str(1:char_n)
! 1/2-saturations
         do jp=1,npmax
           char_n=0
           do jp2=1,npmax
             write(msgBuf,'(e9.3)')kg(jp)
             char_str=char_str(1:char_n)//msgBuf(1:10)
             char_n=char_n+10
           enddo
           write(IniUnit3,'(A)')char_str(1:char_n)
         enddo
         char_n=0
         do jp2=1,npmax
           char_str=char_str(1:char_n)//'----------'
           char_n=char_n+10
         enddo
         write(IniUnit3,'(A)')char_str(1:char_n)
! predator prey-preference
         do jp=1,npmax
           char_n=0
           do jp2=1,npmax
             write(msgBuf,'(e9.3)')graz_pref(jp,jp2)
             char_str=char_str(1:char_n)//msgBuf(1:10)
             char_n=char_n+10
           enddo
           write(IniUnit3,'(A)')char_str(1:char_n)
         enddo
c<><><><><><><><><><><><><><><><><><><><><><><><><><>
         close(IniUnit3)
c-----------------------------------------------
c        write out organic matter remineralisation rates
         CALL MDSFINDUNIT( IniUnit3, mythid )
         open(IniUnit4,file='plankton-orgmat.dat',status='unknown')
! DOM remineralisation rates
         char_n=0
         do io=1,iomax-iChl
           if (io.ne.iSili) then
             write(msgBuf,'(e9.3)')remin(io,1)*pday
             char_str=char_str(1:char_n)//msgBuf(1:10)
           endif
           char_n=char_n+10
         enddo
         write(IniUnit4,'(A)')char_str(1:char_n)
! POM remineralisation rates
         char_n=0
         do io=1,iomax-iChl
           write(msgBuf,'(e9.3)')remin(io,2)*pday
           char_str=char_str(1:char_n)//msgBuf(1:10)
           char_n=char_n+10
         enddo
         write(IniUnit4,'(A)')char_str(1:char_n)
c<><><><><><><><><><><><><><><><><><><><><><><><><><>
         close(IniUnit4)
c-----------------------------------------------
 111     format(1e12.4)
 112     format(2e12.4)
 113     format(3e12.4)
 114     format(4e12.4)
 115     format(5e12.4)
 116     format(6e12.4)
 118     format(8e12.4)
 120     format(10e12.4)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         CALL LEF_ZERO( inputFe,myThid )
         CALL LEF_ZERO( sur_par,myThid )
#ifdef NUT_SUPPLY
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
             DO k=1,nR
               nut_wvel(i,j,k,bi,bj) = 0. _d 0
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
#endif

#ifdef ALLOW_PAR_DAY
      DO iPAR=1,2
       DO bj=myByLo(myThid), myByHi(myThid)
        DO bi=myBxLo(myThid), myBxHi(myThid)
         DO k=1,nR
          DO j=1-Oly,sNy+Oly
           DO i=1-Olx,sNx+Olx
            PARday(i,j,k,bi,bj,iPAR) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      IF ( .NOT. ( startTime .EQ. baseTime .AND.  nIter0 .EQ. 0
     &     .AND. pickupSuff .EQ. ' ') ) THEN
COJ should probably initialize from a file when nIter0 .EQ. 0
        CALL DARWIN_READ_PICKUP( nIter0, myThid )
      ENDIF
#endif
c
#ifdef ALLOW_TIMEAVE
c set arrays to zero if first timestep
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          CALL TIMEAVE_RESET(PARave,   Nr,  bi, bj, myThid)
          CALL TIMEAVE_RESET(PPave,    Nr,  bi, bj, myThid)
#ifdef DAR_DIAG_DIVER
          CALL TIMEAVE_RESET(Diver1ave,   Nr,  bi, bj, myThid)
          CALL TIMEAVE_RESET(Diver2ave,   Nr,  bi, bj, myThid)
          CALL TIMEAVE_RESET(Diver3ave,   Nr,  bi, bj, myThid)
          CALL TIMEAVE_RESET(Diver4ave,   Nr,  bi, bj, myThid)
#endif
c         CALL TIMEAVE_RESET(SURave,   1,  bi, bj, myThid)
          WRITE(msgbuf,'(A)')
     &      'QQ start timeave'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT , mythid)

          do k=1,Nr
            DAR_TimeAve(bi,bj,k)=0. _d 0
          enddo
         ENDDO
        ENDDO
#endif /* ALLOW_TIMEAVE */

#ifdef CHECK_CONS
coj find unused units for darwin_cons output
        CALL MDSFINDUNIT( DAR_cons_unit1, mythid )
        open(DAR_cons_unit1,file='darwin_cons_C.txt',status='unknown')
        CALL MDSFINDUNIT( DAR_cons_unit2, mythid )
        open(DAR_cons_unit2,file='darwin_cons_N.txt',status='unknown')
#endif

c test....................
c       write(6,*)'finishing darwin_init_vari '
c test....................
        WRITE(msgBuf,'(A)')
     &  '// ======================================================='
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') '// Darwin init variables >>> END <<<'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)')
     &  '// ======================================================='
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )

         RETURN
         END
#endif  /*ALLOW_QUOTA*/
#endif  /*ALLOW_DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==========================================================

