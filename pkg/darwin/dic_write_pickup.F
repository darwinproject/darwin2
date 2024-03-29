C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/dic_write_pickup.F,v 1.3 2013/10/09 17:14:37 stephd Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN

#ifdef ALLOW_CARBON

CBOP
C !ROUTINE: DIC_WRITE_PICKUP

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_WRITE_PICKUP( permPickup,
     I                             suff, myTime, myIter, myThid )

C !DESCRIPTION:
C     Writes DIC arrays (needed for a restart) to a pickup file

C !USES: ===============================================================
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DARWIN_FLUX.h"

C !INPUT PARAMETERS: ===================================================
C     permPickup :: write a permanent pickup
C     suff       :: suffix for pickup file (eg. ckptA or 0000000010)
C     myTime     :: Current time in simulation
C     myIter     :: Current iteration number in simulation
C     myThid     :: My Thread Id number
      LOGICAL permPickup
      CHARACTER*(*) suff
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     == Local variables ==
      CHARACTER*(MAX_LEN_FNAM) fn
      LOGICAL glf
      _RL     timList(1)
      INTEGER prec, j, nj
      INTEGER listDim, nWrFlds
      PARAMETER( listDim = 2 )
      CHARACTER*(8) wrFldList(listDim)
      CHARACTER*(MAX_LEN_MBUF) msgBuf

c     IF ( DIC_pickup_write_mdsio ) THEN
        prec = precFloat64
        WRITE(fn,'(A,A)') 'pickup_dic.',suff
        j = 0

C     Firstly, write 3-D fields as consecutive records,

C-    switch to 2-D fields:
        nj = -j*Nr

C       record number < 0 : a hack not to write meta files now:
        j = j + 1
        nj = nj-1
#ifdef pH_3D
        CALL WRITE_REC_3D_RL( fn, prec, Nr, pH, nj, myIter, myThid )
        IF (j.LE.listDim) wrFldList(j) = 'DIC_pH3d'
#else
        CALL WRITE_REC_3D_RL( fn, prec, 1, pH, nj, myIter, myThid )
        IF (j.LE.listDim) wrFldList(j) = 'DIC_pH2d'
#endif

C--------------------------
        nWrFlds = j
        IF ( nWrFlds.GT.listDim ) THEN
          WRITE(msgBuf,'(2A,I5,A)') 'DIC_WRITE_PICKUP: ',
     &     'trying to write ',nWrFlds,' fields'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A,I5,A)') 'DIC_WRITE_PICKUP: ',
     &     'field-list dimension (listDim=',listDim,') too small'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DIC_WRITE_PICKUP (list-size Pb)'
        ENDIF

#ifdef ALLOW_MDSIO
C     uses this specific S/R to write (with more informations) only meta
C     files
        j  = 1
        nj = ABS(nj)
        IF ( nWrFlds*Nr .EQ. nj ) THEN
          j  = Nr
          nj = nWrFlds
        ENDIF
        glf  = globalFiles
        timList(1) = myTime
        CALL MDS_WR_METAFILES( fn, prec, glf, .FALSE.,
     &                         0, 0, j, ' ',
     &                         nWrFlds, wrFldList,
     &                         1, timList, oneRL,
     &                         nj, myIter, myThid )
#endif /* ALLOW_MDSIO */
C--------------------------

c     ENDIF /* DIC_pickup_write_mdsio */


      RETURN
      END
#endif  /*ALLOW_CARBON*/

#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================
