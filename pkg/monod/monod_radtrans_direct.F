C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/monod/monod_radtrans_direct.F,v 1.2 2012/08/24 19:45:36 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: MONOD_RADTRANS_DIRECT

C !INTERFACE: ==========================================================
      subroutine MONOD_RADTRANS_DIRECT(
     I                   H,rmud,Edsf,Essf,a_k,bt_k,bb_k,kmax,
     O                   Edbot,Esbot,Eubot,Estop,Eutop,
     O                   tirrq,tirrwq,
     O                   amp1, amp2,
     I                   myThid)

C !DESCRIPTION:
c
c  Model of irradiance in the water column.  Accounts for three
c  irradiance streams [Ackleson, Balch, Holligan, JGR, 1994],
c
c  Edbot = direct downwelling irradiance in W/m2 per waveband
c  Esbot = diffuse downwelling irradiance in W/m2 per waveband
c  Eubot = diffuse upwelling irradiance in W/m2 per waveband
c
c  all defined at the bottom of each layer.  Also computed are Estop,
c  Eutop at the top of each layer which should be very close to Esbot,
c  Eubot of the layer above.
c
c  The Ed equation is integrated exactly, Es and Eu are computed by
c  solving a set of linear equation for the amplitudes in the exact
c  solution [see, e.g., Kylling, Stamnes, Tsay, JAC, 1995].
c  The boundary condition in the deepest wet layer is
c  downward-decreasing modes only (i.e., zero irradiance at infinite
c  depth, assuming the optical properties of the last layer).
c
c  Also computed are scalar radiance and PAR at the grid cell center
c  (both in uEin/m2/s).
c
C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"             /* Nr */
#include "EEPARAMS.h"
#include "MONOD_SIZE.h"
#include "SPECTRAL_SIZE.h"    /* tlam */
#include "SPECTRAL.h"         /* WtouEin */
#include "WAVEBANDS_PARAMS.h" /* darwin_PAR_ilamLo/Hi
                                 darwin_radmodThresh
                                 darwin_rmus darwin_rmuu */

C !INPUT PARAMETERS: ===================================================
C     H     :: layer thickness (including hFacC!)
C     rmud  :: inv.cosine of direct (underwater solar) zenith angle
C     Edsf  :: direct downwelling irradiance below surface per waveband
C     Essf  :: diffuse downwelling irradiance below surface per waveband
C     a_k   :: absorption coefficient per level and waveband (1/m)
C     bt_k  :: total scattering coefficient per level and waveband (1/m)
C              = forward + back scattering coefficient
C     bb_k  :: backscattering coefficient per level and waveband (1/m)
C     kmax  :: maximum number of layers to compute
      _RL H(Nr)
      _RL rmud
      _RL Edsf(tlam), Essf(tlam)
      _RL a_k(Nr,tlam), bt_k(Nr,tlam), bb_k(Nr,tlam)
      INTEGER kmax
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C     Edbot  :: direct downwelling irradiance at bottom of layer
C     Esbot  :: diffuse downwelling irradiance at bottom of layer
C     Eubot  :: diffuse upwelling irradiance at bottom of layer
C     Estop  :: diffuse downwelling irradiance at top of layer
C     Eutop  :: diffuse upwelling irradiance at top of layer
C     tirrq  :: total scalar irradiance at cell center (uEin/m2/s)
C     tirrwq :: total scalar irradiance at cell center per waveband
C     amp1   :: amplitude of downward increasing mode
C     amp2   :: amplitude of downward decreasing mode
      _RL Edbot(tlam,Nr),Esbot(tlam,Nr),Eubot(tlam,Nr)
      _RL Estop(tlam,Nr),Eutop(tlam,Nr)
      _RL tirrq(Nr)
      _RL tirrwq(tlam,Nr)
      _RL amp1(tlam,Nr), amp2(tlam,Nr)
CEOP

#ifdef DAR_RADTRANS

C !LOCAL VARIABLES: ====================================================
      INTEGER k, nl, kbot
      _RL Edtop(tlam,Nr)
      _RL Etopwq, Ebotwq
      _RL zd
      _RL rmus,rmuu
      _RL cd,au,Bu,Cu
      _RL as,Bs,Cs,Bd,Fd
      _RL bquad,D
      _RL kappa1,kappa2,denom
      _RL c1,c2
      _RL r2(Nr),r1(Nr),x(Nr),y(Nr)
      _RL ed(Nr),e2(Nr),e1(Nr)
      _RL a3d(2*Nr), b3d(2*Nr), c3d(2*Nr), y3d(2*Nr)
      _RL rd, ru
      data rd /1.5 _d 0/   !these are taken from Ackleson, et al. 1994 (JGR)
      data ru /3.0 _d 0/

      rmus = darwin_rmus
      rmuu = darwin_rmuu

c     find deepest wet layer
      kbot = kmax
      DO WHILE (H(kbot).EQ.0 .AND. kbot.GT.0)
        kbot = kbot - 1
      ENDDO

      DO nl = 1,tlam
       DO k=1,Nr
        Edtop(nl,k) = 0.0
        Estop(nl,k) = 0.0
        Eutop(nl,k) = 0.0
        Edbot(nl,k) = 0.0
        Esbot(nl,k) = 0.0
        Eubot(nl,k) = 0.0
        amp1(nl,k) = 0.0
        amp2(nl,k) = 0.0
       ENDDO
      ENDDO
      IF (kbot.GT.0) THEN
       DO nl=1,tlam
        IF (Edsf(nl) .GE. darwin_radmodThresh .OR.
     &      Essf(nl) .GE. darwin_radmodThresh) THEN
         DO k=1,kbot
          zd = H(k)
          cd = (a_k(k,nl)+bt_k(k,nl))*rmud
          au = a_k(k,nl)*rmuu
          Bu = ru*bb_k(k,nl)*rmuu
          Cu = au+Bu
          as = a_k(k,nl)*rmus
          Bs = rd*bb_k(k,nl)*rmus
          Cs = as+Bs
          Bd = bb_k(k,nl)*rmud
          Fd = (bt_k(k,nl)-bb_k(k,nl))*rmud
          bquad = Cs + Cu
          D = 0.5*(bquad + SQRT(bquad*bquad - 4.0*Bs*Bu))
          kappa1 = D - Cs
          kappa2 = Cs - Bs*Bu/D  ! == D - Cu
          r1(k) = Bu/D
          r2(k) = Bs/D
          denom = (cd-Cs)*(cd+Cu) + Bs*Bu
          x(k) = -((cd+Cu)*Fd+Bu*Bd)/denom
          y(k) = (-Bs*Fd+(cd-Cs)*Bd)/denom
          ed(k) = EXP(-cd*zd)
          e1(k) = EXP(-kappa1*zd)
          e2(k) = EXP(-kappa2*zd)
         ENDDO

C integrate Ed equation first
         Edtop(nl,1) = Edsf(nl)
         DO k=1,kbot-1
          Edbot(nl,k) = Edtop(nl,k)*ed(k)
          Edtop(nl,k+1) = Edbot(nl,k)
         ENDDO
         Edbot(nl,kbot) = Edtop(nl,kbot)*ed(kbot)

C setup tridiagonal matrix of continuity/boundary conditions
C variables: c2(1), c1(1), c2(2), ..., c1(kbot)
C a3d,b3d,c3d: lower, main and upper diagonal
C y3d: right-hand side
C
C top b.c.: c2(1) + e1(1)*r1(1)*c1(1) = Essf - x(1)*Edsf
         a3d(1) = 0. _d 0  ! not used
         b3d(1) = 1.           ! A(1,1)*c2(1)
         c3d(1) = e1(1)*r1(1)  ! A(1,2)*c1(1)
         y3d(1) = Essf(nl) - x(1)*Edsf(nl)
C continuity at layer boundaries
         DO k=1, kbot-1
           a3d(2*k) = (1. - r2(k)*r1(k+1))*e2(k)  !   A(2k,2k-1)*c2(k)
           b3d(2*k) = r1(k) - r1(k+1)             ! + A(2k,2k  )*c1(k)
           c3d(2*k) = -1. + r2(k+1)*r1(k+1)       ! + A(2k,2k+1)*c2(k+1)
           y3d(2*k)= (x(k+1) - x(k) - r1(k+1)*(y(k+1)-y(k)))*Edbot(nl,k)
           a3d(2*k+1) = 1 - r1(k)*r2(k)                !   A(2k+1,2k  )*c1(k)
           b3d(2*k+1) = r2(k) - r2(k+1)                ! + A(2k+1,2k+1)*c2(k+1)
           c3d(2*k+1) = (-1. + r1(k+1)*r2(k))*e1(k+1)  ! + A(2k+1,2k+2)*c1(k+1)
           y3d(2*k+1)= (y(k+1) - y(k) - r2(k)*(x(k+1)-x(k)))*Edbot(nl,k)
         ENDDO
c bottom boundary condition: c1 = 0
         a3d(2*kbot) = 0. _d 0  !   A(2*kbot,2*kbot-1)*c2(kbot)
         b3d(2*kbot) = 1. _d 0  ! + A(2*kbot,2*kbot  )*c1(kbot)
         c3d(2*kbot) = 0. _d 0  ! not used
         y3d(2*kbot) = 0. _d 0  ! = 0

         CALL SOLVE_TRIDIAGONAL_PIVOT(a3d,b3d,c3d,y3d,2*kbot,myThid)

C compute irradiances
         DO k=1,kbot
          c2 = y3d(2*k-1)
          c1 = y3d(2*k)
          Estop(nl,k) = c2 + r1(k)*e1(k)*c1 + x(k)*Edtop(nl,k)
          Esbot(nl,k) = e2(k)*c2 + r1(k)*c1 + x(k)*Edbot(nl,k)
          Eutop(nl,k) = r2(k)*c2 + e1(k)*c1 + y(k)*Edtop(nl,k)
          Eubot(nl,k) = r2(k)*e2(k)*c2 + c1 + y(k)*Edbot(nl,k)
          amp1(nl,k) = c1
          amp2(nl,k) = c2
         ENDDO
         IF (kbot .LT. Nr) THEN
          Estop(nl,kbot+1) = Esbot(nl,kbot)
          Eutop(nl,kbot+1) = Eubot(nl,kbot)
         ENDIF

C       endif thresh
        ENDIF

        DO k = 1,Nr
C convert to scalar irradiance in quanta
#ifdef DAR_RADTRANS_RMUS_PAR
C        use rmus for all 3 components (?)
         Etopwq = (Edtop(nl,k)+Estop(nl,k)+Eutop(nl,k))*WtouEins(nl)
         Ebotwq = (Edbot(nl,k)+Esbot(nl,k)+Eubot(nl,k))*WtouEins(nl)
         tirrwq(nl,k) = SQRT(Etopwq*Ebotwq)*rmus
#else
C        use appropriate average cosine for each component
         Etopwq = (rmud*Edtop(nl,k)+rmus*Estop(nl,k)+rmuu*Eutop(nl,k))
     &            *WtouEins(nl)
         Ebotwq = (rmud*Edbot(nl,k)+rmus*Esbot(nl,k)+rmuu*Eubot(nl,k))
     &            *WtouEins(nl)
C        and interpolate
         tirrwq(nl,k) = SQRT(Etopwq*Ebotwq)
#endif
        ENDDO

C      enddo nl
       ENDDO
C     endif kbot.gt.0
      ENDIF

      DO k = 1,Nr
C sum PAR range
       tirrq(k) = 0.0
       DO nl = darwin_PAR_ilamLo,darwin_PAR_ilamHi
        tirrq(k) = tirrq(k) + tirrwq(nl,k)
       ENDDO
      ENDDO
c
#endif /* DAR_RADTRANS */

      return
      end

