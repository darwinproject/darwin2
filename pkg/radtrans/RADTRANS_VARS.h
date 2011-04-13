C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/radtrans/RADTRANS_VARS.h,v 1.1 2011/04/13 18:56:26 jahn Exp $
C $Name:  $

#ifdef ALLOW_RADTRANS

C     nutime :: time of precomputed dpsi,eps
C     dpsi   :: Nutation in longitude (degrees)
C     eps    :: Obliquity of the Ecliptic (degrees)
C               (includes nutation in obliquity)
C
      COMMON /RADTRANS_NUTATE_VARS/dpsi,eps,nutime

      INTEGER nutime
      _RL dpsi, eps

#endif /* ALLOW_RADTRANS */

