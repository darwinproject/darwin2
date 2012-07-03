C $Header$
C $Name$

#ifndef DARWIN_OPTIONS_H
#define DARWIN_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_DARWIN

#include "CPP_OPTIONS.h"

CBOP
C    !ROUTINE: DARWIN_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
C options for darwin package
CEOP

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C common

#define READ_PAR
#undef  USE_QSW
#define MINFE
#undef  NUT_SUPPLY
#undef  CONS_SUPP
#undef  RELAX_NUTS
#undef  FLUX_NUTS
#define PORT_RAND
#undef  OLDSEED
#undef  CHECK_CONS

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#ifdef ALLOW_MONOD

#undef  OLD_GRAZE
#undef  ALLOW_DIAZ
#undef  OLD_NSCHEME
#undef  ALLOW_MUTANTS
#undef  GEIDER
#undef  WAVEBANDS

#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER

c ANNA turn inhbition on/off
c curently only set for Wavebands, LL Pro's
c#define INHIB_GEID
#define INHIB_WAVES

C diagnostic chlorophyll
#undef  DAR_DIAG_CHL

C average PAR daily and store previous day
#undef  ALLOW_PAR_DAY

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C monod dependencies
#ifdef DAR_DIAG_CHL
#define ALLOW_PAR_DAY
#endif

#endif /* ALLOW_MONOD */
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#ifdef ALLOW_QUOTA

c light & dynamic chlorophyll
#define  GEIDER
#define  DYNCHL

cbenw - turn on quota model
#define QUOTA
#define AMMON
#define NITRITE
#undef PQUOTA
#undef SQUOTA
#define FQUOTA

cbenw - switching to enable one zoo debug version
#define SWITCH3

#define QUOTA_DIAG_LIMIT
#undef  UNCERTAINTY

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c quota dependencies
c
#ifdef  DYNCHL
#define GEIDER
#endif
c
#endif /* ALLOW_QUOTA */
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c common dependencies
c
#if defined(WAVEBANDS) || defined(OASIM)
#define ALLOW_SPECTRAL
#endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#endif /* ALLOW_DARWIN */
#endif /* DARWIN_OPTIONS_H */
