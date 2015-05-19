C $Header: /u/gcmpack/MITgcm_contrib/darwin2/verification/quota_1d/code/DARWIN_OPTIONS.h,v 1.6 2015/05/19 14:37:07 benw Exp $
C $Name:  $

#ifndef DARWIN_OPTIONS_H
#define DARWIN_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_DARWIN

#include "CPP_OPTIONS.h"

CBOP
C    !ROUTINE: DARWIN_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
C options for darwin package and sub-packages
CEOP

C ======================================================================
C options for top-level darwin package
C ======================================================================

#define READ_PAR
#undef  USE_QSW
#define MINFE
#undef  NUT_SUPPLY
#undef  CONS_SUPP
#define PORT_RAND
#undef  OLDSEED
#undef  CALC_RATE_TOTALS

#undef  RELAX_NUTS
#undef  FLUX_NUTS

#undef  CHECK_CONS

C ======================================================================
C options for monod package
C ======================================================================
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

C ======================================================================
C options for quota package
C ======================================================================
#ifdef ALLOW_QUOTA
c PFT or pure trait-based model
#define ALLOWPFT
#define SWITCHING
#define ONEGRAZER

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
cbenw - mixotrophy
#define DIFFLIMIT
 
#define IRON_SED_SOURCE
#define IRON_SED_SOURCE_VARIABLE

#define QUOTA_DIAG_LIMIT
#undef  UNCERTAINTY

C ----------------------------------------------------------------------
C quota dependencies
C ----------------------------------------------------------------------
c
#ifdef  DYNCHL
#define GEIDER
#endif
c
#endif /* ALLOW_QUOTA */

C ======================================================================
C common dependencies
C ======================================================================
c

C ======================================================================
C overrides
C ======================================================================
C if you want to override dependent options, do it here:
C
C #define TEMP_RANGE

#endif /* ALLOW_DARWIN */
#endif /* DARWIN_OPTIONS_H */
