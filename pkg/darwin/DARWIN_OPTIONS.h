C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/DARWIN_OPTIONS.h,v 1.6 2012/08/23 21:49:33 jahn Exp $
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

#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER
#undef  DAR_DIAG_GROW
#undef  DAR_DIAG_ACDOM
#undef  DAR_DIAG_ABSORP
#undef  DAR_DIAG_SCATTER
#undef  DAR_DIAG_PART_SCATTER
#undef  DAR_DIAG_IRR
#undef  DAR_DIAG_IRR_AMPS

C ======================================================================
C options for monod package
C ======================================================================
#ifdef ALLOW_MONOD

#undef  OLD_GRAZE
#undef  ALLOW_DIAZ
#undef  ALLOW_DENIT
#undef  DENIT_RELAX
#undef  OLD_NSCHEME
#undef  ALLOW_MUTANTS
#undef  ALLOW_CDOM

#undef NOTEMP
#define TEMP_VERSION 1
#define TEMP_RANGE

#undef TWO_SPECIES_SETUP
#undef NINE_SPECIES_SETUP

#undef SER_GRAZ

#undef  GEIDER
#undef  OASIM
#undef  WAVEBANDS
#undef  DAR_CALC_ACDOM
#undef  DAR_RADTRANS
#undef  DAR_RADTRANS_USE_MODEL_CALENDAR

C truncate to 2 downward decreasing modes a la Aas
C only relevant if old, non-iterative solver is used (darwin_radtrans_niter = -1)
#undef  DAR_RADTRANS_DECREASING

C use rmus for all components to convert to scalar irradiance
C (not recommended)
#undef  DAR_RADTRANS_RMUS_PAR

C define this to turn of reading of phyto backscattering spectra
C and revert to fixed backscat ratios darwin_bbphy(nabp) set in data.darwin
#undef  DAR_NONSPECTRAL_BACKSCATTERING_RATIO

C diagnostic chlorophyll
#undef  DAR_DIAG_CHL

C average PAR daily and store previous day
#undef  ALLOW_PAR_DAY

C ----------------------------------------------------------------------
C monod dependencies
C ----------------------------------------------------------------------
c if two or nine species setup we don't want specific temperature ranges
#ifdef  TWO_SPECIES_SETUP
#undef TEMP_RANGE
#endif
#ifdef  NINE_SPECIES_SETUP
#undef TEMP_RANGE
#endif

c can use either denit_relax or allow_denit but not both
#ifdef ALLOW_DENIT
#undef DENIT_RELAX
#endif

#ifdef DAR_DIAG_CHL
#define ALLOW_PAR_DAY
#endif

#endif /* ALLOW_MONOD */

C ======================================================================
C options for quota package
C ======================================================================
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
