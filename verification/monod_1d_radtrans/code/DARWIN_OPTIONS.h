C $Header: /u/gcmpack/MITgcm_contrib/darwin2/verification/monod_1d_radtrans/code/DARWIN_OPTIONS.h,v 1.2 2012/07/30 15:21:51 jahn Exp $
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
c options for darwin package
CEOP

#define READ_PAR
#define MINFE
#undef  NUT_SUPPLY
#undef  CONS_SUPP
#undef  ALLOW_DIAZ
#undef  ALLOW_MUTANTS
#define PORT_RAND

#undef NOTEMP
#define TEMP_VERSION 1
#define TEMP_RANGE

#undef TWO_SPECIES_SETUP
#undef NINE_SPECIES_SETUP

#define GEIDER
#define WAVEBANDS 
#define OASIM
#define DAR_CALC_ACDOM
#define DAR_RADTRANS
#define DAR_RADTRANS_USE_MODEL_CALENDAR
#define DAR_RADTRANS_ITERATIVE
#undef  DAR_NONSPECTRAL_BACKSCATTERING_RATIO

#undef  RELAX_NUTS
#undef  FLUX_NUTS

#undef  CHECK_CONS
#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER
#undef  DAR_DIAG_GROW
#define DAR_DIAG_ACDOM
#define DAR_DIAG_ABSORP
#define DAR_DIAG_SCATTER
#define DAR_DIAG_PART_SCATTER
#define DAR_DIAG_IRR
#define DAR_DIAG_IRR_AMPS

C diagnostic chlorophyll
#undef  DAR_DIAG_CHL

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C dependencies
c if two or nine species setup we don't want specific temperature ranges
#ifdef  TWO_SPECIES_SETUP
#undef TEMP_RANGE
#endif
#ifdef  NINE_SPECIES_SETUP
#undef TEMP_RANGE
#endif

#ifdef DAR_DIAG_CHL
#define ALLOW_PAR_DAY
#endif
 
#endif /* ALLOW_DARWIN */
#endif /* DARWIN_OPTIONS_H */
