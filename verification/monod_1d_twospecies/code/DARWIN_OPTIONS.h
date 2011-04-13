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
c options for darwin package
CEOP

#define READ_PAR
#define MINFE
#undef  NUT_SUPPLY
#undef  CONS_SUPP
#undef  OLD_GRAZE
#undef  ALLOW_DIAZ
#undef  OLD_NSCHEME
#undef  ALLOW_MUTANTS
#define PORT_RAND
#undef OLDSEED

#undef NOTEMP
#define TEMP_VERSION 1
#define TEMP_RANGE

#define TWO_SPECIES_SETUP
#undef NINE_SPECIES_SETUP

#undef  GEIDER
#undef  ALLOW_CARBON

c ANNA turn wavebands on/off
#undef WAVEBANDS 
c#define WAVEBANDS

#undef  CHECK_CONS
#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER

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
