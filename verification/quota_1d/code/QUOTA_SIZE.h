C $Header$
C $Name$

c QUOTA_SIZE.h
c description: Size specification for generalized ecosystem quota model
c      Ben Ward, Summer/Autumn 2009
c
c adapted from
c DARWIN_SIZE.h 
c description: Size specification for generalized ecosystem model
c      Mick Follows, Scott Grant Fall/Winter   2005
c      Stephanie Dutkiewicz      Spring/Summer 2006
c
c iimax = no of modelled nutrient elements
c npmax = no of "functional groups" of plankton
c komax = no of "functional groups" of organic matter
c
         INTEGER iimax ! inorganic elements
         INTEGER nsize ! number of size classes
         INTEGER npmax ! total plankton types
         INTEGER iomax ! biomass elements (includes C, Si & Chl)
         INTEGER iSil  ! switch for dynamic silicate (1 or 0)
         INTEGER iChl  ! switch for dynamic chlorophyll (1 or 0)
         INTEGER komax ! organic matter types (must be 2)
         PARAMETER(iimax=5,iomax=4,komax=2) ! elements, biomasses and organics
         PARAMETER(npmax=16)                ! size and taxonomic classes
         PARAMETER(iSil=0,iChl=1)           ! flags
c
c inutrient       :: index of first nutrient in Ptracer
c ibiomass        :: index for biomass of first plankton class in Ptracer
c iorgmat         :: index for biomass of first organic matter class in Ptracer
c nDarwin         :: total number of ptracers used by Darwin Quota model
c
c remember to bring the fields in data.ptracers in the right order !
c
      INTEGER inutrient
      INTEGER ibiomass
      INTEGER iorgmat
      INTEGER nDarwin
      INTEGER iDIC
      INTEGER iNO3,iNO2,iNH4
      INTEGER iPO4
      INTEGER iFeT
      INTEGER iSi
      INTEGER iCarb,iNitr,iPhos,iSili,iIron,iChlo 
c
      PARAMETER (inutrient  = 1 )
      PARAMETER (ibiomass   = inutrient + iimax )
      PARAMETER (iorgmat    = ibiomass  + npmax * iomax )
      PARAMETER (nDarwin    = iorgmat   + 2 * (iomax-iChl) - iSil - 1 )
c
c (N.B. ptracers lists biomass of element i for all plankton, then moves on to the next element i+1)
c (N.B. ptracers lists mass of element i for all OM classes, then moves on to the next element i+1)
c
c specify order of nutrient elements
c set indices of unused elements to zero
c mandatory
      PARAMETER (iDIC = 1)
      PARAMETER (iNO3 = 2)
c optional (but preserve order)
      PARAMETER (iNO2 = 3)
      PARAMETER (iNH4 = 4)
      PARAMETER (iPO4 = 0)
      PARAMETER (iFeT = 5)
      PARAMETER (iSi  = 0)

c
c specify order of cellular element reservoirs & organic matter
c set indices of unused elements to zero
c mandatory
      PARAMETER (iCarb = 1)
      PARAMETER (iNitr = 2)
c optional (but preserve order)
      PARAMETER (iPhos = 0)
      PARAMETER (iIron = 3)
      PARAMETER (iSili = 0)
      PARAMETER (iChlo = 4)
c
c organic indices must be in the following order
c    1 C
c    2 N
c    ? P
c    ? Fe
c    ? Si
c Last Chl
c
