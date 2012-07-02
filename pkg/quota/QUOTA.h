C $Header$
C $Name$

c QUOTA.h
c description: Key parameter arrays for generalized ecosystem quota model
c      Ben Ward, Summer/Autumn 2009
c
c adapted from
c DARWIN.h
c description: Key parameter arrays for generalized ecosystem model
c      Mick Follows, Scott Grant Fall/Winter 2005
c      Stephanie Dutkiewicz Spring/Summer 2006
c
c iimax = no of modelled nutrient elements
c npmax = no of "functional groups" of plankton
c komax = no of "functional groups" of organic matter
c
c Parameter arrays
c NOTE: list not complete!!
c
c vmaxi(iimax,npmax)            - maximum uptake rate for inorganic nutrients (and DIC)
c kn(iimax,npmax)               - 1/2 saturation concentration for inorganic nutrient uptake
c biosynth                      - cost of biosynthesis
c amminhib                      - ammonium inhibition
c qcarbon(npmax)                - Mean cell carbon content
c qmin(iomax,npmax)		- minimum quota required to sustain growth
c qmax(iomax,npmax)		- satiated quota (i.e. required for max growth rate)
c alphachl(npmax)               - chl specific initial slope of PI curve
c chl2nmax	               	- maximum cellular chl to nitrogen ratio
c pp_opt(npmax)                 - optimal prd_pry preference ratio
c pp_sig(npmax)                 - std.dev. of prd_pry preference ratio
c graz_pref(npmax,npmax)	- size ratio based lognormal grazing preference matrix
c graz(npmax)			- maximum grazing rate of predator
c kg(npmax)             	- 1/2 saturation prey biomass for grazing
c hollingtype                   - exponent of Holling grazing function
c ass_eff                       - maximum assimilation efficiency
c assim_graz(npmax,npmax) 	- assimilation efficiency of grazing
c beta_graz(iomax-iChl,npmax)   - fraction of sloppy feeding to DOM
c kexc(iimax,npmax)		- biomass specific excretion rate
c kmort(npmax)   		- plankton specific mortality rate
c beta_mort(iomax-iChl,npmax)   - fraction of nortality to DOM
c respiration(npmax)            - respiration rate
c
c biosink(npmax)		- plankton sinking rate
c bioswim(npmax)                - plankton swimming rate
c orgsink(iimax,komax)          - organic matter sinking rate
c remin(iomax-iChl,komax)       - organic matter remineralisation rate
c amm2nrite                     - ammonium to nitrite oxidation rate
c nrite2nrate                   - nitrite to nitrate oxidation rate
c k_w                           - light attenuation by water
c k_chl                         - light attenuation by chlorophyll
c
c use_NO3(npmax)                - nitrate use on/off
c use_Si(npmax)                 - silicate use on/off
c autotrophy(npmax)             - degree of autotrophy (1 = pure auto, 0 = pure hetero)
c pft(npmax)			- phytoplankton functional type
c 1=prochlorococcus, 2=synechococcus, 3=small eukaryotes,
c 4=diatoms, 5=dinoflagellates, 6=generic grazers
c
         COMMON/quota_ecoparam/
     &          vmaxi,kn,
     &          biosynth,
     &          amminhib,
     &          qcarbon,
     &          qmin,qmax,
     &          alphachl,chl2nmax,
     &          pp_opt,pp_sig,graz_pref,
     &          graz,kg,
     &          ass_eff,assim_graz,beta_graz,
     &          kexc,
     &          kmort,beta_mort,
     &          respiration,
     &          biosink,bioswim,orgsink,
     &          remin,
     &          autotrophy,
     &          pft,
     &          amm2nrite,nrite2nrate,
     &          k_w,k_chl,
     &          phytoTempCoeff, phytoTempExp1,
     &          phytoTempExp2, phytoTempOptimum, phytoDecayPower,
     &          tempcoeff1, tempcoeff2_small,
     &          tempcoeff2_big, tempcoeff3,
     &          tempmax, temprange, tempnorm, tempdecay,
     &          ngrowfac, ilight,
     &          PAR0,
#ifdef IRON_SED_SOURCE
     &          depthfesed,fesedflux,fesedflux_pcm,
#endif
#ifdef PART_SCAV
     &          scav_rat, scav_inter, scav_exp,
#endif
     &          nsource,ngroups,
     &          use_NO3,use_Si,
     &          hollingtype

         _RL vmaxi(iimax,npmax)
         _RL kn(iimax,npmax)
         _RL biosynth
         _RL amminhib
         _RL qcarbon(npmax)
         _RL qmin(iomax,npmax)
         _RL qmax(iomax,npmax)
         _RL alphachl(npmax)
         _RL chl2nmax
         _RL pp_opt(npmax),pp_sig(npmax)
         _RL graz_pref(npmax,npmax)
         _RL graz(npmax)
         _RL kg(npmax)
         _RL ass_eff
         _RL assim_graz(npmax,npmax)
         _RL beta_graz(iomax-iChl,npmax)
         _RL kexc(iimax,npmax)
         _RL kmort(npmax)
         _RL beta_mort(iomax-iChl,npmax)
         _RL respiration(npmax)
         _RL biosink(npmax)
         _RL bioswim(npmax)
         _RL orgsink(komax)
         _RL remin(iomax-iChl,komax)
c
         _RL autotrophy(npmax)
c
         _RL amm2nrite
         _RL nrite2nrate
         _RL PAR0
         _RL k_w
         _RL k_chl
         _RL phytoTempCoeff(npmax)
         _RL phytoTempExp1(npmax)
         _RL phytoTempExp2(npmax)
         _RL phytoTempOptimum(npmax)
         _RL phytoDecayPower(npmax)
         _RL tempcoeff1, tempcoeff2_big
         _RL tempcoeff2_small, tempcoeff3
         _RL tempmax, temprange, tempnorm,tempdecay
         _RL ngrowfac,ilight
#ifdef IRON_SED_SOURCE
         _RL depthfesed, fesedflux, fesedflux_pcm
#endif
#ifdef PART_SCAV
         _RL scav_rat, scav_inter, scav_exp
#endif
         INTEGER hollingtype
         INTEGER pft(npmax)
         INTEGER nsource(npmax)
         INTEGER ngroups
         INTEGER use_NO3(npmax)
         INTEGER use_Si(npmax)
c
c        Power law intercepts and exponents, with errors
c        Note, while most errors are relative (times-divide), some are absolute (plus-minus)
         COMMON/quota_allometry/biovol,biovolmin,biovfactor,
     &          a_vmaxi,     b_vmaxi,
     &          a_kn,        b_kn,
     &          a_qcarbon,   b_qcarbon,
     &          a_respir,    b_respir,
     &          a_qmin,      b_qmin,
     &          a_qmax,      b_qmax,
     &          a_alphachl,  b_alphachl,
     &          a_graz,      b_graz,
     &          a_kg,        b_kg,
     &          a_prdpry,    b_prdpry,
     &          a_kexc,      b_kexc,
     &          a_biosink,   b_biosink,
     &          a_bioswim,   b_bioswim,
     &          a_mort,      b_mort,
     &          a_beta_graz, b_beta_graz,
     &          a_beta_mort, b_beta_mort,
     &          ae_vmaxi,    be_vmaxi,
     &          ae_kn,       be_kn,
     &          ae_qcarbon,  be_qcarbon,
     &          ae_respir,   be_respir,
     &          ae_qmin,     be_qmin,
     &          ae_qmax,     be_qmax,
     &          ae_alphachl, be_alphachl,
     &          ae_graz,     be_graz,
     &          ae_kg,       be_kg,
     &          ae_prdpry,   be_prdpry,
     &          ae_kexc,     be_kexc,
     &          ae_biosink,  be_biosink,
     &          ae_bioswim,  be_bioswim,
     &          ae_mort,     be_mort,
     &          ae_beta_graz,be_beta_graz,
     &          ae_beta_mort,be_beta_mort

         _RL biovol(npmax),biovolmin,biovfactor
         _RL a_vmaxi(iimax) ,b_vmaxi(iimax)
         _RL a_kn(iimax)    ,b_kn(iimax)
         _RL a_qcarbon      ,b_qcarbon
         _RL a_respir       ,b_respir
         _RL a_qmin(iomax)  ,b_qmin(iomax)
         _RL a_qmax(iomax)  ,b_qmax(iomax)
         _RL a_alphachl     ,b_alphachl
         _RL a_graz         ,b_graz
         _RL a_kg           ,b_kg
         _RL a_prdpry       ,b_prdpry
         _RL a_kexc(iomax)  ,b_kexc(iomax)
         _RL a_biosink      ,b_biosink
         _RL a_bioswim      ,b_bioswim
         _RL a_mort         ,b_mort
         _RL a_beta_graz(iomax-iChl)
         _RL b_beta_graz(iomax-iChl)
         _RL a_beta_mort(iomax-iChl)
         _RL b_beta_mort(iomax-iChl)
         _RL ae_vmaxi(iimax),be_vmaxi(iimax)
         _RL ae_kn(iimax)   ,be_kn(iimax)
         _RL ae_qcarbon     ,be_qcarbon
         _RL ae_respir      ,be_respir
         _RL ae_qmin(iomax) ,be_qmin(iomax)
         _RL ae_qmax(iomax) ,be_qmax(iomax)
         _RL ae_alphachl    ,be_alphachl
         _RL ae_graz        ,be_graz
         _RL ae_kg          ,be_kg
         _RL ae_prdpry      ,be_prdpry
         _RL ae_kexc(iomax) ,be_kexc(iomax)
         _RL ae_biosink     ,be_biosink
         _RL ae_bioswim     ,be_bioswim
         _RL ae_mort        ,be_mort
         _RL ae_beta_graz(iomax-iChl)
         _RL be_beta_graz(iomax-iChl)
         _RL ae_beta_mort(iomax-iChl)
         _RL be_beta_mort(iomax-iChl)

c ---------------------------------------------------------
c    k0, kc         - light coefficients
         COMMON/quota_lightparam/ k0, kc, parfrac, parconv
         _RL k0, kc, parfrac, parconv

c    chlpmax        - chl:p max ratio
c    chlpmin        - chl:p min ratio
c    istar          - photoadaptation light limit
c
c ---------------------------------------------------------
c
         COMMON/quota_chlparam/ chlpmax, chlpmin, istar
         _RL chlpmax, chlpmin, istar
c
c ---------------------------------------------------------
c
         COMMON/quota_diags_param/
     &          diver_thresh0, diver_thresh1,
     &          diver_thresh2, diver_thresh3,
     &          diver_thresh4

         _RL diver_thresh0, diver_thresh1
         _RL diver_thresh2, diver_thresh3
         _RL diver_thresh4
c
c ---------------------------------------------------------
c
