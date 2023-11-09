#include "rundeck_opts.h"

      MODULE obio_dim

!  Biological constituents and their units
!  P(1) = nitrate (uM or micro-moles/lt=mili-moles/m^3)
!  P(2) = ammonium (uM)
!  P(3) = silica (uM)
!  P(4) = iron (nM)
!  P(5) = diatoms (mg chl m-3)
!  P(6) = chlorophytes (mg chl m-3)
!  P(7) = cyanobacteria (mg chl m-3)
!  P(8) = coccolithophores (mg chl m-3)
!  P(9) = herbivores (mg chl m-3)
!  Detrital components
!  det(1) = N/C detritus (ugC/l)
!  det(2) = silica detritus (uM)
!  det(3) = iron detritus (nM)
!  Carbon components
!  C(1) = DOC (uM)
!  C(2) = DIC (uM)
!  P(11) = alkalinity (uM)
!  pCO2 (uatm)
!  alk  (umolC/kg)
!  Ca_det_calc   Ca in detritus calcite
!  O(1) = oxygen (mM)

      implicit none

      integer, parameter :: nnut=4
     .                     ,nchl=4
     .                     ,nzoo=1
     .                     ,ntyp=nnut+nchl+nzoo
     .                     ,ndet=3
     .                     ,ncar=2
#ifdef TRACERS_Alkalinity
#ifdef TOPAZ_params
     .                     ,nalk=2   ! 1- alk, 2-ca_det_calc
#else
     .                     ,nalk=1   ! 1- alk
#endif
#else
     .                     ,nalk=0   ! alkalinity as a function of salinity
#endif
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     .                     ,no2=1    ! oxygen
#endif
#ifdef TRACERS_abio_O2
     .                     ,nabo2=1 ! abiotic oxygen
#endif
#endif
#ifdef TRACERS_degC
     .                      ,nndegC = 1 !@PL nondegradable carbon tracer. Has to be added at end for warm carbon initialization
#endif
      integer, parameter :: ntrac = nnut+nchl+nzoo+ndet+ncar
#ifdef TRACERS_Alkalinity
     .                            + nalk
#endif
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     .                            + no2
#endif
#ifdef TRACERS_abio_O2
     .                            + nabo2
#endif
#endif
#ifdef TRACERS_degC
     .                            + nndegC !@PL
#endif


      integer, parameter :: ndimc = ntyp+ndet+ncar
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
      integer, parameter :: ndimo2 = ntyp+ndet+ncar+nalk+no2
#endif
#ifdef TRACERS_abio_O2
#ifdef TRACERS_bio_O2
      integer, parameter :: ndimabo2 = ntyp+ndet+ncar+nalk+no2+nabo2
#else
      integer, parameter :: ndimabo2 = ntyp+det+ncar+nalk+nabo2
#endif
#endif
#endif
#ifdef TRACERS_degC
#ifdef TRACERS_abio_O2
      integer,parameter :: ndimndegC = ndimabo2+nndegC
#elif (defined TRACERS_bio_O2)
      integer,parameter :: ndimndegC = ndimo2+nndegC
#else
      integer,parameter :: ndimndegC = ntyp+ndet+ncar+nalk+nndegC
#endif
#endif


      integer, parameter :: 
     .                      nh=200,   !number of depths for mean irradiance
     .                      nch=48,   !number of chl values for mean irrad
     .                      ncd=41    !number of cdom values for mean irrad

      integer, parameter :: npr=15,   !number of spectral values in par
     .                      nhn=12,   !number hourly oasim values per day
     .                      npar=npr  !same as npr


      integer, parameter :: nrg=13    !number of oceanographic basins

#ifdef TRACERS_Alkalinity
      integer, parameter :: ALK_CLIM=2
#else
      integer, parameter :: ALK_CLIM=0    !0-Alk is function of Salinity
                                          !1-Alk is from climatology (GLODAP annmean)
                                          !2-Alk is prognostic
#endif

c --- diagno_bio      output obio-model fields and diagnostic messages
      logical, public:: diagno_bio

      END MODULE obio_dim
