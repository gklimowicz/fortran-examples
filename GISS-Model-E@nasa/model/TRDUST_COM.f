#include "rundeck_opts.h"
      module trdust_mod
!@sum  tracers_dust dust/mineral tracer parameter and variable declarations
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen
!@ver 3.0

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      use constant, only : By6
      use resolution, only : Im,Jm,Lm
      use tracer_com, only: ntm_dust, ntm_clay
      use OldTracer_mod, only: trName, MAX_LEN_NAME

      IMPLICIT NONE

      integer :: i

!@param By8 0.25d0/2d0
      real( kind=8 ), parameter :: By8 = 0.25D0/2D0
!@param By4 1D0/4D0
      real( kind=8 ), parameter :: By4 = 1D0/4D0

!@param dust_names  names of soil dust aerosol tracers
      character(len=MAX_LEN_NAME),allocatable,dimension(:) ::
     &     dust_names(:)

!@param dAridSoils  median particle diameter by volume of arid soils
      real( kind=8 ), parameter :: dAridSoils = 3.4d0 ! Kok, PNAS (2011)
!@param sigmaAridSoils  particle diameter standard deviation by volume of
!@+       arid soils
      real( kind=8 ), parameter :: sigmaAridSoils = 3.0d0 ! Kok, PNAS (2011)
!@param Cv normalization constant for emitted volume size distribution
!@+       according to Brittle Fragmentation Theory (Kok, PNAS 2011)
      real( kind=8 ), parameter :: Cv = 12.63 ! [um]
!@param lambda  crack propagation length for emitted volume size distribution
!@+       according to Brittle Fragmentation Theory (Kok, PNAS 2011)
      real( kind=8 ), parameter :: lambda = 12 ! [um]

!@param nDustBins  number of size classes for soil dust aerosol tracers
      integer, parameter :: nDustBins = 6
!@param dustBounds  particle diameter bounds of soil dust aerosol bins [1 um]
      real(kind=8), parameter, dimension( nDustBins+1 ) :: dustBounds =
     &     (/ 0.1d0, 2.0d0, 4.0d0, 8.0d0, 16.0d0, 32.0d0, 64.0d0 /)

!@param nSubClays  number of sub bins in clay size class of soil dust aerosols
      integer, parameter :: nSubClays = 4
!@param subClayBounds  particle diameter bounds of clay sub bins [1 um]
      real(kind=8), parameter, dimension( nSubClays+1 ) :: subClayBounds
     &     = (/ 0.1d0, 0.2d0, 0.5d0, 1.d0, 2.d0 /)

!@param nDustBinsFull  number of soil dust bins for radiation
      integer, parameter :: nDustBinsFull = nSubClays+nDustBins-1
!@param dustBoundsRadia  particle diameter bounds of soil dust aerosol bins
!@+       for radiation [1 um]
      real( kind=8 ), parameter, dimension( nDustBinsFull + 1 ) ::
     &     dustBoundsRadia = (/ subClayBounds, dustBounds( 3:nDustBins
     &     +1 ) /)

!@var subClayWeights  weights for masses in the sub bins of the clay size
!@+     class for each soil dust tracer
      real( kind=8 ), dimension( nSubClays, ntm_clay ) :: subClayWeights
     &     = reshape ( (/ ( 0.009d0, 0.081d0, 0.234d0, 0.676d0, i = 1
     &     ,ntm_clay ) /), (/ nSubClays, ntm_clay /) )

c**** rundeck parameter to switch between different emission schemes
c****
!@dbparam imDust: 0: scheme using PDF of wind speed (default)
!@+               1: prescribed AeroCom emissions
!@+               2: legacy emission scheme using third power of wind speeds
!@+                  (only works with 72x46 horizontal resolution)
!@+               3: use size distribution of AeroCom for model calculated
!@+                  dust/mineral emission flux (PDF scheme)
!@+               4: use the size distribution of the mineral version of the
!@+                  dust aerosols for the emitted dust flux (PDF scheme)
!@+               5: same as 4, but no dust emission where no AeroCom sources
      integer :: imDust=0

c**** legacy emission code (Tegen, I. and R. Miller, JGR (1998))
c**** declarations for emission scheme using third power of wind speed
c****
!@param CWiCub uplift factor [kg*s**2/m**5] for all size classes
      real( kind=8 ), parameter :: CWiCub = 52.D-9
!@param FClWiCub fraction [1] of uplifted clay for scheme using cubes of
!@+     wind speed
!@param FSiWiCub fractions [1] of uplifted silt for scheme using cubes of
!@+     wind speed
      real( kind=8 ) :: FClWiCub=By6, FSiWiCub=By8
!@var hbaij accumulated precipitation - evaporation balance  [kg/m^2]
!@var ricntd no. of hours with negative precipitation - evaporation balance [1]
      real( kind=8 ), allocatable, dimension(:,:) :: hbaij, ricntd
!@var dryhr  number of hours with evaporation-precipitation greater Zero
!@+          to allow dust emission
!@var frclay fraction of clay
!@var frsilt fraction of silt
!@var vtrsh  threshold wind speed above which dust emis. is allowed [m/s]
      real( kind=8 ), allocatable, dimension(:,:) :: dryhr, frclay,
     &     frsilt,vtrsh

c**** for legacy wet deposition
c**** declaration for simple wet deposition scheme
!@var prelay distributed array with some prec info needed for simple wet
!@+          deposition scheme
      real( kind=8 ), allocatable, dimension(:,:,:) :: prelay

c**** current default emission scheme (Cakmur, R. et al. (2004))
c**** declarations for emission scheme using probability density function of
c**** wind speed

!@param CWiPdf uplift factor [kg*s**2/m**5] for all size classes of soil dust
      real( kind=8 ), parameter :: CWiPdf = 12.068996D-9
!@dparam scaleDustEmission [1] scales emitted dust mass using the same factor
!@+     for all size bins
      real( kind=8 ) :: scaleDustEmission = 1.d0
!@dparam fracClayPDFscheme emission parameter [1] for uplifted clay
!@dparam fracSiltPDFscheme emission parameter [1] for uplifted silt
      real( kind=8 ) :: fracClayPDFscheme = 1.d0
      real( kind=8 ) :: fracSiltPDFscheme = 1.d0
!@dparam vegetationERS switches on/off reading in ERS vegetation proxy data
      integer :: vegetationERS = 1
!@dparam prefDustSources switches on/off reading in preferred dust sources
      integer :: prefDustSources = 1
!@var ers_data field of ERS data
      real( kind=8 ), allocatable, dimension(:,:,:) :: ers_data
!@var dustSourceFunction distribution of preferred dust sources
      real(kind=8),allocatable,dimension(:,:) :: dustSourceFunction

      ! note that Lim,Ljm happen to match kim,kjm in wspdf_mod
      integer, parameter :: Lim = 234, Ljm = 234, Lkm = 22
!@var x1,x2,x3 indices of lock up table for emission
      real( kind=8 ) :: x1(Lim), x2(Ljm), x3(Lkm)
!@var table array of lock up table for emission local to each grid box
      real( kind=8 ) :: table( Lim, Ljm, Lkm )

      !REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: table
      !ALLOCATE(table(Lim,Ljm,Lkm),STAT=ier)

!@var wsubtke_com distributed array of subscale turbulent term
!@var wsubwd_com distributed array of subscale dry convective term
!@var wsubwm_com distributed array of subscale moist convective term
      real( kind=8 ), allocatable, dimension(:,:) :: wsubtke_com
     &     ,wsubwd_com,wsubwm_com

c****
c**** declarations for prescribed daily dust emissions
c****
!@param nAerocomDust Number of AEROCOM dust classes
      integer, parameter :: nAerocomDust = 4
!@var d_dust prescribed daily dust emissions [kg] (e.g. AEROCOM)
      real( kind=8 ), allocatable, dimension(:,:,:,:) :: d_dust

c**** additional declarations for dust tracers with mineralogical composition
!@param nMinerals  number of different minerals
      integer, parameter :: nMinerals = 8
!@param mineralNames  names of Minerals
      character( len=4 ), parameter, dimension( nMinerals ) ::
     &     mineralNames = (/'Illi', 'Kaol', 'Smec', 'Calc', 'Quar',
     &     'Feld', 'Feox', 'Gyps' /)
!@param dustBinNamesIn dust bin names in mineral fraction input file
      character(len=5), parameter, dimension( nDustBins ) ::
     &     dustBinNamesIn = (/ 'Clay ', 'Silt1', 'Silt2', 'Silt3',
     &     'Silt4', 'Silt5' /)
!@var mineralFractions distribution of mineral fractions in soils for each
!+    mineralogical soil dust tracer
      real( kind=8 ), allocatable, dimension(:,:,:) :: mineralFractions

!@dbparam calcEffectiveRadius  flag whether to calculate or prescribe
!@+         effective radius of minerals from particle size distribution for
!@+         radiation calculations (0: prescribed=default; 1:calculated)
      integer :: calcEffectiveRadius = 0

!@param radiusClay  dry radius of clay dust/minerals [um]
      real( kind=8 ), parameter, dimension( nSubClays ) :: radiusClay =
     &     (/ 0.132d0, 0.23d0, 0.416d0, 0.766d0 /)
!@param radiusSil1  dry radius of silt1 dust/minerals [um]
      real( kind=8 ), parameter :: radiusSil1 = 1.386d0
!@param radiusSil2  dry radius of silt2 dust/minerals [um]
      real( kind=8 ), parameter :: radiusSil2 = 2.773d0
!@param radiusSil3  dry radius of silt3 dust/minerals [um]
      real( kind=8 ), parameter :: radiusSil3 = 5.545d0
!@param radiusSil4  dry radius of silt4 dust/minerals [um]
      real( kind=8 ), parameter :: radiusSil4 = 11.090d0
!@param radiusSil5  dry radius of silt5 dust/minerals [um]
      real( kind=8 ), parameter :: radiusSil5 = 22.0d0
!@param radiusMinerals  dry radius of all dust/minerals [um]
      real( kind=8 ), parameter, dimension( nDustBinsFull ) ::
     &     radiusMinerals = (/ radiusClay, radiusSil1, radiusSil2,
     &     radiusSil3, radiusSil4, radiusSil5 /)

!@param dryEffRadMinerals  dry effective radius of minerals for radiation [um]
      real( kind=8 ), parameter, dimension( nDustBinsFull ) ::
     &     dryEffRadMinerals = (/ radiusMinerals /)

#ifdef TRACERS_MINERALS
!@param densityIllite  particle density of Illite [kg/m^3]
!@+            (measured; http://www.mindat.org/min-2011.html)
!      real(kind=8), parameter :: densityIllite = 2.795d3
      real( kind=8 ), parameter :: densityIllite = 2.5725d3 ! average Illi-Smec
!@param densityKaolinite  particle density of Kaolinite [kg/m^3]
!@+            (calculated; http://www.mindat.org/min-2011.html)
      real( kind=8 ), parameter :: densityKaolinite = 2.63d3
!@param densitySmectite  particle density of Smectite [kg/m^3]
!@+            (for Montmorillonite)
!      real(kind=8), parameter :: densitySmectite = 2.35d3
      real( kind=8 ), parameter :: densitySmectite = 2.5725d3! average Illi-Smec
!@param densityCalcite  particle density of Calcite [kg/m^3]
!@+            (measured; http://www.mindat.org/min-859.html)
      real( kind=8 ), parameter :: densityCalcite = 2.71d3
!@param densityQuartz  particle density of Quartz [kg/m^3]
!@+            (measured; http://www.mindat.org/min-3337.html)
      real( kind=8 ), parameter :: densityQuartz = 2.655d3
!@param densityFeldspar  particle density of Feldspar [kg/^3]
!@+            (average Plagioclase)
      real( kind=8 ), parameter :: densityFeldspar = 2.68d3 
!@param densityHematite  particle density of average of Hematite and
!@+       Goethite [kg/m^3]
!@+            (Hematite measured: 5.26d3; http://www.mindat.org/min-1856.html)
!@+            (Goethite measured; 4.28d3; http://www.mindat.org/min-1719.html)
      real( kind=8 ), parameter :: densityHematite = 4.77d3
!@param densityGypsum  particle density of Gypsum [kg/m^3]
!@+            (measured; http://www.mindat.org/min-1784.html)
      real( kind=8 ), parameter :: densityGypsum = 2.312d3

!@dbparam soilDustEmissionDistr  chose type of dust distribution at emission
!@+         0: Kandler 2009 low/medium dust (default); 1: Kandler 2009 high dust
!@+         only used for weighting the sub clay classes and calculating
!@+         effective radii, which are parameters in the radiation code
      integer :: soilDustEmissionDistr = 0
!@dbparam brittleFactor  one-dimension array of size nMinerals with factors to
!@+         weight fractionation of minerals during emission
      real( kind=8 ), dimension( nMinerals ) :: brittleFactor = 1.d0
!@dbparam frIronOxideInAggregate  fraction of Iron(hydr)oxide mineral in other
!@+         mineral/Iron oxide aggregate
      real( kind=8 ) :: frIronOxideInAggregate = 0.05d0 ! arbitrary assumption
!@dbparam noAggregateByTotalFeox  fraction of not aggregated Iron(hydr)oxide
!@+         minerals to all Iron(hydr)oxide in minerals (not aggregated +
!@+         aggregated)
      real( kind=8 ) :: noAggregateByTotalFeox = 0.1d0 ! arbitrary assumption
!@dbparam calcMineralAggrProb  1: calculate mineral aggregation probabilities
!@+         from mineral fractions in soil; 0: use prescribed
!@+         noAggregateByTotalFeox
      integer :: calcMineralAggrProb = 1

!@var mineralIndex  index to map parameters of minerals to mineralogical tracers
      integer, dimension( ntm_dust ) :: mineralIndex
#endif

c**** Parameters for dust/mineral tracer specific diagnostics
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
!@param nDustEmij index of dust emission in ijts_isrc
      integer, parameter :: nDustEmij = 1
!@param nDustEm2ij index of dust emission according to cubic scheme
!@+                in ijts_isrc
      integer, parameter :: nDustEm2ij = 2
!@param nDustTurbij index of dust dry turbulent deposition in ijts_isrc
      integer, parameter :: nDustTurbij = 3  ! not used?
!@param nDustEv1ij index of number of dust events below threshold wind
!@+                in ijts_spec
!@param nDustEv2ij index of number of dust events above threshold wind
!@+                in ijts_spec
!@param nDustWthij index of threshold velocity in ijts_spec
      integer, parameter :: nDustEv1ij = 1, nDustEv2ij = 2, nDustWthij =
     &     3
#endif
!@dbparam to_conc_soildust: For printout of 3D soil dust aerosol concentration
!@+   in kg/m3
!@+   to_conc_soildust = 0: printout is as defined by to_volume_MixRat
!@+   to_conc_soildust = 1: printout is in kg/m3
      integer :: to_conc_soildust = 0

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
!@param nDustEmjl index of dust emission in jls_source
      integer, parameter :: nDustEmjl = 1
!@param nDustEm2jl index of dust emission according to cubic scheme
!@+                in jls_source
      integer, parameter :: nDustEm2jl = 2
!@param nDustTurbjl index of dust dry turbulent deposition in jls_source
      integer, parameter :: nDustTurbjl = 3
!@param nDustEv1jl index of number of dust events below threshold wind
!@+                in jls_spec
!@param nDustEv2jl index of number of dust events above threshold wind
!@+                in jls_spec
!@param nDustWthjl index of threshold velocity in ijts_spec
      integer, parameter :: nDustEv1jl = 1, nDustEv2jl = 2, nDustWthjl =
     &     3
#endif

c**** Variables for specific subdaily soil dust aerosol diagnostics
!@var dustDiagSubdd structured type for specific subdaily dust
!+                  aerosol diagnostics
!@var  %dustEmission   soil dust aerosol emission [kg/m^2/s]
!@var  %dustEmission2  soil dust emission calculated using cubed formula
!@+                    [kg/m^2/s] (only for diagnostic purposes)
!@var  %dustDepoTurb   turbulent soil dust aerosol deposition  [kg/m^2/s]
!@var  %dustDepoGrav   gravitational settling of soil dust aerosols [kg/m^2/s]
!@var  %dustMassInPrec wet deposition of soil dust aerosols [k]
!@var  %dustSurfMixR   surface mixing mixing ratio of soil dust aerosols [kg/kg]
!@var  %dustSurfConc   surface concentration of soil dust aerosols [kg/m^3]
!@var  %dustMass       mass of soil dust aerosol [kg]
!@var  %dustConc       dust concentration [kg/m] (later divided by grid box area)
      type dustDiagSubdd
      real( kind=8 ), allocatable, dimension(:,:,:) :: dustEmission
     &     , dustEmission2, dustDepoTurb, dustDepoGrav, dustMassInPrec
     &     , dustSurfMixR, dustSurfConc
      real( kind=8 ), allocatable, dimension(:,:,:,:) :: dustMass,
     &     dustConc
      end type dustDiagSubdd
!@var dustDiagSubdd_acc structured variable to accumulate data for
!+                      subdaily dust aerosol diagnostics
      type(dustDiagSubdd) :: dustDiagSubdd_acc

#endif /* TRACERS_DUST || TRACERS_MINERALS || TRACERS_AMP || TRACERS_TOMAS */

      end module trdust_mod

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      subroutine alloc_dust(grid)
!@sum  alloc_dust allocates dust/mineral tracer arrays
!@auth Jan Perlwitz

      use domain_decomp_atm, only : dist_grid
      use resolution, only : Lm
      use tracer_com, only : Ntm_dust
      use trdust_mod
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
      use AbstractCalendar_mod
      use model_com, only: calendar

      implicit none

      type(dist_grid),intent(in) :: grid

      integer :: i_0h, i_1h, j_1h, j_0h
      integer :: ier
      logical, save :: qfirst = .true.
      integer :: maxDaysInYear

      if (.not. qfirst) return
      qfirst = .false.

      i_0h=grid%i_strt_halo
      i_1h=grid%i_stop_halo
      j_0h=grid%j_strt_halo
      j_1h=grid%j_stop_halo

      allocate(dust_names(ntm_dust))

      maxDaysInYear = calendar%getMaxDaysInYear()

      allocate(hbaij(i_0h:i_1h,j_0h:j_1h),ricntd(i_0h:i_1h,j_0h:j_1h),
     &     dryhr(i_0h:i_1h,j_0h:j_1h),frclay(i_0h:i_1h,j_0h:j_1h),
     &     frsilt(i_0h:i_1h,j_0h:j_1h),vtrsh(i_0h:i_1h,j_0h:j_1h),
     &     dustSourceFunction(i_0h:i_1h,j_0h:j_1h),
     &     ers_data(i_0h:i_1h,j_0h:j_1h,INT_MONTHS_PER_YEAR),
     &     wsubtke_com(i_0h:i_1h,j_0h:j_1h),
     &     wsubwd_com(i_0h:i_1h,j_0h:j_1h),
     &     wsubwm_com(i_0h:i_1h,j_0h:j_1h),
     &     prelay(i_0h:i_1h,j_0h:j_1h,LM),
     &     d_dust(i_0h:i_1h,j_0h:j_1h,nDustBins,maxDaysInYear),
     &     mineralFractions( i_0h:i_1h, j_0h:j_1h, max( nDustBins,
     &     ntm_dust ) ), STAT=ier)

      dustSourceFunction( i_0h:i_1h, j_0h:j_1h ) = 1.d0

      ers_data( i_0h:i_1h, j_0h:j_1h, : ) = 1.d0

      d_dust(i_0h:i_1h,j_0h:j_1h,:,:)=0.D0

      mineralFractions( i_0h:i_1h, j_0h:j_1h, : ) = 0.d0

      allocate(dustDiagSubdd_acc%dustEmission(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustEmission2(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustDepoTurb(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustDepoGrav(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustMassInPrec(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustSurfMixR(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustSurfConc(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustMass(i_0h:i_1h,j_0h:j_1h,Lm
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustConc(i_0h:i_1h,j_0h:j_1h,Lm
     &     ,Ntm_dust))

      dustDiagSubdd_acc%dustEmission = 0.D0
      dustDiagSubdd_acc%dustEmission2 = 0.D0
      dustDiagSubdd_acc%dustDepoTurb = 0.D0
      dustDiagSubdd_acc%dustDepoGrav = 0.D0
      dustDiagSubdd_acc%dustMassInPrec = 0.D0
      dustDiagSubdd_acc%dustSurfMixR = 0.D0
      dustDiagSubdd_acc%dustSurfConc = 0.D0
      dustDiagSubdd_acc%dustMass = 0.D0
      dustDiagSubdd_acc%dustConc = 0.D0

      return
      end subroutine alloc_dust
#endif /* TRACERS_DUST || TRACERS_MINERALS || TRACERS_AMP || TRACERS_TOMAS */
