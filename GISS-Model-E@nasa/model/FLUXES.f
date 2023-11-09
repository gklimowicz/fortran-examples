#include "rundeck_opts.h"

!@var  itype  1, ocean; 2, ocean ice; 3, land ice; 4, land
      module itype_enum
      integer, parameter :: ITYPE_MIN = 1
      integer, parameter :: ITYPE_OCEAN = 1
      integer, parameter :: ITYPE_OCEANICE = 2
      integer, parameter :: ITYPE_LANDICE = 3
      integer, parameter :: ITYPE_LAND = 4
      integer, parameter :: ITYPE_MAX = 4
      end module itype_enum

      MODULE EXCHANGE_TYPES ! todo: move to another file.
      use dist_grid_mod, only : dist_grid
#ifdef CUBED_SPHERE
      !USE cs2ll_utils, only : cs2llint_type,ll2csint_type
#else
      use domain_decomp_1d, only : band_pack_type
#endif
      use vector_integer_mod, only: vector_integer=>vector
      IMPLICIT NONE

      type simple_bounds_type ! todo: move to another module
         INTEGER :: I_0,I_1, J_0,J_1  ! bounds of domain
         INTEGER :: I_0H,I_1H, J_0H,J_1H ! bounds of arrays
         LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
         INTEGER, DIMENSION(:), POINTER :: IMAXJ
      end type simple_bounds_type

      type, extends(simple_bounds_type) :: atmsrf_xchng_vars

!@var surf_name name of surface type
         character(len=16) :: surf_name

!@var itype4 which of the 4 main surface types this instance is
         integer :: itype4=0

!@var grid a pointer to the grid object whose domain bounds were
!@+   used to allocate an instance of this type
         type(dist_grid), pointer :: grid

!@var xxx_exports are memory spaces holding selected fields for
!@+   en-masse operations (averaging over surface types,
!@+   pointer management, etc.)

!----------------------------------------------------------------
!----------------------------------------------------------------

!@var srfflx_exports contains the flux fields declared below
!@+   that are calculated during SURFACE
         real*8, dimension(:,:,:), pointer ::
     &        srfflx_exports=>null()
         real*8, dimension(:,:), pointer ::
!@var E0 net energy flux at surface [J m-2]
!@var SOLAR absorbed solar radiation [J m-2]
!@var TRHEAT net LW flux accumulation [J m-2]
     &      E0,SOLAR,TRHEAT
!@var DMUA,DMVA momentum flux from atmosphere (kg/m s)
!@+   On atmospheric A grid (tracer point)
     &     ,DMUA, DMVA
!@var EVAPOR evaporation (kg/m^2) 
!@var SENSHT sensible heat flux accumulation [J m-2]
!@var LATHT latent heat flux accumulation [J m-2]
     &     ,EVAPOR,SENSHT,LATHT
     &     ,UFLUX1,VFLUX1 ! (temporary redundancy with dmua, dmva)
     &     ,DTH1  ! (temporary) first layer temp. increment
     &     ,DQ1   ! (temporary) first layer humidity increment
!@+   TODO: have sea ice code refer to these
!@var RUNO runoff [kg m-2]
!@var ERUNO energy of runoff [J m-2]
     &     ,RUNO, ERUNO


!----------------------------------------------------------------
!----------------------------------------------------------------

!@var srfstate_exports contains the fields declared below
!@+   which characterize the state of a surface component
!@+   for flux calculations and diagnostics
         real*8, dimension(:,:,:), pointer ::
     &        srfstate_exports=>null()
         real*8, dimension(:,:), pointer ::
!@var GTEMP temperature of surface [degC]
!@var GTEMP2 "ground" temperature of "second" layer [degC]
!@var GTEMPR radiative ground temperature over surface type [K]
!@var GTEMPS skin temperature over surface type [degC]
!@var SNOW,SNOWFR,SNOWDP snow mass, fraction, depth
     &      GTEMP,GTEMP2,GTEMPR,GTEMPS
     &     ,SNOW,SNOWFR,SNOWDP

!----------------------------------------------------------------
!----------------------------------------------------------------

!@var pbl_exports contains the fields declared below
!@+   which are calculated by the PBL scheme and are
!@+   needed for flux calculations and diagnostics
         real*8, dimension(:,:,:), pointer ::
     &        pbl_exports=>null()
         real*8, dimension(:,:), pointer ::
!@var WSAVG  SURFACE WIND MAGNITUDE (M/S)
     &      WSAVG
!@var USAVG,VSAVG,TSAVG,QSAVG reference-height surf. wind components, temp, humidity
     &     ,USAVG,VSAVG,TSAVG,QSAVG,RSAVG
!@var cmgs drag coefficient (dimensionless surface momentum flux)
!@var chgs Stanton number   (dimensionless surface heat flux)
!@var cqgs Dalton number    (dimensionless surface moisture flux)
     &     ,cmgs,chgs,cqgs
!@var USTAR_pbl friction velocity (sqrt of srfc mom flux) (m/s)
     &     ,ustar_pbl
     &     ,lmonin_pbl
!@var WSAVG     SURFACE WIND MAGNITUDE (M/S)
!@var TSAVG     SURFACE AIR TEMPERATURE (K)
!@var QSAVG     SURFACE AIR SPECIFIC HUMIDITY (1)
!@var RSAVG     SURFACE AIR RELATIVE HUMIDITY (1)
!@var USAVG     SURFACE U WIND
!@var VSAVG     SURFACE V WIND
!@var TAUAVG    SURFACE MOMENTUM TRANSFER (TAU)
!@var TGVAVG    virtual temperature of the ground (K)
!@var gustiwind wind gustiness (m/s)
!@var dblavg    boundary layer height (m)
!@var w2_l1     vertical component of t.k.e. at gcm layer 1
!@var ciaavg    cross-isobar angle
!@var khsavg    vertical diffusivity at "surface" height (m2/s)
!@var wspdf     mean surface wind calculated from PDF of wind speed [m/s]
     &     ,tauavg,tgvavg,qgavg
     &     ,w2_l1,gustiwind,dblavg,rhoavg
     &     ,ciaavg,khsavg,wspdf
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
!@var wsgcm magnitude of the GCM surface wind - ocean currents [m/s]
!@var wsubtke turbulent kinetic energy velocity scale [m/s]
!@var wsubwd dry convective velocity scale [m/s]
!@var wsubwm moist convective velocity scale [m/s]
     &     ,wsgcm,wsubwd,wsubtke,wsubwm
#endif

!----------------------------------------------------------------
!----------------------------------------------------------------

!@var atm_exports_phase1 contains the fields declared below
!@+   which are calculated by the atmospheric model before
!@+   the flux calcuations in SURFACE, and which do not vary
!@+   during the sub-timesteps in SURFACE
!@+   (see atm_exports_phasesrf below).
         real*8, dimension(:,:,:), pointer ::
     &        atm_exports_phase1=>null()
         real*8, dimension(:,:), pointer ::
     &      SRFP  ! SRFP actual surface pressure (hecto-Pascals)
     &     ,SRFPK ! srfp**kapa
     &     ,AM1   ! first-layer air mass (kg/m2)
     &     ,BYAM1 ! 1/AM1
     &     ,P1    ! center pressure of first layer (mb)
!@var PREC precipitation [kg m-2]
!@var EPREC energy of preciptiation [J m-2]
     &     ,PREC,EPREC
!@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
     &     ,COSZ1
!@var FLONG, FSHORT downwelling longwave, shortwave radiation at surface
     &     ,FLONG,FSHORT
     &     ,TRUP_in_rad ! LW emission by surface during rad. timestep.


!----------------------------------------------------------------
!----------------------------------------------------------------

!@var atm_exports_phasesrf contains the atmospheric fields declared below
!@+   which are available before the flux calcuations in SURFACE, but which
!@+   are updated during the sub-timesteps in SURFACE.
         real*8, dimension(:,:,:), pointer ::
     &        atm_exports_phasesrf=>null()
         real*8, dimension(:,:), pointer ::
!@var TEMP1 pot. temp. of first layer w.r.t. 1 mb (K)
!@var Q1 specific humidity of first layer
!@var U1,V1 wind components of first layer (A-grid)
     &      TEMP1
     &     ,Q1
     &     ,U1,V1
! The "layer 1" arrays TEMP1,Q1,U1,V1,P1,AM1 will be redefined
! to correspond to the full boundary layer depth.

!----------------------------------------------------------------
!----------------------------------------------------------------

#ifdef TRACERS_ON

!@var trsrfflx_exports contains the tracer flux fields declared below
!@+   that are calculated during SURFACE
         real*8, dimension(:,:,:,:), pointer ::
     &        trsrfflx_exports=>null()
         real*8, dimension(:,:,:), pointer ::
!@var TRSRFFLX interactive surface sources/sinks for tracers (kg/m2/s)
     &      TRSRFFLX
#ifdef TRACERS_WATER
!@var TREVAPOR tracer evaporation  (kg/m^2) 
!@var TRUNO tracer runoff (kg/m^2)
     &     ,TREVAPOR,TRUNO
#endif
#ifdef TRACERS_DRYDEP
!@var TRDRYDEP tracer dry deposition by type (kg/m^2) (positive down)
     &     ,TRDRYDEP
#endif


!----------------------------------------------------------------
!----------------------------------------------------------------

!@var trsrfstate_exports contains the tracer fields declared below
!@+   which characterize the state of a surface component
!@+   for flux calculations and diagnostics
         real*8, dimension(:,:,:,:), pointer ::
     &        trsrfstate_exports=>null()
         real*8, dimension(:,:,:), pointer ::
!@var GTRACER ground concentration of tracer (kg/kg)
     &        GTRACER

!----------------------------------------------------------------
!----------------------------------------------------------------

!@var trpbl_exports contains the tracer fields declared below
!@+   which are calculated by the PBL scheme and are
!@+   needed for flux calculations and diagnostics
         real*8, dimension(:,:,:,:), pointer ::
     &        trpbl_exports=>null()
         real*8, dimension(:,:,:), pointer ::
!@var travg near-surface tracer mixing ratio
!@var travg_byvol near-surface tracer concentration (kg/m3)
     &        travg,travg_byvol
#ifdef TRACERS_DRYDEP
!@var dep_vel dry deposition velocity
!@var gs_vel gravitational settling velocity
!@var drydflx dry deposition flux
     &       ,dep_vel,gs_vel,drydflx
#ifdef ACCMIP_LIKE_DIAGS
!@var stomatal_dep_vel turbulent deposition velocity via stomata(m/s)
        real*8, dimension(:,:), pointer :: stomatal_dep_vel ! just one tracer
#endif
#endif


!----------------------------------------------------------------
!----------------------------------------------------------------

!@var tratm_exports_phase1 contains the fields declared below
!@+   which are calculated by the atmospheric model before
!@+   the flux calcuations in SURFACE, and which do not
!@+   vary during the sub-timesteps in SURFACE (see
!+    tratm_exports_phasesrf below).
         real*8, dimension(:,:,:,:), pointer ::
     &        tratm_exports_phase1=>null()
         real*8, dimension(:,:,:), pointer ::
     &      TRFLUX_prescr ! prescribed, non-interactive emissions (kg/m2/s)
#ifdef TRACERS_WATER
!@var TRPREC tracers in precip (kg/m^2)
     &     ,TRPREC
#endif

!----------------------------------------------------------------
!----------------------------------------------------------------

!@var tratm_exports_phasesrf contains the atmospheric fields declared below
!@+   which are available before the flux calcuations in SURFACE, but which
!@+   are updated during the sub-timesteps in SURFACE.
         real*8, dimension(:,:,:,:), pointer ::
     &        tratm_exports_phasesrf=>null()
         real*8, dimension(:,:,:), pointer ::
     &      TRM1  ! first-layer tracer mass (kg/m2).  Todo: conc. instead?
! See note above on "layer 1" arrays.

#endif /* TRACERS_ON */

         REAL*8, DIMENSION(:,:), POINTER ::
!@var ftype fraction of the gridcell occupied by this patch = flice*fhc
     &      FTYPE
!@var fhc fraction of ICE-COVERED area occupied by this patch
     &     ,FHC
!@var LAT latitude of gridbox (radians)
     &     ,LAT
!@var WORK1,WORK2 temporary workspace
     &     ,WORK1,WORK2

!@var AIJ a pointer to modelE diag accumulation arrays
         REAL*8, DIMENSION(:,:,:), POINTER :: AIJ

#ifdef TRACERS_ON
!@var TAIJN a pointer to modelE diag accumulation arrays
         REAL*8, DIMENSION(:,:,:,:), POINTER :: TAIJN
#endif

!@var JREG lat/lon array defining regions for modelE AREG diagnostics
         INTEGER, DIMENSION(:,:), POINTER :: JREG

#ifndef CUBED_SPHERE
!@var dlatm latitudinal gridsize in minutes, for certain regrid routines
         real*8 :: dlatm
!@var sini,cosi sin(lon),cos(lon) for vector regrids
         real*8, dimension(:), pointer :: sini,cosi
#endif

#ifdef TRACERS_ON
         integer :: ntm=0
#endif


!@var ipbl flag for whether pbl properties were found at last timestep
      integer, pointer, dimension(:,:) :: ipbl

! The following fields are not fundamental to atmosphere-surface interaction,
! as they are specific to a particular PBL scheme employing semi-prognostic
! subgrid profiles.  However, it is currently convenient to declare them here
! but let them be allocated by the PBL scheme.
!@var uabl boundary layer profile for zonal wind
!@var vabl boundary layer profile for meridional wind
!@var tabl boundary layer profile for temperature
!@var qabl boundary layer profile for humidity
!@var eabl boundary layer profile for turbulent KE (calc. on sec. grid)
      real*8, pointer, dimension(:,:,:) :: uabl,vabl,tabl,qabl,eabl
#ifdef TRACERS_ON
!@var trabl boundary layer profile for tracers
      real*8, pointer, dimension(:,:,:,:) :: trabl
#endif

      end type atmsrf_xchng_vars

      type, extends(atmsrf_xchng_vars) :: atmocn_xchng_vars
         logical :: updated=.false. ! indicates that values have been updated by relevant component
         REAL*8, DIMENSION(:,:), POINTER ::
!@var FOCEAN ocean fraction
     &      FOCEAN
!@var FLOWO,EFLOWO mass, energy from rivers into ocean (kg/m^2, J/m^2)
     &     ,FLOWO, EFLOWO
!@var GMELT,EGMELT mass,energy from glacial melt into ocean (kg/m^2, J/m^2)
     &     ,GMELT, EGMELT
! While not inherently atmospheric, FLOWO/GMELT et al. are here because river
! discharge and icebergs are currently calculated on the atmospheric grid.
!@var UOSURF, VOSURF ocean surface velocity (cell center) (m/s)
!@+   components defined along true N/S and E/W directions
!@+   At the NP, U points from 90E to 90W, V from IDL to GM
!@+   At the SP, U points from 9OW to 90E, V from GM to IDL 
     &     ,UOSURF,VOSURF
!@var OGEOZA ocean surface height geopotential (m^2/s^2)
     &     ,OGEOZA
!@var MLHC ocean mixed layer heat capacity (J/m^2 C) 
     &     ,MLHC
!@var SSS sea surface salinity (ppt)
     &     ,SSS
#ifdef STANDALONE_OCEAN
!@var SSSOBS,SSSRESFAC observed salinity and restoration factor toward it
!@var RSIOBS observed sea ice fraction
     &     ,SSSOBS,SSSRESFAC
     &     ,RSIOBS
#endif
#ifdef TRACERS_WATER
!@var TRFLOWO tracer in river runoff into ocean (kg/m^2)
         REAL*8, DIMENSION(:,:,:), POINTER :: TRFLOWO
#ifdef TRACERS_OCEAN
!@var TRGMELT tracer from glacial melt into ocean (kg/m^2)
         REAL*8, DIMENSION(:,:,:), POINTER :: TRGMELT
#endif
#endif
!@var TRGASEX  tracer gas exchange (mol,CO2/m^2/s)
         REAL*8, DIMENSION(:,:,:), POINTER :: TRGASEX
         REAL*8, DIMENSION(:,:), allocatable ::
     &     DIRVIS,DIFVIS,DIRNIR,DIFNIR
C**** array of Chlorophyll data for use in ocean albedo calculation
!@var CHL Chlorophyll concentration data (mgr/m**3)
         REAL*8, DIMENSION(:,:), POINTER :: CHL
         logical :: chl_defined ! df: temporary until ocalbedo is made more general

!@var eflow_gl global integral of eflowo
         real*8 :: eflow_gl=0.
         logical :: need_eflow_gl=.false.

!@var modd5s,jm_budg,area_of_zone,conserv,nofm
!@+   permit an ocean model to accumulate modelE-style
!@+   conservation diagnostics.  see diag_com.
         integer :: modd5s,jm_budg
         real*8, dimension(:), pointer :: area_of_zone
         real*8, dimension(:,:), pointer :: consrv=>NULL()
         integer, dimension(:,:), pointer :: nofm

! Some atmosphere-declared tracer info for uses within ocean codes.
! See TRACER_COM.f
         real*8, dimension(:), pointer :: trw0
         type(vector_integer) :: gasex_index
         integer :: n_co2n=0
         real*8, dimension(:), allocatable :: vol2mass

#ifdef TRACERS_OCEAN
!@var natmtrcons,tconsrv,nofmt
!@+   permit an ocean model to accumulate modelE-style
!@+   tracer conservation diagnostics.  see diag_com.
         integer :: natmtrcons
         real*8, dimension(:,:,:), pointer :: tconsrv
         integer, dimension(:,:), pointer :: nofmt

#ifndef TRACERS_ON
         integer :: ntm=0
#endif

#endif

      end type atmocn_xchng_vars

      type, extends(atmsrf_xchng_vars) :: atmice_xchng_vars
         REAL*8, DIMENSION(:,:), POINTER ::
!@var FOCEAN ocean fraction
     &      FOCEAN
!@var E1 net energy flux at layer 1
     &     ,E1
!@var UISURF, VISURF dynamic ice surface velocity (Atm A grid) (m/s)
!@+   directions as for UOSURF/VOSURF
     &     ,UISURF,VISURF
!@var FWSIM fresh water sea ice mass (kg/m^2) (used for qflux model)
!@var HSICNV,MSICNV heat and fresh water sea ice mass convergence
!@+   after advsi (kg/m^2) (used for qflux model)
     &     ,FWSIM,MSICNV,HSICNV
!@var RSI fraction of water area covered in ice
!@var SNOWI,ZSNOWI amount, thickness of snow on sea ice (kg/m^2, m)
     &     ,RSI,SNOWI,ZSNOWI
!@var USI,VSI ice velocities (m/s)
     &     ,USI,VSI ! temporary while ice still on atm grid
     &     ,SNOAGE ! really belongs to icestate

! Some arrays and indices for diagnostic purposes.  they are
! placed in this atm-ice interface type because the diagnostics
! are administrated by the atm model.  Many of them were
! introduced to maintain identicality of diagnostics through
! code reorganizations and can be eliminated if small changes
! in diagnostics are acceptable.
         REAL*8, DIMENSION(:,:), POINTER ::
     &     RSIstart,MSIsave,SNOWsave,TICEsave,TI1save,SSI1save,SSI2save
     &     ,SIHC,SNTOSI,SITOPMLT,MSNFLOOD,HSNFLOOD, SNOWsave2,MSIsave2
     &     ,MUSI,MVSI ,HUSI,HVSI ,SUSI,SVSI
         INTEGER :: IJ_RSNW,IJ_SNOW,IJ_RSIT,IJ_ZSNOW
         INTEGER :: IJ_MLTP,IJ_TSICE
     &     ,IJ_SIHC,IJ_SNTOSI,IJ_MSNFLOOD,IJ_HSNFLOOD,IJ_SIBOTMLT
     &     ,IJ_SMFX,IJ_SIGRLT,IJ_FWIO,IJ_HTIO,IJ_STIO,IJ_SIGRFR
     &     ,IJ_SIGRCG,IJ_SSI1,IJ_SSI2,IJ_TSI,IJ_F0OI,IJ_SISNWF
     &     ,IJ_RSOI,IJ_MSI,IJ_SITOPMLT,IJ_SITF,IJ_SIMASS,IJ_SIVOL
!@var IJ_[MHS][UV]SI indices for sea ice mass/heat/salt transport diags
     &     ,IJ_MUSI,IJ_MVSI,IJ_HUSI,IJ_HVSI,IJ_SUSI,IJ_SVSI,IJ_dHSI_Dyn
         INTEGER :: J_IMELT,J_HMELT,J_SMELT
     &     ,j_implm,j_implh
     &     ,j_rsnow,j_rsi,j_ace1,j_ace2,j_snow
         INTEGER :: itoice,itlkice,itocean,itlake
#ifdef TRACERS_WATER
         LOGICAL, DIMENSION(:), ALLOCATABLE :: DO_ACCUM
         REAL*8, DIMENSION(:,:,:), POINTER :: TUSI,TVSI,TRSIsum
         integer :: tij_icocflx, tij_seaice, tij_tusi, tij_tvsi
#endif

      end type atmice_xchng_vars

      ! -----------------------------------------------------
      ! For coupling of atmosphere with glacial ice.
      type, extends(atmsrf_xchng_vars) :: atmgla_xchng_vars


!@var xxx_exports_gla memory spaces holding selected landice-specific fields
!@+   for bulk processing operations (averaging over patches etc.)
         real*8, dimension(:,:,:), pointer ::
     &        srfflx_exports_gla=>null()!,srfstate_exports_gla
#ifdef TRACERS_ON
         real*8, dimension(:,:,:,:), pointer ::
     &        trsrfflx_exports_gla=>null()!,trsrfstate_exports_gla
#endif

         !@var E1 net energy flux at layer 1
         REAL*8, DIMENSION(:,:), POINTER ::
     &     E1
     &    ,TGRND,TGR4  ! temporary temps for surface flux calcs
         !@var IMPLM [kg m-2] implicit mass flux at bottom of domain
         !@var IMPLH [J m-2] implicit energy flux at bottom of domain
         REAL*8, DIMENSION(:,:), POINTER :: IMPLM,IMPLH
#ifdef TRACERS_WATER
!@var IMPLT implicit tracer flux at bottom of domain
         REAL*8, DIMENSION(:,:,:), POINTER :: IMPLT
#endif

         !@var SNOWLI,SNOWFR,SNOWDP snow mass, fraction, depth
         !@+   TODO: move to atmsrf_xchng_vars and unify the names
         !@+   for these in seaice/landice/land components.
c         REAL*8, DIMENSION(:,:), POINTER :: SNOWLI,SNOWFR,SNOWDP

      end type atmgla_xchng_vars
      ! -----------------------------------------------------

      type, extends(atmsrf_xchng_vars) :: atmlnd_xchng_vars
!@var bare_soil_wetness bare_soil_wetness (1)
         REAL*8, DIMENSION(:,:), POINTER :: bare_soil_wetness
!@var snowe snow amount seen by radiation
         REAL*8, DIMENSION(:,:), POINTER :: snowe
ccc FR_SNOW_RAD is snow fraction for albedo computations
ccc actually it should be the same as FR_SNOW_IJ but currently the snow
ccc model can't handle fractional cover for thick snow (will fix later)
         REAL*8, DIMENSION(:,:,:), POINTER :: FR_SNOW_RAD

      end type atmlnd_xchng_vars

      type, extends(simple_bounds_type) :: iceocn_xchng_vars
!@var grid a pointer to the grid object whose domain bounds were
!@+   used to allocate an instance of this type
         type(dist_grid), pointer :: grid

         REAL*8, DIMENSION(:,:), POINTER ::
!@var FWATER water fraction of gridbox
     &      FWATER
!@var RSI fraction of water area covered in ice
     &     ,RSI
!@var CORIOL coriolis parameter (1/s)
     &     ,CORIOL
!@var SOLAR solar radiation penetrating the ice absorbed by ocean [J m-2]
     &     ,SOLAR
!@var APRESS total atmos + sea ice pressure (at base of sea ice) (Pa)
     &     ,APRESS
!@var RUNOSI run off from sea/lake ice after surface (kg/m^2)
!@var ERUNOSI energy of run off from sea/lake ice after surface [J m-2]
!@var SRUNOSI salt in run off from sea/lake ice after surface (kg/m^2)
     &     ,RUNOSI, ERUNOSI, SRUNOSI
!@var MELTI,EMELTI,SMELTI mass,energy,salt from simelt into ocn (kg/m^2,J/m^2)
     &     ,MELTI, EMELTI, SMELTI
!@var RUNPSI run off from sea/lake ice after precip (kg/m^2)
!@var ERUNPSI energy of run off from sea/lake ice after precip [J m-2]
!@var SRUNPSI salt in run off from sea/lake ice after precip (kg/m^2)
     &     ,RUNPSI, ERUNPSI, SRUNPSI
!@var DMUI,DMVI momentum flux from sea ice to ocean (kg/m s)
!@+   On C grid for now
     &     ,DMUI, DMVI
!@var fmsi_io,fhsi_io,fssi_io basal ice-ocean fluxes (kg or J/m^2)
     &     ,fmsi_io,fhsi_io,fssi_io
!@var UI2rho Ustar*2*rho ice-ocean friction velocity on atmospheric grid
     &     ,UI2rho
!@var UOSURF,VOSURF ocean surface velocity (m/s)
     &     ,UOSURF,VOSURF
!@var OGEOZA ocean surface height geopotential (m^2/s^2)
     &     ,OGEOZA
!@var MLHC ocean mixed layer heat capacity (J/m^2 C) 
     &     ,MLHC
!@var MLDLK mixed layer depth in lake (m)
!@var DLAKE depth of lake (m)
!@var GLAKE lake heat content (J/m2)
     &     ,mldlk,dlake,glake ! lakes only

C**** DMSI,DHSI,DSSI are fluxes for ice formation within water column
!@var DMSI mass flux of sea ice 1) open water and 2) under ice (kg/m^2)
!@var DHSI energy flux of sea ice 1) open water and 2) under ice [J m-2]
!@var DSSI salt flux in sea ice 1) open water and 2) under ice (kg/m^2)
         REAL*8, DIMENSION(:,:,:), POINTER ::
     &      DMSI, DHSI, DSSI

#ifdef TRACERS_WATER
         REAL*8, DIMENSION(:,:,:), POINTER ::
!@var TRUNPSI tracer in run off from sea/lake ice after precip (kg/m^2)
!@var TRUNOSI tracer in run off from sea/lake ice after surface (kg/m^2)
!@var TRMELTI tracer from simelt into ocean (kg/m^2)
!@var ftrsi_io ice-ocean tracer fluxes under ice (kg/m^2)
     &     TRUNPSI, TRUNOSI, TRMELTI, ftrsi_io
#endif
#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER) /* huh? */
!@var DTRSI tracer flux in sea ice under ice and on open water (kg/m^2)
         REAL*8, DIMENSION(:,:,:,:), POINTER :: DTRSI
#endif

!@var dlatm latitudinal gridsize in minutes, for certain regrid routines
         real*8 :: dlatm
#ifdef CUBED_SPHERE
#else
!@var pack_a2i,pack_i2a contain info for redistributing data from
!@+   atmos. domains to icedyn domains and vice versa.
         type(band_pack_type), pointer :: pack_a2i,pack_i2a
#endif

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
!@var ntm number of tracers participating in sea ice formation.
         integer :: ntm=0
#endif

      end type iceocn_xchng_vars

      interface alloc_xchng_vars
        module procedure alloc_atmsrf_xchng_vars
        module procedure alloc_atmocn_xchng_vars
        module procedure alloc_atmice_xchng_vars
        module procedure alloc_atmgla_xchng_vars
        module procedure alloc_atmlnd_xchng_vars
        module procedure alloc_iceocn_xchng_vars
      end interface alloc_xchng_vars

      CONTAINS

      subroutine set_simple_bounds_type(grd_dum,bds)
      USE DOMAIN_DECOMP_1D, ONLY : getDomainBounds,DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(simple_bounds_type) :: BDS

      INTEGER :: I_0H,I_1H, J_0H,J_1H, I_0,I_1, J_0,J_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      call getDomainBounds(grd_dum,
     &     I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1,
     &     I_STRT_HALO=I_0H, I_STOP_HALO=I_1H,
     &     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      bds % I_0H = I_0H
      bds % I_1H = I_1H
      bds % J_0H = J_0H
      bds % J_1H = J_1H

      bds % I_0 = I_0
      bds % I_1 = I_1
      bds % J_0 = J_0
      bds % J_1 = J_1

      bds % HAVE_SOUTH_POLE = HAVE_SOUTH_POLE
      bds % HAVE_NORTH_POLE = HAVE_NORTH_POLE

      ALLOCATE( bds % IMAXJ(J_0H:J_1H) )

      bds % IMAXJ(:) = I_1
      IF(HAVE_SOUTH_POLE) bds%IMAXJ(J_0) = 1
      IF(HAVE_NORTH_POLE) bds%IMAXJ(J_1) = 1

      return
      end subroutine set_simple_bounds_type

      subroutine alloc_atmsrf_xchng_vars(grd_dum,this,that)
      USE CONSTANT, only : tf
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmsrf_xchng_vars) :: THIS
      TYPE(atmsrf_xchng_vars), optional :: THAT
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#ifdef TRACERS_ON
      integer :: ntm
#endif

      if(present(that)) then ! only using this routine to set pointers

        ! qtys that are not in any bundle
        that%itype4 = this%itype4
        that%surf_name = this%surf_name
        that%ftype => this%ftype
        that%fhc => this%fhc
        that%ipbl => this%ipbl

        call alloc_atmsrf_xchng_bundles(grd_dum,this,that,.true.)

      else

#ifdef TRACERS_ON
      ntm = this%ntm
#endif

      call set_simple_bounds_type(grd_dum,this%simple_bounds_type)

#ifndef CUBED_SPHERE
      allocate(this%sini(grd_dum%im_world),this%cosi(grd_dum%im_world))
#endif

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      ALLOCATE(
     &          this % FTYPE   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FHC     ( I_0H:I_1H , J_0H:J_1H ),

     &          this % LAT     ( I_0H:I_1H , J_0H:J_1H ),

     &          this % WORK1   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % WORK2   ( I_0H:I_1H , J_0H:J_1H ),

     &   STAT = IER)

      this % ftype = 0.

      allocate(this % ipbl(i_0h:i_1h,j_0h:j_1h))
      this % ipbl = 0

      call alloc_atmsrf_xchng_bundles(grd_dum,this,this,.false.)

! are the inits to zero necessary if make_bundle is doing them?
      this % runo(:,:) = 0.
      this % eruno(:,:) = 0.
      this % GTEMP = 0.    ! initialize at 0 C
      this % GTEMP2 = 0.   ! initialize at 0 C
      this % GTEMPR = TF   ! initialize at 273 K
      this % snow = 0.
      this % snowfr = 0.
      this % snowdp = 0.
#ifdef TRACERS_WATER
      this % truno(:,:,:) = 0.
#endif

      endif

      return
      end subroutine alloc_atmsrf_xchng_vars

      subroutine alloc_atmsrf_xchng_bundles(grd_dum,this_mem,this,
     &     ptrs_only)
      use bundle_maker
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmsrf_xchng_vars) :: THIS_MEM,THIS
      logical :: ptrs_only
      INTEGER :: I_0H, I_1H, J_1H, J_0H
#ifdef TRACERS_ON
      integer :: ntm
#endif

#ifdef TRACERS_ON
      ntm = this%ntm
#endif

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      call make_bundle(this_mem%atm_exports_phase1
     &     ,i_0h,i_1h,j_0h,j_1h
     &     ,this%srfp,this%srfpk,this%p1,this%am1,this%byam1
     &     ,this%prec,this%eprec
     &     ,this%cosz1,this%flong,this%fshort,this%trup_in_rad
     &     )
      if(ptrs_only)
     &     this%atm_exports_phase1 => this_mem%atm_exports_phase1

#ifdef TRACERS_ON
      call make_bundle_lij(this_mem%tratm_exports_phase1
     &     ,ntm,i_0h,i_1h,j_0h,j_1h
     &     ,this%trflux_prescr
#ifdef TRACERS_WATER
     &     ,this%trprec
#endif
     &     )
      if(ptrs_only)
     &     this%tratm_exports_phase1 => this_mem%tratm_exports_phase1
#endif

      call make_bundle(this_mem%atm_exports_phasesrf
     &     ,i_0h,i_1h,j_0h,j_1h
     &     ,this%temp1,this%q1,this%u1,this%v1
     &     )
      if(ptrs_only)
     &     this%atm_exports_phasesrf => this_mem%atm_exports_phasesrf


#ifdef TRACERS_ON
      call make_bundle_lij(this_mem%tratm_exports_phasesrf
     &     ,ntm,i_0h,i_1h,j_0h,j_1h
     &     ,this%trm1
     &     )
      if(ptrs_only)
     &    this%tratm_exports_phasesrf => this_mem%tratm_exports_phasesrf
#endif

      call make_bundle(this_mem%srfflx_exports
     &     ,i_0h,i_1h,j_0h,j_1h
     &     ,this%evapor
     &     ,this%sensht,this%latht,this%trheat,this%solar,this%e0
     &     ,this%dmua,this%dmva
     &     ,this%dTH1,this%dQ1,this%uflux1,this%vflux1
     &     ,this%runo,this%eruno
     &     )
      if(ptrs_only)
     &     this%srfflx_exports => this_mem%srfflx_exports

      call make_bundle(this_mem%srfstate_exports
     &     ,i_0h,i_1h,j_0h,j_1h
     &     ,this%gtemp,this%gtemp2,this%gtempr,this%gtemps
     &     ,this%snow,this%snowfr,this%snowdp
     &     )
      if(ptrs_only)
     &     this%srfstate_exports => this_mem%srfstate_exports

#ifdef TRACERS_ON
      call make_bundle_lij(this_mem%trsrfflx_exports
     &     ,ntm,i_0h,i_1h,j_0h,j_1h
     &     ,this%trsrfflx
#ifdef TRACERS_WATER
     &     ,this%trevapor
     &     ,this%truno
#endif
#ifdef TRACERS_DRYDEP
     &     ,this%trdrydep
#endif
     &     )
      if(ptrs_only)
     &     this%trsrfflx_exports => this_mem%trsrfflx_exports

      call make_bundle_lij(this_mem%trsrfstate_exports
     &     ,ntm,i_0h,i_1h,j_0h,j_1h
     &     ,this%gtracer
     &     )
      if(ptrs_only)
     &     this%trsrfstate_exports => this_mem%trsrfstate_exports
#endif

      call make_bundle(this_mem%pbl_exports
     &     ,i_0h,i_1h,j_0h,j_1h
     &     ,this%tsavg
     &     ,this%qsavg
     &     ,this%rsavg
     &     ,this%usavg
     &     ,this%vsavg
     &     ,this%wsavg
     &     ,this%tauavg
     &     ,this%tgvavg
     &     ,this%qgavg
     &     ,this%w2_l1
     &     ,this%gustiwind
     &     ,this%dblavg
     &     ,this%rhoavg
     &     ,this%ciaavg
     &     ,this%khsavg
     &     ,this%cmgs
     &     ,this%chgs
     &     ,this%cqgs
     &     ,this%ustar_pbl
     &     ,this%lmonin_pbl
     &     ,this%wspdf
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
     &     ,this%wsgcm
     &     ,this%wsubwd
     &     ,this%wsubtke
     &     ,this%wsubwm
#endif
#if (defined TRACERS_DRYDEP ) && (defined ACCMIP_LIKE_DIAGS)
     &     ,this%stomatal_dep_vel
#endif
     &     )
      if(ptrs_only)
     &     this%pbl_exports => this_mem%pbl_exports

#ifdef TRACERS_ON

      call make_bundle_lij(this_mem%trpbl_exports
     &     ,ntm,i_0h,i_1h,j_0h,j_1h
     &     ,this%travg,this%travg_byvol
#ifdef TRACERS_DRYDEP
     &     ,this%dep_vel
     &     ,this%gs_vel
     &     ,this%drydflx
#endif
     &     )
      if(ptrs_only)
     &     this%trpbl_exports => this_mem%trpbl_exports
#endif

      return
      end subroutine alloc_atmsrf_xchng_bundles

      subroutine avg_patches_pbl_exports(grid,patches,avg,rel)
      use domain_decomp_1d, only : dist_grid
      implicit none
      type(dist_grid) :: grid
c gfortran prob. if passed as class() args
c      class(atmsrf_xchng_vars) :: patches(:),avg
      type(atmsrf_xchng_vars) :: patches(:),avg
      logical, intent(in), optional :: rel
c
      integer :: i,j,k,l,np
      real*8 :: ft
      logical :: rel_
c

      np = size(patches)
      if(np == 1) return

      rel_ = .false.
      if(present(rel)) then
        rel_ = rel
      endif

      do l=1,size(avg%pbl_exports,3)
        avg%pbl_exports(:,:,l) = 0.
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          if (rel_) then
            ft = patches(k)%fhc(i,j)
          else
            ft = patches(k)%ftype(i,j)
          end if
          avg%pbl_exports(i,j,l) = avg%pbl_exports(i,j,l) +
     &         patches(k)%pbl_exports(i,j,l)*ft
        enddo
        enddo
        enddo
      enddo

#ifdef TRACERS_ON
      do l=1,size(avg%trpbl_exports,4)
        avg%trpbl_exports(:,:,:,l) = 0.
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          if (rel_) then
            ft = patches(k)%fhc(i,j)
          else
            ft = patches(k)%ftype(i,j)
          end if
          avg%trpbl_exports(:,i,j,l) = avg%trpbl_exports(:,i,j,l) +
     &         patches(k)%trpbl_exports(:,i,j,l)*ft
        enddo
        enddo
        enddo
      enddo
#endif

      return
      end subroutine avg_patches_pbl_exports

      subroutine avg_patches_srfflx_exports(grid,patches,avg,rel)
      use domain_decomp_1d, only : dist_grid
      implicit none
      type(dist_grid) :: grid
c gfortran prob. if passed as class() args
c      class(atmsrf_xchng_vars) :: patches(:),avg
      type(atmsrf_xchng_vars) :: patches(:),avg
      logical, intent(in), optional :: rel
c
      integer :: i,j,k,l,np
      real*8, dimension(:,:,:), allocatable :: ftype
      logical :: rel_
c

      np = size(patches)
      if(np == 1) return

      allocate(
     &     ftype(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,np))

      rel_ = .false.
      if(present(rel)) then
        rel_ = rel
      endif

      if(rel_) then
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          ftype(i,j,k) = patches(k)%fhc(i,j)
        enddo
        enddo
        enddo
      else
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          ftype(i,j,k) = patches(k)%ftype(i,j)
        enddo
        enddo
        enddo
      endif

      do l=1,size(avg%srfflx_exports,3)
        avg%srfflx_exports(:,:,l) = 0.
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          avg%srfflx_exports(i,j,l) = avg%srfflx_exports(i,j,l) +
     &         patches(k)%srfflx_exports(i,j,l)*ftype(i,j,k)
        enddo
        enddo
        enddo
      enddo

#ifdef TRACERS_ON
      do l=1,size(avg%trsrfflx_exports,4)
        avg%trsrfflx_exports(:,:,:,l) = 0.
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          avg%trsrfflx_exports(:,i,j,l) =
     &    avg%trsrfflx_exports(:,i,j,l) +
     &         patches(k)%trsrfflx_exports(:,i,j,l)*ftype(i,j,k)
        enddo
        enddo
        enddo
      enddo
#endif

      deallocate(ftype)

      return
      end subroutine avg_patches_srfflx_exports

      subroutine avg_patches_srfflx_exports_gla(grid,patches,avg,rel)
      use domain_decomp_1d, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(atmgla_xchng_vars) :: patches(:),avg
      logical, intent(in), optional :: rel
c
      integer :: i,j,k,l,np
      real*8, dimension(:,:,:), allocatable :: ftype
      logical :: rel_
c

      np = size(patches)
      if(np == 1) return

      allocate(
     &     ftype(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,np))

      rel_ = .false.
      if(present(rel)) then
        rel_ = rel
      endif

      if(rel_) then
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          ftype(i,j,k) = patches(k)%fhc(i,j)
        enddo
        enddo
        enddo
      else
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          ftype(i,j,k) = patches(k)%ftype(i,j)
        enddo
        enddo
        enddo
      endif

      do l=1,size(avg%srfflx_exports_gla,3)
        avg%srfflx_exports_gla(:,:,l) = 0.
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          avg%srfflx_exports_gla(i,j,l) =
     &    avg%srfflx_exports_gla(i,j,l) +
     &         patches(k)%srfflx_exports_gla(i,j,l)*ftype(i,j,k)
        enddo
        enddo
        enddo
      enddo

#ifdef TRACERS_ON
      do l=1,size(avg%trsrfflx_exports_gla,4)
        avg%trsrfflx_exports_gla(:,:,:,l) = 0.
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          avg%trsrfflx_exports_gla(:,i,j,l) =
     &    avg%trsrfflx_exports_gla(:,i,j,l) +
     &         patches(k)%trsrfflx_exports_gla(:,i,j,l)*ftype(i,j,k)
        enddo
        enddo
        enddo
      enddo
#endif

      deallocate(ftype)

      return
      end subroutine avg_patches_srfflx_exports_gla

      subroutine avg_patches_srfstate_exports(grid,patches,avg,rel)
      use domain_decomp_1d, only : dist_grid
      implicit none
      type(dist_grid) :: grid
c gfortran prob. if passed as class() args
c      class(atmsrf_xchng_vars) :: patches(:),avg
      type(atmsrf_xchng_vars) :: patches(:),avg
      logical, intent(in), optional :: rel
c
      integer :: i,j,k,l,np
      real*8, dimension(:,:,:), allocatable :: ftype
      logical :: rel_
c
      np = size(patches)
      if(np == 1) return

      allocate(
     &     ftype(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,np))

      rel_ = .false.
      if(present(rel)) then
        rel_ = rel
      endif

      if(rel_) then
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          ftype(i,j,k) = patches(k)%fhc(i,j)
        enddo
        enddo
        enddo
      else
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          ftype(i,j,k) = patches(k)%ftype(i,j)
        enddo
        enddo
        enddo
      endif

      do l=1,size(avg%srfstate_exports,3)
        avg%srfstate_exports(:,:,l) = 0.
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          avg%srfstate_exports(i,j,l) = avg%srfstate_exports(i,j,l) +
c    &         patches(k)%srfstate_exports(i,j,l)*ftype(i,j,k)
c workaround for uninitialized patches%srfstate_exports multiply by zero
     &    merge( patches(k)%srfstate_exports(i,j,l)*ftype(i,j,k),
     &           0d0, ftype(i,j,k) > 0 )
        enddo
        enddo
        enddo
      enddo

#ifdef TRACERS_ON
      do l=1,size(avg%trsrfstate_exports,4)
        avg%trsrfstate_exports(:,:,:,l) = 0.
        do k=1,np
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          avg%trsrfstate_exports(:,i,j,l) =
     &    avg%trsrfstate_exports(:,i,j,l) +
     &         patches(k)%trsrfstate_exports(:,i,j,l)*ftype(i,j,k)
        enddo
        enddo
        enddo
      enddo
#endif

      deallocate(ftype)

      return
      end subroutine avg_patches_srfstate_exports

      subroutine alloc_atmocn_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmocn_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#ifdef TRACERS_ON
      integer :: ntm
#endif

      this%itype4 = 1

#ifdef TRACERS_ON
      ntm = this%ntm
#endif

      call alloc_atmsrf_xchng_vars(grd_dum,this%atmsrf_xchng_vars)
      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO
      ALLOCATE( this % FOCEAN  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % GMELT   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % EGMELT  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FLOWO   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % EFLOWO  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % UOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % OGEOZA  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SSS     ( I_0H:I_1H , J_0H:J_1H ),
#ifdef STANDALONE_OCEAN
     &          this % SSSOBS  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SSSRESFAC ( I_0H:I_1H , J_0H:J_1H ),
     &          this % RSIOBS  ( I_0H:I_1H , J_0H:J_1H ),
#endif
     &          this % MLHC    ( I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_WATER
     &          this % TRFLOWO ( NTM , I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_OCEAN
     &          this % TRGMELT ( NTM , I_0H:I_1H , J_0H:J_1H ),
#endif
#endif
     &          this%TRGASEX(this%gasex_index%getsize() ,
     &                                 I_0H:I_1H , J_0H:J_1H ),
     &          this % CHL     ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)

      this % UOSURF = 0.
      this % VOSURF = 0.
      this % OGEOZA = 0.
      this % MLHC = 0.

      this % TRGASEX = 0.
! I dimensioned trgasex by ntm_gasexch to avoid confusion elsewhere
! in the code.  If ntm differs from ntm_gasexch, fluxes need to be
! stored in the appropriate positions in trgasex or this%trgasex.  - MK

      this % CHL = 0.
      this%chl_defined=.false.

#ifndef STANDALONE_OCEAN
        allocate(this % DIRVIS  ( I_0H:I_1H , J_0H:J_1H ),
     &           this % DIFVIS  ( I_0H:I_1H , J_0H:J_1H ),
     &           this % DIRNIR  ( I_0H:I_1H , J_0H:J_1H ),
     &           this % DIFNIR  ( I_0H:I_1H , J_0H:J_1H ),
     &           STAT = IER)
#endif

      this % modd5s = -999

      return
      end subroutine alloc_atmocn_xchng_vars

      subroutine alloc_atmice_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmice_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#ifdef TRACERS_WATER
      integer :: ntm
#endif

      this%itype4 = 2

#ifdef TRACERS_WATER
      ntm = this%ntm
#endif

      call alloc_atmsrf_xchng_vars(grd_dum,this%atmsrf_xchng_vars)
      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO
      ALLOCATE( this % E1 ( I_0H:I_1H , J_0H:J_1H ),
     &          this % UISURF ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VISURF ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FWSIM  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MSICNV ( I_0H:I_1H , J_0H:J_1H ),
     &          this % HSICNV ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FOCEAN ( I_0H:I_1H , J_0H:J_1H ),
     &          this % USI    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VSI    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MUSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MVSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % HUSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % HVSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SUSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SVSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SNOAGE ( I_0H:I_1H , J_0H:J_1H ),
     &          this % ZSNOWI ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
      this % UISURF = 0.
      this % VISURF = 0.
      this % MSICNV = 0.
      this % HSICNV = 0.
      this % SNOAGE = 0.
      this % ZSNOWI = 0.

#ifdef TRACERS_WATER
      ALLOCATE( this % TUSI   (I_0H:I_1H, J_0H:J_1H, NTM),
     &          this % TVSI   (I_0H:I_1H, J_0H:J_1H, NTM),
     &          this % TRSIsum(NTM, I_0H:I_1H, J_0H:J_1H),
     &   STAT = IER)
#endif

      ALLOCATE( this % MSIsave(I_0H:I_1H, J_0H:J_1H),
     &          this % SNOWsave(I_0H:I_1H, J_0H:J_1H),
     &          this % TICEsave(I_0H:I_1H, J_0H:J_1H),
     &          this % TI1save(I_0H:I_1H, J_0H:J_1H),
     &          this % SSI1save(I_0H:I_1H, J_0H:J_1H),
     &          this % SSI2save(I_0H:I_1H, J_0H:J_1H),
     &          this % SIHC(I_0H:I_1H, J_0H:J_1H),
     &          this % RSIstart(I_0H:I_1H, J_0H:J_1H),
     &          this % SNTOSI(I_0H:I_1H, J_0H:J_1H),
     &          this % SITOPMLT(I_0H:I_1H, J_0H:J_1H),
     &          this % MSNFLOOD(I_0H:I_1H, J_0H:J_1H),
     &          this % HSNFLOOD(I_0H:I_1H, J_0H:J_1H),
     &          this % MSIsave2(I_0H:I_1H, J_0H:J_1H),
     &          this % SNOWsave2(I_0H:I_1H, J_0H:J_1H),
     &   STAT = IER)

      this % MSIsave = 0.
      this % SNOWsave = 0.
      this % TICEsave = 0.
      this % TI1save = 0.
      this % SIHC = 0.
      this % RSIstart = 0.
      this % SNTOSI = 0.
      this % SITOPMLT = 0.
      this % MSNFLOOD = 0.
      this % HSNFLOOD = 0.

      return
      end subroutine alloc_atmice_xchng_vars

      subroutine alloc_atmgla_xchng_vars(grd_dum,this)
      use bundle_maker
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmgla_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#ifdef TRACERS_ON
      integer :: ntm
#endif

#ifdef TRACERS_ON
      ntm = this%ntm
#endif

      this%itype4 = 3

      call alloc_atmsrf_xchng_vars(grd_dum,this%atmsrf_xchng_vars)
      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      ALLOCATE( this % TGRND ( I_0H:I_1H , J_0H:J_1H ),
     &          this % TGR4  ( I_0H:I_1H , J_0H:J_1H ),
     &     STAT = IER)

      call make_bundle(this%srfflx_exports_gla
     &     ,i_0h,i_1h,j_0h,j_1h
     &     ,this%e1,this%implm,this%implh
     &     )

#ifdef TRACERS_WATER
      call make_bundle_lij(this%trsrfflx_exports_gla
     &     ,ntm,i_0h,i_1h,j_0h,j_1h
     &     ,this%implt
     &     )
#endif

      return
      end subroutine alloc_atmgla_xchng_vars

      subroutine alloc_atmlnd_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmlnd_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      this%itype4 = 4

      call alloc_atmsrf_xchng_vars(grd_dum,this%atmsrf_xchng_vars)
      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO
      ALLOCATE(
     &     this % bare_soil_wetness ( I_0H:I_1H , J_0H:J_1H ),
     &     this % snowe ( I_0H:I_1H , J_0H:J_1H ),
     &     this % fr_snow_rad(2,i_0h:i_1h,j_0h:j_1h),
     &     STAT = IER)
      this % snowe = 0.
      return
      end subroutine alloc_atmlnd_xchng_vars

      subroutine alloc_iceocn_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(iceocn_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      integer :: ntm
#endif

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      ntm = this%ntm
#endif


      call set_simple_bounds_type(grd_dum,this%simple_bounds_type)

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      ALLOCATE( this % RUNOSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          this % ERUNOSI ( I_0H:I_1H , J_0H:J_1H ), 
     &          this % SRUNOSI ( I_0H:I_1H , J_0H:J_1H ),
     &          this % RUNPSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          this % SRUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          this % ERUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DMUI    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DMVI    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MELTI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % EMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % APRESS  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % RSI     ( I_0H:I_1H , J_0H:J_1H ),
     &          this % fmsi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          this % fhsi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          this % fssi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          this % UI2rho  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SOLAR   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FWATER  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % CORIOL  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % UOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % OGEOZA  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MLHC    ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT=IER)

       ALLOCATE( this % DMSI  (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           this % DHSI  (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           this % DSSI  (  2  , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
       this % DMSI = 0.
       this % DHSI = 0.
       this % DSSI = 0.

#ifdef TRACERS_WATER
      ALLOCATE( this % TRUNPSI (NTM, I_0H:I_1H, J_0H:J_1H),
     &          this % TRUNOSI (NTM, I_0H:I_1H, J_0H:J_1H),
     &          this % TRMELTI (NTM, I_0H:I_1H, J_0H:J_1H),
     &          this % FTRSI_IO(NTM, I_0H:I_1H, J_0H:J_1H),
     &   STAT = IER)
      this % TRUNPSI = 0.
      this % TRUNOSI = 0.
      this % TRMELTI = 0.
      this % FTRSI_IO = 0.
#endif
#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER) /* huh? */
      ALLOCATE( this % DTRSI(NTM, 2, I_0H:I_1H, J_0H:J_1H),
     &   STAT = IER)
      this % DTRSI = 0.
#endif

      return
      end subroutine alloc_iceocn_xchng_vars

      END MODULE EXCHANGE_TYPES

#ifndef STANDALONE_OCEAN

      MODULE FLUXES
!@sum  FLUXES contains the fluxes between various atm-grid components
!@auth Gavin Schmidt
      USE RESOLUTION, only : im,jm,lm
      USE DOMAIN_DECOMP_ATM, ONLY : grid
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,atmice_xchng_vars,
     &     atmsrf_xchng_vars,atmlnd_xchng_vars,atmgla_xchng_vars
#ifdef TRACERS_ON
#ifndef SKIP_TRACER_SRCS
      USE Tracer_mod, only: ntsurfsrcmax,nt3Dsrcmax
#endif
#endif
      use Dictionary_mod, only : sync_param, get_param
      IMPLICIT NONE

!@dbparam NIsurf: DT_Surface  =  DTsrc/NIsurf
      INTEGER :: NIsurf = 2

!@var Fxx fraction of gridbox of type xx (land,ocean,...)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLAND
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FOCEAN
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLICE
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FEARTH0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLAKE0

#ifdef GLINT2
!@var FLICE_ICEMODEL Fraction of gridbox that is landice that
!     comes from a GLINT2-related ice model.
!     NOTE: FLICE_GLINT2 < FLICE
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLICE_GLINT2
#endif

!@param NSTYPE number of surface types for radiation purposes
      !INTEGER, PARAMETER :: NSTYPE=4

!@dbparam UOdrag parameter that decides whether ocean.ice velocities
!@+   feed into drag calculation in surface (default = 0)
      INTEGER :: UOdrag = 0

!@var uflux1 surface turbulent u-flux (=-<uw>) 
!@var vflux1 surface turbulent v-flux (=-<vw>)
!@var tflux1 surface turbulent t-flux (=-<tw>)
!@var qflux1 surface turbulent q-flux (=-<qw>)
      real*8, allocatable, dimension(:,:) :: 
     &        uflux1,vflux1,tflux1,qflux1

C**** The E/FLOWO, E/S/MELTI, E/GMELT arrays are used to flux quantities 
C**** to the ocean that are not tied to the open water/ice covered 
C**** fractions. This is done separately for river flow, complete
C**** sea ice melt and iceberg/glacial melt.

!@var PREC precipitation (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PREC
!@var EPREC energy of preciptiation [J m-2]
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: EPREC
!@var PRECSS precipitation from super-saturation (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PRECSS

#ifdef IRRIGATION_ON
!@var Actual irrigation rate (& energy) [m/s] [W/m2]
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:) :: irrig_water_act
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:) :: irrig_energy_act
#ifdef TRACERS_WATER
!@var Actual irrigation tracer rate [kg/s] 
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: irrig_tracer_act
#endif
#endif

C**** fluxes associated with variable lake fractions
!@var DMWLDF  water deficit over land surface (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DMWLDF
!@var DGML energy associated with DMWLDF (J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DGML

#ifdef TRACERS_ON
!@var TRSOURCE non-interactive surface sources/sinks for tracers (kg/s)
#ifndef SKIP_TRACER_SRCS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trsource
#endif
!@var TRFLUX1 total surface flux for each tracer (kg/m2/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: trflux1
!@var TR3DSOURCE 3D sources/sinks for tracers (kg/s)
#ifndef SKIP_TRACER_SRCS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:):: tr3Dsource
#endif

#ifdef TRACERS_WATER
!@var TRPREC tracers in precip (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: TRPREC

C**** fluxes associated with variable lake fractions
!@var DTRL tracers associate with DMWLDF (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DTRL

#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
!@var trprec_dust dust/mineral tracers in precip [kg]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:):: trprec_dust
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
!@var pprec precipitation at previous time step [kg/m^2]
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: pprec
!@var pevap evaporation at previous time step (land only) [kg/m^2]
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: pevap
!@var dust_flux_glob global array of dust emission flux [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: dust_flux_glob
#endif
!@var dust_flux2_glob global array of cubic dust emission flux (for diags only)
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: dust_flux2_glob

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
#ifdef TRACERS_DRYDEP
!@var depo_turb_glob global array of flux due to dry turb. dep. of tracers
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: depo_turb_glob
!@var depo_grav_glob global array of flux due to gravit. settling of tracers
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: depo_grav_glob
#endif
#endif

#endif

!@param after_atm_phase1, during_srfflx flags to pass to
!@+     routines in SURFACE indicating the calling point
      integer, parameter ::
     &     after_atm_phase1=0, during_srfflx=1

!@var atm{ocn,ice,gla,lnd}s derived-type strucures containing
!@+   variables needed for atmospheric interactions with
!@+   water, floating ice, glacial ice, and the land surface.
      type(atmocn_xchng_vars) :: atmocns(1) ! ocean and lakes
      type(atmice_xchng_vars) :: atmices(1) ! ocean and lakes
#ifdef GLINT2
      type(atmgla_xchng_vars), allocatable, dimension(:) :: atmglas_hp !(1-min(nhc,2)/2:nhc) ! glacial ice
#endif
      type(atmgla_xchng_vars), allocatable, dimension(:) :: atmglas !(1-min(nhc,2)/2:nhc) ! glacial ice
      type(atmlnd_xchng_vars) :: atmlnds(1) ! land surface

!@var atm{ocn,ice,gla,lnd} pointers to the index of atm{ocn,ice,gla,lnd}s
!@+   containing the average over patches (index = either 0 or 1, see notes below)
      type(atmocn_xchng_vars), pointer :: atmocn ! ocean and lakes
      type(atmice_xchng_vars), pointer :: atmice ! ocean and lakes
      type(atmgla_xchng_vars), pointer :: atmgla ! glacial ice
      type(atmlnd_xchng_vars), pointer :: atmlnd ! land surface

      target :: atmocns,atmices,atmglas,atmlnds
#ifdef GLINT2
      target :: atmglas_hp
#endif

!@var atmsrf contains atm-surf interaction quantities averaged over
!@+   all surface types.
      type(atmsrf_xchng_vars) :: atmsrf

!@var asflx4 an array for looping over atmocn,atmice,atmgla,atmlnd
      type(atmsrf_xchng_vars), dimension(4) :: asflx4

!@var asflx an array for looping over atmocns,atmices,atmglas,atmlnds
      integer :: !, parameter ::
     &     nptchs
!     &     =ubound(atmocns,1)
!     &     +ubound(atmices,1)
!     &     +ubound(atmglas,1)
!     &     +ubound(atmlnds,1)
!      type(atmsrf_xchng_vars), dimension(nptchs) :: asflx
      type(atmsrf_xchng_vars), dimension(:), allocatable :: asflx

!@param p[12]xxx lower and upper bounds for a given surface type in
!@+     the asflx array
      integer :: !, parameter ::
     &     p1ocn != 1
     &    ,p2ocn != p1ocn+ubound(atmocns,1)-1
     &    ,p1ice != 1+p2ocn
     &    ,p2ice != p1ice+ubound(atmices,1)-1
     &    ,p1gla != 1+p2ice
     &    ,p2gla != p1gla+ubound(atmglas,1)-1
     &    ,p1lnd != 1+p2gla
     &    ,p2lnd != p1lnd+ubound(atmlnds,1)-1

! Notes on atmxxxs, asflx arrays:
!
! Currently statically sized - will eventually become allocatable
! when the number of possible patches per type becomes dynamic.
!
! If a surface type xxx has only one patch, the lower bound of
! the array atmxxxs is equal to 1, and atmxxx => atmxxxs(1)
! If a surface type has multiple patches, the lower bound of
! the array atmxxxs is set to 0, with the index 0 reserved for the
! average over patches of that type (atmxxx => atmxxxs(0))
!
! Positions 1-4 of the array asflx4 point to the same elements of atmxxxs
! as atmocn,atmice,atmgla,atmlnd respectively.
!
! asflx(p1ocn:p2ocn) => atmocns(1:)
! asflx(p1ice:p2ice) => atmices(1:)
! asflx(p1gla:p2gla) => atmglas(1:)
! asflx(p1lnd:p2lnd) => atmlnds(1:)
!
! When all surface types have multiple patches:
!
! asflx4(1,2,3,4), atm{ocn,ice,gla,lnd}  => atm{ocn,ice,gla,lnd}s(0)
!
! When each surface type has only one patch:
!
! asflx4(1,2,3,4), atm{ocn,ice,gla,lnd}  => atm{ocn,ice,gla,lnd}s(1)
!
! For asflx and asflx4, the "point to" operation for a given index
! happens via a call to alloc_xchng_vars which sets all the pointers
! individually.  This implementation will likely be improved.

      END MODULE FLUXES

      SUBROUTINE ALLOC_FLUXES !(grd_dum)
!@sum   Initializes FLUXES arrays
!@auth  Rosalinda de Fainchtein
      USE EXCHANGE_TYPES, only : alloc_xchng_vars
      USE DOMAIN_DECOMP_ATM, ONLY : GRD_DUM=>GRID,
     &     HALO_UPDATE,HASSOUTHPOLE,HASNORTHPOLE
#ifdef GLINT2
      USE DOMAIN_DECOMP_ATM, ONLY : glint2
      use glint2_modele
#endif
      USE FLUXES
      USE GEOM, only : lat2d

#ifndef CUBED_SPHERE
#ifndef SCM
#define USE_DLATM
#endif
#endif

#ifdef USE_DLATM
      USE GEOM, only : dlatm,sinip,cosip
#endif
#ifdef TRACERS_ON
      USE tracer_com,ONLY : NTM
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
     &     ,Ntm_dust
#endif
      use tracer_com, only : gasex_index, n_co2n
#endif
      USE ATM_COM, only : temperature_istart1
      USE Dictionary_mod
      use pario, only : par_open,par_close,read_dist_data

      IMPLICIT NONE
      !TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      INTEGER :: fid
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: I, J, I_0, I_1, J_1, J_0
      INTEGER :: IER,K
      character(len=2) :: c2
      integer :: NHC_LOCAL = 1

#ifdef GLINT2
      nhc_local = glint2_modele_nhc(glint2)
#else
      call sync_param( "NHC", nhc_local)
#endif
      call sync_param( "NIsurf", NIsurf )
      call sync_param( "UOdrag", UOdrag )

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      I_0 = grd_dum%I_STRT
      I_1 = grd_dum%I_STOP
      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP
!      print *,'ALLOCATE atmglas', 1-min(nhc_local,2)/2, nhc_local
#ifdef GLINT2
      ALLOCATE(atmglas_hp(1-min(nhc_local,2)/2:nhc_local), STAT=IER)
#endif
      ! nhc=1 --> lbound = 1
      ! nhc>1 --> lbound = 0   (zeroth element is for sum)
      ALLOCATE(atmglas(1-min(nhc_local,2)/2:nhc_local), STAT=IER)
      ALLOCATE(FLAND(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(FOCEAN(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(FLICE(I_0H:I_1H,J_0H:J_1H), STAT = IER)
#ifdef GLINT2
      ALLOCATE(FLICE_GLINT2(I_0H:I_1H,J_0H:J_1H), STAT = IER)
#endif
      ALLOCATE(FEARTH0(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(FLAKE0(I_0H:I_1H,J_0H:J_1H), STAT = IER)

C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
C**** Note that FLAKE0 is read in only to provide initial values
C**** Actual array is set from restart file.
      fid = par_open(grid,'TOPO','read')
      call read_dist_data(grid,fid,'focean',FOCEAN)
      call read_dist_data(grid,fid,'flake',FLAKE0)
      call read_dist_data(grid,fid,'fgrnd',FEARTH0)
      call read_dist_data(grid,fid,'fgice',FLICE)
      call par_close(grid,fid)

#ifdef GLINT2
        ! Fix it up with GLINT2
        call glint2_modele_compute_fgice(glint2, .true.,
     &         FLICE_GLINT2,
     &         FLICE, FEARTH0, FOCEAN, FLAKE0,
     &     grid%i_strt_halo, grid%j_strt_halo)
#endif

      CALL HALO_UPDATE(GRD_DUM, FOCEAN)
      CALL HALO_UPDATE(GRD_DUM, FEARTH0)


C**** Deal with single -> double precision problems and potential
C**** ocean/lake inconsistency. Adjust FLAKE0 and FLICE if necessary.
      DO J=J_0,J_1
      DO I=I_0,I_1 !IMAXJ(J)
        IF (FOCEAN(I,J).gt.0) THEN
          FLAND(I,J)=1.-FOCEAN(I,J) ! Land fraction if focean>0
          IF (FLAKE0(I,J).gt.0) THEN
            WRITE(6,*) "Ocean and lake cannot co-exist in same grid box"
     *       ,i,j,FOCEAN(I,J),FLAKE0(I,J)
            FLAKE0(I,J)=0
          END IF
        ELSEIF (FLAKE0(I,J).gt.0) THEN
          FLAND(I,J)=1.-FLAKE0(I,J)  ! for initialization only
        ELSE
          FLAND(I,J)=1.              ! for initialization only
        END IF
C**** Ensure that no round off error effects land with ice and earth
        IF (FLICE(I,J)-FLAND(I,J).gt.-1d-4 .and. FLICE(I,J).gt.0) THEN
          FLICE(I,J)=FLAND(I,J)
        END IF
      END DO
      END DO
      CALL HALO_UPDATE(GRD_DUM, FLAKE0)
      CALL HALO_UPDATE(GRD_DUM, FLICE)

      If (HASSOUTHPOLE(GRD_DUM)) Then
         FLAND(2:IM,1)=FLAND(1,1)
         FLICE(2:IM,1)=FLICE(1,1)
      End If
      If (HASNORTHPOLE(GRD_DUM)) Then
         FLAND(2:IM,JM)=FLAND(1,JM)
         FLICE(2:IM,JM)=FLICE(1,JM)
      End If

      !I-J arrays
      ALLOCATE( uflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          vflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          tflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          qflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          PREC    ( I_0H:I_1H , J_0H:J_1H ),
     &          EPREC   ( I_0H:I_1H , J_0H:J_1H ),
     &          PRECSS  ( I_0H:I_1H , J_0H:J_1H ),
     &          DMWLDF  ( I_0H:I_1H , J_0H:J_1H ),
     &          DGML    ( I_0H:I_1H , J_0H:J_1H ),
#ifdef IRRIGATION_ON
     &          irrig_water_act ( I_0H:I_1H , J_0H:J_1H ),
     &          irrig_energy_act( I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_WATER
     &          irrig_tracer_act(NTM, I_0H:I_1H , J_0H:J_1H),
#endif
#endif
     &   STAT=IER)

!TRACERS_ON**********************

#ifdef TRACERS_ON
      !(I,J,:,:)  array
#ifndef SKIP_TRACER_SRCS
      ALLOCATE(trsource (I_0H:I_1H,J_0H:J_1H,ntsurfsrcmax,NTM)
     &  ,STAT = IER)
      trsource = 0.
#endif

      !(I,J,:) arrays
      ALLOCATE(trflux1 ( I_0H:I_1H , J_0H:J_1H , NTM    ),
     &   STAT = IER)

      !I-J-L-:-: array
#ifndef SKIP_TRACER_SRCS
      ALLOCATE( tr3Dsource(I_0H:I_1H,J_0H:J_1H,LM,nt3Dsrcmax,NTM)
     &  ,STAT = IER)
#endif

#ifdef TRACERS_WATER
                                                    !(:)-(:)-I-J arrays
       !:-I-J arrays
       ALLOCATE( TRPREC  ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           DTRL    ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      ALLOCATE(trprec_dust(Ntm_dust,I_0H:I_1H ,J_0H:J_1H),STAT=ier)
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      ALLOCATE(pprec(I_0H:I_1H,J_0H:J_1H),STAT = IER)
      pprec = 0
      ALLOCATE(pevap(I_0H:I_1H,J_0H:J_1H),STAT = IER)
      pevap = 0
      ALLOCATE(dust_flux_glob(I_0H:I_1H,J_0H:J_1H,Ntm_dust),STAT = IER)
#ifdef TRACERS_DRYDEP
      ALLOCATE(depo_turb_glob(I_0H:I_1H,J_0H:J_1H,Ntm)
     &     ,STAT = IER)
      ALLOCATE(depo_grav_glob(I_0H:I_1H,J_0H:J_1H,Ntm)
     &     ,STAT = IER)
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      ALLOCATE(dust_flux2_glob(I_0H:I_1H,J_0H:J_1H,Ntm_dust),STAT = IER)
#endif

#endif

      atmsrf%surf_name = 'srf'
#ifdef TRACERS_ON
      atmsrf%ntm = ntm
#endif
      call alloc_xchng_vars(grd_dum,atmsrf)
      atmsrf%TSAVG(:,:)=temperature_istart1

      do k=lbound(atmocns,1),ubound(atmocns,1)
        write(c2,'(i2.2)') k
        atmocns(k)%surf_name = 'ocn'//c2
#ifdef TRACERS_ON
        atmocns(k)%ntm = ntm
        atmocns(k)%gasex_index = gasex_index
        atmocns(k)%n_co2n = n_co2n
#endif
        call alloc_xchng_vars(grid,atmocns(k))
        atmocns(k)%grid => grd_dum
        atmocns(k)%lat(:,:) = lat2d(:,:)
        atmocns(k)%focean(:,:) = focean(:,:)
#ifdef USE_DLATM
        atmocns(k)%dlatm = dlatm
        atmocns(k)%sini(:) = sinip(:)
        atmocns(k)%cosi(:) = cosip(:)
#endif
      enddo

      do k=lbound(atmices,1),ubound(atmices,1)
        write(c2,'(i2.2)') k
        atmices(k)%surf_name = 'ice'//c2
#ifdef TRACERS_ON
        atmices(k)%ntm = ntm
#endif
        call alloc_xchng_vars(grid,atmices(k))
        atmices(k)%grid => grd_dum
        atmices(k)%lat(:,:) = lat2d(:,:)
        atmices(k)%focean(:,:) = focean(:,:)
#ifdef USE_DLATM
        atmices(k)%dlatm = dlatm
        atmices(k)%sini(:) = sinip(:)
        atmices(k)%cosi(:) = cosip(:)
#endif
      enddo

      do k=lbound(atmglas,1),ubound(atmglas,1)
        write(c2,'(i2.2)') k
#ifdef GLINT2
        atmglas_hp(k)%surf_name = 'gla'//c2
#endif
        atmglas(k)%surf_name = 'gla'//c2

#ifdef TRACERS_ON
#ifdef GLINT2
        atmglas_hp(k)%ntm = ntm
#endif
        atmglas(k)%ntm = ntm
#endif
#ifdef GLINT2
        call alloc_xchng_vars(grid,atmglas_hp(k))
#endif
        call alloc_xchng_vars(grid,atmglas(k))
#ifdef GLINT2
        atmglas_hp(k)%grid => grd_dum
        atmglas_hp(k)%lat(:,:) = lat2d(:,:)
#endif
        atmglas(k)%grid => grd_dum
        atmglas(k)%lat(:,:) = lat2d(:,:)
      enddo

      do k=lbound(atmlnds,1),ubound(atmlnds,1)
        write(c2,'(i2.2)') k
        atmlnds(k)%surf_name = 'lnd'//c2
#ifdef TRACERS_ON
        atmlnds(k)%ntm = ntm
#endif
        call alloc_xchng_vars(grid,atmlnds(k))
        atmlnds(k)%grid => grd_dum
        atmlnds(k)%lat(:,:) = lat2d(:,:)
      enddo

! set pointers for composite surface types and individual surface patches

      nptchs =
     &      ubound(atmocns,1)
     &     +ubound(atmices,1)
     &     +ubound(atmglas,1)
     &     +ubound(atmlnds,1)
      allocate(asflx(nptchs))

      p1ocn = 1
      p2ocn = p1ocn+ubound(atmocns,1)-1
      p1ice = 1+p2ocn
      p2ice = p1ice+ubound(atmices,1)-1
      p1gla = 1+p2ice
      p2gla = p1gla+ubound(atmglas,1)-1
      p1lnd = 1+p2gla
      p2lnd = p1lnd+ubound(atmlnds,1)-1

      atmocn => atmocns(lbound(atmocns,1))
      atmice => atmices(lbound(atmices,1))
      atmgla => atmglas(lbound(atmglas,1))
      atmlnd => atmlnds(lbound(atmlnds,1))

      call alloc_xchng_vars(grid,atmocn%atmsrf_xchng_vars,asflx4(1))
      call alloc_xchng_vars(grid,atmice%atmsrf_xchng_vars,asflx4(2))
      call alloc_xchng_vars(grid,atmgla%atmsrf_xchng_vars,asflx4(3))
      call alloc_xchng_vars(grid,atmlnd%atmsrf_xchng_vars,asflx4(4))

      do k=p1ocn,p2ocn
        call alloc_xchng_vars(grid,
     &       atmocns(1+k-p1ocn)%atmsrf_xchng_vars,asflx(k))
      enddo
      do k=p1ice,p2ice
        call alloc_xchng_vars(grid,
     &       atmices(1+k-p1ice)%atmsrf_xchng_vars,asflx(k))
      enddo
      do k=p1gla,p2gla
#ifdef GLINT2
        call alloc_xchng_vars(grid,
     &       atmglas_hp(1+k-p1gla)%atmsrf_xchng_vars,asflx(k))
#endif
        call alloc_xchng_vars(grid,
     &       atmglas(1+k-p1gla)%atmsrf_xchng_vars,asflx(k))
      enddo
      do k=p1lnd,p2lnd
        call alloc_xchng_vars(grid,
     &       atmlnds(1+k-p1lnd)%atmsrf_xchng_vars,asflx(k))
      enddo

      END SUBROUTINE ALLOC_FLUXES

      subroutine def_rsf_fluxes(fid)
!@sum  def_rsf_fluxes defines structure of coupler arrays in restart files
!@auth M. Kelley
!@ver  beta
      use fluxes
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      character(len=17) :: ijstr
      character(len=64) :: dimstr,vname
      integer :: ipatch

      ijstr='(dist_im,dist_jm)'
      !call defvar(agrid,fid,atmocn%gtemp,'gtemp'//ijstr)
      !call defvar(agrid,fid,atmocn%gtempr,'gtempr'//ijstr)
      call defvar(grid,fid,atmocn%gtemp,'asst'//ijstr)
      call defvar(grid,fid,atmocn%gtempr,'atempr'//ijstr)
      call defvar(grid,fid,atmocn%sss,'sss'//ijstr)
      call defvar(grid,fid,atmocn%ogeoza,'ogeoza'//ijstr)
      call defvar(grid,fid,atmocn%uosurf,'uosurf'//ijstr)
      call defvar(grid,fid,atmocn%vosurf,'vosurf'//ijstr)
      call defvar(grid,fid,atmocn%mlhc,'mlhc'//ijstr)
#ifdef TRACERS_ON
      call defvar(grid,fid,atmocn%gtracer,
     &                     'gtracer(ntm,dist_im,dist_jm)')
#endif

      do ipatch = 1,size(asflx)
        dimstr='_'//trim(asflx(ipatch)%surf_name)// 
     &       '(npbl,dist_im,dist_jm)'
        vname = 'uabl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%uabl,vname)
        vname = 'vabl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%vabl,vname)
        vname = 'tabl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%tabl,vname)
        vname = 'qabl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%qabl,vname)
        vname = 'eabl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%eabl,vname)
        dimstr='_'//trim(asflx(ipatch)%surf_name)//ijstr
        vname = 'cmgs'//dimstr
        call defvar(grid,fid,asflx(ipatch)%cmgs,vname)
        vname = 'chgs'//dimstr
        call defvar(grid,fid,asflx(ipatch)%chgs,vname)
        vname = 'cqgs'//dimstr
        call defvar(grid,fid,asflx(ipatch)%cqgs,vname)
        vname = 'ipbl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%ipbl,vname)
        vname = 'ustar_pbl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%ustar_pbl,vname)
        vname = 'lmonin_pbl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%lmonin_pbl,vname)
#ifdef TRACERS_ON
        dimstr='_'//trim(asflx(ipatch)%surf_name)// 
     &       '(npbl,ntm,dist_im,dist_jm)'
        vname = 'trabl'//dimstr
        call defvar(grid,fid,asflx(ipatch)%trabl,vname)
#endif
      enddo

      call defvar(grid,fid,atmsrf%tsavg,'tsavg(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%qsavg,'qsavg(dist_im,dist_jm)')
c      call defvar(grid,fid,atmsrf%dclev,'dclev(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%usavg,'usavg(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%vsavg,'vsavg(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%wsavg,'wsavg(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%tauavg,'tauavg(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%tgvavg,'tgvavg(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%qgavg,'qgavg(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%ustar_pbl
     &   ,'ustar_pbl(dist_im,dist_jm)')
      call defvar(grid,fid,atmsrf%lmonin_pbl
     &   ,'lmonin_pbl(dist_im,dist_jm)')

      return
      end subroutine def_rsf_fluxes

      subroutine new_io_fluxes(fid,iaction)
!@sum  new_io_fluxes read/write coupler arrays from/to restart files
!@auth M. Kelley
!@ver  beta
      use model_com, only : iowrite,ioread
      use fluxes
      use domain_decomp_atm, only: grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      character(len=16) :: suffix
      character(len=64) :: vname
      integer :: ipatch

      select case (iaction)
      case (iowrite)            ! output to restart file
        !call write_dist_data(grid,fid,'gtemp',atmocn%gtemp)
        !call write_dist_data(grid,fid,'gtempr',atmocn%gtempr)
        call write_dist_data(grid,fid,'asst',atmocn%gtemp)
        call write_dist_data(grid,fid,'atempr',atmocn%gtempr)
        call write_dist_data(grid,fid,'sss',atmocn%sss)
        call write_dist_data(grid,fid,'ogeoza',atmocn%ogeoza)
        call write_dist_data(grid,fid,'uosurf',atmocn%uosurf)
        call write_dist_data(grid,fid,'vosurf',atmocn%vosurf)
        call write_dist_data(grid,fid,'mlhc',atmocn%mlhc)
#ifdef TRACERS_ON
        call write_dist_data(grid,fid,'gtracer',atmocn%gtracer, jdim=3)
#endif

        do ipatch = 1,size(asflx)
          suffix = '_'//asflx(ipatch)%surf_name
          vname = 'uabl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%uabl, jdim=3)
          vname = 'vabl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%vabl, jdim=3)
          vname = 'tabl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%tabl, jdim=3)
          vname = 'qabl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%qabl, jdim=3)
          vname = 'eabl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%eabl, jdim=3)
          vname = 'cmgs'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%cmgs)
          vname = 'chgs'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%chgs)
          vname = 'cqgs'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%cqgs)
          vname = 'ipbl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%ipbl)
          vname = 'ustar_pbl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%ustar_pbl)
          vname = 'lmonin_pbl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%lmonin_pbl)
#ifdef TRACERS_ON
          vname = 'trabl'//suffix
          call write_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%trabl, jdim=4)
#endif
        enddo

        call write_dist_data(grid,fid,'tsavg',atmsrf%tsavg)
        call write_dist_data(grid,fid,'qsavg',atmsrf%qsavg)
c        call write_dist_data(grid,fid,'dclev',atmsrf%dclev)
        call write_dist_data(grid,fid,'usavg',atmsrf%usavg)
        call write_dist_data(grid,fid,'vsavg',atmsrf%vsavg)
        call write_dist_data(grid,fid,'wsavg',atmsrf%wsavg)
        call write_dist_data(grid,fid,'tauavg',atmsrf%tauavg)
        call write_dist_data(grid,fid,'tgvavg',atmsrf%tgvavg)
        call write_dist_data(grid,fid,'qgavg',atmsrf%qgavg)
        call write_dist_data(grid,fid,'ustar_pbl',atmsrf%ustar_pbl)
        call write_dist_data(grid,fid,'lmonin_pbl',atmsrf%lmonin_pbl)

      case (ioread)             ! input from restart file
        !call read_dist_data(grid,fid,'gtemp',atmocn%gtemp)
        !call read_dist_data(grid,fid,'gtempr',atmocn%gtempr)
        call read_dist_data(grid,fid,'asst',atmocn%gtemp)
        call read_dist_data(grid,fid,'atempr',atmocn%gtempr)
        call read_dist_data(grid,fid,'sss',atmocn%sss)
        call read_dist_data(grid,fid,'ogeoza',atmocn%ogeoza)
        call read_dist_data(grid,fid,'uosurf',atmocn%uosurf)
        call read_dist_data(grid,fid,'vosurf',atmocn%vosurf)
        call read_dist_data(grid,fid,'mlhc',atmocn%mlhc)
#ifdef TRACERS_ON
        call read_dist_data(grid,fid,'gtracer',atmocn%gtracer, jdim=3)
#endif
        do ipatch = 1,size(asflx)
          suffix = '_'//asflx(ipatch)%surf_name
          vname = 'uabl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%uabl, jdim=3)
          vname = 'vabl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%vabl, jdim=3)
          vname = 'tabl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%tabl, jdim=3)
          vname = 'qabl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%qabl, jdim=3)
          vname = 'eabl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%eabl, jdim=3)
          vname = 'cmgs'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%cmgs)
          vname = 'chgs'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%chgs)
          vname = 'cqgs'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%cqgs)
          vname = 'ipbl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%ipbl)
          vname = 'ustar_pbl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%ustar_pbl)
          vname = 'lmonin_pbl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%lmonin_pbl)
#ifdef TRACERS_ON
          vname = 'trabl'//suffix
          call read_dist_data(grid, fid, trim(vname),
     &         asflx(ipatch)%trabl, jdim=4)
#endif
        enddo

        call read_dist_data(grid,fid,'tsavg',atmsrf%tsavg)
        call read_dist_data(grid,fid,'qsavg',atmsrf%qsavg)
c        call read_dist_data(grid,fid,'dclev',atmsrf%dclev)
        call read_dist_data(grid,fid,'usavg',atmsrf%usavg)
        call read_dist_data(grid,fid,'vsavg',atmsrf%vsavg)
        call read_dist_data(grid,fid,'wsavg',atmsrf%wsavg)
        call read_dist_data(grid,fid,'tauavg',atmsrf%tauavg)
        call read_dist_data(grid,fid,'tgvavg',atmsrf%tgvavg)
        call read_dist_data(grid,fid,'qgavg',atmsrf%qgavg)
        call read_dist_data(grid,fid,'ustar_pbl',atmsrf%ustar_pbl)
        call read_dist_data(grid,fid,'lmonin_pbl',atmsrf%lmonin_pbl)

      end select
      return
      end subroutine new_io_fluxes

#endif /* not STANDALONE_OCEAN */
