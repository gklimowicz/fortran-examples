#include "rundeck_opts.h"
#ifdef SKIP_TRACERS_RAD
#undef TRACERS_ON
#endif
      MODULE RAD_COM
!@sum  RAD_COM Model radiation arrays and parameters
!@auth Original Development Team
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : lm_req
      USE RADPAR, only : S0,nraero_aod=>NTRACE
      use AbstractOrbit_mod, only: AbstractOrbit
#ifdef TRACERS_AMP
      USE AERO_CONFIG, ONLY: NMODES
#endif
#ifdef TRACERS_TOMAS
      USE TOMAS_AEROSOL, only: icomp
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      use trdust_mod, only: nSubClays
      use tracer_com, only: ntm_dust, ntm_clay, ntm_sil1, ntm_sil2,
     &     ntm_sil3, ntm_sil4, ntm_sil5
#endif
!@var S0 solar 'constant' needs to be saved between calls to radiation
      IMPLICIT NONE
      SAVE

!@dbparam NRad:   DT_Rad      =  NRad*DTsrc
      INTEGER :: NRad = 5
!@var MODRD: if MODRD=0 do radiation, else skip
      INTEGER :: MODRD

C**** DEFAULT ORBITAL PARAMETERS FOR EARTH
C**** Note PMIP runs had specified values that do not necesarily
C**** coincide with those used as the default, or the output of ORBPAR.
C****                    OMEGT          OBLIQ        ECCEN
C**** DEFAULT (2000 AD): 282.9          23.44        0.0167
C**** PMIP CONTROL:      282.04         23.446       0.016724
C**** PMIP 6kyr BP:      180.87         24.105       0.018682
C**** PMIP LGM (21k):    294.42         22.949       0.018994
!@param OMEGT_def precession angle (degrees from vernal equinox)
      real*8, parameter :: omegt_def = 282.9d0
!@param OBLIQ_def obliquity angle  (degrees)
      real*8, parameter :: obliq_def = 23.44d0
!@param ECCN_def eccentricity
      real*8, parameter :: eccn_def  = .0167d0
!@var OMEGT,OBLIQ,ECCN actual orbital parameters used
      real*8 OMEGT,OBLIQ,ECCN

C**** Database parameters to control orbital parameter calculation
C**** Note: setting variable_orb_par=0, orb_par_year_bp=-50 (=year 2000)
C**** does not produce exactly the same as the default values.
!@dbparam variable_orb_par 1 if orbital parameters are time dependent
!@+       1 : use orb par from year "JYEAR - orb_par_year_bp"
!@+       0 : use orb par from year orb_par_year_bp (BP=before 1950)
!@+      -1 : set eccn/obliq/omegt to orb_par(1:3)
!@+    else : set eccn/obliq/omegt to defaults of orb_par
      integer :: variable_orb_par = -2
!@dbparam orb_par_year_bp = offset from model_year or 1950 (fixed case)
      integer :: orb_par_year_bp = 0
!@dbparam orb_par :: directly specifies orbital parameters
      real*8, dimension(3) :: orb_par = (/ eccn_def, obliq_def,
     *     omegt_def /)

!@var dimrad_sv dimension sum of input fields saved for radia_only runs
      INTEGER, PARAMETER :: dimrad_sv=IM*JM*(7*LM+3*LM_REQ+24)
!@var RQT Radiative equilibrium temperatures above model top
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RQT
!@var Tchg Total temperature change in adjusted forcing runs
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: Tchg
!@var SRHR(0) Solar   raditive net flux into the ground          (W/m^2)
!@var TRHR(0) Thermal raditive downward flux into ground(W/O -StB*T^4)(W/m^2)
!@*   Note: -StB*T^4 is added in SURFACE, since T varies betw. rad. calls
!@var SRHR(1->LM) Solar   raditive heating rate (W/m^2)  (short wave)
!@var TRHR(1->LM) Thermal raditive heating rate (W/m^2)  (long wave)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SRHR,TRHR
!@var TRSURF upward thermal radiation at the surface from rad step W/m2
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRSURF
!@var FSF Solar Forcing over each type (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: FSF
!@var FSRDIR Solar incident at surface, direct fraction (1)
!@var DIRVIS Direct beam solar incident at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FSRDIR,DIRVIS
!@var SRVISSURF Incident solar direct+diffuse visible at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SRVISSURF
!@var SRDN Total incident solar at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SRDN  ! saved in rsf
!@var FSRDIF diffuse visible incident solar at surface
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FSRDIF
!@var DIRNIR direct  nir     incident solar at surface
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DIRNIR
!@var DIFNIR diffuse nir     incident solar at surface
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DIFNIR
!@var srnflb_save  Net solar radiation (W/m^2)
!@var trnflb_save  Net thermal radiation (W/m^2)
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: srnflb_save,trnflb_save
!@var TAUSUMW,TAUSUMI column-sum water,ice cloud opt. depths (for diags)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: TAUSUMW,TAUSUMI
#ifdef mjo_subdd
!@var OLR_acc, OLR_cnt --  Net thermal radiation at TOA (W/m^2) for SUBDD
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: OLR_acc
      REAL*8 :: OLR_cnt = 0.d0
!@var SWHR,LWHR,SWHR_cnt,LWHR_cnt -- shortwave/longwave heating rates for SUBDD (C/d)
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: SWHR,LWHR
      REAL*8 :: SWHR_cnt = 0.d0
      REAL*8 :: LWHR_cnt = 0.d0
!@var swu_avg,swu_cnt -- upward shortwave fluxes at srf for SUBDD (C/d)
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: swu_avg
      REAL*8 :: swu_cnt = 0.d0
#endif
#ifdef TRACERS_ON

!@var DIAG_FC Controls the number of radiation calls for the calculation of
!@+           aerosol radiative forcing. One call if =1, multiple calls if
!@+           =2, with their number depending on the aerosol scheme used.
!@+           Use =2 sparingly, it is s l o w. Default is 1. No calls if zero.
      integer :: diag_fc=1
! nraero_xxxx are the aerosol-specific nraero_aod (old ntrace) components of
! aerosol-active species in radiation. nraero_aod=sum(nraero_xxxx)
!@var nraero_aod Number of aerosol types in optical depth calculations
!@var nraero_rf Number of aerosol types in forcing calculations, which is
!@+             different from nraero_aod when DIAG_FC=1 (default)
      integer :: nraero_rf=0
#ifdef TRACERS_AEROSOLS_Koch
#ifdef SULF_ONLY_AEROSOLS
      integer, parameter :: nraero_koch=1
#else
#ifdef TRACERS_AEROSOLS_VBS
#ifdef TRACERS_AEROSOLS_SOA
      integer, parameter :: nraero_koch=5
#else
      integer, parameter :: nraero_koch=4
#endif  /* TRACERS_AEROSOLS_SOA */
#else
#ifdef TRACERS_AEROSOLS_SOA
      integer, parameter :: nraero_koch=6
#else
      integer, parameter :: nraero_koch=5
#endif  /* TRACERS_AEROSOLS_SOA */
#endif  /* TRACERS_AEROSOLS_VBS */
#endif  /* SULF_ONLY_AEROSOLS */
#else
      integer, parameter :: nraero_koch=0
#endif  /* TRACERS_AEROSOLS_Koch */

#ifdef TRACERS_NITRATE
      integer, parameter :: nraero_nitrate=1
#else
      integer, parameter :: nraero_nitrate=0
#endif  /* TRACERS_NITRATE */

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      integer, parameter :: nraero_clay = nSubClays * ntm_clay
      integer, parameter :: nraero_dust = nraero_clay + ntm_sil1 +
     &     ntm_sil2 + ntm_sil3 + ntm_sil4 + ntm_sil5
!@var nr_soildust First index of dust tracers in radiation (nraero_aod)
      integer :: nr_soildust = 0
#else
      integer, parameter :: nraero_clay = 0
      integer, parameter :: nraero_dust = 0
#endif  /* TRACERS_DUST */

!@var nraero_OMA Number of OMA tracers that have an AOD value
!@var nraero_AMP Number of AMP tracers that have an AOD value
!@var nraero_TOMAS Number of TOMAS tracers that have an AOD value
      integer :: nraero_OMA=0
      integer :: nraero_AMP=0
      integer :: nraero_TOMAS=0

#ifdef TRACERS_AEROSOLS_SEASALT
      integer, parameter :: nraero_seasalt=2
#else
      integer, parameter :: nraero_seasalt=0
#endif  /* TRACERS_AEROSOLS_SEASALT */

#ifdef TRACERS_ON
!@var njaero max expected rad code tracers passed to photolysis
!@var nraero_aod_rsf value of nraero_aod found in the rsf file
!@var nraero_rf_rsf value of nraero_rf found in the rsf file
!@var save_dry_aod_rsf value of save_dry_aod found in the rsf file
!@var tau_as All-sky aerosol optical saved 1:nraero_aod not 1:ntm
!@+   This is so clays are separate. Now also used for old parameter
!@+   mxfastj: Number of aerosol/cloud types currently active in the model
!@var tau_cs Same as tau_as for clear-sky
!@var tau_dry Same as tau_as for dry aerosol (RH=0%)
      integer :: njaero ! nraero_aod+2 cloud types (water/ice)
      integer :: nraero_aod_rsf=0
      integer :: nraero_rf_rsf=0
      integer :: save_dry_aod_rsf=0
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: tau_as
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: tau_cs
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: tau_dry
#ifdef CACHED_SUBDD
!@var abstau_as Same as tau_as for absorption
!@var abstau_cs Same as tau_cs for absorption
!@var abstau_dry Same as tau_dry for absorption
!@var swfrc Shortwave aerosol radiative forcing
!@var lwfrc Shortwave aerosol radiative forcing
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: abstau_as
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: abstau_cs
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: abstau_dry
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: swfrc,lwfrc
#endif  /* CACHED_SUBDD */
#endif
#endif
!@var CFRAC Total cloud fraction as seen be radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: CFRAC ! saved in rsf
!@var RCLD Total cloud optical depth as seen be radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RCLD ! saved in rsf
!@var chem_tracer_save 3D O3, CH4 saved elsewhere for use in radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: chem_tracer_save!saved rsf
#ifdef GCC_COUPLE_RAD
!@var GCCco2_tracer_save 3D CO2 saved elsewhere for use in radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GCCco2_tracer_save!saved rsf
#endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
!@var stratO3_tracer_save 3D stratOx saved elsewhere for use in rad code
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::stratO3_tracer_save!saved rsf
#endif
!@var rad_to_chem save 3D quantities from radiation code for use in
!@+   chemistry (or rest of model). 1=Ozone, 2=aerosol ext, 3=N2O, 4=CH4,
!@+   5=CFC11+CFC12 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: rad_to_chem !saved in rsf
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: rad_to_file
#ifdef GCC_COUPLE_RAD
!@var GCCco2rad_to_chem save 3D quantities from radiation code 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GCCco2rad_to_chem !saved in rsf
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GCCco2rad_to_file
#endif
!@var KLIQ Flag indicating dry(0)/wet(1) atmosphere (memory feature)
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: KLIQ ! saved in rsf
!@dbparam Ikliq 0,1,-1 initialize kliq as dry,equil,current model state
      INTEGER :: Ikliq = -1  !  get kliq-array from restart file
!@dbparam RHfix const.rel.humidity passed to radiation for aeros. tests
      REAL*8 :: RHfix = -1.  !  pass the current model rel.humidity
!@dbparam dalbsnX global coeff for snow alb change by black carbon depos
      REAL*8 ::  dalbsnX = 0.
!@dbparam albsn_yr year of blk carb depos used for snow alb. reduction
      INTEGER ::  albsn_yr = 1951

!     variables related to aerosol indirect effects:
!     (CDNC=cloud droplet number concentration)
!@dbparam CC_CDNCx scaling factor relating cld cvr change and CDNC change
      REAL*8 :: CC_CDNCX = .0000d0  ! .0036d0
!@dbparam OC_CDNCx scaling factor relating cld opt depth and CDNC change
      REAL*8 :: OD_CDNCX = .0000d0  ! .007d0
!@var pcdnc,vcdnc pressure,vertical profile for cld.cvr change
      real*8, parameter, dimension(7) ::
     * pcdnc=(/984.d0, 964.d0, 934.d0, 884.d0, 810.d0, 710.d0, 550.d0/)
     *,vcdnc=(/ .35d0,  .20d0,  .10d0,  .17d0,  .10d0,  .08d0,   0.d0/)
!@var cdncl = vcdnc interpolated to current vertical resolution
      real*8 cdncl(LM)

!@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: COSZ1
!@var COSZ_day Mean Solar Zenith angle for current day
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: COSZ_day
!@var SUNSET Time of sunset for current day (radians from local noon)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SUNSET
!@dbparam S0X solar constant multiplication factor
      REAL*8 :: S0X = 1.
!@dbparam S0_yr,S0_day obs.date of solar constant (if 0: time var)
      INTEGER :: S0_yr = 1951 , S0_day = 182
!@dbparam CO2X,... scaling factors for CO2 N2O CH4 CFC11 CFC12 XGHG
      REAL*8 :: CO2X=1.,N2OX=1.,CH4X=1., CFC11X=1.,CFC12X=1.,XGHGX=1.
     *         ,O2X=1.,NO2X=1.,N2CX=1.,YGHGX=2.,SO2X=0.
     *         ,CH4X_RADoverCHEM=1.d0
!@dbparm ref_mult factor to control REFDRY from rundeck
      REAL*8 :: ref_mult = 1.
!@dbparam GHG_yr,GHG_day obs.date of well-mixed GHgases (if 0: time var)
      INTEGER :: GHG_yr = 1951 , GHG_day = 182
!@dbparam Volc_yr,Volc_day obs.date of Volc.Aerosols (if 0: time var)
!@+   special cases: Volc_yr=-1    : 150-yr mean 1850-1999
!@+                  Volc_yr=-2010 : current year up to 2010 then
!@+                                  repeat volcanos from 100 yrs ago
!@+                  Volc_yr=-2000 : older way of creating future volc
      INTEGER :: Volc_yr = 1951 , Volc_day = 182
!@dbparam Aero_yr obs.year of troposph.Aerosols (if 0: use current yr)
      INTEGER :: Aero_yr = 1951    ! always use annual cycle
!@dbparam dust_yr nominal year for prescribed dust climatology (if 0: use current yr)
      INTEGER :: dust_yr = 1951    ! always use annual cycle
!@dbparam O3_yr obs.year of Ozone (if 0: use current year)
      INTEGER :: O3_yr = 1951      ! always use annual cycle
!@dbparam H2OstratX strat_water_vapor, cloud, Ozone scaling factor
      REAL*8 :: H2OstratX = 1. , cldX = 1. , O3X = 1.
!@dbparam H2ObyCH4 if not 0: add CH4 produced H2O into layers 1->LM
      REAL*8 :: H2ObyCH4 = 1.
!@var dH2O  zonal H2O-prod.rate in kg/m^2/ppm_CH4/second in layer L
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: dH2O
!@var RSDIST,SIND,COSD orbit related variables computed once a day
      REAL*8 :: RSDIST,SIND,COSD
!@var ALB is SRNFLB(1)/(SRDFLB(1)+1.D-20),PLAVIS,PLANIR,ALBVIS,ALBNIR,
!@+       SRRVIS,SRRNIR,SRAVIS,SRANIR (see RADIATION)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: ALB

!@var SALB (1.-broadband surface albedo) - saved in rsf
      REAL*8, POINTER, DIMENSION(:,:) :: SALB   ! = ALB(:,:,1)
!      EQUIVALENCE (SALB,ALB)

#ifdef ALTER_RADF_BY_LAT
!@var FULGAS_lat multiplicative factors for altering FULGAS by latitude
!@+ for non-transient runs. (greenhouse gas regional forcing)
      real*8, dimension(13,46):: FULGAS_lat !rad not model grid, 13 gasses
!@var FS8OPX_lat multiplicative factors for altering FS8OPX by latitude
!@+ for non-transient runs. (aerosol regional forcing) SOLAR
!@var FT8OPX_lat multiplicative factors for altering FT8OPX by latitude
!@+ for non-transient runs. (aerosol regional forcing) THERMAL
      real*8, dimension(8,46):: FS8OPX_lat,FT8OPX_lat !rad not model grid
                                                !8 groups of aerosols
#endif
!@dbparam rad_interact_aer =1 for radiatively active non-chem tracers
      INTEGER :: rad_interact_aer = 0  ! defaults to 0
!@dbparam clim_interact_chem=1 for radiatively active chem tracers
!@+ also affects chemisty to humidity feedback
      INTEGER :: clim_interact_chem= 0  ! defaults to 0

!@dbparam nradfrc sets frequency of inst. rad. forcing calculations
      INTEGER :: nradfrc=1 ! do them every nrad*nradfrc physics steps
!                nradfrc=0: skip all, no repeated radiation calculations
C**** the radiative forcing level for instantaneous forcing calcs is set
C**** using the rad_forc_lev parameter.
!@dbparam rad_forc_lev = 0 for TOA, 1 for LTROPO (default=0)
      INTEGER :: rad_forc_lev = 0

!@dbparam cloud_rad_forc = 1 for calculation of cloud radiative forcing
      INTEGER :: cloud_rad_forc = 0

!@dbparam TAero_aod_diag = 1 outputs offline aerosol optical properties,
!@+   = 2 outputs band 6 only. Note this only works for background aerosols,
!@+   not tracers.
      INTEGER :: TAero_aod_diag = 0
!@dbparam aer_rad_forc = 1 for calculation of aerosol radiative forcing
!@+   note this only works for background aerosols, not tracers
      INTEGER :: aer_rad_forc = 0

!@var co2ppm Current CO2 level as seen by radiation
      REAL*8 :: co2ppm = 280.    ! set a reasonable default value
      
C**** Local variables initialised in init_RAD
!@var PLB0,QL0 global parts of local arrays (to avoid OMP-copyin)
      REAL*8, DIMENSION(LM_REQ)       :: PLB0,SHL0
!@var ntrix_aod Indexing array for aerosol optical depth tracer names
!@var ntrix_rf Indexing array for aerosol radiative forcing tracer names
      INTEGER, allocatable, DIMENSION(:) :: ntrix_aod,ntrix_rf
!@var WTTR weighting array for optional aerosol-ratiation interactions
      REAL*8, allocatable, DIMENSION(:) :: WTTR

#ifdef CUBED_SPHERE
!@var JM_DH2O number of latitudes in CH4->H2O input file
!@var LAT_DH2O latitudes in CH4->H2O input file (converted to radians)
      integer, parameter :: jm_dh2o=18
      real*8 :: lat_dh2o(jm_dh2o)
#endif

!@dbparam snoage_def determines how snowage is calculated:
!@+       = 0     independent of temperature
!@+       = 1     only when max daily local temp. over type > 0
      integer :: snoage_def = 0
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SNOAGE
      class (AbstractOrbit), allocatable :: orbit

!@dbparam chl_from_obio =1 to use chl from obio when computing ocean albedo
      INTEGER :: chl_from_obio = 0
!@dbparam chl_from_seawifs =1 to use chl from SeaWIFs when computing ocn albedo
      INTEGER :: chl_from_seawifs = 0

      contains

      subroutine radiationSetOrbit(anOrbit)
      class (AbstractOrbit), intent(in) :: anOrbit
      allocate(orbit, source=anOrbit)
      end subroutine radiationSetOrbit

      END MODULE RAD_COM

      SUBROUTINE ALLOC_RAD_COM(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Rodger Abel

      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE DOMAIN_DECOMP_ATM, ONLY : getDomainBounds
      USE RESOLUTION, ONLY : IM, JM, LM
      USE ATM_COM, ONLY : LM_REQ
#ifdef TRACERS_ON
      USE tracer_com,ONLY : NTM
#endif
      USE RAD_COM, ONLY : RQT,Tchg,SRHR,TRHR,FSF,FSRDIR,SRVISSURF,TRSURF
     *     ,SRDN, CFRAC, RCLD, chem_tracer_save,rad_to_chem,rad_to_file
     *     ,KLIQ, COSZ1, COSZ_day, SUNSET, dH2O, ALB, SALB, SNOAGE
     *     ,srnflb_save, trnflb_save
     *     ,FSRDIF,DIRNIR,DIFNIR,TAUSUMW,TAUSUMI,DIRVIS
#ifdef GCC_COUPLE_RAD
     *     ,GCCco2_tracer_save,GCCco2rad_to_chem,GCCco2rad_to_file
#endif
#ifdef mjo_subdd
     *     ,SWHR_cnt,LWHR_cnt,SWHR,LWHR,OLR_acc,OLR_cnt
     *     ,swu_avg,swu_cnt
#endif
#ifdef CUBED_SPHERE
     &     ,JM_DH2O
#endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &     ,stratO3_tracer_save
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_0H, J_1H
      INTEGER :: IER

      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE( RQT(LM_REQ, I_0H:I_1H, J_0H:J_1H),
     *     Tchg(LM+LM_REQ, I_0H:I_1H, J_0H:J_1H),
     *     SRHR(0:LM, I_0H:I_1H, J_0H:J_1H),
     *     TRHR(0:LM, I_0H:I_1H, J_0H:J_1H),
     *     TRSURF(4, I_0H:I_1H, J_0H:J_1H),
     *     FSF(4,I_0H:I_1H, J_0H:J_1H),
     *     FSRDIR(I_0H:I_1H, J_0H:J_1H),
     *     DIRVIS(I_0H:I_1H, J_0H:J_1H),
     *     SRVISSURF(I_0H:I_1H, J_0H:J_1H),
     *     FSRDIF(I_0H:I_1H, J_0H:J_1H),
     *     DIRNIR(I_0H:I_1H, J_0H:J_1H),
     *     DIFNIR(I_0H:I_1H, J_0H:J_1H),
     *     TAUSUMW(I_0H:I_1H, J_0H:J_1H),
     *     TAUSUMI(I_0H:I_1H, J_0H:J_1H),
     *     SRDN(I_0H:I_1H, J_0H:J_1H),
     *     CFRAC(I_0H:I_1H, J_0H:J_1H),
     *     RCLD(LM, I_0H:I_1H, J_0H:J_1H),
#ifdef GCC_COUPLE_RAD
     *     GCCco2_tracer_save(LM, I_0H:I_1H, J_0H:J_1H),
     *     GCCco2rad_to_chem(LM, I_0H:I_1H, J_0H:J_1H),
     *     GCCco2rad_to_file(LM, I_0H:I_1H, J_0H:J_1H),
#endif
     *     chem_tracer_save(2,LM, I_0H:I_1H, J_0H:J_1H),
     *     rad_to_chem(5, LM, I_0H:I_1H, J_0H:J_1H),
     *     rad_to_file(5, LM, I_0H:I_1H, J_0H:J_1H),
     *     SNOAGE(3,I_0H:I_1H,J_0H:J_1H),
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     *     stratO3_tracer_save(LM, I_0H:I_1H, J_0H:J_1H),
#endif
     *     KLIQ(LM,4, I_0H:I_1H, J_0H:J_1H),
     *     COSZ1   (I_0H:I_1H, J_0H:J_1H),
     *     COSZ_day(I_0H:I_1H, J_0H:J_1H),
     *     SUNSET  (I_0H:I_1H, J_0H:J_1H),
#ifdef CUBED_SPHERE
     *     dH2O(JM_DH2O, LM, 12),
#else
     *     dH2O(J_0H:J_1H, LM, 12),
#endif
     *     ALB(I_0H:I_1H, J_0H:J_1H, 9),
     &     srnflb_save(I_0H:I_1H,J_0H:J_1H,Lm),
     &     trnflb_save(I_0H:I_1H,J_0H:J_1H,Lm),
#ifdef mjo_subdd
     *     OLR_acc(I_0H:I_1H,J_0H:J_1H),
     *     SWHR(I_0H:I_1H,J_0H:J_1H,Lm),
     *     LWHR(I_0H:I_1H,J_0H:J_1H,Lm),
     *     swu_avg(I_0H:I_1H,J_0H:J_1H),
#endif
#ifdef TRACERS_SPECIAL_Shindell
#endif
     *     STAT=IER)

#ifdef mjo_subdd
      OLR_acc=0.
      OLR_cnt=0.
      SWHR=0.
      LWHR=0.
      SWHR_cnt=0.
      LWHR_cnt=0.
      swu_avg=0.
      swu_cnt=0.
#endif
      KLIQ = 1
      dH2O = 0.
      SALB => ALB(:,:,1)
      SRVISSURF = 0
      FSF=0
      TRSURF=0
      RETURN
      END SUBROUTINE ALLOC_RAD_COM

      subroutine def_rsf_rad(fid)
!@sum  def_rsf_rad defines radiation array structure in restart files
!@auth M. Kelley
!@ver  beta
      use rad_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
#ifdef TRACERS_ON
      use TRDIAG_COM, only: save_dry_aod
#endif
      implicit none
      integer fid   !@var fid file id

      call defvar(grid,fid,s0,'s0')
      call defvar(grid,fid,rqt,'rqt(lm_req,dist_im,dist_jm)')
      call defvar(grid,fid,kliq,'kliq(lm,four,dist_im,dist_jm)')
      call defvar(grid,fid,srhr,'srhr(zero_to_lm,dist_im,dist_jm)')
      call defvar(grid,fid,trhr,'trhr(zero_to_lm,dist_im,dist_jm)')
      call defvar(grid,fid,trsurf,'trsurf(nstype,dist_im,dist_jm)')
      call defvar(grid,fid,fsf,'fsf(nstype,dist_im,dist_jm)')
      call defvar(grid,fid,fsrdir,'fsrdir(dist_im,dist_jm)')
      call defvar(grid,fid,srvissurf,'srvissurf(dist_im,dist_jm)')
      call defvar(grid,fid,srdn,'srdn(dist_im,dist_jm)')
      call defvar(grid,fid,cfrac,'cfrac(dist_im,dist_jm)')
      call defvar(grid,fid,salb,'salb(dist_im,dist_jm)')
      call defvar(grid,fid,fsrdif,'fsrdif(dist_im,dist_jm)')
      call defvar(grid,fid,dirnir,'dirnir(dist_im,dist_jm)')
      call defvar(grid,fid,difnir,'difnir(dist_im,dist_jm)')
      call defvar(grid,fid,rcld,'rcld(lm,dist_im,dist_jm)')
      call defvar(grid,fid,snoage,'snoage(d3,dist_im,dist_jm)')

#ifdef TRACERS_ON
      if (nraero_aod > 0) then
        call defvar(grid,fid,nraero_aod,'nraero_aod')
        call defvar(grid,fid,save_dry_aod,'save_dry_aod')
        call defvar(grid,fid,tau_as,
     &       'tau_as(dist_im,dist_jm,lm,nraero_aod)')
        call defvar(grid,fid,tau_cs,
     &       'tau_cs(dist_im,dist_jm,lm,nraero_aod)')
        if (save_dry_aod>0) then
          call defvar(grid,fid,tau_dry,
     &         'tau_dry(dist_im,dist_jm,lm,nraero_aod)')
        endif
#ifdef CACHED_SUBDD
        call defvar(grid,fid,abstau_as,
     &       'abstau_as(dist_im,dist_jm,lm,nraero_aod)')
        call defvar(grid,fid,abstau_cs,
     &       'abstau_cs(dist_im,dist_jm,lm,nraero_aod)')
        if (save_dry_aod>0) then
          call defvar(grid,fid,abstau_dry,
     &         'abstau_dry(dist_im,dist_jm,lm,nraero_aod)')
        endif
        call defvar(grid,fid,nraero_rf,'nraero_rf')
        if (nraero_rf>0) then
          call defvar(grid,fid,swfrc,'swfrc(dist_im,dist_jm,nraero_rf)')
          call defvar(grid,fid,lwfrc,'lwfrc(dist_im,dist_jm,nraero_rf)')
        endif
#endif  /* CACHED_SUBDD */
      endif
#ifdef GCC_COUPLE_RAD
      call defvar(grid,fid,GCCco2_tracer_save,
     &     'GCCco2_tracer_save(lm,dist_im,dist_jm)')
      call defvar(grid,fid,GCCco2rad_to_chem,
     &     'GCCco2rad_to_chem(lm,dist_im,dist_jm)')
#endif

#if (defined TRACERS_SPECIAL_Shindell) 
      call defvar(grid,fid,chem_tracer_save,
     &     'chem_tracer_save(two,lm,dist_im,dist_jm)')
      call defvar(grid,fid,rad_to_chem,
     &     'rad_to_chem(five,lm,dist_im,dist_jm)')
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      call defvar(grid,fid,strato3_tracer_save,
     &     'strato3_tracer_save(lm,dist_im,dist_jm)')
#endif
#endif
#ifdef TRACERS_DUST
      call defvar(grid,fid,srnflb_save,
     &     'srnflb_save(dist_im,dist_jm,lm)')
      call defvar(grid,fid,trnflb_save,
     &     'trnflb_save(dist_im,dist_jm,lm)')
#endif
#endif  /* TRACERS_ON */
      return
      end subroutine def_rsf_rad

      subroutine new_io_rad(fid,iaction)
!@sum  new_io_rad read/write radiation arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
#ifdef TRACERS_ON
      USE tracer_com , only : NTM
      use TRDIAG_COM, only: save_dry_aod
#endif
      use rad_com
      use domain_decomp_atm, only : grid, getDomainBounds
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file

      integer :: I_0H, I_1H
      integer :: J_0H, J_1H

      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_data(grid, fid,'s0', s0)
        call write_dist_data(grid, fid,'rqt',  rqt, jdim=3)
        call write_dist_data(grid, fid,'kliq', kliq, jdim=4)
        call write_dist_data(grid, fid,'srhr', srhr, jdim=3)
        call write_dist_data(grid, fid,'trhr', trhr, jdim=3)
        call write_dist_data(grid, fid,'trsurf', trsurf, jdim=3)
        call write_dist_data(grid, fid,'fsf',  fsf, jdim=3)
        call write_dist_data(grid, fid,'salb', salb)
        call write_dist_data(grid, fid,'fsrdir', fsrdir)
        call write_dist_data(grid, fid,'srvissurf', srvissurf)
        call write_dist_data(grid, fid,'fsrdif', fsrdif)
        call write_dist_data(grid, fid,'dirnir', dirnir)
        call write_dist_data(grid, fid,'difnir', difnir)
        call write_dist_data(grid, fid,'srdn',   srdn)
        call write_dist_data(grid, fid,'cfrac',  cfrac)
        call write_dist_data(grid, fid,'rcld', rcld, jdim=3)
        call write_dist_data(grid, fid,'snoage', snoage,jdim=3)
#if (defined GCC_COUPLE_RAD)
        call write_dist_data(grid,fid,
     &       'GCCco2_tracer_save', GCCco2_tracer_save, jdim=3)
        call write_dist_data(grid,fid,
     &       'GCCco2rad_to_chem', GCCco2rad_to_chem, jdim=3)
#endif
#if (defined TRACERS_SPECIAL_Shindell)
        call write_dist_data(grid,fid,
     &       'chem_tracer_save', chem_tracer_save, jdim=4)
        call write_dist_data(grid,fid,
     &       'rad_to_chem', rad_to_chem, jdim=4)
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        call write_dist_data(grid,fid,
     &       'strato3_tracer_save', strato3_tracer_save, jdim=3)
#endif
#endif
#ifdef TRACERS_DUST
        call write_dist_data(grid,fid,'srnflb_save',srnflb_save)
        call write_dist_data(grid,fid,'trnflb_save',trnflb_save)
#endif
#ifdef TRACERS_ON
        if (nraero_aod > 0) then
          call write_data(grid, fid,'nraero_aod', nraero_aod)
          call write_data(grid, fid,'save_dry_aod', save_dry_aod)
          call write_dist_data(grid,fid,'tau_as',tau_as)
          call write_dist_data(grid,fid,'tau_cs',tau_cs)
          if (save_dry_aod>0) then
            call write_dist_data(grid,fid,'tau_dry',tau_dry)
          endif
#ifdef CACHED_SUBDD
          call write_dist_data(grid,fid,'abstau_as',abstau_as)
          call write_dist_data(grid,fid,'abstau_cs',abstau_cs)
          if (save_dry_aod>0) then
            call write_dist_data(grid,fid,'abstau_dry',abstau_dry)
          endif
          call write_data(grid, fid,'nraero_rf', nraero_rf)
          if (nraero_rf>0) then
            call write_dist_data(grid,fid,'swfrc',swfrc)
            call write_dist_data(grid,fid,'lwfrc',lwfrc)
          endif
#endif  /* CACHED_SUBDD */
        endif
#endif  /* TRACERS_ON */
      case (ioread)
        call read_data(grid, fid,'s0', s0, bcast_all=.true.)
        call read_dist_data(grid, fid,'rqt',  rqt, jdim=3)
        call read_dist_data(grid, fid,'kliq', kliq, jdim=4)
        call read_dist_data(grid, fid,'srhr', srhr, jdim=3)
        call read_dist_data(grid, fid,'trhr', trhr, jdim=3)
        call read_dist_data(grid, fid,'trsurf', trsurf, jdim=3)
        call read_dist_data(grid, fid,'fsf',  fsf, jdim=3)
        call read_dist_data(grid, fid,'salb', salb)
        fsrdir = 0.; srvissurf = 0.
        call read_dist_data(grid, fid,'fsrdir', fsrdir)
        call read_dist_data(grid, fid,'srvissurf', srvissurf)
        dirvis = fsrdir*srvissurf ! reconstruct when restarting.
        call read_dist_data(grid, fid,'fsrdif', fsrdif)
        call read_dist_data(grid, fid,'dirnir', dirnir)
        call read_dist_data(grid, fid,'difnir', difnir)
        call read_dist_data(grid, fid,'srdn',   srdn)
        call read_dist_data(grid, fid,'cfrac',  cfrac)
        call read_dist_data(grid, fid,'rcld', rcld, jdim=3)
        call read_dist_data(grid, fid,'snoage', snoage,jdim=3)
#ifdef GCC_COUPLE_RAD
        call read_dist_data(grid,fid,
     &       'GCCco2_tracer_save', GCCco2_tracer_save, jdim=3)
        call read_dist_data(grid,fid,
     &       'GCCco2rad_to_chem', GCCco2rad_to_chem, jdim=3)
#endif
#if (defined TRACERS_SPECIAL_Shindell) 
        call read_dist_data(grid,fid,
     &       'chem_tracer_save', chem_tracer_save, jdim=4)
        call read_dist_data(grid,fid,
     &       'rad_to_chem', rad_to_chem, jdim=4)
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        call read_dist_data(grid,fid,
     &       'strato3_tracer_save', strato3_tracer_save, jdim=3)
#endif
#endif
#ifdef TRACERS_DUST
        call read_dist_data(grid,fid,'srnflb_save',srnflb_save)
        call read_dist_data(grid,fid,'trnflb_save',trnflb_save)
#endif
#ifdef TRACERS_ON
        if (.not.allocated(tau_as)) then
          call read_data(grid,fid,'nraero_aod',nraero_aod_rsf,
     &                   bcast_all=.true.)
          call read_data(grid,fid,'save_dry_aod',save_dry_aod_rsf,
     &                   bcast_all=.true.)
          if (nraero_aod_rsf /= 0) then
            allocate(tau_as(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod_rsf))
            allocate(tau_cs(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod_rsf))
            if (save_dry_aod_rsf>0) then
              allocate(tau_dry(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod_rsf))
            endif
#ifdef CACHED_SUBDD
            allocate(abstau_as(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod_rsf))
            allocate(abstau_cs(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod_rsf))
            if (save_dry_aod_rsf>0) then
             allocate(abstau_dry(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod_rsf))
            endif
            call read_data(grid,fid,'nraero_rf',nraero_rf_rsf,
     &                     bcast_all=.true.)
            if (nraero_rf_rsf>0) then
              allocate(swfrc(I_0H:I_1H,J_0H:J_1H,nraero_rf_rsf))
              allocate(lwfrc(I_0H:I_1H,J_0H:J_1H,nraero_rf_rsf))
            endif
#endif  /* CACHED_SUBDD */
          endif
        endif
        if (allocated(tau_as)) then ! needs to be separate from previous if
          call read_dist_data(grid,fid,'tau_as',tau_as)
          call read_dist_data(grid,fid,'tau_cs',tau_cs)
          if (save_dry_aod_rsf>0) then
            call read_dist_data(grid,fid,'tau_dry',tau_dry)
          endif
#ifdef CACHED_SUBDD
          call read_dist_data(grid,fid,'abstau_as',abstau_as)
          call read_dist_data(grid,fid,'abstau_cs',abstau_cs)
          if (save_dry_aod_rsf>0) then
            call read_dist_data(grid,fid,'abstau_dry',abstau_dry)
          endif
          if (nraero_rf_rsf>0) then
            call read_dist_data(grid,fid,'swfrc',swfrc)
            call read_dist_data(grid,fid,'lwfrc',lwfrc)
          endif
#endif  /* CACHED_SUBDD */
        endif
#endif  /* TRACERS_ON */
      end select
      return
      end subroutine new_io_rad

      subroutine read_rad_ic
!@sum   read_rad_ic read radiation coldstart initial conditions file.
      use rad_com, only : snoage
      use domain_decomp_atm, only : grid
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      implicit none
      integer fid   !@var fid unit number of read/write

      if(file_exists('GIC')) then
        ! Read snow age using old-style IC (from rsf)
        fid = par_open(grid,'GIC','read')
        call read_dist_data(grid, fid, 'snoage', snoage,jdim=3)
        call par_close(grid,fid)
      else
        ! Newer cold-start IC files contain only the fundamental state variables.
        ! Set snow age to zero (Initial snow albedo irrelevant for cold starts).
        snoage = 0d0
      endif
      return
      end subroutine read_rad_ic

      MODULE DIAG_COM_RAD
      implicit none

      integer ::
     &      j_h2och4=1
     &     ,j_pcldss=1
     &     ,j_pcldmc=1
     &     ,j_clddep=1
     &     ,j_pcld=1
     &     ,j_srincp0=1
     &     ,j_srnfp0=1
     &     ,j_srnfp1=1
     &     ,j_srincg=1
     &     ,j_srnfg=1
     &     ,j_brtemp=1
     &     ,j_trincg=1
     &     ,j_hsurf=1
     &     ,j_hatm=1
     &     ,j_plavis=1
     &     ,j_planir=1
     &     ,j_albvis=1
     &     ,j_albnir=1
     &     ,j_srrvis=1
     &     ,j_srrnir=1
     &     ,j_sravis=1
     &     ,j_sranir=1
     &     ,j_trnfp0=1
     &     ,j_trnfp1=1
     &     ,j_clrtoa=1
     &     ,j_clrtrp=1
     &     ,j_tottrp=1
#ifdef HEALY_LM_DIAGS
     *     ,j_vtau=1
     *     ,j_ghg=1
#endif

      integer ::
     &      jl_srhr=1
     &     ,jl_trcr=1
     &     ,jl_totcld=1
     &     ,jl_sscld=1
     &     ,jl_mccld=1
     &     ,jl_wcld=1
     &     ,jl_icld=1
     &     ,jl_wcod=1
     &     ,jl_icod=1
     &     ,jl_wcsiz=1
     &     ,jl_icsiz=1
     &     ,jl_wcldwt=1
     &     ,jl_icldwt=1

      integer ::
     &      ij_pmccld=1
     &     ,ij_trnfp0=1
     &     ,ij_cldcv=1
     &     ,ij_pcldl=1
     &     ,ij_pcldm=1
     &     ,ij_pcldh=1
     &     ,ij_pcldl_ss=1
     &     ,ij_cldtppr=1
     &     ,ij_srvis=1
     &     ,ij_rnfp1=1
     &     ,ij_srnfp0=1
     &     ,ij_srincp0=1
     &     ,ij_srnfg=1
     &     ,ij_srincg=1
     &     ,ij_btmpw=1
     &     ,ij_srref=1
     &     ,ij_frmp=1
     &     ,ij_clr_srincg=1
     &     ,ij_CLDTPT=1
     &     ,ij_cldt1t=1
     &     ,ij_cldt1p=1
     &     ,ij_cldcv1=1
     &     ,ij_wtrcld=1
     &     ,ij_icecld=1
     &     ,ij_optdw=1
     &     ,ij_optdi=1
     &     ,ij_swcrf=1
     &     ,ij_lwcrf=1
     &     ,ij_srntp=1
     &     ,ij_trntp=1
     &     ,ij_clr_srntp=1
     &     ,ij_clr_trntp=1
     &     ,ij_clr_srnfg=1
     &     ,ij_clr_trdng=1
     &     ,ij_clr_sruptoa=1
     &     ,ij_clr_truptoa=1
     &     ,ij_swdcls=1
     &     ,ij_swncls=1
     &     ,ij_lwdcls=1
     &     ,ij_swnclt=1
     &     ,ij_lwnclt=1
     &     ,ij_srvdir=1
     &     ,ij_srvissurf=1
     &     ,ij_chl=-1
     &     ,ij_swaerrf=1
     &     ,ij_lwaerrf=1
     &     ,ij_swaersrf=1
     &     ,ij_lwaersrf=1
     &     ,ij_swaerrfnt=1
     &     ,ij_lwaerrfnt=1
     &     ,ij_swaersrfnt=1
     &     ,ij_lwaersrfnt=1
     &     ,ij_swcrf2=1
     &     ,ij_lwcrf2=1
     &     ,ij_siswd=1
     &     ,ij_siswu=1
     &     ,ij_lwprad=1
     &     ,ij_iwprad=1
     &     ,ij_h2och4 = 1
     &     ,ij_sw_cs_noa=1
     &     ,ij_lw_cs_noa=1
     &     ,ij_sw_as_noa=1
     &     ,ij_lw_as_noa=1

#ifdef ACCMIP_LIKE_DIAGS
!@var IJ_fcghg GHG forcing diagnostics (2=LW,SW, 4=CH4,N2O,CFC11,CFC12)
      integer, dimension(2,4) :: ij_fcghg
#endif
      
      integer ::
     &      ijl_rc=1
     &     ,ijl_cf=1
     &     ,ijl_QLrad=1
     &     ,ijl_QIrad=1
     &     ,ijl_wtrtau=1
     &     ,ijl_icetau=1

      integer ::
     &      idd_cl7=1
     &     ,idd_ccv=1
     &     ,idd_isw=1
     &     ,idd_palb=1
     &     ,idd_galb=1
     &     ,idd_aot=1
     &     ,idd_aot2=1
     &     ,idd_absa=1


      END MODULE DIAG_COM_RAD
