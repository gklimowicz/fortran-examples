R6TomaSSP585F40oQ40.R GISS ModelE Lat-Lon Atmosphere Model, transient ocn/atm OMA tracers

WARNING: NON-PRODUCTION VERSION. I.E.:
R6TomaSSP585F40oQ40 = E6TomaSSP585F40oQ40 (production version) but change for regression
                      testing. (ISTART=2, INPUTZ section, initial_GHG_setup, etc.)

E6TomaSSP585F40oQ40 = adds coupled ocean based on CMIP6 run E212Tomaf10aF40oQ40_2.R
                and set up for transition between historical and SSP run (ISTART=9)
                This example uses climate-interactive CH4 wetlands+tundra
                emissions and SSP585 (IAMC-REMIND-MAGPIE-ssp585-1-1)

E6TomaF40: Same as E6F40, with OMA tracers and computed aerosol
           indirect effect, including Shindell chemistry

Lat-lon: 2x2.5 degree horizontal resolution
F40: 40 vertical layers with standard hybrid coordinate, top at .1 mb
Atmospheric composition for year 1850
Ocean climatology prescribed from years 1876-1885, CMIP6
Uses turbulence scheme (no dry conv), grav.wave drag
Time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
Filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define STDHYB                   ! standard hybrid vertical coordinate
#define ATM_LAYERING L40         ! 40 layers, top at .1 mb
#define NEW_IO                   ! new I/O (netcdf) on
#define NEW_IO_SUBDD
#define CACHED_SUBDD
#define IRRIGATION_ON
#define SWFIX_20151201
#define NO_HDIURN                ! exclude hdiurn diagnostics
#define MODIS_LAI
#define CHECK_OCEAN                  ! needed to compile aux/file CMPE002
#define SIMPLE_MESODIFF
#define OCN_LAYERING L40_5008m
#define ODIFF_FIXES_2017
#define EXPEL_COASTAL_ICEXS
#define NEW_BCdalbsn
!---> generic tracers code start
#define TRAC_ADV_CPU             ! timing index for tracer advection on
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
!  OFF #define CALCULATE_LIGHTNING ! Calculate lightning flash rates when NOx is not needed
!  OFF #define AUTOTUNE_LIGHTNING  ! Automatically generate lightning tuning parameters (present-day only)
!<--- generic tracers code end
!---> chemistry start
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
!  OFF #define AUXILIARY_OX_RADF ! radf diags for climatology or tracer Ozone
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
#define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define DYNAMIC_BIOMASS_BURNING  ! alter biomas burning my flammability
#define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
#define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
#define ACCMIP_LIKE_DIAGS  ! adds many diags as defined by ACCMIP project
!<--- chemistry end
!---> OMA start
#define TRACERS_DUST             ! include dust tracers
#define TRACERS_DUST_Silt4       ! include 4th silt size class of dust
#define TRACERS_AEROSOLS_SEASALT ! seasalt
#define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
#define TRACERS_AEROSOLS_SOA     ! Secondary Organic Aerosols
!  OFF #define SOA_DIAGS                ! Additional diagnostics for SOA
#define TRACERS_NITRATE
#define TRACERS_HETCHEM
!<--- OMA end
#define BC_ALB                    !optional tracer BC affects snow albedo
#define CLD_AER_CDNC              !aerosol-cloud interactions
#define BLK_2MOM                  !aerosol-cloud interactions
!  OFF #define NUDGE_ON                 ! nudge the meteorology
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm144x90                           ! horizontal resolution is 144x90 -> 2x2.5deg
AtmLayering                         ! vertical resolution
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform
ORES_1Qx1 OFFT288E                  ! ocean horiz res 1.25x1deg

IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

#include "latlon_source_files"
#include "modelE4_source_files"
#include "dynamic_ocn_source_files_CMIP6"
OCN_Int_LATLON                      ! atm-ocn regrid routines
#include "tracer_shared_source_files"
#include "tracer_shindell_source_files"
#include "tracer_OMA_source_files"
TRDIAG                              ! new i/o
SUBDD
CLD_AEROSOLS_Menon_MBLK_MAT_E29q BLK_DRV ! aerosol-cloud interactions
CLD_AER_CDNC                        ! aerosol-cloud interactions wrapper
! flammability_drv flammability       ! Olga's fire model

Components:
#include "E4_components_nc"    /* without "Ent" */
tracers
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT /* needed for "Ent" only */
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
#include "IC_144x90_input_files"
! Next line for OIC for these regression tests only: I traced back E212Tomaf10aF40oQ40_2.R to the
! last parent run that used ISTART=2: E198F40oQ40.R from Larissa. Which used this OIC.
OIC=altocnbc288x180_20170717/tempsalt_dec1_288x180_phc3.0_ext.nc     ! Levitus ocean intial conditions
#include "dynamic_ocn_288x180_input_files_CMIP6_istart8or9"
TOPO=Z2HX2fromZ1QX1N.BS1.nc        ! surface fractions and topography (1 cell Bering Strait)
ICEDYN_MASKFAC=iceflowmask_144x90.nc

TDISS=altocnbc288x180_20170717/TIDAL_e_v2_1QX1.HB.nc
TDISS_N=tdiss/Jayne2009_288x180.nc
POROS=altocnbc288x180_20170717/poros.nc

RVR=RD_Fd.nc             ! river direction file
NAMERVR=RD_Fd.names.txt  ! named river outlets

#include "land144x90_input_files_SSP585"
#include "rad_input_files_SSP585"
#include "rad_144x90_input_files_CMIP6"
#include "chemistry_input_files"
#include "chemistry_144x90_input_files"
#include "dust_tracer_144x90_input_files"
#include "dry_depos_144x90_input_files"
#include "chem_emiss_144x90_input_files_SSP585"
#include "ch4_interactive_wetlands_files"
#include "aerosol_OMA_input_files_SSP585"

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
R6TomaSSP585F40oQ40 (SSP case dynamic ocean atm tracer model with OMA and Shindell chemistry)

&&PARAMETERS

#include "dynamic_ocn_params"
ocean_use_qus=1     ! Advection uses the quadratic upstream scheme
DTO=112.5
ocean_use_tdmix=1  ! tdmix scheme for meso mixing
ocean_use_gmscz=1  ! vertically variation of meso diffusivity, option 1
ocean_kvismult=2.  ! mult. factor for meso diffusivity
ocean_enhance_shallow_kmeso=1 ! stronger meso mixing in shallow water
ocean_use_tdiss=1  ! simple tidally induced diapycnal diffusivity

#include "sdragF40_params"
#include "gwdragF40_params"

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! The following two lines are only used when aerosol/radiation interactions are off
FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.3,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.625  ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 2.
use_vmp=1
radius_multiplier=1.1

PTLISO=0.        ! pressure(mb) above which radiation assumes isothermal layers
H2ObyCH4=0.      ! if =1. activates stratospheric H2O generated by CH4 without interactive chemistry
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

initial_GHG_setup = 1 ! Set to 0 after initial setup.

#include "atmCompos_transient_params_SSP"
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use abs(o3_yr) instead!
!!!!!!!!!!!!!!!!!!!!!!!
#include "aerosol_OMA_params_CMIP6"
#include "dust_params_vmp_oma"
#include "common_tracer_params_SSP"
#include "chemistry_params_CMIP6"
#include "ch4_params_with_emissions_CMIP6"
! Lines like these next two may be needed in some cases (like ISTART=8 starts) to
! prevent tracers reinitialization:
! itime_tr0=-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0
! allowSomeChemReinit = 0 ! If set to 0, disallows certain sections of chemistry
!   code from reinitializing (including stratospheric model Q( ) variable).
!   Also prevents fractional application of chemistry terms for first X-timesteps
!   from ItimeI.
! ---- for interactive wetlands -----
nn_or_zon=0     ! int dist method 1=zonal avg, 0=nearest neighbor
int_wet_dist=1  ! turn on(1)/off(0) interacive SPATIAL wetlands
ice_age=0.      ! if not 0 no wetl emis for lats poleward of +/- this in deg
ns_wet=12       ! index of CH4 source that is the wetlands (dumb, I know)
exclude_us_eu=0 ! to exclude (=1) the U.S. and E.U. from inter wetl dist
topo_lim=205.d0 ! upper limit of topographic variation for new wetlands
sat_lim=-9.d0   ! lower limit on surf air temp for new wetlants
gw_ulim=100.d0  ! upper limit on ground wetness for new wetlands
gw_llim=18.d0   ! lower limit on ground wetness for new wetlands
SW_lim=27.d0    ! lower limit on SW downward flux for new wetlands
! -----------------------------------

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! surface interaction computed NIsurf times per source time step
NRAD=5           ! radiation computed NRAD times per source time step
#include "diag_params"
! save3dAOD=1      ! needed if 3D AOD (itAOD or ictAOD) SUBDDs are on and adiurn_dust=0

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960        ! write fort.1.nc or fort.2.nc every NDISK source time step
&&END_PARAMETERS

&INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
/
