E6TmatrixF40clim2000.R GISS ModelE Lat-Lon Atmosphere Model, climatological ocn/atm MATRIX tracers

E6TmatrixF40clim2000: E6TmatrixF40clim but year 2000 settings
               (e.g. 9-year averages centered around nominal date)

Lat-lon: 2x2.5 degree horizontal resolution
F40: 40 vertical layers with standard hybrid coordinate, top at .1 mb
Atmospheric composition for year 2000
Ocean climatology prescribed from years 1996-2005, CMIP6
Uses turbulence scheme (no dry conv), grav.wave drag
Time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
Filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define STDHYB                   ! standard hybrid vertical coordinate
#define ATM_LAYERING L40         ! 40 layers, top at .1 mb
#define NEW_IO                   ! new I/O (netcdf) on
#define IRRIGATION_ON
#define MODIS_LAI
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
!  OFF #define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define DYNAMIC_BIOMASS_BURNING  ! alter biomas burning my flammability
#define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
#define ACCMIP_LIKE_DIAGS  ! adds many diags as defined by ACCMIP project
!<--- chemistry end
!---> MATRIX start
#define TRACERS_AMP
#define TRACERS_AMP_M1
!<--- MATRIX end
#define BC_ALB                    !optional tracer BC affects snow albedo
#define CLD_AER_CDNC              !aerosol-cloud interactions
#define BLK_2MOM                  !aerosol-cloud interactions
!  OFF #define NUDGE_ON                 ! nudge the meteorology
#define CACHED_SUBDD
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm144x90                           ! horizontal resolution is 144x90 -> 2x2.5deg
AtmLayering                         ! vertical resolution
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform

IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

#include "latlon_source_files"
#include "modelE4_source_files"
#include "static_ocn_source_files"

#include "tracer_shared_source_files"
#include "tracer_shindell_source_files"
#include "tracer_AMP_source_files"
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
#include "static_ocn_2000_144x90_input_files"
RVR=RD_Fd.nc             ! river direction file
NAMERVR=RD_Fd.names.txt  ! named river outlets

#include "land144x90_input_files"
#include "rad_input_files"
#include "rad_144x90_input_files_CMIP6clim"
#include "chemistry_input_files_nosoa"
#include "chemistry_144x90_input_files"
#include "dust_tracer_144x90_input_files"
#include "dry_depos_144x90_input_files"
#include "chem_emiss_144x90_input_files_CMIP6clim"
#include "aerosol_MATRIX_input_files_CMIP6clim"

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E6TmatrixF40clim2000 (climatological prescribed ocean atmospheric tracer model with MATRIX and Shindell chemistry)

&&PARAMETERS
#include "static_ocn_params"
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

#include "atmCompos_2000_params"
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use abs(o3_yr) instead!
!!!!!!!!!!!!!!!!!!!!!!!
#include "aerosol_MATRIX_params_CMIP6"
#include "dust_params_vmp_matrix"
#include "common_tracer_params_CMIP6"
#include "chemistry_params_CMIP6"
#include "ch4_params_CMIP6_2000"

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
