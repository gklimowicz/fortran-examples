nonProduction_E4TcadC12  GISS Model E Tom Clune 02/28/2011

nonProduction_E4TcadC12.R  = combination of template EC12 and E4TcadF40 decks, then 
altered a bit.

This is for faster testing of tracer code, not scrutinized for "science" purposes.

Preprocessor Options
#define NEW_IO                   ! new I/O (netcdf) on
#define NO_HDIURN                ! exclude hdiurn diagnostics
#define OLD_BCdalbsn
! OFF  #define MODIS_LAI
!---> generic tracers code start
#define TRAC_ADV_CPU             ! timing index for tracer advection on
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
#define TRACERS_DRYDEP           ! default dry deposition
#define ALLOW_MORE_DRYDEP_NTYPE  ! larger dimension needed for 8x10 regridded VEGTYPE
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
!<--- generic tracers code end
!---> chemistry start
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
! OFF  #define RAD_O3_GCM_HRES     ! Use GCM horiz resl to input rad code clim Ozone
#define RAD_O3_2010              ! 2010 ozone dataset
!  OFF #define AUXILIARY_OX_RADF ! radf diags for climatology or tracer Ozone
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
!  OFF #define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define DYNAMIC_BIOMASS_BURNING  ! alter biomas burning my flammability
!  OFF #define CALCULATE_LIGHTNING ! turn on Colin Price lightning when TRACERS_SPECIAL_Shindell off
!  OFF #define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!  OFF #define ACCMIP_LIKE_DIAGS  ! adds many diags as defined by ACCMIP project
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
#define BIN_OLSON   ! use binary versions of vegtype and LAI files for tracers
!  OFF #define NUDGE_ON                 ! nudge the meteorology
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm36x24                     ! horizontal resolution is 36x24 -> 8x10
AtmL12 STRAT_DUM             ! vertical resolution is 12 layers -> 10mb
DIAG_RES_M                   ! diagnostics
FFT36                        ! Fast Fourier Transform

IO_DRV                              ! new i/o

    ! GISS dynamics
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
!   STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

#include "tracer_shared_source_files"
#include "tracer_shindell_source_files"
#include "tracer_OMA_source_files"
TRDIAG                              ! new i/o

#include "latlon_source_files"
#include "modelE4_source_files"
CLD_AEROSOLS_Menon_MBLK_MAT_E29q BLK_DRV ! aerosol-cloud interactions
CLD_AER_CDNC            ! aerosol-cloud interactions wrapper
! flammability_drv flammability       ! Olga's fire model

#include "static_ocn_source_files"

Components:
#include "E4_components_nc"    /* without "Ent" */
tracers
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB /* needed for "Ent" only */
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
#include "IC_36x24_input_files"
#include "static_ocn_1950_36x24_input_files"  /* 1880 for 8x10 needs to be made */
RVR=RD8X10.nc            ! river direction file
NAMERVR=RD8X10.names.txt ! named river outlets

#include "land36x24_input_files"
#include "rad_input_files"      ! CO2profile is disregarded since LM<27
#include "rad_36x24_input_files"

#include "chemistry_input_files"
#include "chemistry_36x24_input_files"

#include "dust_tracer_36x24_input_files"
#include "dry_depos_36x24_input_files"

#include "chem_emiss_36x24_input_files"

#include "aerosol_OMA_36x24_input_files"
Ox_ref=o3_zeros_36x24x49.nc

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG8X10                      ! special regions-diag

Label and Namelist:  (next 2 lines)
TCC12_E4Tcad (ModelE4 8x10, 12 layers, with Tcad tracers for testing only)

&&PARAMETERS
#include "static_ocn_params"
#include "sdragC12_params"
#include "gwdragC12_params"

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! The following two lines are only used when aerosol/radiation interactions are off
FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.54  ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00  ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 1.
use_vmp=1
radius_multiplier=1.1

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=0.      ! activates strat.H2O generated by CH4
                 ! [Turn off when Shindell tracer CH4 on & clim_interact_chem=1]
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

initial_GHG_setup = 1 ! Set to 0 after initial setup.

#include "atmCompos_36x24_params"
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols
#include "aerosol_OMA_36x24_params"
#include "dust_params_vmp_oma"
#include "chemistry_36x24_params"
! The following 2 lines OVERWRITE the include chemistry_params values!!
! ch4_init_sh=1.750      ! init cond/fixed conditions SH CH4 ppmv
! ch4_init_nh=1.855      ! init cond/fixed conditions NH CH4 ppmv

DTsrc=1800.      ! cannot be changed after a run has been started
DT=450.          ! could be 900 if DTsrc=3600.  do not forget nda5 et al.
! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)
#include "diag_params"
! save3dAOD=1      ! needed if 3D AOD (itAOD or ictAOD) SUBDDs are on and adiurn_dust=0

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=480
&&END_PARAMETERS

 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
/
