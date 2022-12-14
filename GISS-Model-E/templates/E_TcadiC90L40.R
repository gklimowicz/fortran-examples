E_TcadiC90L40.R                                     jan perlwitz 09/2011
E_TcadiC90L40: E4TcadC90L40 + computed aerosol 1st indirect effect

E4C90L40 = E4F40 on cubed sphere grid 6x90x90 grid boxes
E4F40 = modelE as frozen in April 2010:
Cubed Sphere grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics: finite Volume ; physics 30 min.; radiation 2.5 hrs
filters: none

Preprocessor Options
#define TRAC_ADV_CPU             ! timing index for tracer advection on
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
#define TRACERS_DUST             ! include dust tracers
#define TRACERS_DUST_Silt4       ! include 4th silt size class of dust
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
!  OFF #define AUXILIARY_OX_RADF ! radf diags for climatology or tracer Ozone
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
#define TRACERS_AEROSOLS_SEASALT ! seasalt
#define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
#define TRACERS_AEROSOLS_SOA     ! Secondary Organic Aerosols
!  OFF #define SOA_DIAGS                ! Additional diagnostics for SOA
#define TRACERS_NITRATE
#define TRACERS_HETCHEM
#define BC_ALB                      ! optional tracer BC affects snow albedo
#define CLD_AER_CDNC                ! aerosol-cloud interactions
#define BLK_2MOM                    ! aerosol-cloud interactions
!  OFF #define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define DYNAMIC_BIOMASS_BURNING  ! alter biomas burning by flammability
!  OFF #define CALCULATE_LIGHTNING ! turn on Colin Price lightning when TRACERS_SPECIAL_Shindell off
!  OFF #define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!  OFF #define NUDGE_ON                 ! nudge the meteorology
!  OFF #define ACCMIP_LIKE_DIAGS  ! adds many diags as defined by ACCMIP project
#define NEW_IO
#define CALC_GWDRAG
#define SET_SOILCARBON_GLOBAL_TO_ZERO
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
AtmCS90                           ! 90 Cube-Sphere Grid
AtmL40                             ! vertical resolution is 40 layers -> 0.1mb
     ! Codes used by the cubed-atmosphere configuration (FV dynamics)
#include "cubed_sphere_source_files"

TRDUST_COM TRDUST TRDUST_DRV        ! dust tracer specific code
#include "tracer_shared_source_files"
TRDIAG                              ! for offline postprocessing
#include "tracer_shindell_source_files"
#include "tracer_OMA_source_files"

STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)
#include "modelE4_source_files"
CLD_AEROSOLS_Menon_MBLK_MAT_E29q BLK_DRV ! aerosol-cloud interactions
CLD_AER_CDNC            ! aerosol-cloud interactions wrapper
! flammability_drv flammability       ! Olga's fire model

#include "static_ocn_source_files"

Components:
#include "E4_components"    /* without "Ent" */
tracers
Ent
dd2d

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB  /* needed for "Ent" only */
OPTS_dd2d = NC_IO=PNETCDF
FVCUBED = YES

Data input files:
#include "IC_CS90_input_files"
#include "static_ocn_1880_CS90_input_files"
!#include "static_ocn_2000_CS90_input_files"

RVR=RDdistocean_CS90_EM.nc             ! river direction file
NAMERVR=RDdistocean_CS90_EM.names.txt  ! named river outlets

! OFF VEG_DENSE=veg_dense_C90_from_2x2.5

! OFF ! ----- for interactive wetlands -----
! OFF PREC_NCEP=ncep_prec_w_2wk_lag_2x2.5_C90
! OFF TEMP_NCEP=ncep_g1temp_2x2.5_C90
! OFF BETA_NCEP=beta_p_ch4_4x5_2x2.5gf_C90
! OFF ALPHA_NCEP=alpha_t_ch4_4x5_2x2.5gf_C90

#include "landCS90_input_files"
#include "rad_input_files"
#include "rad_C90_input_files"

#include "chemistry_input_files"
#include "chemistry_C90_input_files"

#include "dust_tracer_C90_input_files"
#include "dry_depos_C90_input_files"

#include "chem_emiss_C90_input_files"

#include "aerosol_OMA_C90_input_files"

MSU_wts=MSU_SSU_RSS_weights.txt     ! MSU-diag
REG=REG.txt                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E_TcadiC90L40 (E4TcadC90L40 + computed aerosol 1st indirect effect)


&&PARAMETERS
#include "static_ocn_params"
#include "sdragCS90_params"
#include "gwdragCS90_params"

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.65 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.20 ! below 850mb and MC regions; tune this last  to get rad.balance

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=0.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

initial_GHG_setup = 1 ! Set to 0 after initial setup.

#include "atmCompos_1850_params"
!#include "atmCompos_2000_params"
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use abs(o3_yr) instead!
!!!!!!!!!!!!!!!!!!!!!!!
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols
#include "aerosol_OMA_params"
#include "dust_params_oma"
#include "common_tracer_params"
#include "chemistry_params"

DTsrc=1800.      ! cannot be changed after a run has been started
DT=1800.         ! for FV dynamics, set same as DTsrc

NIsurf=1         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)

#include "diag_params"

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960
&&END_PARAMETERS

 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
/
