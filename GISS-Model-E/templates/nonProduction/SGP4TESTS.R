SGP4TESTS.R GISS Model E      M. Kelley 10/2013

A temporary copy of SCM.R for regression testing using pre-extracted
location-dependent input data for the SGP site

Initial framework for truly single-column mode for Model E.

Preprocessor Options
#define SCM
#define NEW_IO
End Preprocessor Options

Object modules: (in order of decreasing priority)

AtmL40
AtmRes

SCM_COM SCM
ATMDYN_SCM
ATMDYN_SCM_EXT

CLOUDS_COM CLOUDS2 CLOUDS2_DRV

SURFACE
PBL_COM PBL_DRV PBL ATURB

LANDICE LANDICE_COM SURFACE_LANDICE LANDICE_DRV

GHY_COM GHY_DRV

VEG_DRV

ENT_DRV  ENT_COM   ! + Ent          ! new vegetation

LAKES_COM LAKES

OCN_DRV OCEAN OCNML

SEAICE SEAICE_DRV ICEDYN_DUM

RAD_COM RADIATION RAD_DRV
COSZ_2D RAD_UTILS ALBEDO READ_AERO ocalbedo

DIAG_COM DIAG DEFACC QUICKPRT

ATM_DRV ATMDYN_COM ATM_UTILS ATM_COM

MODEL_COM
IO_DRV
MODELE
MODELE_DRV

QUS_COM QUSDEF

FLUXES

STRAT_DUM

Components:
shared MPI_Support solvers giss_LSM dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB

Data input files:

! SCM input files
SCM_NML=SCM_ARM.nml                                 ! input variable namelist with units
SCM_PS=sgp60varanarucC1.c1.20050101.000000.cdf      ! surface pressure
SCM_SFLUX=sgp60varanarucC1.c1.20050101.000000.cdf   ! surface heat fluxes
SCM_TSKIN=sgp60varanarucC1.c1.20050101.000000.cdf   ! skin temperature
! if horizontal and geostrophic wind profiles are specified, horizontal are initial
SCM_WIND=sgp60varanarucC1.c1.20050101.000000.cdf    ! horizontal wind profiles
!SCM_GEO=sgp60varanarucC1.c1.20050101.000000.cdf     ! geostrophic wind profiles
SCM_TEMP=sgp60varanarucC1.c1.20050101.000000.cdf    ! temperature profile(s)
SCM_WVMR=sgp60varanarucC1.c1.20050101.000000.cdf    ! water vapor mixing ratio profile(s)
SCM_OMEGA=sgp60varanarucC1.c1.20050101.000000.cdf   ! large-scale vertical wind
!SCM_LS_V=sgp60varanarucC1.c1.20050101.000000.cdf    ! large-scale vert adv flux div profile(s)
SCM_LS_H=sgp60varanarucC1.c1.20050101.000000.cdf    ! large-scale horiz adv flux div profile(s)

! Topography, area fractions of surface types
TOPO=SGP_extractions/Z2HX2fromZ1QX1N.nc

! Atm. initial conditions (temperature, wind, humidity, surface pressure)
! on the model's vertical grid and consistent with the
! orography from TOPO.  While some SCM modes
! of operation may subsequently overwrite the values from
! this file, the model has not yet been programmed to
! skip reading this file if it is absent.
AIC=SGP_extractions/AIC.RES_F40.D771201.nc

! Optional: if absent, ozone is set to zero.
! Zero stratospheric ozone is usually a bad idea though.
O3file=SGP_extractions/o3_2005_shindelltrop_144x90x49_1850-1997_ple.nc

! Optional: if absent, dust is set to zero
DUSTaer=SGP_extractions/dust_mass_CakmurMillerJGR06_144x90x20x7x12_unlim.nc

! Optional: if absent and MADAER flag not set, aerosols are zero.
! If these files are omitted, rundeck parameters od_cdncx and cc_cdncx
! must not be set (or if set, set to zero.)
! Currently, these files must be used as a group (will change in future).
TAero_SUL=SGP_extractions/SUL_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_SSA=SGP_extractions/SSA_Koch2008_kg_m2_144x90x20h.nc
TAero_NIT=SGP_extractions/NIT_Bauer2008_kg_m2_144x90x20_1890-2000h.nc
TAero_OCA=SGP_extractions/OCA_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_BCA=SGP_extractions/BCA_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_BCB=SGP_extractions/BCB_Koch2008_kg_m2_144x90x20_1890-2000h.nc

! These ocean files are only needed if the TOPO file for the
! SCM case contains a nonzero ocean fraction.  Certain SCM
! modes of operation may also prescribe ocean surface conditions
! via mechanisms other than these files. SICE and ZSIFAC can be
! omitted if the simulation location is free of sea ice.
OSST=SGP_extractions/OST_144x90.B.1975-1984avg.Hadl1.nc  ! SST
! rsi var. in SICE is sea ice fraction, ZSIFAC var. dm is used to get ice thickness
SICE=SGP_extractions/SICE_144x90.B.1975-1984avg.Hadl1.nc
ZSIFAC=SGP_extractions/SICE_144x90.B.1975-1984avg.Hadl1.nc
!OCNML=XXX ! see comment above regarding the OCNML source code.

! These land-surface files are only needed if the TOPO file for
! the SCM case contains a nonzero land fraction.  But note
! that the GIC contains initial conditions for surface types other
! than the land surface and that the SOILIC method is an optional
! replacement for the land-surface GIC.
VEG=SGP_extractions/V144X90_no_crops.ext.nc
CROPS=SGP_extractions/CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
SOIL=SGP_extractions/S144X900098M.ext.nc
TOP_INDEX=SGP_extractions/top_index_144x90_a.ij.ext.nc


! Optional land-surface file which prescribes soil initial
! conditions in intuitive intensive units (temperature,
! relative wetness, snow depth) rather than the extensive
! units (total heat and water per layer) of the arrays in
! the file GIC.
!SOILIC=/path/to/SOILIC.nc

! Optional file to specify a "background" roughness
! length for the land surface that is not tied to
! the scale-dependent topographic standard deviation in
! TOP_INDEX as per the procedure developed for 8x10 Model II.
! Note that the actual roughness length for the
! land surface is taken as the maximum of the vegetation-derived
! value and the value from either ROUGHL or TOP_INDEX.
!ROUGHL=/path/to/ROUGHL.nc

! Initial conditions for surface components from an arbitrary restart file
GIC=SGP_extractions/GIC.144X90.DEC01.1.ext_2.nc

! All input files below this line are location-independent.

GHG=GHG.Mar2009.txt
RADN1=sgpgxg.table8
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800
RADN4=LWCorrTables33k
RADN5=H2Ocont_MT_CKD
RADN3=miescatpar.abcdv2
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.SATO.1850-1999.Apr02_hdr
RADN9=solar.CMIP6official.ann1850-2299.nc ! need KSOLAR=2
RADNE=topcld.trscat8

! optional files for optional diagnostics
ISCCP=ISCCP.tautables
MSU_wts=MSU_SSU_RSS_weights.txt

Label and Namelist:
SGP4TESTS (SCM at SGP site for regression tests)


&&PARAMETERS

! SCM parameters
SCM_lon=-96.25             ! Southern Great Plains site longitude (deg)
SCM_lat=37.                ! Southern Great Plains site latitude (deg)
SCM_area=49370385348.1287  ! nominal grid box area (m2) from 144x90 grid
SCM_sfc=1                  ! 1:land,2:ocean
SCM_z0m=0.0005             ! surface roughness height (m)
SCM_alb=0.3                ! mid-visible surface albedo (-)
SCM_tau=10800.             ! nudging time constant

! Southern Great Plains Target Coordinates (degrees E/N).
! lon_targ and lat_targ only affect the solar zenith angle and
! corolis parameter (the latter is used in PBL parameterizations).
lon_targ=-96.25
lat_targ=37. 
! nominal gridbox area (m2) from a 144x90 lon-lat grid at lat_targ
nominal_area=49370385348.1287

DTsrc=1800.     ! Atm. physics timestep.
NIsurf=1        ! Number of surface physics timesteps per atm. physics timestep.
NRAD=1          ! Full radiation calculation every NRAD physics timesteps.

! cloud tuning parameters
U00a=.60
U00b=1.00
wmui_multiplier=2.0
entrainment_cont1=.4

! parameters that control temporally varying inputs:
! if set to 0, the current (day/) year is used: transient run
master_yr=1979
!crops_yr=1850  ! if -1, crops in VEG-file is used
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
volc_yr=-1
!volc_day=182
!aero_yr=1850
!o3_yr=-1850
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

! radiation flags
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
PTLISO=15.      ! press(mb) above which rad. assumes isothermal layers
madaer=3        ! indicates use of TAero_XXX aerosol files by radiation.


! parameters affecting diagn. output
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

! KCOPY=1: save acc and alternating checkpoint files only.
! KCOPY=2: save end-of-month rsf also (probably not useful for SCM).
KCOPY=1

! save alternating checkpoint files every Ndisk timesteps
! (checkpoint also saved when model reaches end time).
! Useful for perusing instantaneous states or checking
! checking whether model states are identical after coding
! rearrangements.
Ndisk=1440

! restart state is saveable every nssw timesteps.
Nssw=2

! SCM-useful GCM-native subdaily diagnostics system not yet imported to master branch
!SUBDD=' '        ! no sub-daily frequency diags
!NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)

! KOCEAN=0 means prescribed surface ocean conditions.  This parameter is currently
! mandatory even if ocean is absent at the SCM location.
KOCEAN=0

! variable lakes.  Probably unimportant for typical SCM simulation lengths.
! note that lakes can only evolve in response to local precip, evap, runoff.
variable_lk=1

&&END_PARAMETERS

 &INPUTZ
 YEARI=2000,MONTHI=1,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=2000,MONTHE=1,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=2000,MONTHE=1,DATEE=1,HOURE=1,
/

