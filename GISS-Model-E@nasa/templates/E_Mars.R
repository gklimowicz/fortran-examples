likeMars.R M20 Desert World with selected Mars-like parameters, inputs    M. Kelley 7/2013

Preprocessor Options
#define PLANET_PARAMS likeMars ! parameter set selector in shared/PlanetParams_mod.F90
#define NEW_IO
#define DO_CO2_CONDENSATION
#define DO_CO2_CONDENSATION_APPLY
End Preprocessor Options

Object modules: (in order of decreasing priority)
     ! resolution-specific source codes
Atm72x46                   ! horizontal resolution is 72x46 -> 4x5deg
AtmL20                     ! vertical resolution is 20 layers -> 0.1mb
DIAG_RES_M
STRAT_DUM
FFT72                              ! Fast Fourier Transform
IO_DRV                             ! new i/o

     ! GISS dynamics w/o gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers

    ! lat-lon grid specific source codes
AtmRes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO ocalbedo ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics
OCN_DRV                             ! driver for ocean-grid components
OCEAN OCNML                         ! ocean modules

Components:
shared MPI_Support solvers giss_LSM 
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB    

Data input files:

! Initial atmospheric temperature, specific humidity, wind, surface pressure.
! If rundeck parameter initial_psurf_from_topo=1, initial surface pressure is
! reset to be consistent with orography.
AIC=planet/Mars/AIC.coldmars.nc

! Initial soil temperature, relative wetness, and snow depth.
SOILIC=planet/desert_world/soilic_drysoil.nc

! This file contains Martian orography and defines the surface to be 100% land.
! Nominal longitudes might be shifted 180 degrees.
!TOPO=planet/Mars/marstopo.nc
TOPO=planet/Mars/marstopo_cap13K.nc

! If file ROUGHL is present, the variable top_dev in TOP_INDEX will only affect
! snow masking and its precise value is probably not crucial.
TOP_INDEX=planet/desert_world/stdev_72x46_desertworld.nc

! The presence of this file causes the model to read roughness length from
! it rather than computing roughness from topographic standard deviation
! as in Model II.  This allows users to control this field directly.
! Units: meters.
ROUGHL=planet/desert_world/z0m_72x46_desertworld.nc

! Soil textures file (array q).  Sand (imt=1) is a reasonable texture to start with.
! Arrays sl and qk do not matter for dry conditions.  dz is layer thickness.
SOIL=planet/desert_world/soil_allsand.nc

! In modelE, the percentages of bright and dark soil are read from
! the vegetation file.  The ratio of the two soils is chosen to give a
! surface albedo of roughly 15% for dry conditions.
VEG=planet/Mars/veg_allbare_alb15.nc

! N2O, CH4, and CFCs are set to zero in this file, which also has the
! preindustrial CO2 mixing ratio serving as the basis for adjusting
! CO2x.  Note that the absence of ozone and aerosol files in this
! rundeck implies that those species are also zero.
GHG=planet/desert_world/GHG.CO2only.txt

! radiation input files
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
RADN3=miescatpar.abcdv2
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADNE=topcld.trscat8

!RADN4=LWCorrTables33k ! correction factors for Earth conditions are not applicable

!Ent needs this
SOILCARB_global=soilcarb_top30cm_4x5.nc

Label and Namelist:
likeMars (M20 Desert World with selected Mars-like parameters, inputs)

&&PARAMETERS

! calculate initial surface pressure for hydrostatic consistency with topography
initial_psurf_from_topo=1

! linear damping timescales (sec) for winds in the top 3 layers.
! These probably need to be tuned.  Change the number of values here
! to change the number of layers in which damping is applied.
rtau=300000.,200000.,100000.

! minimum and maximum allowed column mass (kg/m2) for error checking
mincolmass=0.
maxcolmass=500.

planetName = 'Mars'
! source: http://en.wikipedia.org/wiki/Mars
eccentricity = 0.093
obliquity = 25.19d0 ! degrees
longitudeAtPeriapsis = 251.0! degrees
siderealOrbitalPeriod  = 59354294.4 ! seconds
siderealRotationPeriod = 88642.6848 ! seconds
meanDistance = 1.52366231 ! AU


! scaling factor for solar brightness is now from mean distance.
s0x=1.00

! scaling factor for O2 amounts. 0.7% of present-day Earth gives 0.146%.
o2x=0.007

! scaling factor for CO2 amounts.  Increase by the ratio of 95% molar CO2 to
! whatever concentration is in the GHG file (which has ppmv units)
!co2x=1.
co2x=3333.

! -1 uses Thekaekara solar spectrum, so no RADN9 input file required
KSOLAR=-1

! A dry planet will have no lakes, but we set this in case someone
! eventually decides to start with an initially wet soil.
variable_lk=1 ! lets lake fraction vary with time

! Physics timestep (sec). Cannot be changed after a run has been started.
DTsrc=1850.

! Dynamics timestep.  This is half of the value for 4x5 Earth runs.
! Note that DTsrc/DT should be an even integer.
DT=230.

! parameters that control the Shapiro filter
DT_XUfilter=460. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=460. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! Number of surface timesteps per physics timestep.  Probably does not
! need to be changed for M20 Mars.
NIsurf=1        ! increase as layer 1 gets thinner

! Number of physics timesteps per radiation timestep.  Default is 5.
nrad=1

! save alternating checkpoint files every Ndisk physics timesteps
Ndisk=1440

! KCOPY=1: save acc and alternating checkpoint files only.  KCOPY=2: save rsf also
KCOPY=2

nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc

! restart state is saveable every nssw timesteps.
nssw=2          ! until diurnal diagn. are fixed, nssw should be even

! the model requires the following two parameters to be set, though they
! do not affect Mars runs
master_yr=1850
KOCEAN=0

&&END_PARAMETERS

 &INPUTZ
 YEARI=0001,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=0001,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=0001,MONTHE=12,DATEE=1,HOURE=1,
/

