SCM.R GISS Model E      M. Kelley 10/2013

Initial framework for truly single-column mode for Model E.

This is a functioning rundeck, not a template (excepting the lines
containing  /path/to/user/directory/extractions - see notes below)

SCM-irrelevant codes and input files are excluded.
Template #includes should be refactored so that this exclusion happens automatically.

SCM case: SGP Jan 2005
For other cases, change one or more of the following as necessary:
(1) SCM input variable namelist (SCM_NML) and SCM input files (SCM_PS, SCM_SFLUX, etc.)
(2) SCM parameters SCM_lon, SCM_lat, and files extracted from gridded data at SCM_lon, SCM_lat
    See notes in the input files section on a pre-scripted extraction procedure etc.
    SCM_area and other SCM parameters may also be changed.

Preprocessor Options
#define CACHED_SUBDD
#define SCM
#define NEW_IO
End Preprocessor Options

Object modules: (in order of decreasing priority)

SUBDD

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
! VEG_COM VEGETATION 
ENT_DRV ENT_COM ! + Ent

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
shared MPI_Support solvers giss_LSM 
dd2d
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

! The set of forcings for a particular SCM test case typically does not include
! all of the data required to run Model E.  Each line below of the form
!   SHORTNAME=/path/to/user/directory/extractions/filename.nc
! corresponds to a location-dependent input dataset for which one of the following
! two options must be chosen.  The second option is provided for convenience.
! (1) Replace /path/to/user/directory/extractions/filename.nc with the
!     path to a single-column file that has already been created somehow.
!     Note that the GCM file-reading infrastructure considers single-column
!     files to be on a horizontal grid of size 1, and therefore the
!     two horizontal dimensions must be retained in netcdf variables
!     (which must have the netcdf names that model E expects).
!     For any files previously extracted from global data via option (2),
!     make sure that the path does not contain substring "/extractions/"
!     if option (2) will be used create other files.
! (2) Use exec/extract_scm.sh to sample gridded files at location
!     lon_targ, lat_targ.
!     Firstly,
!       replace /path/to/user/directory with a real path, preferably which
!        (a) contains a string denoting the SCM location/case being run
!        (b) is unlikely to be chosen by any other users on the system
!       Habits (a) and (b) will prevent clutter and accidental overwrites.
!       Note that
!        (a) substring "/extractions/" must be retained in each path
!        (b) the rundeck paths indicate the resulting single-column files
!            to be read by the model.  GCMSEARCHPATH (from your modelErc)
!            is the location of the gridded file from which to extract the column
!        (c) There is no requirement that all files from which data are extracted are
!            on the same grid - dimension and coordinate information is scanned per-file.
!        (d) As currently programmed (10/2013), extract_scm.sh does require that
!            each file possesses 1D coordinate variables named lon(lon) and lat(lat).
!            Extraction from arbitrary grids (cubed-sphere etc.) will soon be enabled.
!        (e) NCO must be installed on your system and in your $PATH.
!     Secondly, execute
!        extract_scm.sh THISRUNDECK.R


! Topography, area fractions of surface types
TOPO=/path/to/user/directory/extractions/Z2HX2fromZ1QX1N.nc

! Atm. initial conditions (temperature, wind, humidity, surface pressure)
! on the model's vertical grid and consistent with the
! orography from TOPO.  While some SCM modes
! of operation may subsequently overwrite the values from
! this file, the model has not yet been programmed to
! skip reading this file if it is absent.
AIC=/path/to/user/directory/extractions/NCARIC.144x90.D7712010_ext.nc

! Optional: if absent, ozone is set to zero.
! Zero stratospheric ozone is usually a bad idea though.
O3file=/path/to/user/directory/extractions/o3_2005_shindelltrop_144x90x49_1850-1997_ple.nc

! Optional: if absent, dust is set to zero
DUSTaer=/path/to/user/directory/extractions/dust_mass_CakmurMillerJGR06_144x90x20x7x12_unlim.nc

! Optional: if absent and MADAER flag not set, aerosols are zero.
! If these files are omitted, rundeck parameters od_cdncx and cc_cdncx
! must not be set (or if set, set to zero.)
! Currently, these files must be used as a group (will change in future).
TAero_SUL=/path/to/user/directory/extractions/SUL_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_SSA=/path/to/user/directory/extractions/SSA_Koch2008_kg_m2_144x90x20h.nc
TAero_NIT=/path/to/user/directory/extractions/NIT_Bauer2008_kg_m2_144x90x20_1890-2000h.nc
TAero_OCA=/path/to/user/directory/extractions/OCA_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_BCA=/path/to/user/directory/extractions/BCA_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_BCB=/path/to/user/directory/extractions/BCB_Koch2008_kg_m2_144x90x20_1890-2000h.nc

! These ocean files are only needed if the TOPO file for the
! SCM case contains a nonzero ocean fraction.  Certain SCM
! modes of operation may also prescribe ocean surface conditions
! via mechanisms other than these files. SICE and ZSIFAC can be
! omitted if the simulation location is free of sea ice.
!OSST=/path/to/user/directory/extractions/OST_144x90.1996-2005avg.HadISST1.1.nc  ! SST
! rsi var. in SICE is sea ice fraction, ZSIFAC var. dm is used to get ice thickness
!SICE=/path/to/user/directory/extractions/SICE_144x90.1996-2005avg.HadISST1.1.nc
!ZSIFAC=/path/to/user/directory/extractions/SICE_144x90.1996-2005avg.HadISST1.1.nc
!OCNML=XXX ! see comment above regarding the OCNML source code.

! These land-surface files are only needed if the TOPO file for
! the SCM case contains a nonzero land fraction.  But note
! that the GIC contains initial conditions for surface types other
! than the land surface and that the SOILIC method is an optional
! replacement for the land-surface GIC.
VEG=/path/to/user/directory/extractions/V144X90_no_crops.ext.nc
CROPS=/path/to/user/directory/extractions/CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
SOIL=/path/to/user/directory/extractions/S144X900098M.ext.nc
SOILCARB_global=/path/to/user/directory/extractions/soilcarb_top30cm_2x2.5.nc
TOP_INDEX=/path/to/user/directory/extractions/top_index_144x90_a.ij.ext.nc


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
GIC=/path/to/user/directory/extractions/GIC.144X90.DEC01.1.ext_2.nc

! All input files below this line are location-independent.

GHG=GHG.Mar2004.txt
RADN1=sgpgxg.table8
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800
RADN4=LWCorrTables33k
RADN5=H2Ocont_MT_CKD
RADN3=miescatpar.abcdv2
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.SATO.1850-1999.Apr02_hdr
RADN9=solar.DBbglean.ann850-2000.uvflux.nc
RADNE=topcld.trscat8

! optional files for optional diagnostics
ISCCP=ISCCP.tautables
MSU_wts=MSU_SSU_RSS_weights.txt

Label and Namelist:
SCM (documenting the Single Column Model)


&&PARAMETERS

! SCM parameters
SCM_lon=-97.49             ! Southern Great Plains site longitude (deg)
SCM_lat=36.61              ! Southern Great Plains site latitude (deg)
SCM_area=49370385348.1287  ! nominal grid box area (m2) from 144x90 grid
SCM_sfc=1                  ! 1:land,2:ocean
SCM_z0m=0.0005             ! surface roughness height (m)
SCM_alb=0.3                ! mid-visible surface albedo (-)
SCM_tau=10800.             ! nudging time constant

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

SUBDD='u v t q rh z p_3d p_surf prec mcp ssp snowfall snowdp qcl qci'
SUBDD1='cldss cldmc cldss_2d totcld totcld_diag'
SUBDD2='gtempr shflx lhflx ustar pblht pwv lwp iwp tau_ss tau_mc'
SUBDD3='olrrad olrcs lwds lwdscs lwus swds swus swdf egcm'
SUBDD4='dq_turb dth_turb dq_mc dth_mc dq_ss dth_ss dth_sw dth_lw dth_rad'
SUBDD5='dq_ls dth_ls dq_nudge dth_nudge'
SUBDD6='isccp_sunlit isccp_ctp isccp_tau isccp_lcld isccp_hcld'
NSUBDD=1         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
SCM_PlumeDiag=0  !to save Plume diagnostics set SCM_PlumeDiag=1
WRITE_ONE_FILE=1 ! all outputs to a single file

! KOCEAN=0 means prescribed surface ocean conditions.  This parameter is currently
! mandatory even if ocean is absent at the SCM location.
KOCEAN=0

! variable lakes.  Probably unimportant for typical SCM simulation lengths.
! note that lakes can only evolve in response to local precip, evap, runoff.
variable_lk=1

&&END_PARAMETERS

 &INPUTZ
 YEARI=2005,MONTHI=1,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=2005,MONTHE=1,DATEE=31,HOURE=23,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=2005,MONTHE=1,DATEE=31,HOURE=23,
/
