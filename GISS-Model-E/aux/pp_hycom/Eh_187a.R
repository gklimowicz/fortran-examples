Eh_187a.R GISS Model E  coupled version     tnl 09/22/16

TNL the first attempt to setup hycom master branch

Preprocessor Options
#define NEW_IO                   ! new I/O (netcdf) on
#define SWFIX_20151201
#define NO_HDIURN                ! exclude hdiurn diagnostics
#define MODIS_LAI
#define CHECK_OCEAN                  ! needed to compile aux/file CMPE002
! #define TRACERS_GASEXCH_Natassa    ! special tracers to be passed to ocean
! #define TRACERS_GASEXCH_Natassa    ! special tracers to be passed to ocean
! #define TRACERS_HYCOM_Ventilation
#define ATM2x2h                      ! 2x2.5 40 layer atm
#define HYCOM1degRefined             ! 1deg   refined hycom (387x360)
! #define HYCOM1degUnrefined         ! 1deg unrefined hycom (359x360)
! #define HYCOM26layers              ! use 26 layers in hycom
#define HYCOM32layers                ! use 32 layers in hycom
#define EXPEL_COASTAL_ICEXS          ! Attempt decrease the coastal ice thicknes
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm144x90                           ! horizontal resolution is 144x90 -> 2x2.5deg
AtmL40                              ! vertical resolution is 40 layers -> 0.1mb
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform

IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

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
SparseCommunicator_mod              ! sparse gather/scatter module
hycom_arrays|-r8| hycom_dim|-r8| kprf_arrays|-r8|
kprf_arrays_loc_renamer|-r8| hycom_atm|-r8|
hycom_arrays_glob|-r8| hycom_arrays_glob_renamer|-r8|
hycom_scalars|-r8| hycom_dim_glob|-r8|
hycom |-r8| OCEAN_hycom|-r8|        ! ocean model - driver
advfct|-r8|                         ! advection
archyb|-r8|                         ! archiving in netcdf
barotp|-r8|                         ! barotropic eqn.
bigrid|-r8|                         ! basin grid
blkprf|-r8|                         ! block data
cdf_io|-r8|                         ! netcdf utility
cnuity|-r8|                         ! continuity eqn.
convec|-r8|                         ! convection
cpler |-r8|                         ! coupler
diapfx|-r8|                         ! diapycnal diffusion
dpthuv|-r8| dpudpv|-r8|             ! off-center depth
eice  |-r8|                         ! ice forming
geopar|-r8|                         ! geography related parameters
hybgn1|-r8|                         ! grid generator
inicon|-r8| inigis|-r8| inikpp|-r8| ! initial conditions
matinv|-r8| mxkprf|-r8| mxlayr|-r8| ! mixing scheme
momtum|-r8|                         ! momemtum Eqn.
prtetc|-r8|                         ! print routines, etc.
reflux|-r8|                         ! flux conversion
sigetc|-r8|                         ! eqn.of state, etc.
thermf|-r8|                         ! thermal forcing
trcadv|-r8|                         ! tracer advection
tsadvc|-r8| advem|-r8|              ! advecting t/s

Components:
shared MPI_Support solvers giss_LSM
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m  ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
AIC=AIC.RES_F40.D771201.nc      ! observed init cond (atm. only)
GIC=GIC.144X90.DEC01.1.ext_1.nc ! initial ground conditions
temp_ini=temp387x360x32jan_sig1a_topo2009_oct2015.txt  ! 3-d temperature as initial condition
salt_ini=saln387x360x32jan_sig1a_topo2009_oct2015.txt  ! 3-d salinity as initial condition
pout_ini=pout387x360x32jan_sig1a_topo2009_oct2015.txt  ! 3-d layer pressure as initial condition
latlonij=latlon387x360_oct2015.4bin                    ! lat & lon at 4 positions in each ocean grid box
ibasin=ibasin387x360_topo2009_col3_oct2015.txt         ! ocean basin mask
wgt_a2o=wgt_a2o_atm2x2h_h387x360_topo2009_may2016.8bin ! coupler weights for flux/scalar from atm to ocean, both on A grid
wgt_o2a=wgt_o2a_atm2x2h_h387x360_topo2009_may2016.8bin ! coupler weights for flux/scalar from ocean to atm, both on A grid
cososino=cososino387x360_oct2015.8bin                  ! cos/sin of i,j axis angle on ocean grid
kpar=seawifs_kpar_h387x360_topo2009_oct2015.t4big      ! monthly/annual seawifs_kpar data
hycomtopo=depth387x360_2009.4big  	                   ! topography used in ocean model with Baltic Sea
TOPO=Z144X90N.h387x360_topo2009_may2016.nc             ! surface fractions and topography
RVR=RD_modelE_Fa.RVR_h387x360_topo2009_oct2015.nc      ! river direction file
NAMERVR=RD_Fb.names.txt  ! named river outlets

ICEDYN_MASKFAC=iceflowmask_144x90.nc

CDN=CD144X90.ext.nc
VEG=V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.nc
LAIMAX=V144x90_EntMM16_lai_max_trimmed_scaled_ext1.nc
HITEent=V144x90_EntMM16_height_trimmed_scaled_ext1.nc
LAI=V144x90_EntMM16_lai_trimmed_scaled_ext1.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
SOIL=S144X900098M.ext.nc
TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
GLMELT=GLMELT_144X90_gas.OCN.nc
RADN1=sgpgxg.table8                           ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800  ! rad.tables and history files
RADN4=LWCorrTables33k                         ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-2012.May13_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean2015.ann1610-2014_hdr ! need KSOLAR=2
RADNE=topcld.trscat8

ISCCP=ISCCP.tautables
GHG=GHG_RCP45.txt ! use GHG.Jul2009.txt for runs that start before 1850
dH2O=dH2O_by_CH4_monthly
DUSTaer=TcadiAR5_aero/144x90/DUST_Tcadi2012_Bauer_kg_m2_144x90x40.nc
BC_dep=TcadiAR5_aero/144x90/BC_dep_Tcadi2012_Bauer_kg_m2_s_144x90_1850-2100.nc
! updated aerosols need MADAER=3
TAero_SUL=TcadiAR5_aero/144x90/SUL_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100.nc
TAero_SSA=TcadiAR5_aero/144x90/SSA_Tcadi2012_Bauer_kg_m2_144x90x40.nc
TAero_NIT=NIT_Bauer2008_kg_m2_144x90x40_1850-2100h.nc
TAero_OCA=TcadiAR5_aero/144x90/OCA_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100.nc
TAero_BCA=TcadiAR5_aero/144x90/BCA_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100.nc
TAero_BCB=TcadiAR5_aero/144x90/BCB_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100.nc


! O3file soon to be replaced by one from latest chemistry code
O3file=jan2012_o3_shindell_144x90x49x12_1850-2010_ple.nc
Ox_ref=o3_2010_shindell_144x90x49_April1850.nc

MSU_wts=MSU.RSS.weights.data      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
Eh_187a (E187 master branch the first attempt)

&&PARAMETERS

! parameters set for coupled ocean runs:
KOCEAN=1            ! ocn is prognostic
OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
OTIDE = 0           !  Ocean tides are not used
variable_lk=1
init_flake=1
! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.00010  ! threshold (1/s) for triggering deformation waves
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.1       ! default is 0.5
CDEF=1.6       ! tuning factor for deformation -> momentum flux
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF


! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! The following two lines are only used when aerosol/radiation interactions are off
FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.64 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.   ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 2.
use_vmp=1
radius_multiplier=1.1

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
!crops_yr=1850  ! if -1, crops in VEG-file is used
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
volc_yr=-1
!volc_day=182
!aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=1850
dalbsnX=.024
!o3_yr=-1850
!aer_int_yr=1850    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)

! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=2          ! saving acc + rsf
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=48          ! until diurnal diags are fixed, Nssw has to be even
Ndisk=480

itest=-1         ! default is -1
jtest=-1         ! default is -1
iocnmx=2         ! default is 2
brntop=50.       ! default is 50.
brnbot=200.      ! default is 200.
diapyn=2.e-7     ! default is 2.e-7
diapyc=.2e-4     ! default is .2e-4
jerlv0=1         ! default is 1; 0 => read dataset
bolus_biharm_constant=0   ! bolus_biharm_constant=1 uses thkdff=0.05 or 0.10 m/s
bolus_laplc_constant=1       ! bolus_laplc_constant=1  uses thkdff=0.01 or 0.02 m/s
bolus_laplc_exponential=0 ! bolus_laplc_exponential=1 uses thkdff=0.03 m/s; thkdff_bkgd is hard-wired to 0.003m/s
thkdff=.01
&&END_PARAMETERS

 &INPUTZ
 YEARI=1900,MONTHI=01,DATEI=01,HOURI=00, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1900,MONTHE=01,DATEE=01,HOURE=00,   KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1900,MONTHE=01,DATEE=02,HOURE=00,
/
