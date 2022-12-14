E3hyc01_mpi.R.R GISS Model E  coupled version      rar   8/31/2005

E3hyc01_mpi.R: E3AoM20A + hycom kpp refinement + topo_20w

modelE1 (3.0) 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1880
ocean: coupled to HYCOM Version 0.9.2 KPP
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag)
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W direction (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm72x46                          ! horizontal resolution is 72x46 -> 4x5deg
AtmL20 STRAT_DUM                 ! vertical resolution is 20 layers -> 0.1mb
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
              ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_M                          ! ESMF
      FFT72                         ! utilities
POUT                                ! post-processing output
hycom |$(R8) -O2 | OCEAN_hycom|$(R8) -O2 | ! ocean driver
advfct|$(R8) -O2 |                         ! advection
archyb|$(R8) -O2 |                         ! continuity eqn. 
barotp|$(R8) -O2 |                         ! barotropic eqn. 
bigrid|$(R8) -O2 |                         ! basin grid
blkd10|$(R8) -O2 | blkpp2|$(R8) -O2 |       ! block data
cnuitb|$(R8) -O2 |                         ! continuity eqn.
cpler |$(R8) -O2 |                         ! coupler
dpthuv|$(R8) -O2 | dpudpv|$(R8) -O2 | ! off-center depth  
eic8  |$(R8) -O2 |                         ! ice forming
geopar|$(R8) -O2 |                         ! geography related parameters
hybg05|$(R8) -O2 |                         ! grid generator 
inirfn|$(R8) -O2 | inigis|$(R8) -O2 | inikpp|$(R8) -O2 | ! initial conditions
matinv|$(R8) -O2 | mxkprf|$(R8) -O2 | ! KPP mixing scheme
momtum|$(R8) -O2 |                         ! momemtum Eqn.
prtetc|$(R8) -O2 |                         ! print routines, etc.
reflux|$(R8) -O1 |                         ! flux conversion
sigetc|$(R8) -O2 |                         ! eqn.of state, etc.
thermf|$(R8) -O2 |                         ! thermal forcing
trcadv|$(R8) -O2 |                         ! tracer advection 
tsadvc|$(R8) -O2 |                         ! T/S advection 
hybrid_mpi_omp_renamer|-O2 |            ! ESMF
hybrid_mpi_omp_coupler|-O2 |            ! ESMF

Components:
MPI_Support shared

Data input files:
AIC=AIC.RES_M20A.D771201.nc    !initial conditions (atm.) needs GIC,OIC ISTART=2
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc ! initial conditions (ground) and 300 year spin-up
CDN=CD4X500S.ext.nc
  ! VEG=V72X46.1.cor2.ext
VEG=V72X46.1.cor2_no_crops.ext.nc
CROPS=CROPS2007_72X46N.cor4_nocasp.nc  ! veg. fractions, crops history
SOIL=S4X50093.ext.nc
TOPO=Z72X46N.2deg_rfn_20w              !!! hycom
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR.2deghycom_20w.bin         !!! hycom
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800              ! rad.tables and history files
RADN4=LWCorrTables33k              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_Ma_2008
RADN3=miescatpar.abcdv2
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr     ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
#include "rad_72x46_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46_a.ij.ext.nc
MSU_wts=MSU_SSU_RSS_weights.txt
GLMELT=GLMELT_4X5.OCN.nc   ! glacial melt distribution

Label and Namelist:
E3hyc01_mpi.R (bihar=.1,hybgn1_5m,pump[20:200],full kpp)

DTFIX=300
&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1        ! ocn is prognostic
!glmelt_fac_nh=0.
!glmelt_fac_sh=1.34
!river_fac=1.05

! parameters usually not changed when switching to coupled ocean:

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_SDRAG=1         ! conserve ang. mom.

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

 
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.59      ! U00ice up => nethtz0 down (alb down); goals: nethtz0=0,plan.alb=30%
U00wtrX=1.40    ! U00wtrX+.01=>nethtz0+.7                for global annual mean
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1850 ! if -1, crops in VEG-file is used  ? 1979
s0_yr=1850
s0_day=182
ghg_yr=1850
ghg_day=182
!volc_yr=1850
volc_yr=-1
volc_day=182
aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1850
dalbsnX=.015
o3_yr=1850
!!!!   SNP SBP SSP ANP ONP OBP BBP  SUI ANI  OCI BCI OCB BCB
!! dfl:1.0,1.0,1.0,1.0,2.5,2.5,1.9, 1.0,1.0, 2.5,1.9,2.5,1.9
aermix=1.0,1.0,1.0,1.0,2.5,2.5,1.9, 1.0,1.0, 2.5,1.9,2.5,1.9

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
DTsrc=1800.
DT=450.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480       ! use =48 except on halem
Nssw=48         ! check flag once a day
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time if isccp-diags are not essential
! cloud_rad_forc=0 ! use =1 to activate this diagnostic (doubles radiation calls !)
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
&&END_PARAMETERS

 &INPUTZ
   YEARI=1800,MONTHI=01,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=1800,MONTHE=01,DATEE=3,HOURE=0, KDIAG=13*9,
   ISTART=2,IRANDI=0,YEARE=1800,MONTHE=01,DATEE=3,HOURE=0, KDIAG=13*9, 
/
link_20w
