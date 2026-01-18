E1oCS90L40.R GISS Model E CS atmosphere coupled ocean     Denis   09/24/2009

E1oCS90L40:  CS90 cubed sphere with 1Qx1 32 layers dynamic ocean

E1oCS90L40: replace this section by a description of what distinguishes this run ?
          Use as many lines as you need. Look carefully at all the possible    ?
          choices, particularly the lines containing '?'. In some cases, you   ?
          will have to pick the appropriate choice to make this rundeck work   ?
          The final rundeck should contain no '?'
          Check and modify the rest of the description below:                  ?
ocean: coupled to GISS ocean model (Russell - Schmidt)                      ?
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag) ?
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs  ?
filters: U,V in E-W direction (after every dynamics time step)              ?
         sea level pressure (after every physics time step)                 ?

Run Options
STACKSIZE=524288

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
#define OCN_LAYERING L32
#define AG2OG_PRECIP_BUNDLE
#define OG2AG_TOC2SST_BUNDLE
#define AG2OG_OCEANS_BUNDLE
#define OG2AG_OCEANS_BUNDLE
#define BUNDLE_INTERP
#define NEW_IO
End Preprocessor Options

Object modules: (in order of decreasing priority)
AtmCS90                           ! 90 Cube-Sphere Grid
AtmL40                             ! vertical resolution is 40 layers -> 0.1mb
FFTW_COM          
ORES_1Qx1 OFFT288E                  ! ocean horiz res 1.25x1deg
MODEL_COM GNOM_CS IO_DRV
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
dd2d_utils pario_nc pario_fbsa
regrid regrid_com 
ALLOC_DRV			    ! allocate global distributed arrays
ATMDYN_COM ATM_DUM
ATM_UTILS                           ! utilities for some atmospheric quantities
FV_UTILS FV_CS_Mod FV_INTERFACE     ! FV dynamical core wrapper
QUS_COM QUSDEF              ! advection of tracers
QUScubed
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H                 ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN ICEDYN_DRV                       ! ice dynamics modules
ADVSIcs
ODIAG_COM OCEAN_COM OSTRAITS_COM OGEOM  ! dynamic ocean modules
OCNDYN  OCNDYN2  OTIDELL                     ! dynamic ocean routines
OCN_Interp OSTRAITS OCNGM OCNKPP             ! dynamic ocean routines
OCEANR_DIM AFLUXES OFLUXES
ODIAG_PRT                              ! ocean diagnostic print out
OCNFUNTAB                              ! ocean function look up table
SNOW_DRV SNOW                          ! snow model
RAD_COM RAD_DRV RADIATION COSZ_2D   ! radiation modules
RAD_UTILS ALBEDO READ_AERO             ! radiation and albedo
DIAG_COM DIAG DEFACC QUICKPRT          ! diagnostics
DIAG_RES_F
DIAG_ZONALcs
GCDIAGcs cs2ll_utils
POUT                                   ! post-processing output
SparseCommunicator_mod                 ! sparse gather/scatter module

Components:
MPI_Support shared dd2d CS_Support

Data input files:
AIC=AIC_CS90.nc      ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC_CS90_1.nc      ! initial conditions (ground)
OIC=OIC288X180.D1201.nc      ! Levitus ocean intial conditions
TOPO_OC=OZ288X180N.nc        ! ocean fraction and topography
OFTAB=OFTABLE_NEW               ! ocean function table
KBASIN=KB288X180.modelE.nc   ! ocean basin designations
OSTRAITS=OSTRAITS_288x180.nml ! parameterized straits info
TOPO=Z_C90fromZ1QX1N.nc
SOIL=SOIL_CS90.nc ! soil/topography bdy.conds
VEG=V_CS90_144X90_5percent.nc
CROPS=CROPS_CS90.nc
CDN=CD_CS90.nc                     ! surf.drag coefficient
REG=REG.txt                        ! special regions-diag
RVR=RDdistocean_CS90_EM.nc             ! river direction file
NAMERVR=RDdistocean_CS90_EM.names.txt  ! named river outlets
TOP_INDEX=top_index_CS90.nc      ! only used if #define DO_TOPMODEL_RUNOFF
GLMELT=GLMELT_CS90.nc  ! glacial melt distribution
REMAP=remap288-180C90-90.nc
RADN1=sgpgxg.table8             ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800            ! rad.tables and history files
RADN4=LWCorrTables33k            ! rad.tables and history files
RADN5=H2Ocont_Ma_2000           ! rad.tables and history files
RADN3=miescatpar.abcdv2
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h
TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr     ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
#include "rad_C90_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
MSU_wts=MSU_SSU_RSS_weights.txt

Label and Namelist:
E1oCS90L40 (32 ocean layers; 1850 atm.,the current modelE version)

DTFIX=180
&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1        ! ocn is prognostic
variable_lk=0
init_flake=1

! parameters usually not changed when switching to coupled ocean:

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.004,.0004  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0004       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_SDRAG=1         ! conserve ang. mom.

OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
OTIDE = 0           !  Ocean tides are not used

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

 
U00a=0.72   ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.40   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.60       ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.47     ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
madaer=3    ! updated aerosols

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1850  ! if -1, crops in VEG-file is used
s0_yr=1850
s0_day=182
ghg_yr=1850
ghg_day=182
volc_yr=-1
volc_day=182
aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1850
dalbsnX=.015
o3_yr=-1850

! parameters that control the Shapiro filter
DT_XUfilter=112.5 ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=112.5 ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
DTsrc=1800.
DT=1800.
DTO=112.5
NIsurf=2        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=960       ! use =48 except on halem
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time if isccp-diags are not essential
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
&&END_PARAMETERS

 &INPUTZ
   YEARI=1900,MONTHI=12,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=1951,MONTHE=12,DATEE=1,HOURE=0, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1900,MONTHE=12,DATEE=1,HOURE=1,
/
