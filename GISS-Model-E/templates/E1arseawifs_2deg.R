test3a.R GISS Model E  coupled version          aromanou  02/14/2009

test3a: obio in GISS ocean based on Larissa's E1F40o32.R
   2x2.5x40 layers modelE version, 1850 atm.; 32 layers in the ocean
          NOTE: new ocean initial condition OIC=OIC.WOA98.2HX2.L32.D1201
modelE1 (3.0) 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)       ?
atmospheric composition from year 1880 ? 1979                               ?
ocean: coupled to GISS ocean model (Russell - Schmidt)                      ?
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag) ?
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs  ?
filters: U,V in E-W direction (after every dynamics time step)              ?
         sea level pressure (after every physics time step)                 ?

Run Options
STACKSIZE=524288

Preprocessor Options
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
!#define TRACERS_OCEAN               ! GISS Ocean tracers activated
!#define TRACERS_OCEAN_INDEP         ! independently defined ocn tracers
!#define TRACERS_OceanBiology
!#define OBIO_ON_GISSocean
!#define pCO2_ONLINE
!!!!#define CHL_from_SeaWIFs - replaced by run-time parameter
#define OCN_LAYERING L32
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm144x90                  ! horizontal resolution is 144x90 -> 2x2.5deg
AtmL40                      ! vertical resolution is 40 layers -> 0.1mb
ORES_2Hx2 OFFT144E                  ! ocean horiz res 2x2.5deg
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
ODIAG_COM OCEAN_COM OSTRAITS_COM OGEOM ! dynamic ocean modules
OCNDYN  OCNDYN2  OTIDELL                 ! dynamic ocean routines
OCN_Interp OCN_Int_LATLON                ! dynamic ocean routines
OSTRAITS OCNGM OCNKPP                    ! dynamic ocean routines
OCEANR_DIM AFLUXES OFLUXES
ODIAG_PRT                              ! ocean diagnostic print out
OCNFUNTAB                              ! ocean function look up table
SNOW_DRV SNOW                          ! snow model
RAD_COM RAD_DRV RADIATION              ! radiation modules
RAD_UTILS ALBEDO READ_AERO             ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT          ! diagnostics
DIAG_ZONAL GCDIAGb                     ! grid dependent code for lat-circle dia
DIAG_RES_F                             ! diagnostics (resolution dependent)
      FFT144                           ! utilities
POUT                                   ! post-processing output
SparseCommunicator_mod                 ! sparse gather/scatter module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OCN_TRACER_COM
!OCN_TRACER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!   OCEAN BIOLOGY      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!obio_dim         |$(R8)|
!obio_incom       |$(R8)|
!obio_com         |$(R8)|
!obio_forc        |$(R8)|
!obio_init        |$(R8)|
!obio_bioinit     |$(R8)|
!obio_model       |$(R8)|
!obio_daysetrad   |$(R8)|
!obio_daysetbio   |$(R8)|
obio_ocalbedo    |$(R8)|
obio_reflectance |$(R8)|
!obio_sfcirr      |$(R8)|
!obio_edeu        |$(R8)|
!obio_ptend       |$(R8)|
!obio_carbon      |$(R8)|
!obio_update      |$(R8)|
!obio_sinksettl   |$(R8)|

!!!ar!!!obio_oasimhr     |$(R8)|

Components:
MPI_Support shared

Data input files:
!!! AIC=AIC.RES_F40.D771201.nc         ! observed init cond (atm. only) ISTART=2
AIC=1JAN2501.rsfE3F40o32        ! Larissa's restart              ISTART=8
GIC=GIC.144X90.DEC01.1.ext_1.nc      ! initial ground conditions      ISTART=2
OIC=OIC.E2HX2.L32.D1201.nc      ! Levitus ocean intial conditions
OFTAB=OFTABLE_NEW               ! ocean function table
KBASIN=KB144X90.modelE.nc       ! ocean basin designations
TOPO_OC=OZ2HX2fromZ1QX1N.nc     ! ocean fraction and topography
OSTRAITS=OSTRAITS_144x90.nml    ! parameterized straits info
CDN=CD144X90.ext.nc             ! neutral drag coefficient
VEG=V144X90_no_crops.ext.nc     ! vegatation file
CROPS=CROPS2007_144X90N_nocasp.nc  ! crops
SOIL=S144X900098M.ext.nc        ! soil properties
TOPO=Z2HX2fromZ1QX1N.nc         ! surface fractions and topography
REG=REG2X2.5                    ! special regions-diag
RVR=RD_modelE_F.nc             ! river direction file
NAMERVR=RD_modelE_F.names.txt  ! named river outlets
RADN1=sgpgxg.table8             ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800            ! rad.tables and history files
RADN4=LWCorrTables33k            ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
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
#include "rad_144x90_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_144x90_a.ij.ext.nc
MSU_wts=MSU_SSU_RSS_weights.txt
GLMELT=GLMELT_144X90_gas.OCN.nc   ! glacial melt distribution
!!!!!!!!!!!!!!!!!!! obio  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfle1=abw25b.dat                         ! seawater spectral absorp. and
scatt. coefs
cfle2=acbc25b.dat                        ! phytoplankton spectrl absorp. and
scatt. coefs
!!!pco2table=pco2.tbl.asc                   ! table to compute pco2 vals from
sst,sss,dic,alk
                                            ! if not defined pCO2_ONLINE
nitrates_inicond=no3_nodc_annmean_90x144.nc    ! initial cond for nitrates (NODC)
silicate_inicond=sio2_nodc_annmean_90x144.nc   ! initial cond for silicate (NODC)
dic_inicond=dic_glodap_annmean_90x144.nc       ! initial cond for dic (GLODAP)
alk_inicond=alk_glodap_annmean_90x144.nc       ! initial cond/forcing for alk
(GLODAP)
!!!oasimdirect=oasimdirect_20w_new          ! spectral light components
                                            ! if not defined
atmFe_inicond=iron_gocart_1x1mon_90x144.nc     ! GOCART iron flux
atmFedirect1=iron_ron_195x180_20w.asc    ! Ron Miller's dust fluxes
facirr=facirr.asc                        ! factors for mean irradiance w/in
water
eda_esa_ratios=eda_esa_ratios.asc        ! ratios of radiation spectral
components
!!!!!!!!!!!!!!!!!!! obio_rad  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHL_DATA=CHL_WG_2x2.5                    !CHL_WG_4x5 or CHL_WG_2x2.5
                                         !in GISS ocean grid
                                         !to be used with CHL_from_SeaWIFs

Label and Namelist:
test3a (32 ocean layers; 1850 atm.,the current modelE version)

DTFIX=180
&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1        ! ocn is prognostic
variable_lk=1
init_flake=1

! parameters usually not changed when switching to coupled ocean:

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38  39  40
vsdragl=0.021,0.041,0.077,0.125,0.22,0.275,0.276,0.447,0.96,0.92,  0.91,1.22,1.53,0.3,0.6,0.83, 1.
ANG_SDRAG=1         ! conserve ang. mom.

OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
OTIDE = 0           !  Ocean tides are not used

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)


U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.60       ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.47     ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

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
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
DTsrc=1800.
DT=225.
DTO=225.
NIsurf=2        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480       ! use =48 except on halem
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time if isccp-diags are not essential
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
nssw=48         ! obio needs that in order to always restart from hour 0
                ! then we need to do setups for a whole day

!parameters that affect CO2 gas exchange
atmCO2=368.6      !uatm for year 2000

chl_from_seawifs = 1
&&END_PARAMETERS

 &INPUTZ
   YEARI=2501,MONTHI=1,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=2521,MONTHE=1,DATEE=1,HOURE=0, KDIAG=13*0,
   ISTART=8,IRANDI=0, YEARE=2501,MONTHE=1,DATEE=2,HOURE=0,IWRITE=1,JWRITE=1,
/
