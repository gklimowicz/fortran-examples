E1oM20agnat.R GISS Model E  coupled version          larissa   01/22/2009

E1oM20agnat: replace this section by a description of what distinguishes this run ?
       Use as many lines as you need. Look carefully at all the possible    ?
       choices, particularly the lines containing '?'. In some cases, you   ?
       will have to pick the appropriate choice to make this rundeck work   ?
       The final rundeck should contain no '?'
       Check and modify the rest of the description below:                  ?
modelE1 (3.0) 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)       ?
atmospheric composition from year 1880 ? 1979                               ?
ocean: coupled to GISS ocean model (Russell - Schmidt)                      ?
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag) ?
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs  ?
filters: U,V in E-W direction (after every dynamics time step)              ?
         sea level pressure (after every physics time step)                 ?

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
#define TRACERS_OCEAN               ! GISS Ocean tracers activated
#define TRACERS_AGE_OCEAN           ! Natassa's Ocean Age tracer code
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm72x46                          ! horizontal resolution is 72x46 -> 4x5deg
AtmL20 STRAT_DUM                 ! vertical resolution is 20 layers -> 0.1mb
RES_5x4_L13                         ! ocean horiz res 4x5deg, 13 vert layers  
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
              ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM        ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
ODIAG_COM OCEAN_COM OSTRAITS_COM OGEOM   ! dynamic ocean modules
OCNDYN  OCNDYN2  OTIDELL                 ! dynamic ocean routines
OCN_Interp OCN_Int_LATLON                ! dynamic ocean routines
OSTRAITS OCNGM OCNKPP                    ! dynamic ocean routines
OCEANR_DIM AFLUXES OFLUXES
ODIAG_PRT                              ! ocean diagnostic print out
OCNFUNTAB                           ! ocean function look up table
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_M                          ! diagnostics (resolution dependent)
      FFT72 OFFT72E                 ! utilities
POUT                                ! post-processing output
SparseCommunicator_mod              ! sparse gather/scatter module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  tracer part  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TRACER_COM TRACERS_DRV              ! configurable tracer code
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
OCN_TRACER
OCN_TRACER_COM
!!!!!!!!!!!!!!!!!!!!!!!!!!!  tracer part  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Components:
tracers MPI_Support shared

Data input files:
! AIC=AIC.RES_M20A.D771201           ! initial conditions (atm.)     needs GIC,OIC ISTART=2
! GIC=GIC.E046D3M20A.1DEC1955        ! initial conditions (ground)         and 300 year spin-up
! OIC=OIC4X5LD.Z12.gas1.CLEV94.DEC01 ! ocean initial conditions
! AIC=1JAN2012.rsfE051oM20A            ! full IC (GIC,OIC not needed) ISTART=8 (spun up 380 yrs)
AIC=AIC/1JAN2600.rsfE13AoM20           ! istart=5
OFTAB=OFTABLE_NEW                    ! ocean function table
KBASIN=KB4X513.OCN.gas1.nc           ! ocean basin designations
TOPO_OC=OZ72X46N_gas.1_nocasp.nc     ! ocean bdy.cond
OSTRAITS=OSTRAITS_72x46.nml          ! parameterized straits info
CDN=CD4X500S.ext.nc
  ! VEG=V72X46.1.cor2.ext
VEG=V72X46.1.cor2_no_crops.ext.nc
CROPS=CROPS2007_72X46N.cor4_nocasp.nc  ! veg. fractions, crops history
SOIL=S4X50093.ext.nc
TOPO=Z72X46N_gas.1_nocasp.nc ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets
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
E1oM20agnat (coupled version, 4x5, 20 lyrs, 1880 atm; use up to 72 (or 80) columns and ??
up to 60 (or 52) columns here to describe your run)?<- col 53  to  72 ->   80 ->
DTFIX=300
&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1        ! ocn is prognostic
OTIDE = 0       !  Ocean tides are not used

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
U00wtrX=1.39    ! U00wtrX+.01=>nethtz0+.7                      for global annual mean
!?1979 U00wtrX=1.38
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1880 ! if -1, crops in VEG-file is used  ? 1979
s0_yr=1880    !? 1979
s0_day=182
ghg_yr=1880   !? 1979
ghg_day=182
volc_yr=1880  !? or -1 to get mean volc.aerosols
volc_day=182
aero_yr=1880  !? 1979
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1880 !? 1979
dalbsnX=.015       ! should be .024
o3_yr=-1880    !? 1979

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
   YEARI=1901,MONTHI=1,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=2000,MONTHE=1,DATEE=2,HOURE=0, KDIAG=13*0,
   ISTART=5,IRANDI=0, YEARE=1901,MONTHE=1,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
/

! Instructions for related rundeck types
! ======================================
! the "frozen (or 'slush') version" of 2006 paper E1oM20 -> EofzM20
!     -------------------------------------------           =======
!     replace in "Object modules"     the 4 files
! CLOUDS2    PBL    ATURB    RADIATION    RAD_DRV         by:
! CLOUDS2_E1 PBL_E1 ATURB_E1 RADIATION_E1 RAD_DRV_E1
!     replace in "Data input files:" RADN2            by:
! RADN2=radfil33k  ! RADN4 and RADN5 are not used
!     set in &&PARAMETERS : variable_lk=0 ! lake fractions are fixed in time
!                           river_fac=1.04
!                           wsn_max=0.    ! do not restrict snow depth
!                           glmelt_on=2   ! skip annual adjustment of glacial melt
!                           glmelt_fac_nh=2.91
!                           glmelt_fac_sh=1.98
!                           oBottom_drag=0
!                           oCoastal_drag=0
