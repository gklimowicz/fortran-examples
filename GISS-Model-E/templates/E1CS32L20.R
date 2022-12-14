E1CS32L20.R GISS Model E  2004 modelE             your_name  1/13/07

E1CS32L20: description/motivation    (any number of lines) ?
   (delete)  Lines you may want to inspect in any case contain a "?"
   (delete)  At the end you'll find instructions on how to modify this rundeck
   (delete)  for other simple ocean models and change in vert/hor resolution
modelE C32 cubed-sphere grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1880 (or 1979)      (look below for "_yr")  ?
ocean data: prescribed, 1876-1885 (or 1975-1984) climatology  (see OSST/SICE) ?
uses turbulence scheme, simple strat.drag (not grav.wave drag)
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
filters:    U,V in E-W direction (after every dynamics time step)
            sea level pressure (after every physics time step)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define NEW_IO
End Preprocessor Options

Run Options
STACKSIZE=524288

Object modules: (in order of decreasing priority)
AtmCS32                           ! 32 Cube-Sphere Grid
AtmL20 STRAT_DUM                 ! vertical resolution is 20 layers -> 0.1mb
DIAG_RES_M FFTW_COM       
MODEL_COM GEOM_CS IO_DRV            ! model variables and geometry
!GNOM_CS                            ! GNOMONIC cubed sphere geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
pario_fbsa
regrid regrid_com
ALLOC_DRV                           ! allocate global distributed arrays
ATMDYN_COM ATM_DUM !MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
FV_UTILS FV_CS_Mod FV_INTERFACE     ! FV dynamical core wrapper
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
QUScubed                            ! cubed-sphere adaptation of QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN_DUM
!ICEDYN_DRV ICEDYN                   ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION COSZ_2D   ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC QUICKPRT       ! diagnostics
DIAG_ZONALcs GCDIAGcs cs2ll_utils   ! grid-dependent code for lat-circle diags

Components:
MPI_Support shared dd2d CS_Support

Data input files:
    ! start up from restart file of earlier run
! AIC=1DECxxxx.rsfEyyyy           ! initial conditions (atm./ground), no GIC, ISTART=8
! AICfv=1DECxxxx.fvEyyyy          ! initial conditions (fv internal) only for ISTART=9
! AICdfv=1DECxxxx.dfvEyyyy        ! tendencies                       only for ISTART=9
    ! or start up from observed conditions
AIC=AIC_CS32.nc      ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC_CS32_1.nc      ! initial conditions (ground)
    ! ocean data for "prescribed ocean" runs : climatological ocean
! prescr. climatological ocean (1 yr of data)
OSST=OST_CS32.nc
OSST_eom=OST_CS32.nc
! prescr. climatological sea ice
SICE=SICE_CS32.nc
SICE_eom=SICE_CS32.nc
ZSIFAC=SICE_CS32.nc
OCNML=Z1O.B4X5.cor                ! mixed layer depth (needed for post processing)
!                                             (end of section 1 of data input files)
    ! resolution dependent files
TOPO=Z_CS32_4X5.nc
SOIL=SOIL_CS32.nc ! soil/topography bdy.conds
! VEG=V72X46.1.cor2   ! or:      ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V_CS32_144X90_1percent.nc
CROPS=CROPS_CS32 
CDN=CD_CS32.nc                   ! surf.drag coefficient
REG=REG.txt                      ! special regions-diag
RVR=RDdistocean_CS32.nc             ! river direction file
NAMERVR=RDdistocean_CS32.names.txt  ! named river outlets
TOP_INDEX=top_index_CS32.nc      ! only used if #define DO_TOPMODEL_RUNOFF
GLMELT=GLMELT_CS32.nc   ! glacial melt distribution
!                                             (end of section 2 of data input files)
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
RADN9=solar.lean02.ann.uvflux_hdr      ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
#include "rad_C32_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
MSU_wts=MSU_SSU_RSS_weights.txt

Label and Namelist:
E1CS32L20 (ModelE1 Cubed Sphere C32, 20 lyrs, 1880 atm/ocn;
up to 60 (or 52) columns here to describe your run)?<- col 53  to  72 ->   80 ->
DTFIX=300

&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0 ! 0 or 1 , use =0 if ocn is prescribed, use =1 if ocn is predicted
Kvflxo=0 ! use 1 ONLY to save VFLXO daily to prepare for q-flux run ?

variable_lk=1 ! let lakes grow or shrink in horizontal extent
init_flake=0

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

! tuning param.: this setting works for both 1880 and 1979
U00ice=.63      ! U00ice+.01 =>dBal=1.5,dPl.alb=-.9%   goals:Bal=0,plan.alb=30%
U00wtrX=1.34    ! U00wtrX+.01=>dBal=0.7,dPl.alb=-.25%  Bal=glb.ann NetHt at z0
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=-1 no land ice fixup, 1 Lacis' scheme)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1880 ! if -1, crops in VEG-file is used   ! =1979 , also change OSST,SICE
s0_yr=1880                                         ! =1979 , also change OSST,SICE
s0_day=182
ghg_yr=1880                                        ! =1979 , also change OSST,SICE
ghg_day=182
volc_yr=1880                                       ! =1979 , also change OSST,SICE
volc_day=182
aero_yr=1880                                       ! =1979 , also change OSST,SICE
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1880                                      ! =1979 , also change OSST,SICE
dalbsnX=.024
o3_yr=-1880                                        ! =1979 , also change OSST,SICE

! parameters that control the Shapiro filter
DT_XUfilter=0. ! Shapiro filter on U in E-W direction; (not used)
DT_XVfilter=0. ! Shapiro filter on V in E-W direction; (not used)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

MFILTR=0 ! no slp-filter needed for fv dynamics

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=1800.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf  ? =3 to also save "oda"-files
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
cloud_rad_forc=0 ! use =1 to activate this diagnostic (doubles radiation calls !)
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
nssw=2          ! until diurnal diagn. are fixed, nssw should be even
&&END_PARAMETERS

 &INPUTZ
   YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, IYEAR1=1949 ! or earlier
   YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
 /
