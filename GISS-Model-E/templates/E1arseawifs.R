E1arseawifs.R GISS Model E  2004 modelE           aromanou     8/07/08

E1arseawifs: spin-off E1M20.R
             modelE + prescribed SST + seawifs chlorophyl
             test albedo-chlorophyl interaction

modelE 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850 (or 1979)      (look below for "_yr")  ?
ocean data: prescribed, 1876-1885 (or 1975-1984) climatology  (see OSST/SICE) ?
uses turbulence scheme, simple strat.drag (not grav.wave drag)
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
filters:    U,V in E-W direction (after every dynamics time step)
            sea level pressure (after every physics time step)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
!!!!!#define CHL_from_SeaWIFs            ! read in SeaWIFs - repaled by run-time parameer
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm72x46                          ! horizontal resolution is 72x46 -> 4x5deg
AtmL20 STRAT_DUM                 ! vertical resolution is 20 layers -> 0.1mb
DIAG_RES_M FFT72 
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
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
                                    ! utilities
POUT                                ! post-processing output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
obio_ocalbedo|$(R8)|                  ! changes to ocean albedo due to chlorophyl

Components:
MPI_Support shared

Data input files:
    ! start up from restart file of earlier run
! AIC=1DECxxxx.rsfEyyyy           ! initial conditions (atm./ground), no GIC, ISTART=8
    ! or start up from observed conditions
AIC=AIC.RES_M20A.D771201.nc          ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc   ! initial conditions (ground)
    ! ocean data for "prescribed ocean" runs : climatological ocean
! prescr. climatological ocean (1 yr of data)
OSST=OST4X5.B.1876-85avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1876-85avg.Hadl1.1.nc
! prescr. climatological sea ice
SICE=SICE4X5.B.1876-85avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1876-85avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1876-85avg.Hadl1.1.nc
!? for 1979 OSST=OST4X5.B.1975-84avg.Hadl1.1
!? for 1979 SICE=SICE4X5.B.1975-84avg.Hadl1.1
OCNML=Z1O.B4X5.cor.nc                ! mixed layer depth (needed for post processing)
!                                             (end of section 1 of data input files)
    ! resolution dependent files
TOPO=Z72X46N.cor4_nocasp.nc
SOIL=S4X50093.ext.nc ! soil/topography bdy.conds
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops.ext.nc
CROPS=CROPS2007_72X46N.cor4_nocasp.nc  ! veg. fractions, crops history
CDN=CD4X500S.ext.nc               ! surf.drag coefficient
REG=REG4X5                        ! special regions-diag
RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets
TOP_INDEX=top_index_72x46_a.ij.ext.nc  ! only used if #define DO_TOPMODEL_RUNOFF
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
#include "rad_72x46_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
MSU_wts=MSU_SSU_RSS_weights.txt
GLMELT=GLMELT_4X5.OCN.nc   ! glacial melt distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHL_DATA=CHL_WG_4x5                      !CHL_WG_4x5 in GISS ocean grid
                                         !to be used with CHL_from_SeaWIFs

Label and Namelist:
E1arseawifs (ModelE1 4x5, 20 lyrs, 1850 atm/ocn; use up to 72 (or 80) columns and ??
up to 60 (or 52) columns here to describe your run)?<- col 53  to  72 ->   80 ->
DTFIX=300

&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0 ! 0 or 1 , use =0 if ocn is prescribed, use =1 if ocn is predicted
Kvflxo=0 ! use 1 ONLY to save VFLXO daily to prepare for q-flux run ?

variable_lk=1 ! let lakes grow or shrink in horizontal extent
wsn_max=2.   ! restrict snow depth to 2 m-h2o (if 0. snow depth is NOT restricted)

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

! tuning param.: this setting works for 1850; use U00wtrX=1.28 for 1979
 
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.62      ! U00ice+.01 =>dBal=1.5,dPl.alb=-.9%   goals:Bal=0,plan.alb=30%
U00wtrX=1.29    ! U00wtrX+.01=>dBal=0.7,dPl.alb=-.25%  Bal=glb.ann NetHt at z0
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=-1 no land ice fixup, 1 Lacis' scheme)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1850 ! if -1, crops in VEG-file is used   ! =1979 , also change OSST,SICE
s0_yr=1850                                         ! =1979 , also change OSST,SICE
s0_day=182
ghg_yr=1850                                        ! =1979 , also change OSST,SICE
ghg_day=182
volc_yr=-1  ! 1850-1999 mean strat.aeros           ! =1979 , also change OSST,SICE
volc_day=182
aero_yr=-1850                                       ! =-1979 , also change OSST,SICE
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1850                                      ! =1979 , also change OSST,SICE
dalbsnX=.024
o3_yr=-1850                                        ! =-1979 , also change OSST,SICE

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=450.
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

chl_from_seawifs=1
&&END_PARAMETERS

 &INPUTZ
   YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, IYEAR1=1949 ! or earlier
   YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
 /

