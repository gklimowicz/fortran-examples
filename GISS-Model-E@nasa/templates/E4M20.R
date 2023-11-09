E4M20.R GISS Model E  2004 modelE                     rar     07/15/2009
!! E4qsM20.R GISS Model E  2004 modelE 65m q-flux ocn      rar     07/15/2009

!! E4qsM20: E4M20 with 65m q-flux ocean ("!!"-lines are for E4qsM20.R)
E4M20: modelE as frozen in July 2009 without gravity wave drag
modelE 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology  (see OSST/SICE)
uses turbulence scheme, simple strat.drag (not grav.wave drag)
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define NEW_IO
#define OLD_BCdalbsn
End Preprocessor Options

Object modules: (in order of decreasing priority)
     ! resolution-specific source codes
Atm72x46                   ! horizontal resolution is 72x46 -> 4x5deg
AtmL20 STRAT_DUM          ! vertical resolution is 20 layers -> 0.1mb
DIAG_RES_M FFT72   

IO_DRV                              ! new i/o

    ! GISS dynamics
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                    ! advection of Q/tracers

#include "latlon_source_files"
#include "modelE4_source_files"
#include "static_ocn_source_files"

Components:
Ent shared MPI_Support solvers giss_LSM dd2d

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB

Data input files:
    ! resolution dependent files
    ! start up from restart file of earlier run
! AIC=1DECxxxx.rsfEyyyy           ! initial conditions (atm./ground), no GIC, ISTART=8
    ! or start up from observed conditions
AIC=AIC.RES_M20A.D771201.nc          ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc ! initial conditions (ground)
! prescr. climatological ocean (1 yr of data)
OSST=OST4X5.B.1876-85avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1876-85avg.Hadl1.1.nc
! prescr. climatological sea ice
SICE=SICE4X5.B.1876-85avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1876-85avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1876-85avg.Hadl1.1.nc
! For q-flux ocean, replace lines above by the next 2 lines & set KOCEAN=1, ISTART=8
!! AIC=1JAN1961.rsfE4M20.MXL65m   ! = end of preliminary run with KOCEAN=0,Kvflxo=1
!! OHT=OTSPEC.E4M20.MXL65m.1956-1960 ! ocean horizontal heat transport
OCNML=Z1O.B4X5.cor.nc                ! mixed layer depth (needed for post processing)
TOPO=Z72X46N.cor4_nocasp.nc       ! topography
SOIL=S4X50093.ext.nc              ! soil bdy.conds
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops.ext.nc ! veg. fractions
CROPS=CROPS2007_72X46N.cor4_nocasp.nc       ! crops history
soil_textures=soil_textures_top30cm
SOILCARB_global=soilcarb_top30cm_4x5.nc
CDN=CD4X500S.ext.nc               ! surf.drag coefficient
REG=REG4X5                        ! special regions-diag
RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets
TOP_INDEX=top_index_72x46_a.ij.ext.nc  ! only used if #define DO_TOPMODEL_RUNOFF
GLMELT=GLMELT_4X5.OCN.nc   ! glacial melt distribution
    ! resolution independent files
#include "rad_input_files"
#include "rad_72x46_input_files"
MSU_wts=MSU_SSU_RSS_weights.txt

Label and Namelist:
E4M20 (ModelE1 4x5, 20 lyrs, 1850 atm/ocn)

&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0 ! ocn is prescribed
!! KOCEAN=1 ! ocn is computed
Kvflxo=0 ! set =1 after spinup to prepare for q-flux run (edit "I")

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

! U00a and U00b below were tuned by Mike Way using:
! net_rad_planet_hemis = 0.0222975
! trnf_toa = -237.4189 
! srnf_toa = 237.4412
! tsurf = 286.001
! plan_alb = 30.44657

U00a=0.52   ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
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

madaer=3        ! updated aerosols
#include "atmCompos_1850_params"

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
cloud_rad_forc=1 ! use =0 to reduce CPU time (diagnostics doubles radiation calls !)
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
nssw=2          ! until diurnal diagn. are fixed, nssw should be even
&&END_PARAMETERS

 &INPUTZ
   YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1 may be = IYEARI (default) or earlier
   YEARE=1961,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
!  suggested settings for E4qsM20
!! YEARI=1901,MONTHI=1,DATEI=1,HOURI=0,
!! YEARE=1931,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=13*0,
!! ISTART=8,IRANDI=0, YEARE=1901,MONTHE=1,DATEE=1,HOURE=1,
 /
