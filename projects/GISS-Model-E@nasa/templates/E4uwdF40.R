E4uwdF40.R GISS Model E  1850 ocn/atm                     tzhou 03/03/2010
E4uwdF40 : modelE with alternative gravity wave drag: Unresolved Wave Drag
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)

     ! resolution-specific source codes
Atm144x90                  ! horizontal resolution is 144x90 -> 2x2.5deg
AtmL40                      ! vertical resolution is 40 layers -> 0.1mb
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform

IORSF                               ! old i/o

     ! GISS dynamics with alt gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV TQUS_DRV                    ! advection of Q/tracers
! UNRDRAG_COM UNRDRAG UNRDRAG_DRV     ! unresolved wave drag (alternative gravity wave drag)

#include "latlon_source_files"
#include "modelE4_source_files"
#include "static_ocn_source_files"

Components:
Ent shared MPI_Support solvers giss_LSM

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB

Data input files:
    ! resolution dependent files
    ! start up from the restart file of an earlier run ...
! AIC=1....rsfE... ! initial conditions, no GIC needed, use ISTART=8
    ! ... or from observed conditions AIC and model ground data GIC
AIC=AIC.RES_F40.D771201.nc  ! observed init cond (atm. only) ISTART=2
GIC=GIC.144X90.DEC01.1.ext_1.nc   ! initial ground conditions      ISTART=2
! prescr. climatological ocean (1 yr of data)
OSST=OST_144x90.1876-1885avg.HadISST1.1.nc
OSST_eom=OST_144x90.1876-1885avg.HadISST1.1.nc
! prescr. climatological sea ice
SICE=SICE_144x90.1876-1885avg.HadISST1.1.nc
SICE_eom=SICE_144x90.1876-1885avg.HadISST1.1.nc
ZSIFAC=SICE_144x90.1876-1885avg.HadISST1.1.nc
! For q-flux ocean, replace all the above by the next 2 lines, set KOCEAN=1, ISTART=8
!! AIC=1JAN1961.rsfE4F40.MXL65m        ! end of preliminary run with KOCEAN=0
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960   ! ocean horizontal heat transports
OCNML=Z1O.B144x90.nc                    ! mixed layer depth (not used if KOCEAN=0)
CDN=CD144X90.ext.nc
VEG=V144X90_no_crops.ext.nc
CROPS=CROPS2007_144X90N_nocasp.nc
SOIL=S144X900098M.ext.nc
TOPO=Z144X90N_nocasp.nc              ! bdy.cond
REG=REG2X2.5                      ! special regions-diag
RVR=RD_modelE_Fa.nc             ! river direction file
NAMERVR=RD_modelE_Fa.names.txt  ! named river outlets
TOP_INDEX=top_index_144x90_a.ij.ext.nc
Z4var=Z4var144x89                 ! topographic variances multiplied by four
GLMELT=GLMELT_144X90_gas.OCN.nc   ! glacial melt distribution
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
    ! resolution independent files
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800              ! rad.tables and history files
RADN4=LWCorrTables33k              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_SSA=SSA_Koch2008_kg_m2_144x90x20h.nc
TAero_NIT=NIT_Bauer2008_kg_m2_144x90x20_1890-2000h.nc
TAero_OCA=OCA_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_BCA=BCA_Koch2008_kg_m2_144x90x20_1890-2000h.nc
TAero_BCB=BCB_Koch2008_kg_m2_144x90x20_1890-2000h.nc
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-1999.Apr02_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr       ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
#include "rad_144x90_input_files"
GHG=GHG.Mar2009.txt ! use GHG.Jul2009.txt for runs that start before 1850
dH2O=dH2O_by_CH4_monthly
MSU_wts=MSU_SSU_RSS_weights.txt

Label and Namelist:  (next 2 lines)
E4uwdF40 (E4F40 with alternative gravity wave drag: unresolved wave drag)

&&PARAMETERS
! parameters set for choice of ocean model:
KOCEAN=0        ! ocean is prescribed
!! KOCEAN=1        ! ocean is computed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
!  Kvflxo=1     ! saves VFLXO files to prepare for q-flux runs (mkOTSPEC)

variable_lk=1   ! variable lakes

USE_UNR_DRAG=1      !if 1 => SDRAG is turned off and unresolved drag is applied.
                    !if 0 => SDRAG is intact and alternative gwd is not employed.

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.73 ! affects clouds above 850mb w/o MC region;  tune this first to get about 30% high cloud
U00b=1.68 ! affects clouds below 850mb and MC regions; tune this last  to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.57      ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.46    ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
madaer=3    ! updated aerosols
aer_rad_forc=0
cloud_rad_forc=1

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
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
albsn_yr=1850
dalbsnX=.024
o3_yr=-1850
CO2X=1.

variable_orb_par=0
orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD

! parameters that control the Shapiro filter
DT_XUfilter=180. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=180. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=180.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags 0hrly
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
Nssw=2   ! until diurnal diags are fixed, Nssw has to be even
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs
&&END_PARAMETERS

 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default) or earlier
 YEARE=1971,MONTHE=1,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
/
