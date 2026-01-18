E001M23.R GISS Model E                                 gas 06/00

E001M23: modelE1 (3.0) (strat. version with gravity wave drag)
Balanced for 1880 SST + forcings

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm72x46                   ! horizontal resolution is 72x46 -> 4x5deg
AtmL23                      ! vertical resolution is 23 layers -> 0.1mb
DIAG_RES_M FFT72 
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
              ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
ATURB                               ! turbulence in whole atmosphere
! DRYCNV                            ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN  ! or: ICEDYN_DUM ! dynamic ice modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
                                    ! utilities
POUT                                ! post-processing output

Components:
MPI_Support shared

Data input files:
! AIC=AIC.RES_M23.D771201.nc
! GIC=GIC.E046D3M20A.1DEC1955.ext
AIC=1JAN1881.rsfE001M23
OSST=OST4X5.B.1876-85avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1876-85avg.Hadl1.1.nc
SICE=SICE4X5.B.1876-85avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1876-85avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1876-85avg.Hadl1.1.nc
OCNML=Z1O.B4X5.cor.nc   ! needed only for postprocessing
!                                             (end of section 1 of data input files)
TOPO=Z72X46N.cor4_nocasp.nc
SOIL=S4X50093.ext.nc              ! bdy.cond
VEG=V72X46.1.cor2_no_crops.ext.nc
CROPS=CROPS2007_72X46N.cor4_nocasp.nc
CDN=CD4X500S.ext.nc
REG=REG4X5           ! special regions-diag
RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets
TOP_INDEX=top_index_72x46_a.ij.ext.nc
ZVAR=ZVAR4X5.nc         ! topographic variation for gwdrag
!                                             (end of section 2 of data input files)
RADN1=sgpgxg.table8    ! rad.tables
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800              ! rad.tables and history files
RADN4=LWCorrTables33k              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_Ma_2008
RADN3=miescatpar.abcdv2
! RADNA,RADNB are no longer used
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr    ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
#include "rad_72x46_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
MSU_wts=MSU_SSU_RSS_weights.txt
GLMELT=GLMELT_4X5.OCN.nc   ! glacial melt distribution

Label and Namelist:
E001M23 (ModelE1 (3.0) based on B402A - strat. version - 1880 conditions)
R=00BG/B
DTFIX=300
&&PARAMETERS
X_SDRAG=.0005,.00005  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.01     ! lin. sdrag above p_sdrag mb (top layer for M23) except near poles
PP_SDRAG=4.6   ! lin. sdrag above  pp_sdrag mb near poles (top 5 layers for M23)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
WMAX=1000.     ! maximum wind velocity in sdrag; default=200 when GW drag not used

PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000030 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=1.      ! Shear drag coefficient
CMTN=0.25      ! default is 0.5
CDEF=1.5       ! deformation drag coefficient

KOCEAN=0
Kvflxo=0        ! saving VFLXO (daily)

variable_lk=1
wsn_max=2.

 
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice  = .62  ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
u00wtrx = 1.29  ! 1880 conditions

cond_scheme=2  ! more elaborate conduction scheme (GHY, Nancy Kiang)
H2ObyCH4=1.    ! activates strat.H2O generated by CH4
KSIALB=0       ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DTsrc = 1800.        ! half-hour physics time step (default: DTsrc=3600.)
DT=450.,             ! dynamic time step
NIsurf=1,            ! number of surface time steps

! The next 4 parameters correspond to the half-hour physics time step
NDISK = 480,
ndaa = 13,
nda5k = 13,
nda4  = 48,

NSUBDD=0        ! saving sub-daily diags
KCOPY=2         ! saving acc + rsf
isccp_diags=1

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1880 ! if -1, crops in VEG-file is used
s0_yr=1880
ghg_yr=1880
ghg_day=182
s0_day=182
volc_yr=1880
volc_day=182
aero_yr=1880
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1880
dalbsnX=.024
o3_yr=-1880
&&END_PARAMETERS

 &INPUTZ
   YEARI=1880,MONTHI=1,DATEI=1,HOURI=0,  !  from default: IYEAR1=YEARI
   YEARE=1890,MONTHE=1,DATEE=1,HOURE=0,    KDIAG=0,2,2,9*0,9,
   ISTART=8,IRANDI=0, YEARE=1880,MONTHE=1,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
/

! Alternate simple ocean parameterizations  (all require a preliminary run with
! ========================================  specified ocean data, here: E001M23)
! q-flux run based on E001M23 with 65m ocn (sensitivity runs) E001M23 -> E001qsM23
! --------------------------------------                                 =========
!     need last 5-10 yrs of VFLXO-files and final rsf from E001M23
!     replace section 1 of " Data input files" by the 3 lines:
! AIC=1JAN1961.rsfE001M23.MXL65m      ! made by aux/mkOTSPEC  (65m)
! OHT=OTSPEC.E001M23.MXL65m.1956-1960 ! made by aux/mkOTSPEC  (65m)
! OCNML=Z1O.B4X5.cor              ! mixed layer depth (now needed)
!     set in &&PARAMETERS : KOCEAN=1
!     replace the the namelist &INPUTZ by (e.g.)
! &INPUTZ
!   YEARI=1901,MONTHI=1,DATEI=1,HOURI=0,
!   YEARE=1931,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=13*0,
!   ISTART=8,IRANDI=0, YEARE=1901,MONTHE=1,DATEE=1,HOURE=1,
! /

! q-flux run based on E001M23 with 250m ocean                 E001M23 -> E001qM23
! -----------------------------------------                              ========
!     need last 5-10 yrs of VFLXO-files and final rsf from E001M23
!     replace section 1 of " Data input files" by the 3 lines:
! AIC=1JAN1961.rsfE001M23.MXL250m      ! made by aux/mkOTSPEC  (250m)
! OHT=OTSPEC.E001M23.MXL250m.1956-1960 ! made by aux/mkOTSPEC  (250m)
! OCNML=Z1O.B4X5.cor              ! mixed layer depth (now needed)
!     set in &&PARAMETERS : KOCEAN=1 , KCOPY=3  (to create *.odaE001qM23 files)
!     replace the the namelist &INPUTZ by
! &INPUTZ
!   YEARI=1901,MONTHI=1,DATEI=1,HOURI=0,
!   YEARE=2001,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=13*0,
!   ISTART=8,IRANDI=0, YEARE=1901,MONTHE=1,HOURE=1,
! /

! q-flux run with diffusion into deep ocean (needs E001qM23)  E001qM23 -> E001qdM23
! -----------------------------------------                               =========
!     need last 10 yrs of oda-files and final rsf from E001qM23 (KCOPY=3)
!     replace in "Object modules" :                          OCNML -> ODEEP
!     replace section 1 of "Data input files" by the 5 lines:
! AIC=1JAN2001.rsfE001qM23
! OHT=OTSPEC.E001M23.MXL250m.1951-1960
! OCNML=Z1O.B4X5.cor              ! mixed layer depth (now needed)
! TG3M=TG3M.E001qM23  ! made by E001M23_bin/mkdeep.exe (gmake auxdeep ...)
! EDDY=ED4X5 ! eddy diffusion for mixing into deep ocean
!     set in &&PARAMETERS : KOCEAN=1 , KCOPY=2
!     replace the the namelist &INPUTZ by
! &INPUTZ
!   YEARI=1901,MONTHI=1,DATEI=1,HOURI=0,
!   YEARE=2001,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=13*0,
!   ISTART=4,IRANDI=0, YEARE=1901,MONTHE=1,HOURE=1,
! /
! Note: may start also from an earlier qdM23-model run with ISTART=8
