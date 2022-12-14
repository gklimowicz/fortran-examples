E001C12.R GISS Model E                                 gas 06/00

WARNING: The boundary conditions used here may not be what you want
         and only preliminary tuning has been done.
  Please check and see before running
E001C12: modelE1 (3.0) 8x10

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm36x24                          ! horizontal resolution is 36x24 -> 8x10
AtmL12 STRAT_DUM                   ! vertical resolution is 12 layers -> 10mb
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
! pick exactly one of the next 2 choices: ATURB or DRYCNV
ATURB                               ! turbulence in whole atmosphere
! DRYCNV                            ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
OCEAN OCNML                         ! ocean modules
ICEDYN_DUM                          ! dummy ice dynamics
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_C                          ! diagnostics (resolution dependent)
      FFT36                         ! utilities
POUT                                ! post-processing output

Components:
MPI_Support shared

Data input files:
! AIC=DEC1958.rsfB394M12.modelE.16 ! model init cond (atm. and ground) ISTART=7
AIC=AIC.RES_C12.D771201.nc   ! observed init cond   (atm. only)       ISTART=2
GIC=GIC.8X10.modelE_1.nc     ! initial ground conditions              ISTART=2
! OHT=OTSPEC.RunIDM12.M250D ! hor.heat transp.  for q-flux ocean only
! OCNML=Z1O.B4X5.cor.nc        ! mixed layer depth,needed for post-processing only
! prescr. climatological ocean (1 yr of data)
OSST=OST8X10.B.1946-55avg.Hadl1.1.nc
OSST_eom=OST8X10.B.1946-55avg.Hadl1.1.nc
! prescr. climatological sea ice
SICE=SICE8X10.B.1946-55avg.Hadl1.1.nc
SICE_eom=SICE8X10.B.1946-55avg.Hadl1.1.nc
ZSIFAC=SICE8X10.B.1946-55avg.Hadl1.1.nc
!
CDN=CD8X10.modelE.nc
VEG=V8X10.modelE.nc
SOIL=S8X10.modelE.nc
TOPO=Z8X10.modelE.nc ! bdy.cond
REG=REG8X10          ! special regions-diag
RVR=RD8X10.nc                       ! river direction file
NAMERVR=RD8X10.names.txt            ! named river outlets
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
RADN9=solar.lean02.ann.uvflux_hdr  ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
#include "rad_36x24_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_8x10.ij.nc
MSU_wts=MSU_SSU_RSS_weights.txt

Label and Namelist:
E001C12 (modelE1 (3.0) 8x10)

DTFIX=450
&&PARAMETERS
X_SDRAG=.00025,.000025
C_SDRAG=0.
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
do_polefix=0    ! polefix enhancements not yet valid for real 8x10 resolution

KOCEAN=0
 
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.74   ! tune this first to get reas.alb/cldcvr
U00wtrX=1.

! parameters that control the Shapiro filter
DT_XUfilter=900. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=900. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DT=900.,        ! from default: DTsrc=3600.,
H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

SUBDD=' '     ! save SLP at sub-daily frequency
NSUBDD=0       ! saving sub-daily diags 12hrly
Kvflxo=0        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=1
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=12,DATEI=1,HOURI=0, ! to be used with ISTART=2
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0, KDIAG=13*0,
   YEARE=1951,MONTHE=2,
   ISTART=2,IRANDI=0, YEARE=1950,MONTHE=12,DATEE=1,HOURE=1,
/
