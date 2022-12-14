E001M53.R GISS Model E                                 gas 06/00

WARNING: The boundary conditions used here may not be what you want
         and no tuning has yet been done.
  Please check and see before running
E001M53: 53 layer 4x5 model

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm72x46                   ! horizontal resolution is 72x46 -> 4x4.5deg
AtmL53                      ! vertical resolution is 53 layers -> 0.1mb
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
CLOUDS2 CLOUDS2_DRV CLOUDS_COM        ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                              ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
! ATURB                             ! turbulence in whole atmosphere
DRYCNV                              ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
OCEAN OCNML                         ! ocean modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_M                          ! diagnostics (resolution dependent)
      FFT72                         ! utilities
POUT                                ! post-processing output

Components:
MPI_Support shared

Data input files:
AIC=AIC.RES_M53.D771201.nc
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc       ! initial conditions (ground)
! OHT=OTSPEC.RB399AM12.M250D ! not needed if KOCEAN=0
OCNML=Z1O.B4X5.cor.nc   ! needed only for postprocessing
OSST=OST4X5.B.1946-55avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1946-55avg.Hadl1.1.nc
SICE=SICE4X5.B.1946-55avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1946-55avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1946-55avg.Hadl1.1.nc
CDN=CD4X500S.ext.nc
VEG=V72X46.1.cor2_no_crops.ext.nc
SOIL=S4X50093.ext.nc
TOPO=Z72X46N.cor4_nocasp.nc ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets
ZVAR=ZVAR4X5.nc         ! topographic variation for gwdrag
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
E001M53 (4x5, 53 layer model)
R=00BG/B

&&PARAMETERS
X_SDRAG=.0075,.00075  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.01    ! lin. sdrag above .01mb (top 2 layers) except near poles
PP_SDRAG=1.1   ! lin. sdrag above 1.mb near poles
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
WMAX=1000.     ! maximum wind velocity in sdrag; default=200 when GW drag not used
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000037 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=15.     ! Shear drag coefficient
CMTN=0.25      ! default is 0.5
CDEF=1.5       ! deformation drag coefficient

KOCEAN=0
 
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.85   ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
HRMAX=400.   ! tune this to get rad.equilibrium (range: 100.-1500. meters)

H2ObyCH4=1.  ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the Shapiro filter
DT_XUfilter=180. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=180. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

LMCM=26              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DT=180.,             ! from default: DTsrc=3600.,
NIsurf=4,            ! number of surface time steps

NSUBDD=0        ! saving sub-daily diags
Kvflxo=0        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=1

&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0,  !  from default: IYEAR1=YEARI
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,    KDIAG=13*0,
   YEARE=1950,MONTHE=2,
   ISTART=2,IRANDI=0, YEARE=1950,MONTHE=1,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
/
