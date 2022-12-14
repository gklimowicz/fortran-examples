E5TgasM20.R GISS Model E  2002 modelE                 rar 12/01/03
 
E025TgasM20: water isotope tracers T1=water ; T2=O18 ; T3=HDO
now uses supsatfac=0.004
modelE3 with 20 lyrs, top at .1 mb - 1980 atmosphere/ocean
no gravity wave drag;     uses turbulence (not dry convection)
Sdrag: weak linear strat. drag in top layer, near poles down to 20 mb
       ang.mom loss is added in below 150 mb
sea level pressure filter applied every hour, UV-filter used
6-band oice albedo; Hogstrom(1984) pbl drag
Note: Some of these choices may be changed using the PARAMETERs below.
 
Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_WATER               ! include water tracers code
! #define TRACERS_DRYDEP              ! include tracer dry deposition
! #define TRACERS_COSMO
#define TRACERS_SPECIAL_O18         ! include water isotope code
#define NUDGE_ON
End Preprocessor Options
 
Run Options
STACKSIZE=262144
                                                                                                                                                                
Object modules: (in order of decreasing priority)
Atm72x46                          ! horizontal resolution is 72x46 -> 4x5deg
AtmL20 STRAT_DUM                 ! vertical resolution is 20 layers -> 0.1mb
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
              ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
NUDGE
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
TRACER_COM TRACERS_DRV              ! configurable tracer code
TRACERS                             ! generic tracer code
TRACERS_O18                         ! special tracer code for water isotopes
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
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
ICEDYN_DRV ICEDYN  ! or: ICEDYN_DUM ! dynamic sea ice modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_M                          ! diagnostics (resolution dependent)
      FFT72                         ! utilities
POUT                                ! post-processing output
 
Components:
tracers MPI_Support shared

Data input files:
    ! the first group of files is specific to prescribed ocean runs
AIC=AIC.RES_M20A.D771201.nc           ! initial conditions (atm.)
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc  ! initial conditions (ground)
! AIC=1DEC1951.rsfE000   ! or:    ! initial conditions (atm./ground), no GIC, ISTART=8
!AIC=AIC.RES_M20A.D771201          ! initial conditions (atm.)     needs GIC, ISTART=2
!GIC=GIC.E046D3M20A.1DEC1955       ! initial conditions (ground)
!AIC=1JAN1956.rsfE021TgasCM20
! 1880 OSST=OST4X5.B.1876-85avg.Hadl1.1  ! prescr. climatological ocean (1 yr of data)
! 1880 SICE=SICE4X5.B.1876-85avg.Hadl1.1 ! prescr. climatological sea ice
! prescr. climatological ocean (1 yr of data)
OSST=OST4X5.B.1975-84avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1975-84avg.Hadl1.1.nc
! prescr. climatological sea ice
SICE=SICE4X5.B.1975-84avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1975-84avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1975-84avg.Hadl1.1.nc
    ! if the prescr. ocean varies from year to year use instead:
! OSST=OST4X5.B.1871.M02.Hadl1.1  ! ocean data   Feb 1871 - 2002
! SICE=SICE4X5.B.1871.M02.Hadl1.1 ! ocean data   Feb 1871 - 2002
    ! the next group files are specific to q-flux ocean runs
! AIC=E001M20A/1JAN1960.rsfE001M20A.MXL65m   ! AIC/OHT made by aux/mkOTSPEC
! OHT=E001M20A/OTSPEC.E001M20A.MXL65m.1951-1960 ! horizontal ocean heat transport
OCNML=Z1O.B4X5.cor.nc                ! mixed layer depth (use for post processing)
    ! files needed for all models
CDN=CD4X500S.ext.nc               ! surf.drag coefficient
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops.ext.nc
CROPS=CROPS2007_72X46N.cor4_nocasp.nc  ! veg. fractions, crops history
SOIL=S4X50093.ext.nc
TOPO=Z72X46N.cor4_nocasp.nc   ! soil/topography bdy.conds
REG=REG4X5                        ! special regions-diag
RVR=RD4X525.RVR.1.bin                   ! river direction file
RADN1=sgpgxg.table8               ! rad.tables and history files
! RADN2=radfil33k                   !     8/2003 version
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800              ! rad.tables and history files
RADN4=LWCorrTables33k
RADN3=miescatpar.abcdv2
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! RADNA,RADNB are no longer used
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop.  aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr      ! need KSOLAR=2
RADNE=topcld.trscat8
#include "rad_72x46_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46_a.ij.ext.nc
MSU_wts=MSU_SSU_RSS_weights.txt
GLMELT=GLMELT_4X5.OCN.nc   ! glacial melt distribution
ISCCP=ISCCP.tautables
!         NUDGING INPUT
u2000.nc=uwnd.2000.GISS4x5_MANIP.nc
v2000.nc=vwnd.2000.GISS4x5_MANIP.nc
  
 
Label and Namelist:
E5TgasM20 (ModelE1 1980 atm/ocn, 20 layers, Water isotopes)
 
DTFIX=300
 
&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0        ! ocn is prescribed
Kvflxo=0        ! do NOT save VFLXO (daily) (use 1 to prepare for q-flux run)
 
! parameters usually not changed when switching to q-flux ocean:
 
! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.        ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin.  drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
 
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
 
PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers
 
xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)
 
 
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.59      ! U00ice up => nethtz0 down (alb down); goals: nethtz0=0,plan.alb=30%
U00wtrX=1.37    ! U00wtrX+.01=>nethtz0+.5   (alb down);        for global annual mean
!        U00wtrX=1.39    ! use with 1880 atmosphere/ocean
!1979    U00wtrX=1.37    ! use with 1979 atmosphere/ocean
 
RWCLDOX=1.   !  wtr cld particle size *RWCLDx over ocean
RICLDX=1.    !  ice cld particle size * 1(at 0mb)->RICLDx (at 1000mb)
 
anudgeu=0.01            !for nudged winds
anudgev=0.01
 
CO2X=1.
H2OstratX=1.
 
H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
                                                                                                                                                                          
! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1979 ! if -1, crops in VEG-file is used
s0_yr=1979
s0_day=182
ghg_yr=1979
ghg_day=182
volc_yr=1979
volc_day=182
aero_yr=1979
o3_yr=1979
 
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
Ndisk=48        ! use =480 on halem
SUBDD=' ' !'PS TALL RALL DALL CALL'       ! sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags 0hrly
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=1         ! use =13 to save cpu time, but energy cons. diag may not be accurate
nda5s=1         ! use =13 to save cpu time, but energy cons. diag may not be accurate
 
supsatfac=0.004
to_volume_MixRat=1,0,0,0   ! for tracer printout
to_per_mil=0,1,1,1   ! for tracer printout
&&END_PARAMETERS
 
 &INPUTZ
   YEARI=2000,MONTHI=1,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=2000,MONTHE=1,DATEE=1,HOURE=0, KDIAG=0,2,2,10*0,
   ISTART=2,IRANDI=0,
   YEARE=2000,MONTHE=1,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
 / 
