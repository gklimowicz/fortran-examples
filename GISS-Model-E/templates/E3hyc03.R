E3hyc03.R GISS Model E  coupled version      ssun   2/22/2007

E3hyc03: E3AoM20A + hycom kpp refinement + topo_20w + tsadts + partial kpp
       replace this section by a description of what distinguishes this run ?
       Use as many lines as you need. Look carefully at all the possible    ?
       choices, particularly the lines containing '?'. In some cases, you   ?
       will have to pick the appropriate choice to make this rundeck work   ?
       The final rundeck should contain no '?'
       Check and modify the rest of the description below:                  ?
modelE1 (3.0) 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)       ?
atmospheric composition from year 1880 ? 1979                               ?
ocean: coupled to HYCOM Version 0.9.2 KPP
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag) ?
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs  ?
filters: U,V in E-W direction (after every dynamics time step)              ?
         sea level pressure (after every physics time step)                 ?

Preprocessor Options
! #define TRACERS_ON                  ! include tracers code
! #define TRACERS_GASEXCH_Natassa     ! special tracers to be passed to ocean
End Preprocessor Options

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
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_M                          ! diagnostics (resolution dependent)
      FFT72                         ! utilities
POUT                                ! post-processing output
hycom2|$(R8)| OCEAN_hycom|$(R8)|        ! ocean model - driver
advfct|$(R8)|                         ! advection
archyb|$(R8)|                         ! continuity eqn. 
barotp|$(R8)|                         ! barotropic eqn. 
bigrid|$(R8)|                         ! basin grid
blkd2a|$(R8)| blkpp2|$(R8)|             ! block data
cnuitb|$(R8)|                         ! continuity eqn.
cpler |$(R8)|                         ! coupler
diapfc|$(R8)|                         ! diapycnal diffusion
dpthuv|$(R8)| dpudpv|$(R8)|             ! off-center depth  
eic8  |$(R8)|                         ! ice forming
geopar|$(R8)|                         ! geography related parameters
hybgn1a|$(R8)|                        ! grid generator 
inicon|$(R8)| inigis|$(R8)| inikpp|$(R8)| ! initial conditions
matinv|$(R8)| mxkpp2|$(R8)|             ! partial KPP mixing scheme
momtum|$(R8)|                         ! momemtum Eqn.
prtetc|$(R8)|                         ! print routines, etc.
reflux|$(R8)|                         ! flux conversion
sigetc|$(R8)|                         ! eqn.of state, etc.
thermf|$(R8)|                         ! thermal forcing
trcadv|$(R8)|                         ! tracer advection 
tsadts|$(R8)| advem|$(R8)|              ! T/S advection advecting t/s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!  tracer part  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TRACER_COM TRACERS_DRV              ! configurable tracer code
! TRACERS                             ! generic tracer code
! TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
! OCN_TRACER
!!!!!!!use part below only for Natassa's gas exchange experiments !!!!!!
! TRACER_GASEXCH_Natassa              ! tracer functions needed for gas exch expts
!!!!!!!!!!!!!!!!!!!!!! end  tracer part  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Components:
MPI_Support shared

Data input files:
AIC=AIC.RES_M20A.D771201.nc    !initial conditions (atm.) needs GIC,OIC ISTART=2
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc ! initial conditions (ground) and 300 year spin-up
CDN=CD4X500S.ext.nc
  ! VEG=V72X46.1.cor2.ext
VEG=V72X46.1.cor2_no_crops.ext.nc 
CROPS=CROPS2007_72X46N.cor4_nocasp.nc  ! veg. fractions, crops history
SOIL=S4X50093.ext.nc
TOPO=Z72X46N.2deg_rfn_20w              !!! hycom
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR.2deghycom_20w.bin         !!! hycom
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
OCMIP_cfc=OCMIP_cfc.dat
!!!!!!!!!!!!!!!!!!! HYCOM input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
latlonij=latlon195x180_20w.4bin    ! lat & lon at each i,j
hycomtopo=depth195x180_20w.4bin    ! topography used in ocean model
temp_ini=temp181x180x20jan_lt.txt  ! sea surface temperature as initial condition
salt_ini=salt181x180x20jan_lt.txt  ! salinity as initial condition
pout_ini=pout181x180x20jan_lt.txt  ! layer pressure as initial condition
ibasin=ibasin195x180_20w.txt       ! basin mask
flxa2o=flxa2o195x180.8bin          ! coupler weights for flux from atm to ocean
taua2o=taua2o195x180.8bin          ! weights for vector from atm to ocean
ssto2a=ssto2a195x180.8bin          ! weights for sst from ocean to atm
e_o2a=e_o2a195x180.8bin            ! weights for eastward vel from ocean to atm
n_o2a=n_o2a195x180.8bin            ! weights for northward vel from ocean to atm
cososino=cososino195x180.8bin      ! cos/sin of i,j axis angle on ocean grid

Label and Namelist:
E3hyc03 (hycoma+tsadts: bihar=.1,hybgn1,pump[20:200],part.kppx3,enh.den,no conv)

DTFIX=300
&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1        ! ocn is prognostic
variable_lk=0
init_flake=0

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
U00ice=.6        ! U00ice up => nethtz0 down (alb down); goals: nethtz0=0,plan.alb=30%
U00wtrX=1.3      ! U00wtrX+.01=>nethtz0+.7                      for global annual mean

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
volc_yr=1880  !? 1979
volc_day=182
aero_yr=1880  !? 1979
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1880 !? 1979
dalbsnX=.015
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
isccp_diags=0   ! use =0 to save cpu time if isccp-diags are not essential
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
nssw=48         ! do ssw at the end of day
&&END_PARAMETERS

 &INPUTZ
   YEARI=1800,MONTHI=01,DATEI=01,HOURI=00, !  from default: IYEAR1=YEARI
   YEARE=1800,MONTHE=01,DATEE=01,HOURE=04, KDIAG=13*0,
   ISTART=2,IRANDI=0,YEARE=1800,MONTHE=01,DATEE=01,HOURE=01
/
link_ha_feb07
