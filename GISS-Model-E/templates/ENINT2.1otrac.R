ENINT2.1otrac.R      GISS Model E atm-ocean-otracers  June 2018

Earobio_AR6: E6F40 + E200F40oQ40 + new obio :
             new qus tracer moments
             new ABIO calculation
             drop ihra calculations for better avgq
             new pp (rmus,remin nitr/sili)
E200F40oQ40: E199F40oQ40 + 5m maximum for lake ice (LAKE_ICE_MAX=5.) 

Preprocessor Options
#define STDHYB                   ! standard hybrid vertical coordinate
#define ATM_LAYERING L40         ! 40 layers, top at .1 mb
#define NEW_IO                   ! new I/O (netcdf) on
#define NEW_IO_SUBDD
#define CACHED_SUBDD
#define IRRIGATION_ON
#define MODIS_LAI
#define CHECK_OCEAN                  ! needed to compile aux/file CMPE002
#define SIMPLE_MESODIFF
#define OCN_LAYERING L40_5008m

#define ODIFF_FIXES_2017
#define EXPEL_COASTAL_ICEXS
#define NEW_BCdalbsn

#define TRACERS_OCEAN               ! RUSSELL's Ocean tracers activated 
#define TRACERS_OCEAN_INDEP         ! independently defined ocn tracers -- ocean tracers indept of atm tracers
#define OBIO_ON_GISSocean           ! obio on Russell ocean
#define TRACERS_ON                  ! include tracers code
#define TRACERS_OceanBiology
#define pCO2_ONLINE
#define TRACERS_GASEXCH_ocean       ! ANY ocean: special tracers to be passed to ocean
#define TRACERS_GASEXCH_ocean_CO2   ! ANY ocean: special tracers to be passed to ocean
#define OCN_CFC
#define constCO2
#define obio_rhsdiags

End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm144x90                           ! horizontal resolution is 144x90 -> 2x2.5deg
AtmLayering                         ! vertical resolution
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform
ORES_1Qx1 OFFT288E                  ! ocean horiz res 1.25x1deg

IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

    ! lat-lon grid specific source codes
AtmRes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
IRRIGMOD                            ! irrigation module
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO ocalbedo ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics
OCN_DRV                             ! driver for ocean-grid components
OLAYERS                             ! ocean layering options
ODIAG_COM  OCEAN_COM  OGEOM         ! dynamic ocean modules
OCNDYN  OCNDYN2  OTIDELL            ! dynamic ocean routines
OCN_Interp                          ! dynamic ocean routines
OSTRAITS_COM  OSTRAITS              ! straits parameterization
OCNKPP                              ! ocean vertical mixing
OCNMESO_DRV OCNTDMIX OCNGM          ! ocean mesoscale mixing
OCEANR_DIM  OFLUXES
ODIAG_PRT                           ! ocean diagnostic print out
OCNFUNTAB                           ! ocean function look up table
SparseCommunicator_mod              ! sparse gather/scatter module
OCNQUS                              ! QUS advection scheme
OCNGISS_TURB                        ! GISS Turbulence vertical mixing scheme
OCNGISS_SM                          ! GISS Sub-mesoscale mixing scheme
ODIAG_ZONAL
SUBDD

OCN_Int_LATLON                      ! atm-ocn regrid routines

#include "tracer_shared_source_files"
     ! atmospheric tracers
TRDIAG
TRACER_GASEXCH_CO2                  ! tracer functions needed for gas exch expts
TRACER_GASEXCH_CFC                 ! tracer functions needed for gas exch expts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!   OCEAN TRACERS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OCN_TRACER_COM
OCN_TRACER

     ! ocean carbon cycle
obio_dim         |$(R8)|
obio_incom       |$(R8)|
obio_com_R       |$(R8)|
obio_forc        |$(R8)|
obio_init        |$(R8)|
obio_bioinit     |$(R8)|
obio_model       |$(R8)|
obio_daysetrad   |$(R8)|
obio_daysetbio   |$(R8)|
obio_sfcirr      |$(R8)|
obio_edeu        |$(R8)|
obio_ptend       |$(R8)|
obio_carbon      |$(R8)|
obio_update      |$(R8)|
obio_sinksettl_R |$(R8)|
!!!obio_alkalinity  |$(R8)|
!!!obio_diffmod|$(R8)|
obio_chkbalances |$(R8)|
obio_conserv_R |$(R8)|

Components:
tracers
shared MPI_Support solvers giss_LSM 
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT 
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
    ! start from the restart file of an earlier run ...                 ISTART=8
!AIC=1JAN4101.rsfE199F40oQ40.nc    !initial conditions, no GIC needed
AIC=/discover/nobackup/aromanou/TEST/1JAN4190.rsfE200F40oQ40.nc


OFTAB=OFTABLE_NEW                   ! ocean function table
KBASIN=KB288X180.modelE.BS1.nc      ! ocean basin designations (1 cell Bering Strait)
TOPO_OC=altocnbc288x180_20170717/OZ1QX1N.BS1.BAB.PG.HB.GIB.SIC.nc    ! ocean fraction and topography
! TOPO_OC=altocnbc288x180_20170717/OZ1QX1N.BS1.BAB.PG.HB.GIB.nc        ! ocean fraction and topography
TIDES=TIDAL_e_v2_1QX1               ! ocean bottom tidal energy and velocity squared
OSTRAITS=altocnbc288x180_20170717/OSTRAITS_288x180_zspec_BAB.T4.nml  ! parameterized straits info
TOPO=Z2HX2fromZ1QX1N.BS1.nc        ! surface fractions and topography (1 cell Bering Strait)
ICEDYN_MASKFAC=iceflowmask_144x90.nc

TDISS=altocnbc288x180_20170717/TIDAL_e_v2_1QX1.HB.nc
TDISS_N=tdiss/Jayne2009_288x180.nc
POROS=altocnbc288x180_20170717/poros.nc

RVR=RD_Fd.nc             ! river direction file
NAMERVR=RD_Fd.names.txt  ! named river outlets

CDN=CD144X90.ext.nc
VEG=V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ext.nc
LAIMAX=V144x90_EntMM16_lai_max_trimmed_scaled_ext.nc
HITEent=V144x90_EntMM16_height_trimmed_scaled_ext.nc
LAI=V144x90_EntMM16_lai_trimmed_scaled_ext.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
IRRIG=Irrig144x90_1848to2100_FixedFuture_v3.nc
SOIL=S144X900098M.ext.nc
TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
GLMELT=GLMELT_144X90_gas.OCN.nc
RADN1=sgpgxg.table8                           ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800  ! rad.tables and history files
RADN4=LWCorrTables33k                         ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-2014_CMIP6_hdr  ! needs MADVOL=2
RADN8=cloud.epsilon4.72x46
!RADN9=solar.lean2015.ann1610-2014.nc ! need KSOLAR=2
RADN9=solar.CMIP6official.ann1850-2299.nc ! need KSOLAR=2
RADNE=topcld.trscat8

ISCCP=ISCCP.tautables
GHG=GHG.CMIP6.1-2014.txt  !  GreenHouse Gases for CMIP6 runs up to 2014
CO2profile=CO2profile.Jul16-2017.txt ! scaling of CO2 in stratosphere
dH2O=dH2O_by_CH4_monthly

! Begin NINT E2.1 input files

BCdalbsn=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCdalbsn
DUSTaer=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/DUST
TAero_SUL=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/SUL
TAero_SSA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/SSA
TAero_NIT=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/NIT
TAero_OCA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/OCA
TAero_BCA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCA
TAero_BCB=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCB

O3file=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/O3
Ox_ref=o3_2010_shindell_144x90x49_April1850.nc

! End NINT E2.1 input files

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

!!!!!!!!!!!!!!!!!!! obio  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfle1=abw25b.dat                         ! seawater spectral absorp.
                                         ! and scatt. coefs
cfle2=acbc25b.dat                        ! phytoplankton spectrl absorp.
                                         ! and scatt. coefs
!!!!pco2table=pco2.tbl.asc               ! table to compute pco2 values
                                         ! from sst,sss,dic,alk
                                         ! if not defined pCO2_ONLINE
nitrates_inicond=no3_nodc_annmean_180x288.nc    ! initial cond for nitrates (NODC)
silicate_inicond=sio2_nodc_annmean_180x288.nc   ! initial cond for silicate (NODC)
dic_inicond=dic_glodap_annmean_180x288.nc       ! initial cond for dic (GLODAP)
alk_inicond=alk_glodap_annmean_180x288.nc       ! initial cond/forc for alk(GLODAP)
!!!oasimdirect=oasimdirect_20w_new       ! spectral light components
                                         ! if not def OBIO_RAD_coupling
atmFe_inicond=iron_gocart_1x1mon_180x288.nc     ! GOCART iron flux
atmFedirect1=iron_ron_195x180_20w.asc    ! Ron Miller dust fluxes
facirr=facirr.asc                        ! factors for mean irrad w/in water
eda_esa_ratios=eda_esa_ratios.asc        ! ratios of rad spectrl components
!!!!!!!!!!!!!!!!!!! obio_rad  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHL_DATA=CHL_WG_2x2.5zavg                !CHL_WG_4x5 or CHL_WG_2x2.5
!! CHL_DATA=CHL_WG_2x2.5                    !CHL_WG_4x5 or CHL_WG_2x2.5
                                         !in RUSSELL ocean grid
                                         !to be used with CHL_from_SeaWIFs
cfcatm_data=CFC_ATM_Hist_2014.txt   !CFC concentrations in atmosphere


Label and Namelist:  (next 2 lines)
xxtestmaster (E200F40oQ40 + new abio + new avgq)

&&PARAMETERS
ocean_trname = 'OceanAge  Ventilatn WatrMass1  WatrMass2  WatrMass3 aoCFC aoCFC12 aoSF6 abioDIC'
! parameters set for coupled ocean runs:
KOCEAN=1            ! ocn is prognostic
OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
OTIDE = 0           !  Ocean tides are not used
variable_lk=1   ! variable lakes
init_flake=1
ocean_use_qus=1     ! Advection uses the quadratic upstream scheme
DTO=112.5
ocean_use_tdmix=1  ! tdmix scheme for meso mixing
ocean_use_gmscz=1  ! vertically variation of meso diffusivity, option 1
ocean_kvismult=2.  ! mult. factor for meso diffusivity
ocean_enhance_shallow_kmeso=1 ! stronger meso mixing in shallow water
ocean_use_tdiss=1  ! simple tidally induced diapycnal diffusivity

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000055  ! threshold (1/s) for triggering deformation waves
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.1       ! default is 0.5
CDEF=1.6       ! tuning factor for deformation -> momentum flux
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF


! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! The following two lines are only used when aerosol/radiation interactions are off
FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.3,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.655  ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 2.
use_vmp=1
radius_multiplier=1.1

PTLISO=0.        ! pressure(mb) above which radiation assumes isothermal layers
H2ObyCH4=1.      ! if =1. activates stratospheric H2O generated by CH4 without interactive chemistry
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
!crops_yr=1850  ! if -1, crops in VEG-file is used
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
!irrig_yr=1850
volc_yr=-1
!volc_day=182
!aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=1850
dalbsnX=1.
!o3_yr=-1850
!aer_int_yr=1850    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

MADVOL=2

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! surface interaction computed NIsurf times per source time step
NRAD=5           ! radiation computed NRAD times per source time step
! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=1          ! 0: no output; 1: save .acc; 2: unused; 3: include ocean data
KRSF=12          ! 0: no output; X: save rsf at the beginning of every X month
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=48           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960        ! write fort.1.nc or fort.2.nc every NDISK source time step

! parameters that affect CO2 gas exchange
!!! atmCO2=368.6          !uatm for year 2000
!!! atmCO2=0.             !prognostic atmCO2
!!!atmCO2=285.226         !uatm for new preindustrial runs AR5 runs
atmCO2=284.65             !uatm for OCMIP6 preindustrial runs
to_volume_MixRat=1    ! for tracer printout
solFe=0.02            ! default iron solubility
!!!solFe=0.05            ! enhanced iron solubility

&&END_PARAMETERS

&INPUTZ
YEARI=1850,MONTHI=1,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
YEARE=1850,MONTHE=1,DATEE=3,HOURE=0,     KDIAG=12*0,9,
ISTART=8,IRANDI=0, YEARE=1850,MONTHE=1,DATEE=2,HOURE=0,
/

