from rundeck_section import *


# real definitions here:

Grid_M20 = RundeckSection( 
    name("Grid" ),
    obj_files_text("""
     ! resolution-specific source codes
RES_M20AT DIAG_RES_M          ! horiz/vert resolution, 4x5deg, 20 layers -> 0.1mb
FFT72                              ! Fast Fourier Transform
""" )
    )

Grid_F40 = RundeckSection( 
    name("Grid" ),
    obj_files_text("""
RES_stratF40                        ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform
""" )
    )



Main = RundeckSection( 
    name("Main" ),
    components("shared", "MPI_Support", "solvers", "dd2d"),
    obj_files_text("""
MODELE                              ! ModelE cap - initialization and main loop
MODEL_COM                           ! calendar, timing variables
IO_DRV                             ! new i/o
DIAG_COM DIAG DEFACC                ! diagnostics
DIAG_PRT POUT                       ! diagn/post-processing output
""" )
    )

TrAdv_QUS3D = RundeckSection( 
    name("TrAdv_QUS3D" ),
    obj_files_text("""
QUS_DRV QUS3D                       ! advection of Q/tracers
""" )
    )

TrAdv_TQUS = RundeckSection( 
    name("TrAdv_TQUS" ),
    obj_files_text("""
QUS_DRV TQUS_DRV                    ! advection of Q/tracers
""" )
    )

StratDyn = RundeckSection( 
    name("StratDyn" ),
    obj_files_text("""
     ! GISS dynamics with gravity wave drag
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)
""" )
    )

Atm_M20 = RundeckSection(
    input_files_text("""
    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m  ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
AIC=AIC.RES_M20A.D771201.nc          ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc ! initial conditions (ground)
! prescr. climatological ocean (1 yr of data)
OSST=OST4X5.B.1876-85avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1876-85avg.Hadl1.1.nc
! prescr. climatological sea ice
SICE=SICE4X5.B.1876-85avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1876-85avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1876-85avg.Hadl1.1.nc
!! q-flux ocean: use the next line instead,       set KOCEAN=1
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960            ! ocean horizontal heat transports
!! OCNML is not used if KOCEAN=0, but needed in and to prepare for q-flux model
OCNML=Z1O.B4X5.cor.nc             ! mixed layer depth (needed for post processing)
TOPO=Z72X46N.cor4_nocasp.nc       ! topography

RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets

CDN=CD4X500S.ext.nc                 ! surf.drag coefficient
VEG=V72X46.1.cor2_no_crops.ext.nc   ! veg. fractions
CROPS=CROPS2007_72X46N.cor4_nocasp.nc
SOIL=S4X50093.ext.nc                ! soil bdy.conds
TOP_INDEX=top_index_72x46_a.ij.ext.nc  ! only used if #define DO_TOPMODEL_RUNOFF
ZVAR=ZVAR4X5.nc                      ! topographic variation for gwdrag
! probably need these (should convert to 72x46)
!soil_textures=soil_textures_top30cm.nc
SOILCARB_global=soilcarb_top30cm_4x5.nc
GLMELT=GLMELT_4X5.OCN.nc

BC_dep=BC.Dry+Wet.depositions.ann.nc
"""),
    parameters_text("""
xCDpbl=1.
cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.54      ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00      ! below 850mb and MC regions; then tune this to get rad.balance

WMUI_multiplier = 2.

KSIALB=0         ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)

aer_rad_forc=0
cloud_rad_forc=1

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
""")
    )

Atm_F40 = RundeckSection(
    input_files_text("""
    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m  ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
AIC=AIC.RES_F40.D771201.nc      ! observed init cond (atm. only)
GIC=GIC.144X90.DEC01.1.ext_1.nc ! initial ground conditions
! prescr. climatological ocean (1 yr of data)
OSST=OST_144x90.1876-1885avg.HadISST1.1.nc
OSST_eom=OST_144x90.1876-1885avg.HadISST1.1.nc
! prescr. climatological sea ice
SICE=SICE_144x90.1876-1885avg.HadISST1.1.nc
SICE_eom=SICE_144x90.1876-1885avg.HadISST1.1.nc
ZSIFAC=SICE_144x90.1876-1885avg.HadISST1.1.nc
!! q-flux ocean: use the next line instead,       set KOCEAN=1
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960            ! ocean horizontal heat transports
!! OCNML is not used if KOCEAN=0, but needed in and to prepare for q-flux model
OCNML=Z1O.B144x90.nc                               ! mixed layer depth
TOPO=Z2HX2fromZ1QX1N.nc

RVR=RD_modelE_Fa.nc             ! river direction file
NAMERVR=RD_modelE_Fa.names.txt  ! named river outlets

CDN=CD144X90.ext.nc
VEG=V144X90_no_crops.ext.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
SOIL=S144X900098M.ext.nc
TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
GLMELT=GLMELT_144X90_gas.OCN.nc

BC_dep=BC.Dry+Wet.depositions.ann_144x90.nc
"""),
    parameters_text("""
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000045 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.2       ! default is 0.5
CDEF=1.5       ! deformation drag coefficient
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.72 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.68 ! below 850mb and MC regions; tune this last  to get rad.balance

NIsurf=1         ! (surf.interaction NIsurf times per physics time step)
""")
    )



Atm = RundeckSection( 
    name("Atm" ),
    obj_files_text("""
     ! GISS dynamics w/o gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
    ! lat-lon grid specific source codes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB_E1                            ! turbulence in whole atmosphere
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
""" ),
    input_files_text("""
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800            ! rad.tables and history files
RADN4=LWCorrTables33k              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.SATO.1850-1999.Apr02_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean2015.ann1610-2014_hdr ! need KSOLAR=2
RADNE=topcld.trscat8

ISCCP=ISCCP.tautables
GHG=GHG.Mar2009.txt ! use GHG.Jul2009.txt for runs that start before 1850
dH2O=dH2O_by_CH4_monthly
DUSTaer=dust_mass_CakmurMillerJGR06_$resijx20x7x12_unlim.nc
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_$resijx20_1890-2000h.nc
TAero_SSA=SSA_Koch2008_kg_m2_$resijx20h.nc
TAero_NIT=NIT_Bauer2008_kg_m2_$resijx20_1890-2000h.nc
TAero_OCA=OCA_Koch2008_kg_m2_$resijx20_1890-2000h.nc
TAero_BCA=BCA_Koch2008_kg_m2_$resijx20_1890-2000h.nc
TAero_BCB=BCB_Koch2008_kg_m2_$resijx20_1890-2000h.nc
O3file=o3_2005_shindelltrop_$resijx49_1850-1997_ple.nc

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG$resDEG                      ! special regions-diag
"""),
parameters_text("""
! parameters set for choice of ocean model:
KOCEAN=0        ! ocean is prescribed
!! KOCEAN=1        ! ocean is computed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
!  Kvflxo=1     ! saves VFLXO files to prepare for q-flux runs (mkOTSPEC)

variable_lk=1   ! variable lakes

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
!#include "gwdragM20_params"

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
!crops_yr=1850  ! if -1, crops in VEG-file is used
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
volc_yr=-1
!volc_day=182
!aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=1850
dalbsnX=.024
!o3_yr=-1850
!aer_int_yr=1850    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NRAD=5           ! radiation (every NRAD'th physics time step)
! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=2          ! saving acc + rsf
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=480
"""),
    variants(M20=Atm_M20, F40=Atm_F40)
    )


Surface = RundeckSection( 
    name("Surface" ),
    components("giss_LSM"),
    obj_files_text("""
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
LAKES_COM LAKES                     ! lake modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
""" )
    )


OceanPrescribed = RundeckSection( 
    name("Ocean" ),
    obj_files_text("""
OCN_DRV                             ! driver for ocean-grid components
OCEAN OCNML                         ! ocean modules
SEAICE SEAICE_DRV                   ! seaice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
""" )
    )


VegetationEnt = RundeckSection( 
    name("VegetationEnt" ),
    components("Ent"),
    component_options("OPTS_Ent = ONLINE=YES PS_MODEL=FBB"),
    obj_files_text("""
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
!ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
""" ),
    obj_files("ENT_DRV",  "ENT_COM")
    )
