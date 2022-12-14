oceanR1deg.R ocean R standalone version          M. Kelley 6/17/2011

Resolution: 288x180, 32 layers
Prescribed Atmosphere: CORE CNYF_2
Includes sea ice physics and transport

Preprocessor Options
#define NEW_IO
#define STANDALONE_OCEAN
#define OCN_LAYERING L32
!#define OCN_GISS_TURB
!#define OCN_GISS_SM
End Preprocessor Options

Object modules:

OCN_DRV                             ! driver for ocean-grid components

SEAICE SEAICE_DRV                   ! seaice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules

ORES_1Qx1 OFFT288E                  ! ocean horiz res 1.25x1deg
#include "dynamic_ocn_source_files"
OCN_Int_LATLON                      ! atm-ocn regrid routines

MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
FLUXES                              ! atm-ocean exchange variables
SOATM_COM                           ! modelE atmosphere modules
SOATM_DRV                           ! data-driven atm. with modelE procedure names
IO_DRV                              ! i/o
COSZ_2D                             ! solar zenith angle

Components:
shared MPI_Support solvers
dd2d

Component Options:
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
AIC=CORE/empty_AIC
GIC=GIC.288X180.DEC01.1.ext.nc

IAF=CORE/R1/IAF
PREC=CORE/R1/CNYF_2/precip.nc
RAD=CORE/R1/CNYF_2/rad.nc
SLP=CORE/R1/CNYF_2/slp.nc
T10=CORE/R1/CNYF_2/t.nc
Q10=CORE/R1/CNYF_2/q.nc
U10=CORE/R1/CNYF_2/u.nc
V10=CORE/R1/CNYF_2/v.nc
RUNOFF=CORE/R1/runoff.nc
SSS=CORE/R1/srfsal.nc
RSI=CORE/R1/SICE_288x180x365.1975-1984avg.HadISST1.1.nc

#include "dynamic_ocn_288x180_input_files_AR5"

GTAU=sgpgxg.table8

Label and Namelist:  (next 2 lines)
oceanR1deg (standalone 288x180 32-layer ocean R)

&&PARAMETERS

DTsrc=1800.      ! cannot be changed after a run has been started

KOCEAN=1            ! ocn is prognostic

sss_restore_dt=300. ! timescale (days) for relaxing surf salinity back to obs
sss_restore_dtice=10. ! surf. sal. relax. timescale under sea ice

calc_orb_par=1
paleo_orb_yr=100.  !  BP i.e. 1950-paleo_orb_yr AD = 1850 AD

OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
OTIDE = 0           !  Ocean tides are not used
DTO=112.5
!ocean_use_qus=1     ! Advection uses the quadratic upstream scheme

!interannual_forcing=1 ! atm. forcing history varies with year
!iaf_year_start=1948   ! first year of forcing history
!iaf_year_end=2007     ! last year of forcing history

Ndisk=960
KCOPY=2          ! saving acc + rsf
Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
master_yr=1850

&&END_PARAMETERS

 &INPUTZ
   YEARI=1900,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
   YEARE=1901,MONTHE=1,DATEE=1,HOURE=3, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1900,MONTHE=12,DATEE=2,HOURE=1
/
