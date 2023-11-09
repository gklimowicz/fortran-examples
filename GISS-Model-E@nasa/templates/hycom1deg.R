hycom1deg.R  standalone HYCOM run in model E framework      M. Kelley 6/17/2011

Resolution: 387x360, 26 layers
Prescribed Atmosphere: CORE CNYF_2
Includes sea ice physics, but no transport yet

Preprocessor Options
#define NEW_IO
#define HYCOM1deg                    ! 1deg 26 layer hycom (387x360x26)
#define STANDALONE_HYCOM
#define STANDALONE_OCEAN
End Preprocessor Options

Object modules:

OCN_DRV                             ! driver for ocean-grid components

SEAICE SEAICE_DRV                   ! seaice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules

#include "hycom_source_files"

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
GIC=zeroSIhycom1deg.nc

PREC=CORE/H1/CNYF_2/precip.nc
RAD=CORE/H1/CNYF_2/rad.nc
SLP=CORE/H1/CNYF_2/slp.nc
T10=CORE/H1/CNYF_2/t.nc
Q10=CORE/H1/CNYF_2/q.nc
U10=CORE/H1/CNYF_2/u.nc
V10=CORE/H1/CNYF_2/v.nc
RUNOFF=CORE/H1/runoff.nc
SSS=CORE/H1/srfsal.nc
RSI=CORE/H1/SICE_387x360x365.1975-1984avg.HadISST1.1.nc

temp_ini=temp387x360x26dec_hv_z1.txt ! 3-d temperature as initial condition
salt_ini=salt387x360x26dec_hv_z1.txt ! 3-d salinity as initial condition
pout_ini=pout387x360x26dec_hv_z1.txt ! 3-d layer pressure as initial condition
latlonij=latlon387x360.4bin          ! lat & lon at 4 positions in each ocean grid box
ibasin=ibasin387x360.txt             ! ocean basin mask
cososino=cososino387x360.8bin        ! cos/sin of i,j axis angle on ocean grid
kpar=seawifs_kpar_387x360.tbin       ! monthly/annual seawifs_kpar data
hycomtopo=depth387x360.4bin          ! topography used in ocean model with Baltic Sea

GTAU=sgpgxg.table8

Label and Namelist:  (next 2 lines)
hycom1deg (standalone 1-degree 26-layer hycom)

&&PARAMETERS

DTsrc=1800.      ! cannot be changed after a run has been started

KOCEAN=1            ! ocn is prognostic

sss_restore_dt=300. ! timescale (days) for relaxing surf salinity back to obs
sss_restore_dtice=10. ! surf. sal. relax. timescale under sea ice

calc_orb_par=1
paleo_orb_yr=100.  !  BP i.e. 1950-paleo_orb_yr AD = 1850 AD

Ndisk=960
KCOPY=2          ! saving acc + rsf
Nssw=48          ! until diurnal diags are fixed, Nssw has to be even

itest=-1         ! default is -1
jtest=-1         ! default is -1
iocnmx=2         ! default is 2
brntop=50.       ! default is 50.
brnbot=200.      ! default is 200.
diapyn=1.e-7     ! default is 1.e-7
diapyc=.1e-4     ! default is .1e-4
jerlv0=1         ! default is 1
&&END_PARAMETERS

 &INPUTZ
 YEARI=1899,MONTHI=12,DATEI=01,HOURI=00, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1900,MONTHE=12,DATEE=02,HOURE=00, KDIAG=13*0,
 ISTART=2,IRANDI=0, YEARE=1899,MONTHE=12,DATEE=1,HOURE=1
/
