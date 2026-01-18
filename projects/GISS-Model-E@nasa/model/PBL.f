#include "rundeck_opts.h"

#ifdef CUBED_SPHERE
#define PBL_USES_GCM_TENDENCIES
#endif

      MODULE SOCPBL
!@sum module SOCPBL defines subroutines and variables associated with
!@+  the boundary layer physics.
!@+  It sets up npbl(=8) sublayers between the surface (sublayer 1)
!@+  and the middle of the first GCM layer (sublayer npbl),
!@+  and integrates, over these sublayers, the dynamic equations
!@+  for the mean turbulent variables using turbulence models,
!@+  to find the surface values of these variables and
!@+  related fluxes.
!@+  t_pbl_args is a derived type structure which contains all
!@+  input/output arguments for PBL.
!@+  SOCPBL contains the following subroutines:
!@+  advanc,stars,getl,dflux,simil,griddr,tfix
!@+  ccoeff0,getk,e_eqn,t_eqn,q_eqn,uv_eqn,
!@+  t_eqn_sta,q_eqn_sta,uv_eqn_sta,
!@+  inits,tcheck,ucheck,check1,output,rtsafe.
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)


      USE CONSTANT, only : grav,pi,radian,bygrav,teeny,deltx,tf
     &     ,by3,lhe,rgas,rhows,mair,byrhows,sha,shv,shw,stbo,visc_air
      USE SEAICE, only : tfrez
      USE LANDICE, only : snmin
#ifdef TRACERS_ON
      use OldTracer_mod, only: trName, nWATER, tr_wd_TYPE
      USE TRACER_COM, only: NTM
#ifdef TRACERS_SPECIAL_O18
      use TRACER_COM, only: n_water
#endif
#ifdef TRACERS_DRYDEP
      use OldTracer_mod, only: dodrydep
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      use TRACER_COM, only: Ntm_dust, n_soildust
#endif
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: NBINS, xk
#endif
#endif /* TRACERS_DUST || TRACERS_MINERALS || TRACERS_AMP || TRACERS_TOMAS */

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      use trdust_mod,only : nDustBins
#endif
      USE TRIDIAG_MOD, only :  TRIDIAG
      IMPLICIT NONE

      private

      ! public parameters
      public n,zgs,XCDpbl,kappa,emax,skin_effect,ustar_min
     &   ,lmonin_min,lmonin_max,xdelt,calc_wspdf
     &   ,sigma,gamahu,gamams,gamamu,zet1,slope1,zeth
     &   ,se,k_max,kmmin,khmin,kqmin,kemin,emin
     &   ,find_phim0,find_phih


      ! public interfaces
      public advanc,inits,ccoeff0
      public alloc_pbl_args, dealloc_pbl_args

      ! model coefficients (actually a hack, but leave it for now)
      public rimax,ghmin,ghmax,gmmax0,gm_at_rimax,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,s7,s8,c1,c2,c3,c4,c5,b1,b123,b2,prt
     &     ,g0,d1_3,d2_3,d3_3,d4_3,d5_3
     *     ,s0_3,s1_3,s2_3,s3_3,s4_3,s5_3,s6_3
     *     ,g1,g2,g3,g4,g5,g6,g7,g8

      public t_pbl_args

      integer, parameter :: n=8  !@param n  no of pbl. layers
      integer, parameter :: npbl=n
      integer, public :: MAXNTM

c**** t_pbl_args is a derived type structure which contains all
c**** input/output arguments for PBL
c**** Please, use this structure to pass all your arguments to PBL
c**** Do not use global variables for that purpose !
      type t_pbl_args
        ! input:
        real*8 dtsurf,zs1,tgv,tkv
        !@var qg_sat Saturation vapor mixing ratio (kg vapor / kg air in a given volume)
        real*8 qg_sat
        real*8 qg_aver,hemi,tr4
        real*8 evap_max,fr_sat,uocean,vocean,psurf,trhr0
        real*8 tg,elhx,qsol,sss_loc
        logical :: ocean,ddml_eq_1
        ! inout:
        real*8 gusti,tdns,qdns,tprime,qprime,snow
        ! output:
        real*8 us,vs,ws,tsv,qsrf,cm,ch,cq,dskin,ws0
        ! the following args needed for diagnostics
        real*8 psi,dbl,khs,ug,vg,wg,ustar,lmonin,zgs
        real*8 canopy_temperature
!@var wsgcm magnitude of the GCM surface wind - ocean currents [m/s]
!@var wspdf mean surface wind calculated from PDF of wind speed [m/s]
!@var wsubtke turbulent kinetic energy velocity scale [m/s]
!@var wsubwd dry convective velocity scale [m/s]
!@var wsubwm moist convective velocity scale [m/s]
        real*8 :: wsgcm,wspdf,wsubtke,wsubwd,wsubwm
!@var mcfrac fractional area with moist convection in grid cell
        REAL(KIND=8) :: mcfrac

#ifdef TRACERS_ON
c**** Attention: tracer arrays in this structure have dim 1:ntm
c**** while local arrays in PBL code have dim 1:ntx
c**** Tracer input/output
!@var trtop,trs tracer mass ratio in level 1/surface
!@var trsfac, trconstflx factors in surface flux boundary cond.
!@var trgrnd2 factor for correction to landice evap
!@var ntx number of tracers that need pbl calculation
!@var ntix index array to map local tracer number to global
!@var trprime anomalous tracer concentration in downdraft

        real*8, allocatable, dimension(:) :: trtop,trs,trsfac,trconstflx
        real*8, allocatable, dimension(:) :: trdn1,trprime,trgrnd2
        integer ntx
        integer, allocatable, dimension(:) :: ntix
#ifdef TRACERS_SPECIAL_O18
        real*8, allocatable, dimension(:) :: frack
#endif
#ifdef BIOGENIC_EMISSIONS
        real*8 :: emisop
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
c**** input
!@var pbl_args%wearth earth water of first layer [kg/m^2]
!@var pbl_args%aiearth earth ice of first layer [kg/m^2]
!@var pbl_args%wfcs water field capacity of first ground layer [kg/m^2]
        REAL*8 :: wearth,aiearth,wfcs
!@var pbl_args%ers_data ERS data
!@var pbl_args%dustSourceFunction distribution of preferred sources
        real(kind=8) :: ers_data,dustSourceFunction
!@var pbl_args%frclay fraction of clay
!@var pbl_args%frsilt fraction of silt
!@var pbl_args%dryhr number of hours with evaporation-precipitation greater
!@+                  Zero to allow dust emission
!@var pbl_args%vtrsh thresh. wind speed above which dust emis. is allowed [m/s]
        REAL*8 :: frclay,frsilt,dryhr,vtrsh
!@var pbl_args%pprec precipitation at previous time step [kg/m^2]
!@var pbl_args%pevap evaporation at previous time step [kg/m^2]
        REAL*8 :: pprec,pevap
!@var pbl_args%d_dust prescribed daily dust emissions [kg/m^2/s] (e.g. AEROCOM)
        real( kind=8 ) :: d_dust( nDustBins )
!@var pbl_args%mineralFractions  mineral fractions of emitted dust aerosols [1]
        real(kind=8) :: mineralFractions( max( nDustBins, ntm_dust ) )
c**** output
!@var pbl_args%pdfint integral of dust emission probability density function
        REAL*8 :: pdfint
!@var pbl_args%wtrsh velocity threshold for dust emission (depends on soil
!@+                  moisture) [m/s]
!@var pbl_args%dust_event1 number of dust events [1]
!@var pbl_args%dust_event2 number of dust events above velocity threshold
!@+                        of cubic emission scheme (diagnositcs only) [1]
!@var pbl_args%dust_flux1 dust flux [kg/m^2/s]
!@var pbl_args%dust_flux2 dust flux from cubic emission scheme (diagnostics
!@+                       only) [kg/m^2/s]
        REAL*8 :: dust_flux(Ntm_dust),dust_flux2(Ntm_dust)
     *       ,dust_event1,dust_event2,wtrsh
        REAL*8 :: z(npbl),km(npbl-1),gh(npbl-1),gm(npbl-1),zhat(npbl-1)
!@var pbl_args%hbaij accumulated precipitation - evaporation balance  [kg/m^2]
!@var pbl_args%ricntd no. of hours with negative precipitation - evaporation
!@+                   balance [1]
        REAL*8 :: hbaij,ricntd
!@var pbl_args%qdust flag whether conditions for dust emission are fulfilled
        LOGICAL :: qdust
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
        integer :: moddd,ih,ihm
#endif

#if (defined TRACERS_AEROSOLS_SEASALT) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
        real*8 :: ss1_flux,ss2_flux
#endif  /* TRACERS_AEROSOLS_SEASALT || TRACERS_AMP || TRACERS_TOMAS */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)
        real*8 :: DMS_flux
#ifdef TRACERS_TOMAS
!@var tomas_ss_flux : sea-salt emission fraction at each size bin
        real*8,dimension(nbins)::  tomas_ss_flux
#endif /* TRACERS_TOMAS */ 
#endif
#ifdef TRACERS_AEROSOLS_OCEAN
        real*8 :: OCocean_flux
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_DRYDEP
!@var dep_vel turbulent deposition velocity = 1/bulk sfc. res. (m/s)
!@var gs_vel gravitational settling velocity (m/s)
!@var stomatal_dep_vel turbulent deposition velocity via stomata(m/s)
        real*8, pointer, dimension(:) :: dep_vel=>null()
        real*8, pointer, dimension(:) ::  gs_vel=>null()
        real*8 :: stomatal_dep_vel
#endif

#ifdef TRACERS_WATER
!@var tr_evap_max maximum amount of tracer available in ground reservoir
        real*8, pointer, dimension(:) :: tr_evap_max
#endif


        real*8, dimension(:), allocatable  :: Kw_gas,alpha_gas,beta_gas
#endif


      end type t_pbl_args


!@dbparam XCDpbl factor for momentum drag coefficient
      real*8 :: XCDpbl=1d0

      real*8 :: rimax,ghmin,ghmax,gmmax0,gm_at_rimax,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,s7,s8,c1,c2,c3,c4,c5,b1,b123,b2,prt

      !for level 3 model only:
      real*8 :: g0,d1_3,d2_3,d3_3,d4_3,d5_3
     *         ,s0_3,s1_3,s2_3,s3_3,s4_3,s5_3,s6_3
     *         ,g1,g2,g3,g4,g5,g6,g7,g8

C**** boundary layer parameters
      real*8, parameter :: kappa=0.40d0 !@var kappa  Von Karman constant
!@var  zgs  height of the surface layer (nominally 10 m)
      real*8, parameter :: zgs=10. !@var zgs height of surface layer (m)

C**** parameters for surface fluxes
      !Hogstrom 1988,1996:
      real*8, parameter :: sigma=0.95d0,sigma1=1.-sigma,sigma2=1.+sigma
      real*8, parameter :: gamamu=19.0d0,gamahu=11.6d0,gamams=5.3d0,
     *     gamahs=8.d0/sigma

      real*8, parameter :: zet1=0.5d0   !var zet1 critical value of zet=z/lmonin
      real*8, parameter :: slope1=0.1d0 !var slope1 slope of PHI's of zet>zet1

      !Zeng et al 1998:
      real*8, parameter :: zetm=-1.464d0
      real*8, parameter :: zeth=-1.072d0
      real*8, parameter :: se=0.1d0,k_max=500.d0
     &  ,kmmin=1.5d-5,khmin=2.5d-5,kqmin=2.5d-5,kemin=1.5d-5,emin=1d-6

CCC !@var bgrid log-linear gridding parameter
CCC      real*8 :: bgrid

!@var smax,smin,cmax,cmin limits on drag coeffs.
!@var emax limit on turbulent kinetic energy
!@var ustar_min limit on surface friction speed 
      real*8, parameter :: smax=0.25d0,smin=0.005d0,cmax=smax*smax,
     *     cmin=smin*smin,emax=1.d5,ustar_min=1d-2
     *    ,lmonin_min=1d-6,lmonin_max=1d6

!@param xdelt When used in place of deltx in expressions involving
!@+     virtual temperature T, xdelt=0 switches off virtual T effects.
!@+     This allows the PBL module to control the definition of its
!@+     input/output T; internally it always uses virtual T when
!@+     calculating stratification but may use actual T as a
!@+     prognostic variable.
      real*8, parameter :: xdelt=0d0    ! input/output is actual T
!                          xdelt=deltx  !                 virtual T

!@param twoby3 2/3 constant
      real*8, parameter :: twoby3 = 2d0/3d0

!@dbparam skin_effect sets whether skin effects are used or not
      integer :: skin_effect=1  ! Used by default

!@dbparam calc_wspdf if calc_wspdf==1, calculate mean surface
!@+       wind speed from PDF of velocities
      integer :: calc_wspdf=0

#ifdef WATER_PROPORTIONAL
      logical :: force_limit
#endif
      CONTAINS

      subroutine advanc(pbl_args,coriol,utop,vtop,qtop,ztop,mdf
     &     ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &     ,dtdt_gcm,utop_old,vtop_old
     &     ,ilong,jlat,ihc,itype
     &     ,kms,kqs,z0m,z0h,z0q,w2_1,ufluxs,vfluxs,tfluxs,qfluxs
     &     ,u,v,t,q,e
#if defined(TRACERS_ON)
     &     ,tr,trnradius,trndens,trnmm
#endif
     &     )
!@sum time steps the solutions for the boundary layer variables.
!@+  It is called from within the subroutine pbl (in PBL_DRV.f).
!@+  All its outputs are contained in the structure pbl_args
!@+  (an instance of t_pbl_args).
!@auth  Ye Cheng/G. Hartke
c    input:
!@var  coriol  2.*omega*sin(latitude), the coriolis factor
!@var  utop  x component of wind at the top of the layer
!@var  vtop  y component of wind at the top of the layer
!@var  ttop  virtual potential temperature at the top of the layer
!@+          (if xdelt=0, ttop is the actual temperature)
!@var  qtop  moisture at the top of the layer
!@var  tgrnd0 bulk virt. pot. temp. of the ground, at the roughness height
!@+          (if xdelt=0, tgrnd0 is the actual temperature)
!@var  tg  actual ground temp
!@var  elhx relevant latent heat for qg calc
!@var  qgrnd0 initial moisture at the ground, at the roughness height
!@var  qgrnd_sat saturated moisture at the ground, at the roughness height
!@var  evap_max maximal evaporation from unsaturated soil
!@var  fr_sat fraction of saturated soil
!@var  ztop height of the first model layer, approx 200 m if lm=9
!@var  dtime  time step
!@var  ilong  longitude identifier
!@var  jlat   latitude identifier
!@var  ihc    Height class (for debugging printout)
!@var  itype  1, ocean; 2, ocean ice; 3, land ice; 4, land
!@var  psurf surface pressure
!@var  trhr0 incident long wave radiation
!@var  sss_loc SSS at i,j
!@var  dtdt_gcm ttop tendency from processes other than turbulence (K/s)
!@var  utop_old  utop after turb. diffusion during previous timestep
!@var  vtop_old  vtop after turb. diffusion during previous timestep
c   output:
!@var  us  x component of the surface wind (i.e., due east)
!@var  vs  y component of the surface wind (i.e., due north)
!@var  tsv  virtual potential surface temperature
!@+          (if xdelt=0, tsv is the actual temperature)
!@var  qsrf  surface specific moisture
!@var  kms  surface value of km
!@var  khs  surface value of kh
!@var  kqs  surface value of kq
!@var  ws  magnitude of surface wind modified by buoyancy flux (m/s)
!@var  ustar  friction speed
!@var  cm  dimensionless momentum flux at surface (drag coeff.)
!@var  ch  dimensionless heat flux at surface (stanton number)
!@var  cq  dimensionless moisture flux at surface (dalton number)
!@var  z0m  roughness height for momentum (if itype=1 or 2)
!@var  z0h  roughness height for heat
!@var  z0q  roughness height for moisture
!@var  dskin skin-bulk SST or skin-bulk snow temp difference
!  more output (may be duplicated)
!@var US     = x component of surface wind, positive eastward (m/s)
!@var VS     = y component of surface wind, positive northward (m/s)
!@var WSGCM  = magnitude of the GCM surface wind - ocean currents (m/s)
!@var WSPDF  = mean surface wind calculated from PDF of wind speed (m/s)
!@var WS     = magn of GCM surf wind - ocean curr + buoyancy + gust (m/s)
!@var TSV    = virtual potential temperature of the surface (K)
!@+            (if xdelt=0, tsv is the actual temperature)
!@var QS     = surface value of the specific moisture
!@var DBL    = boundary layer height (m)
!@var KMS    = momentum transport coefficient at ZGS (m**2/s)
!@var KHS    = heat transport coefficient at ZGS (m**2/s)
!@var KHQ    = moist transport coefficient at ZGS (m**2/s)
!@var PPBL   = pressure at DBL (mb)
!@var USTAR  = friction speed (square root of momentum flux) (m/s)
!@var CM     = drag coefficient (dimensionless surface momentum flux)
!@var CH     = Stanton number   (dimensionless surface heat flux)
!@var CQ     = Dalton number    (dimensionless surface moisture flux)
!@var z0m   = roughness length for momentum,
!@+           prescribed for itype=3,4 but computed for itype=1,2 (m)
!@var z0h   = roughness length for temperature (m)
!@var z0q   = roughness length for water vapor (m)
!@var UG     = eastward component of the geostrophic wind (m/s)
!@var VG     = northward component of the geostrophic wind (m/s)
!@var MDF    = downdraft mass flux (m/s)
!@var WINT   = integrated surface wind speed over sgs wind distribution
!@var  u  local due east component of wind
!@var  v  local due north component of wind
!@var  t  local virtual potential temperature
!@+          (if xdelt=0, t is the actual temperature)
!@var  q  local specific humidity (a passive scalar)
!@var  e  local turbulent kinetic energy
#if defined(TRACERS_ON)
!@var  trtop  tracer conc. at the top of the layer
!@var  trs  surface tracer conc.
!@var  trsfac  factor for pbl_args%trs in surface boundary condition
!@var  trconstflx  constant component of surface tracer flux
!@var  ntx  number of tracers to loop over
!@var  ntix index of tracers used in pbl
#endif
#if defined(TRACERS_ON) && defined(TRACERS_GASEXCH_ocean)
!@var  Kw_gas  gas exchange transfer velocity at i,j only over ocean
!@var  alpha_gas  solubility of gas
!@var  beta_gas  conversion term  that includes solubility
#endif
#if defined(TRACERS_ON) && defined(TRACERS_WATER)
!@var  tr_evap_max max amount of possible tracer evaporation
#endif
c  internals:
!@var  n     number of the local, vertical grid points
!@var  lscale turbulence length scale. computed on secondary grid.
!@var  z     altitude of primary vertical grid points
!@var  zhat  altitude of secondary vertical grid points
!@var  dzh   dz evaluated at zhat(i)
!@var  dz    dz evaluated at z(i)
!@var  dxi   (z(n)-z(1))/(n-1)
!@var  km    turbulent momentum tranport coefficient.
!@var  kh    turbulent thermometric conductivity. computed
!@var  ke    transport coefficient for the turbulent kinetic energy.
!@var  tv    local virtual potential temperature

#if (defined TRACERS_AEROSOLS_SEASALT) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      use tracers_seasalt, only: read_seasalt_sources
#endif  /* TRACERS_AEROSOLS_SEASALT || TRACERS_AMP || TRACERS_TOMAS */
#ifdef TRACERS_TOMAS
      USE TOMAS_EMIS
#endif 
#if defined(TRACERS_ON)
      use tracer_com, only: tr_mm, gasex_index, n_co2n, n_cfcn
#endif

!@var tdns downdraft temperature in K, (i,j)
!@var qdns downdraft humidity in kg/kg, (i,j)

      use runtimecontrols_mod, only: obio_wspdf

      implicit none

      !-- in/out structure
      type (t_pbl_args), intent(inout) :: pbl_args
      !-- input:
      real*8, intent(in) :: coriol,utop,vtop,qtop,ztop
      real*8, intent(in) :: mdf
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0
      integer, intent(in) :: ilong,jlat,ihc,itype
      real*8, intent(in) :: dtdt_gcm,utop_old,vtop_old
      !-- output:
      real*8, intent(out) :: kms,kqs,z0h,z0q,w2_1
      real*8, intent(out) :: ufluxs,vfluxs,tfluxs,qfluxs
      !-- inout:
      real*8, intent(inout) :: z0m
      real*8, dimension(n),   intent(inout) :: u,v,t,q
      real*8, dimension(n-1), intent(inout) :: e
#if defined(TRACERS_ON)
!@var  tr local tracer profile (passive scalars)
      real*8, dimension(NTM), intent(in) :: trnradius,trndens,trnmm
      real*8, dimension(n,NTM), intent(inout) :: tr
#endif

c**** local vars for input from pbl_args
      real*8 :: evap_max,fr_sat,uocean,vocean,psurf,trhr0,tg,elhx,qsol
      real*8 :: dtime,sss_loc,dbl,ug,vg,tgrnd0,ttop,qgrnd_sat,qgrnd0
      real*8 :: tdns,qdns,tprime,qprime,tr4
      logical :: ocean,ddml_eq_1
      real*8 :: gusti
c**** local vars for output to pbl_args
      real*8 :: us,vs,ws,tsv,qsrf,khs,dskin,ustar,cm,ch,cq,wsgcm,wspdf
      real*8 :: ws0,lmonin
      real*8 :: ws_select
c**** other local vars
      real*8 :: qsat,deltaSST,tgskin,qnet,ts,rhosrf,qgrnd,delt,snow
      real*8 :: deltaSnowT,dqnetdtg,deltatg0
      real*8 :: tstar,qstar,ustar0,test,wstar3,wstar2h,tgrnd,ustar_oc
      real*8 :: bgrid,an2,as2,dudz,dvdz,tau,tgr4skin
      real*8 :: ws02,dm
      real*8, parameter ::  tol=1d-3,w=.5d0
      integer, parameter ::  itmax=5
      integer, parameter :: iprint=0,jprint=41  ! set iprint>0 to debug
      real*8, dimension(n) :: dz,xi,usave,vsave,tsave,qsave
     *       ,usave1,vsave1,tsave1,qsave1,tv
      real*8, dimension(n-1) :: lscale,dzh,xihat,kh,kq,ke,esave,esave1
      integer :: i,iter,ierr  !@var i,iter loop variable
C****
      REAL*8,DIMENSION(n) :: z
      REAL*8,DIMENSION(n-1) :: zhat,km,gm,gh
      REAL*8 :: lmonin_dry
#ifdef TRACERS_ON
      real*8, dimension(n,NTM) :: trsave
      real*8 trcnst,trsf,cqsave,byrho,rh1,evap,visc
      real*8, dimension(n-1) :: kqsave
      integer itr, ngx
#ifdef TRACERS_WATER
      real*8 :: trc2         ! could be passed out....
#ifdef TRACERS_SPECIAL_O18
      real*8 :: trc1,trs1    ! could be passed out....
      real*8 :: fac_cq_tr(NTM),fk
#endif
#endif
#ifndef TRACERS_TOMAS 
#ifdef TRACERS_DRYDEP
      real*8 vgs
      real*8 :: tr_dens, tr_radius ! variable tracer density and size
      logical hydrate
#endif
#endif 
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      INTEGER :: n1
      REAL*8 :: dsrcflx,dsrcflx2
#endif
#ifdef TRACERS_TOMAS
      INTEGER ss_bin,num_bin,du_bin,num_bin2,k,nx
      real*8 ss_num(nbins),dust_num(nbins),tot_dust,tot_seasalt
      real*8 ss_emis
#endif
      real*8 tg1
      if(xdelt /= 0d0) call stop_model(
     &     'PBL.f is not yet compatible with xdelt==deltx',255)

c**** get input from pbl_args structure
      dtime = pbl_args%dtsurf
      tgrnd0 = pbl_args%tgv
      ttop = pbl_args%tkv
      qgrnd_sat = pbl_args%qg_sat
      qgrnd0 = pbl_args%qg_aver
      evap_max = pbl_args%evap_max
      fr_sat = pbl_args%fr_sat
      uocean = pbl_args%uocean
      vocean = pbl_args%vocean
      psurf = pbl_args%psurf
      trhr0 = pbl_args%trhr0
      tg = pbl_args%tg
      tr4 =  pbl_args%tr4
      elhx = pbl_args%elhx
      qsol = pbl_args%qsol
      sss_loc = pbl_args%sss_loc
      ocean = pbl_args%ocean
      cm = pbl_args%cm
      ch = pbl_args%ch
      cq = pbl_args%cq
      dskin = pbl_args%dskin
      dbl = pbl_args%dbl
      khs = pbl_args%khs
      ug = pbl_args%ug
      vg = pbl_args%vg
      ddml_eq_1 = pbl_args%ddml_eq_1
      snow = pbl_args%snow
      
      call griddr(z,zhat,xi,xihat,dz,dzh,zgs,ztop,bgrid,n,ierr)
      if (ierr.gt.0) then
        print*,"advanc: i,j,ihc=",ilong,jlat,ihc
        print*,"advanc: itype,ztop,bgrid=",itype,ztop,bgrid
        call stop_model("PBL error in advanc",255)
      end if

      usave(:)=u(:)
      vsave(:)=v(:)
      tsave(:)=t(:)
      qsave(:)=q(:)
      esave(:)=e(:)

#ifdef TRACERS_ON
      trsave(:,1:pbl_args%ntx)=tr(:,1:pbl_args%ntx)
#endif
      do i=1,n-1
        usave1(i)=usave(i)
        vsave1(i)=vsave(i)
        tsave1(i)=tsave(i)
        qsave1(i)=qsave(i)
        esave1(i)=esave(i)
      end do
      ustar0=0.

      tgrnd=tgrnd0              ! use initial bulk ground temp
      qgrnd=qgrnd0              ! use initial sat humidity
      tgskin=tg                 ! initially assume no skin/bulk difference
      tgr4skin=tr4              ! initially assume no skin/bulk difference
      dskin=0
      ts=t(1)/(1+q(1)*xdelt)

      call getl1(e,zhat,dzh,lscale,n)

      if(ddml_eq_1) then
        tdns=pbl_args%tdns
        qdns=pbl_args%qdns
c       tprime=tdns-t(1)/(1.+xdelt*q(1))
c       qprime=qdns-q(1)
        tprime=tdns-ttop/(1.+xdelt*qtop)
        qprime=qdns-qtop
#ifdef TRACERS_ON
        pbl_args%trprime(1:pbl_args%ntx) =
     &     pbl_args%trdn1(1:pbl_args%ntx)-pbl_args%trtop(1:pbl_args%ntx)
#endif
      else ! either ddml(ilong,jlat).ne.1 or USE_PBL_E1
        tdns=0.d0
        qdns=0.d0
        tprime=0.d0
        qprime=0.d0
#ifdef TRACERS_ON
        pbl_args%trprime(:)=0
#endif
      endif
#ifdef TRACERS_TOMAS
      ss_bin=0
      num_bin=0
      du_bin=0
      num_bin2=0
#endif
      do iter=1,itmax

        call get_tv(t,q,tv,n)

        if(iter.gt.1) then
          call getl(e,u,v,tv,zhat,dzh,lmonin,ustar,lscale,dbl,n)
C**** adjust tgrnd/qgrnd for skin effects over the ocean & lakes & snow
          if ((itype.eq.1 .or. (itype.eq.2 .and. snow.gt.0)) .and.
     &         skin_effect.gt.0) then
c estimate net flux and ustar_oc from current tg,qg etc.
            ts=t(1)/(1+q(1)*xdelt)
            rhosrf=100.*psurf/(rgas*t(1)) ! surface air density
            Qnet= (lhe+tgskin*shv)*cq*rhosrf*(ws*(q(1)-qgrnd)
     &           +gusti*qprime)        ! Latent
     &           + sha*ch*rhosrf*(ws*(ts-tgskin)+gusti*tprime) ! Sensible
     &           +trhr0-stbo*tgr4skin ! LW
            dQnetdtg=shv*cq*rhosrf*(ws*(q(1)-qgrnd)+gusti*qprime)
     &           - sha*ch*rhosrf*ws-4d0*stbo*tgskin*tgskin*tgskin
            deltatg0=tgskin-tg
            
            if (itype.eq.1) then
               ustar_oc=ustar*sqrt(rhosrf*byrhows)
               dskin=deltaSST(Qnet,Qsol,ustar_oc)
               tgskin=0.5*(tgskin+(tg+dskin)) ! smooth changes in iteration
               tgskin=max(tgskin,tf+tfrez(sss_loc)) ! prevent unphysical values
#ifdef SNOW_SKIN_TEMP
            elseif (itype.eq.2 .and. snow.gt.0) then
               dskin=deltaSnowT(Qnet,Qsol,snow,tgskin,dQnetdtg,deltatg0)
               tgskin=0.5*(tgskin+(tg+dskin)) ! smooth changes in iteration
               tgskin=min(tgskin,tf) ! prevent unphysical values
#endif
            end if

            dskin=tgskin-tg     ! net dskin diagnostic
            tgr4skin=(sqrt(sqrt(tr4))+dskin)**4
            qgrnd=qsat(tgskin,elhx,psurf)
            if (ocean) qgrnd=0.98d0*qgrnd  ! use ocean adjustment
            tgrnd=tgskin*(1.+qgrnd*xdelt)
          endif
       endif

        call getk(km,kh,kq,ke,gm,gh,u,v,tv,e,lscale,dzh,n)
        call stars(ustar,tstar,qstar,lmonin,lmonin_dry,tgrnd,qgrnd,ts,
     2             u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,dm,
#ifdef TRACERS_SPECIAL_O18
     *             fac_cq_tr,
#endif
     3             km,kh,kq,dzh,itype,n)
#ifdef TRACERS_ON
        kqsave=kq
        cqsave=cq
#endif

        !@var wstar the convection-induced wind according to
        !@+ M.J.Miller et al. 1992, J. Climate, 5(5), 418-434, Eqs(6-7),
        !@+ for heat and mositure

        if(tv(2).lt.tv(1)) then !convective
          wstar3=-dbl*grav*2.*(tv(2)-tv(1))*kh(1)/((tv(2)+tv(1))*dzh(1))
          wstar2h = wstar3**twoby3
        else
          wstar3=0.
          wstar2h=0.
        endif

        call e_eqn(esave,e,u,v,tv,km,kh,ke,lscale,dz,dzh,
     2                 ustar,dtime,n)

ccc     call e_les(tstar,ustar,wstar3,dbl,lmonin,zhat,lscale,e,n)

        ! Inclusion of gustiness in surface fluxes
        ! Redelsperger et al. 2000, eqn(13), J. Climate, 13, 402-421
        gusti=pbl_args%gusti ! defined in PBL_DRV.f
        ws02=(u(1)-uocean)**2+(v(1)-vocean)**2+wstar2h
        ws0=sqrt(ws02)
        ws=sqrt(ws02+gusti*gusti)

        call q_eqn(qsave,q,kq,dz,dzh,cq,ws,qgrnd_sat,qtop,dtime,n
     &       ,evap_max,fr_sat,ws0,gusti,qprime,qdns,ddml_eq_1)

        call t_eqn(u,v,tsave,t,q,z,kh,kq,dz,dzh
     &       ,ch,ws,tgrnd,ttop,qtop,dtdt_gcm,dtime
     &       ,n,dpdxr,dpdyr,dpdxr0,dpdyr0,ws0,gusti,tprime,tdns
     &       ,qdns,ddml_eq_1)

        call uv_eqn(usave,vsave,u,v,z,km,dz,dzh
     &       ,ustar,cm,z0m,utop,vtop,utop_old,vtop_old
     &       ,dtime,coriol,ug,vg,uocean,vocean,n,dpdxr,dpdyr,dpdxr0
     &       ,dpdyr0)

        if ( ((ttop.ge.tgrnd).and.(lmonin_dry.le.0.)).or.
     &       ((ttop.le.tgrnd).and.(lmonin_dry.ge.0.)) )
     &     call tfix(t,z,ttop,tgrnd
     &              ,lmonin_dry,tstar,ustar,kh(1),n)

c       if(ddml_eq_1) then
c         tprime=tdns-t(1)/(1.+xdelt*q(1))
c         qprime=qdns-q(1)
c       endif

        test=abs(2.*(ustar-ustar0)/(ustar+ustar0))
        if (test.lt.tol) exit

        if (iter.lt.itmax) then
        do i=1,n-1
          u(i)=w*usave1(i)+(1.-w)*u(i)
          v(i)=w*vsave1(i)+(1.-w)*v(i)
          t(i)=w*tsave1(i)+(1.-w)*t(i)
          q(i)=w*qsave1(i)+(1.-w)*q(i)
          e(i)=w*esave1(i)+(1.-w)*e(i)
          usave1(i)=u(i)
          vsave1(i)=v(i)
          tsave1(i)=t(i)
          qsave1(i)=q(i)
          esave1(i)=e(i)
        end do
        end if
        ustar0=ustar
cc      WRITE(0,*) 'advanc iter t(1),q(1) ',iter,t(1),q(1)

      end do ! end of iter loop

      call get_tv(t,q,tv,n)

c**** cannot update ws without taking care that ws used for tracers is
c**** the same as that used for q
      wsgcm=sqrt((u(1)-uocean)**2+(v(1)-vocean)**2)

C**** Preliminary coding for use of sub-gridscale wind distribution
C**** generic calculations for all tracers
C**** To use, uncomment next two lines and adapt the next chunk for
C**** your tracers. The integrated wind value is passed back to SURFACE
C**** and GHY_DRV. This may need to be tracer dependent?
      if(calc_wspdf == 1) then
        delt = t(1)/(1.+q(1)*xdelt) - tgrnd/(1.+qgrnd*xdelt)
        CALL sig(e(1),mdf,dbl,delt,ch,wsgcm,t(1),pbl_args%wsubtke
     &       ,pbl_args%wsubwd,pbl_args%wsubwm)
        CALL get_wspdf(pbl_args%wsubtke,pbl_args%wsubwd,pbl_args%wsubwm
     &       ,pbl_args%mcfrac,wsgcm,wspdf)
      else
        wspdf = 0.
      endif
csgs      sig0 = sig(e(1),mdf,dbl,delt,ch,ws,t(1))
csgsC**** possibly tracer specific coding
csgs      wt = 3.                 ! threshold velocity
csgs      wmin = wt               ! minimum wind velocity (usually wt?)
csgs      wmax = 50.              ! maxmimum wind velocity
csgs      icase=3                 ! icase=3 ==> w^3 dependency
csgs      call integrate_sgswind(sig0,wt,wmin,wmax,ws,icase,wint)

#ifdef TRACERS_ON
C**** tracer calculations are passive and therefore do not need to
C**** be inside the iteration. Use moisture diffusivity.
C**** First, define some useful quantities
      ts=t(1)/(1.+q(1)*xdelt)   ! surface air temp (K)
      visc=visc_air(ts)         ! viscosity
      rhosrf=100.*psurf/(rgas*t(1)) ! surface air density
      byrho=1d0/rhosrf
      tg1 = tgskin-tf ! re-calculate ground T (C)
      rh1=q(1)/qsat(ts,lhe,psurf) ! rel. hum. at surface (wrt water)
      evap=cqsave*rhosrf*(ws*(qgrnd-q(1))-gusti*qprime)  ! net evap

#ifdef TRACERS_DRYDEP
C**** Get tracer deposition velocity (= 1 / bulk sfc resistance)
C**** for all dry deposited tracers
      call get_dep_vel(ilong,jlat,itype,lmonin,dbl,ustar,ts
     &     ,pbl_args%dep_vel,pbl_args%stomatal_dep_vel,trnmm
#ifdef TRACERS_TOMAS
     &     ,pbl_args%gs_vel,pbl_args%psurf
#endif
     &     )
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      CALL dust_emission_constraints(itype,wsgcm,pbl_args)
#endif

#ifdef WATER_PROPORTIONAL
      force_limit = .false.
#endif
C**** loop over tracers
      do itr=1,pbl_args%ntx
C**** Define boundary conditions

C****   1) default air mass tracers
        trcnst=pbl_args%trconstflx(itr)*byrho   ! convert to (conc * m/s)
        trsf=pbl_args%trsfac(itr)

#ifdef TRACERS_WATER
C****   2) water mass tracers
C**** Water tracers need to multiply trsfac/trconstflx by cq*Usurf
C**** and qgrnd_sat (moved from driver routines to deal with skin effects)
        if (tr_wd_TYPE(pbl_args%ntix(itr)).eq.nWATER) then
          trcnst=cqsave*(pbl_args%trconstflx(itr)*ws*qgrnd_sat
     *         -gusti*pbl_args%trprime(itr))
          trc2=pbl_args%trconstflx(itr)
!          if (ddml_eq_1) then   ! hmmm... gusti>0 even if ddml_eq_1=F
!            trsf=pbl_args%trsfac(itr)*cqsave*ws0
!          else
            trsf=pbl_args%trsfac(itr)*cqsave*ws
!          end if
          if (itype.eq.3) then  ! possible correction for large E over LI
            if (evap.gt.pbl_args%snow .and. pbl_args%snow.gt.snmin) then
              trc2 = (pbl_args%snow*pbl_args%trconstflx(itr)+(evap ! weighted mean tracer conc
     *             -pbl_args%snow)*pbl_args%trgrnd2(itr))/evap 
              trcnst=cqsave*(trc2*ws*qgrnd_sat
     *             -gusti*pbl_args%trprime(itr))
            end if
          end if
#ifdef TRACERS_SPECIAL_O18
C**** get fractionation for isotopes
          call get_frac(itype,ws,tg1,trcnst,trsf,evap*byrho
     $         ,trc2,q(1),pbl_args%ntix(itr)
     $         ,fac_cq_tr(itr),trc1,trs1,fk)
          pbl_args%frack(itr)=fk  ! save for output
          trcnst=trc1
          trsf  =trs1
#endif
        end if
#endif

#ifdef TRACERS_DRYDEP
C****   3) dry deposited tracers (including gravitational settling)
C**** Tracer Dry Deposition boundary condition for dry dep tracers:
        if(dodrydep(pbl_args%ntix(itr))) then
C****   get settling velocity
#ifndef TRACERS_TOMAS  
c****   set tracer size and density (not constant anymore)
        tr_radius = trnradius(pbl_args%ntix(itr))
        tr_dens =   trndens(pbl_args%ntix(itr))

C**** need to hydrate the sea salt before determining settling
          hydrate = (trname(pbl_args%ntix(itr)).eq.'seasalt1'.or.
     *         trname(pbl_args%ntix(itr)).eq.'seasalt2'.or.
     *         trname(pbl_args%ntix(itr)).eq.'M_SSA_SS'.or.
     *         trname(pbl_args%ntix(itr)).eq.'M_SSC_SS'.or.
     *         trname(pbl_args%ntix(itr)).eq.'M_SSS_SS')

          if (trnradius(pbl_args%ntix(itr)).gt.0.) then
            pbl_args%gs_vel(pbl_args%ntix(itr))=vgs(rhosrf,rh1
     &           ,tr_radius,tr_dens,visc,hydrate)
          else
            pbl_args%gs_vel(pbl_args%ntix(itr))=0.
          end if
#endif 
!TOMAS - gs_vel for TOMAS is computed in TRDRYDEP.f. 
          trsf=pbl_args%trsfac(itr)*(pbl_args%dep_vel(pbl_args%ntix(itr)
     &         )+pbl_args%gs_vel(pbl_args%ntix(itr)))
        end if
#endif

C****   4) tracers with interactive sources
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)
        select case (trname(pbl_args%ntix(itr)))

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS)
        case ('DMS')
          call read_DMS_sources(ws,itype,ilong,jlat,pbl_args%DMS_flux)
          trcnst=pbl_args%DMS_flux*byrho
#endif

#if (defined TRACERS_AEROSOLS_SEASALT) || (defined TRACERS_AMP)
        case ('seasalt1', 'M_SSA_SS')
          call read_seasalt_sources(ws,itype,1,ilong,jlat
     &         ,pbl_args%ss1_flux,trname(pbl_args%ntix(itr)))
          trcnst=pbl_args%ss1_flux*byrho
        case ('seasalt2', 'M_SSC_SS')
          call read_seasalt_sources(ws,itype,2,ilong,jlat
     &         ,pbl_args%ss2_flux,trname(pbl_args%ntix(itr)))
          trcnst=pbl_args%ss2_flux *byrho
        case ('M_SSS_SS')
          call read_seasalt_sources(ws,itype,1,ilong,jlat
     &         ,pbl_args%ss1_flux,trname(pbl_args%ntix(itr)))
          call read_seasalt_sources(ws,itype,2,ilong,jlat
     &         ,pbl_args%ss2_flux,trname(pbl_args%ntix(itr)))
          trcnst=(pbl_args%ss2_flux + pbl_args%ss1_flux)*byrho
#endif  /* TRACERS_AEROSOLS_SEASALT || TRACERS_AMP */

#ifdef TRACERS_AEROSOLS_OCEAN
        case ('OCocean')
          call read_seasalt_sources(ws,itype,1,ilong,jlat
     &         ,pbl_args%OCocean_flux,trname(pbl_args%ntix(itr)))
          trcnst=pbl_args%OCocean_flux*byrho
#endif  /* TRACERS_AEROSOLS_OCEAN */

#ifdef TRACERS_TOMAS
        case ('ANACL_01','ANACL_02','ANACL_03','ANACL_04', 
     &         'ANACL_05','ANACL_06','ANACL_07','ANACL_08',
     &         'ANACL_09','ANACL_10','ANACL_11','ANACL_12'
#ifdef TOMAS_12_3NM 
     *         ,'ANACL_13','ANACL_14','ANACL_15'
#endif
     &         )
        ss_bin=ss_bin+1
       
        call read_seasalt_sources(ws,itype,ss_bin,ilong,jlat
     &       ,ss_emis,trname(pbl_args%ntix(itr)))
        
        pbl_args%tomas_ss_flux(ss_bin)=ss_emis
        trcnst=pbl_args%tomas_ss_flux(ss_bin)*byrho
        ss_num(ss_bin)=trcnst/sqrt(xk(ss_bin)*xk(ss_bin+1)) 
#endif  /* TRACERS_TOMAS */
        end select
#endif

#ifdef BIOGENIC_EMISSIONS
! Nadine Unger test code for biogenic emissions
! emisop is in units kg C/m2/s
        select case (trname(pbl_args%ntix(itr)))
        case ('Isoprene')
          call isoprene_emission(ilong,jlat,itype,pbl_args%emisop
     &         ,pbl_args%canopy_temperature)
          trcnst=pbl_args%emisop*byrho
        end select
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
ccc dust emission from earth
        SELECT CASE (trname(pbl_args%ntix(itr)))
        case('Clay','Silt1','Silt2','Silt3','Silt4','Silt5','ClayIlli'
     &         ,'ClayKaol','ClaySmec','ClayCalc','ClayQuar','ClayFeld'
     &         ,'ClayHema','ClayGyps','ClayIlHe','ClayKaHe','ClaySmHe'
     &         ,'ClayCaHe','ClayQuHe','ClayFeHe','ClayGyHe','Sil1Quar'
     &         ,'Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps','Sil1Illi'
     &         ,'Sil1Kaol','Sil1Smec','Sil1QuHe','Sil1FeHe','Sil1CaHe'
     &         ,'Sil1GyHe','Sil1IlHe','Sil1KaHe','Sil1SmHe','Sil2Quar'
     &         ,'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps','Sil2Illi'
     &         ,'Sil2Kaol','Sil2Smec','Sil2QuHe','Sil2FeHe','Sil2CaHe'
     &         ,'Sil2GyHe','Sil2IlHe','Sil2KaHe','Sil2SmHe','Sil3Quar'
     &         ,'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps','Sil3Illi'
     &         ,'Sil3Kaol','Sil3Smec','Sil3QuHe','Sil3FeHe','Sil3CaHe'
     &         ,'Sil3GyHe','Sil3IlHe','Sil3KaHe','Sil3SmHe','Sil4Quar'
     &         ,'Sil4Feld','Sil4Calc','Sil4Hema','Sil4Gyps','Sil4Illi'
     &         ,'Sil4Kaol','Sil4Smec','Sil4QuHe','Sil4FeHe','Sil4CaHe'
     &         ,'Sil4GyHe','Sil4IlHe','Sil4KaHe','Sil4SmHe','Sil5Quar'
     &         ,'Sil5Feld','Sil5Calc','Sil5Hema','Sil5Gyps','Sil5Illi'
     &         ,'Sil5Kaol','Sil5Smec','Sil5QuHe','Sil5FeHe','Sil5CaHe'
     &         ,'Sil5GyHe','Sil5IlHe','Sil5KaHe','Sil5SmHe')
          n1=pbl_args%ntix(itr)-n_soildust+1
          CALL local_dust_emission( pbl_args%ntix(itr), wsgcm, pbl_args,
     &         dsrcflx, dsrcflx2 )
          trcnst=dsrcflx*byrho
          pbl_args%dust_flux(n1)=dsrcflx
          pbl_args%dust_flux2(n1)=dsrcflx2
        END SELECT
#endif

#ifdef TRACERS_AMP
        SELECT CASE (trname(pbl_args%ntix(itr)))
        CASE ('M_DD1_DU')
          DO n1=1,2
            CALL local_dust_emission(n1,wsgcm,pbl_args,dsrcflx,
     &           dsrcflx2)
            trcnst=dsrcflx*byrho
            pbl_args%dust_flux(n1)=dsrcflx
            pbl_args%dust_flux2(n1)=dsrcflx2
          END DO
        CASE ('M_DD2_DU')
          DO n1=3,4
            CALL local_dust_emission(n1,wsgcm,pbl_args,dsrcflx,
     &           dsrcflx2)
            trcnst=dsrcflx*byrho
            pbl_args%dust_flux(n1)=dsrcflx
            pbl_args%dust_flux2(n1)=dsrcflx2
          END DO
        CASE ('M_DDD_DU')
          DO n1=1,4
            CALL local_dust_emission(n1,wsgcm,pbl_args,dsrcflx,
     &           dsrcflx2)
            trcnst=dsrcflx*byrho
            pbl_args%dust_flux(n1)=dsrcflx
            pbl_args%dust_flux2(n1)=dsrcflx2
          END DO
        END SELECT
#endif
#ifdef TRACERS_TOMAS  
        SELECT CASE (trname(pbl_args%ntix(itr)))
        CASE ('ADUST_01','ADUST_02','ADUST_03','ADUST_04'
     &          ,'ADUST_05','ADUST_06','ADUST_07','ADUST_08'
     &          ,'ADUST_09','ADUST_10','ADUST_11','ADUST_12'
#ifdef TOMAS_12_3NM 
     *         ,'ADUST_13','ADUST_14','ADUST_15'
#endif
     &         )
           du_bin=du_bin+1
          if(du_bin.eq.1)then
             n1=1
                CALL local_dust_emission(n1,wsgcm,pbl_args,
     &              dsrcflx,dsrcflx2)
                pbl_args%dust_flux(n1)=dsrcflx
                pbl_args%dust_flux2(n1)=dsrcflx2
             n1=2             !silt1
                CALL local_dust_emission(n1,wsgcm,pbl_args,
     &              dsrcflx, dsrcflx2)
                pbl_args%dust_flux(n1)=dsrcflx
                pbl_args%dust_flux2(n1)=dsrcflx2
             n1=3               !silt2
                CALL local_dust_emission(n1,wsgcm,pbl_args,
     &              dsrcflx, dsrcflx2)
                pbl_args%dust_flux(n1)=dsrcflx
                pbl_args%dust_flux2(n1)=dsrcflx2 
            n1=4              !silt3
                CALL local_dust_emission(n1,wsgcm,pbl_args,
     &              dsrcflx, dsrcflx2)
                pbl_args%dust_flux(n1)=dsrcflx
                pbl_args%dust_flux2(n1)=dsrcflx2                 
          endif

          trcnst=pbl_args%dust_flux(1)*byrho*scalesizeclay(du_bin)
     &       +sum(pbl_args%dust_flux(2:4))*byrho*scalesizesilt(du_bin)

          dust_num(du_bin)=trcnst/sqrt(xk(du_bin)*xk(du_bin+1))
          
!TOMAS - Silt3 is out of size range for TOMAS. 
        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04',
     &         'ANUM__05','ANUM__06','ANUM__07','ANUM__08',
     &         'ANUM__09','ANUM__10','ANUM__11','ANUM__12'
#ifdef TOMAS_12_3NM 
     *         ,'ANUM__13','ANUM__14','ANUM__15'
#endif
     &         )
           num_bin2=num_bin2+1
           trcnst=dust_num(num_bin2)+ss_num(num_bin2)  
        END SELECT
#endif

        ngx=gasex_index%getindex(pbl_args%ntix(itr))
        if (pbl_args%ntix(itr)==n_co2n) then
          IF (ocean) THEN  ! OCEAN only
          !note trcnst is already multiplied by byrho in TRACERS_GASEXCH_ocean_CO2_PBL

            if (obio_wspdf) then
              ws_select=wspdf
            else
              ws_select=ws
            endif
#ifdef TRACERS_GASEXCH_ocean
#ifdef TRACERS_GASEXCH_ocean_CO2
            call TRACERS_GASEXCH_ocean_CO2_PBL(tg1,ws_select,
     .          pbl_args%sss_loc,psurf,tr_mm(pbl_args%ntix(itr)),
     .          pbl_args%trconstflx(itr),
     .          byrho,pbl_args%Kw_gas(ngx),pbl_args%alpha_gas(ngx),
     .          pbl_args%beta_gas(ngx),trsf,trcnst,ilong,jlat)
#endif                         
#else
            !call stop_model('gas exchange code missing', 255)
            ! no ocean gas exchange code - assuming trivial BCs
            trsf = 0.d0
            trcnst = 0.d0
#endif                         

!     write(*,'(a,2i5,4e12.4,i5,7e12.4)')'PBL:', 
!    .  ilong,jlat,tg1,ws,pbl_args%sss_loc,psurf,itr,
!    .             pbl_args%trconstflx(itr),
!    .          byrho,pbl_args%Kw_gas(ngx),pbl_args%alpha_gas(ngx),
!    .          pbl_args%beta_gas(ngx),trsf,trcnst
          ELSE
            trsf = 0.d0
            trcnst = 0.d0
          ENDIF
        else if (pbl_args%ntix(itr)==n_cfcn) then
          IF (ocean) THEN  ! OCEAN only
#ifdef TRACERS_GASEXCH_ocean
#ifdef TRACERS_GASEXCH_ocean_CFC
            call TRACERS_GASEXCH_ocean_CFC_PBL(tg1,ws,
     .          pbl_args%sss_loc,psurf,tr_mm(pbl_args%ntix(itr)),
     .          pbl_args%trconstflx(itr),
     .          byrho,pbl_args%Kw_gas(ngx),pbl_args%alpha_gas(ngx),
     .          pbl_args%beta_gas(ngx),trsf,trcnst,ilong,jlat)
#endif
#else
            call stop_model('gas exchange code missing', 255)
#endif
          ENDIF
        endif

C**** solve tracer transport equation
        call tr_eqn(trsave(1,itr),tr(1,itr),kqsave,dz,dzh,trsf
     *       ,trcnst,pbl_args%trtop(itr),
#ifdef TRACERS_WATER
     *       pbl_args%tr_evap_max(itr),fr_sat,
#endif
     *       dtime,n)

#ifdef TRACERS_DRYDEP
C**** put in a check to prevent unphysical solutions. If too much
C**** tracer is being taken out, replace profile with linear one
C**** with maximum allowed flux.
        if (dodrydep(pbl_args%ntix(itr))) then
          if ((trsf*tr(1,itr)-trcnst)*dtime
     &         .gt.pbl_args%trtop(itr)*ztop) then
            do i=1,n
              tr(i,itr)=(pbl_args%trtop(itr)*ztop/dtime+trcnst)/trsf
     &             +(i-1) *pbl_args%trtop(itr)/float(n-1)
            end do
          end if
        end if
#endif
      end do
#endif

      us  = u(1)
      vs  = v(1)
      tsv = t(1)
      qsrf  = q(1)

      kms = km(1)
      khs = kh(1)
      kqs = kq(1)

      ufluxs=km(1)*(u(2)-u(1))/dzh(1)
      vfluxs=km(1)*(v(2)-v(1))/dzh(1)
      tfluxs=kh(1)*(t(2)-t(1))/dzh(1)
      qfluxs=kq(1)*(q(2)-q(1))/dzh(1)

      an2=2.*grav*(tv(n)-tv(n-1))/((tv(n)+tv(n-1))*dzh(n-1))
      dudz=(u(n)-u(n-1))/dzh(n-1)
      dvdz=(v(n)-v(n-1))/dzh(n-1)
      as2=dudz*dudz+dvdz*dvdz
      tau=B1*lscale(n-1)/max(sqrt(2.*e(n-1)),teeny)
      !@var w2_1 the vertical component of 2*e at GCM layer 1
      w2_1=twoby3*e(n-1)-tau*by3*(s7*km(n-1)*as2+s8*kh(n-1)*an2)
      w2_1=max(0.24d0*e(n-1),w2_1) ! Mellor-Yamada 1982

c     call check1(ustar,1,ilong,jlat,2)

c Diagnostics printed at a selected point:

      if ((ilong.eq.iprint).and.(jlat.eq.jprint)) then
        call output(u,v,t,q,e,lscale,z,zhat,dzh,
     2              km,kh,kq,ke,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3              ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4              utop,vtop,ttop,qtop,
     5              dtime,bgrid,ilong,jlat,iter,itype,n)
      endif

c**** copy output to pbl_args
      pbl_args%us = us
      pbl_args%vs = vs
      pbl_args%ws = ws    ! wind speed including buoyancy + gusti
      pbl_args%ws0 = ws0  ! wind speed including buoyancy
      pbl_args%tsv = tsv
      pbl_args%qsrf = qsrf
      pbl_args%cm = cm
      pbl_args%ch = ch
      pbl_args%cq = cq
      pbl_args%dskin = dskin
      !!pbl_args%psi = psi   ! maybe compute it here ?
      pbl_args%dbl = dbl
      pbl_args%khs = khs
      pbl_args%ustar = ustar
      pbl_args%lmonin = lmonin
      pbl_args%zgs = zgs
      pbl_args%gusti = gusti
      pbl_args%tprime=tprime
      pbl_args%qprime=qprime

      pbl_args%wsgcm=wsgcm
      pbl_args%wspdf=wspdf

C**** tracer code output
#ifdef TRACERS_ON
      pbl_args%trs(1:pbl_args%ntx) = tr(1,1:pbl_args%ntx)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      pbl_args%z(:) = z(:)
      pbl_args%zhat(:) = zhat(:)
      pbl_args%km(:) = km(:)
      pbl_args%gm(:) = gm(:)
      pbl_args%gh(:) = gh(:)
#endif
#endif

      return
      end subroutine advanc

      subroutine stars(
     &     ustar,tstar,qstar,lmonin,lmonin_dry,tgrnd,qgrnd,ts,
     2                 u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,dm,
#ifdef TRACERS_SPECIAL_O18
     *                 fac_cq_tr,
#endif
     3                 km,kh,kq,dzh,itype,n)
!@sum computes the friction speed, ustar, the virtual potential
!@+  temperature scale, tstar, and the specific humidity scale,
!@+  qstar. Note that
!@+  surface momentum flux = ustar*ustar
!@+  surface heat flux     = ustar*tstar
!@+  surface moisture flux = ustar*qstar
!@+  It also calculates and outputs the Monin-Obukov length, lmonin, 
!@+  the roughness lengths (z0m,z0h,z0q), the drag coefficient (cm),
!@+  the Stanton number (ch) and the Dalton number (cq)
!@+  by calling subroutine dflux.
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@var USTAR the friction speed
!@var TSTAR the virtual potential temperature scale
!@var QSTAR the moisture scale
!@var LMONIN the Monin-Obukhov length scale
      USE CONSTANT, only : teeny
      implicit none

      integer, intent(in) :: itype,n
      real*8, dimension(n), intent(in) :: u,v,t,q,z
      real*8, dimension(n-1), intent(in) :: km,kh,kq,dzh
      real*8, intent(in) :: tgrnd,qgrnd,ts
      real*8, intent(inout) :: z0m
      real*8, intent(out) :: ustar,tstar,qstar,lmonin,lmonin_dry
      real*8, intent(out) :: z0h,z0q,cm,ch,cq,dm
#ifdef TRACERS_SPECIAL_O18
      real*8, intent(out) :: fac_cq_tr(NTM)
#endif


      real*8 dz,vel1,du1,dv1,dudz,dtdz,dqdz,zgs
      real*8 tflx,qflx,tvflx,tgrndv,tv(2),tstarv,dtv1

      tgrndv = tgrnd*(1.+deltx*qgrnd)
      call get_tv(t,q,tv,2)

      dz     = dzh(1)
      vel1   = sqrt(u(1)*u(1)+v(1)*v(1))
      du1=u(2)-u(1)
      dv1=v(2)-v(1)
      dudz=sqrt(du1*du1+dv1*dv1)/dz
      tflx   = kh(1)*(t(2)-t(1))/dz
      qflx   = kq(1)*(q(2)-q(1))/dz
      tvflx = tflx*(1.+deltx*q(1)) + t(1)*deltx*qflx
      ustar  = sqrt(km(1)*dudz)
      ustar  = max(ustar,ustar_min)
      tstar  = tflx/ustar
      if (abs(tstar).gt.smax*abs(t(1)-tgrnd)) tstar=smax*(t(1)-tgrnd)
      if (abs(tstar).lt.smin*abs(t(1)-tgrnd)) tstar=smin*(t(1)-tgrnd)
      qstar  = qflx/ustar
      if (abs(qstar).gt.smax*abs(q(1)-qgrnd)) qstar=smax*(q(1)-qgrnd)
      if (abs(qstar).lt.smin*abs(q(1)-qgrnd)) qstar=smin*(q(1)-qgrnd)
      tstarv  = tvflx/ustar
      dtv1 = tv(1)-tgrndv
      if (abs(tstarv).gt.smax*abs(dtv1)) tstarv=smax*dtv1
      if (abs(tstarv).lt.smin*abs(dtv1)) tstarv=smin*dtv1

      zgs    = z(1)

      if (ustar.gt.smax*vel1) ustar=smax*vel1
      if (ustar.lt.smin*vel1) ustar=smin*vel1
      if (tstar.eq.0.) tstar=teeny
      if (tstarv.eq.0.) tstarv=teeny
      ustar  = max(ustar,ustar_min)

      lmonin_dry = ustar*ustar*tgrnd/(kappa*grav*tstar)
      if(abs(lmonin_dry).lt.lmonin_min)
     &     lmonin_dry=sign(lmonin_min,lmonin_dry)
      if(abs(lmonin_dry).gt.lmonin_max)
     &     lmonin_dry=sign(lmonin_max,lmonin_dry)
      
      lmonin = ustar*ustar*tgrndv/(kappa*grav*tstarv)
      if(abs(lmonin).lt.lmonin_min) lmonin=sign(lmonin_min,lmonin)
      if(abs(lmonin).gt.lmonin_max) lmonin=sign(lmonin_max,lmonin)

c**** To compute the drag coefficient,Stanton number and Dalton number
      call dflux(lmonin,ustar,vel1,ts,z0m,z0h,z0q,zgs,cm,ch,cq,dm,
#ifdef TRACERS_SPECIAL_O18
     *     fac_cq_tr,
#endif
     *     itype)

      return
      end subroutine stars

      subroutine getl1(e,zhat,dzh,lscale,n)
!@sum getl1 estimates the master length scale, lscale, of the
!@+  turbulence model on the secondary grid, zhat.
!@auth Ye Cheng/G. Hartke
      implicit none

      integer, intent(in) :: n     !@var n  array dimension
      real*8, dimension(n-1), intent(in) :: e,zhat,dzh
      real*8, dimension(n-1), intent(out) :: lscale
      real*8, parameter :: alpha=0.2d0

      real*8 :: sum1,sum2,l0,l1
      integer :: j

      sum1=0.
      sum2=0.
      do j=1,n-1
        sum1=sum1+sqrt(e(j))*zhat(j)*dzh(j)
        sum2=sum2+sqrt(e(j))*dzh(j)
      end do
      l0=alpha*sum1/sum2
      if (l0.lt.zhat(1)) l0=zhat(1)

      do j=1,n-1
        l1=kappa*zhat(j)
        lscale(j)=l0*l1/(l0+l1)
      end do

      return
      end subroutine getl1

      subroutine getl(e,u,v,t,zhat,dzh,lmonin,ustar,lscale,dbl,n)
!@sum getl computes the master length scale, lscale, of the
!@+  turbulence model on the secondary grid, zhat,
!@+  using the formulas by Nakanishi (2001) from the LES data.
!@auth Ye Cheng/G. Hartke
!@var e z-profle of turbulent kinetic energy
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var lscale z-profile of the turbulent dissipation length scale
!@var z vertical grids (main, meter)
!@var zhat vertical grids (secondary, meter)
!@var dzh(j)  z(j+1)-z(j)
!@var dbl PBL height (m)
!@var n number of vertical subgrid main layers

      USE CONSTANT, only : by3
      implicit none

      integer, intent(in) :: n   !@var n  array dimension
      real*8, dimension(n-1), intent(in) :: e,zhat,dzh
      real*8, dimension(n), intent(in) :: u,v,t
      real*8, dimension(n-1), intent(out) :: lscale
      real*8, intent(in) :: dbl,lmonin,ustar

      integer :: i   !@var i  array dimension
      real*8 kz,l0,ls,lb,an2,an,qty,qturb,zeta

      l0=.3d0*dbl ! Moeng and Sullivan 1994

      if (l0.lt.zhat(1)) l0=zhat(1)

      kz=kappa*zhat(1)
      lscale(1)=l0*kz/(l0+kz)

      do i=2,n-1
        kz=kappa*zhat(i)
        zeta=zhat(i)/lmonin
        ! Nakanishi (2001)
        if(zeta.ge.1.) then
          ls=kz/3.7d0
        elseif(zeta.ge.0.) then
          ls=kz/(1.+2.7d0*zeta)
        else
          ls=kz*(1.-100.*zeta)**0.2d0
        endif
        if (t(i+1).gt.t(i)) then
          an2=2.*grav*(t(i+1)-t(i))/((t(i+1)+t(i))*dzh(i))
          an=sqrt(an2)
          qturb=sqrt(2*e(i))
          if(zeta.ge.0.) then
             lb=qturb/an
          else
             qty=(ustar/((-kappa*lmonin*l0*l0)**by3*an))**0.5d0
             lb=qturb*(1.+5.*qty)/an
          endif
        else
          lb=1.d30
        endif
        lscale(i)=l0*ls*lb/(l0*ls+l0*lb+ls*lb)
      end do

      return
      end subroutine getl

      subroutine dflux(lmonin,ustar0,vsurf,ts,z0m,z0h,z0q,zgs,
     *                 cm,ch,cq,dm,
#ifdef TRACERS_SPECIAL_O18
     *                 fac_cq_tr,
#endif
     *                 itype)
!@sum dflux computes the dimensionless surface fluxes of momentum,
!@+  heat and moisture (drag coefficient Cm , Stanton number Ch,
!@+  and Dalton number Cq), with explicit Schmidt number (Sc) and
!@+  Prandtl number (Pr) dependence
!@+  and flexibility for water isotopes.
!@+  It also computes the roughness lengths for momentum, z0m
!@+  (for itype=1 or 2, i.e., surface type ocean or seaice),
!@+  for temperature, z0h, and for water vapor, z0q, all in meters
!@+  (Hartke and Rind, 1997).
!@+  It is called from within subroutine stars.
!@auth  Ye Cheng/G. Hartke (mods by G. Schmidt)
!@var lmonin = Monin-Obukhov length (m)
!@var ustar  = friction speed (sqrt of surface momentum flux) (m/sec)
!@var vsurf  = total surface wind speed, used to limit ustar -> cm
!@var ts     = surface air temperature (K)
!@var zgs    = height of the surface layer (m)
!@var itype  = integer identifying surface type
!@var z0m   = momentum roughness length, prescribed (itype=3,4) (m)
!@var z0m   = roughness length for momentum, computed (itype=1,2)
!@var cm    = drag coefficient for momentum
!@var ch    = Stanton number
!@var cq    = Dalton number
!@var z0h   = roughness length for temperature (m)
!@var z0q   = roughness length for water vapor (m)
!@var Sc    = Schmidt (no relation) number (visc_air_kin/diff)
!@var Pr    = Prandtl number (visc_air_kin/therm_diff)
!@var fac_cq_tr = ratio of cq for water isotopes = f(Sc_tr)

      USE CONSTANT, only : visc_air_kin
#ifdef SCM
      use SCM_COM, only : SCMopt,SCMin
#endif
      implicit none

      real*8,  intent(in) :: lmonin,ustar0,vsurf,zgs,ts
      integer,  intent(in) :: itype
      real*8,  intent(inout) :: z0m
      real*8,  intent(out) :: cm,ch,cq,dm,z0h,z0q
#ifdef TRACERS_SPECIAL_O18
      real*8, intent(out) :: fac_cq_tr(NTM)
      real*8 :: cq_tr(NTM),z0q_tr(NTM),Sc_tr,get_diff_rel
      integer :: itr
#endif
      real*8 :: nu,num,nuh,nuq,dpsim,dpsih,dpsiq
      real*8 ustar,dum
      real*8, parameter :: Sc=0.595d0, Pr=0.71d0

C**** Kinematic viscosity
      nu=visc_air_kin(ts)       ! temperature dependent

      if ((itype.eq.1).or.(itype.eq.2)) then
c *********************************************************************
        ustar = max(ustar0,1.0125d-5)  ! make sure not too small
c Compute roughness lengths using smooth/rough surface formulation:

c       z0m=0.11d0*nu/ustar+0.011d0*ustar*ustar*bygrav ! COARE algorithm
        z0m=0.135d0*nu/ustar+0.018d0*ustar*ustar*bygrav ! Hartke and Rind (1996)

#ifdef SCM
        if( SCMopt%z0m ) z0m=SCMin%z0m ! specified roughness length
#endif

        call getzhq(ustar,z0m,Pr,nu,1.4d-5,z0h)  ! heat
        call getzhq(ustar,z0m,Sc,nu,1.3d-4,z0q)  ! vapour

#ifdef TRACERS_SPECIAL_O18
C**** calculate different z0q for different diffusivities
          do itr=1,NTM
            if (tr_wd_TYPE(itr).eq.nWater) then
              Sc_tr=Sc*get_diff_rel(itr)
              call getzhq(ustar,z0m,Sc_tr,nu,1.3d-4,z0q_tr(itr))
            end if
          end do
#endif

      else
c *********************************************************************
c  For land and land ice, z0m is specified. For z0h and z0q,
c    empirical evidence suggests:
        z0h=z0m*.13533528d0    ! = exp(-2.)
        z0q=z0h
#ifdef TRACERS_SPECIAL_O18
        z0q_tr(:)=z0q
#endif
c *********************************************************************
      endif

      call  getcm(zgs,z0m,lmonin,dm,dpsim,cm)
      call getchq(zgs,z0h,lmonin,dm,dpsih,ch)
      call getchq(zgs,z0q,lmonin,dm,dpsiq,cq)

#ifdef TRACERS_SPECIAL_O18
      do itr=1,NTM
        if (tr_wd_TYPE(itr).eq.nWater)
!     *       call getchq(zgs,z0m,lmonin,dm,z0q_tr(itr),cq_tr(itr),dum)
     *       call getchq(zgs,z0m,lmonin,dm,dpsim,cq_tr(itr))
      end do
      do itr=1,NTM
        if (tr_wd_TYPE(itr).eq.nWater)
     *       fac_cq_tr(itr)=cq_tr(itr)/cq_tr(n_Water)
      end do
#endif

      return
      end subroutine dflux

      subroutine getzhq(ustar,z0m,ScPr,nu,z0min,z0hq)
!@sum getzhq calculates the roughness lengths for heat (z0h)
!@+  and for humidity (z0q), 
!@+  modified from Eqs 5.24, 5.27 and 5.35 in Brutsaert (1982).
!@+  It is called from within subroutine dflux.
!!*** remove z0min for original HR97 code

      implicit none
      real*8, intent(in) :: ustar,z0m,ScPr,z0min,nu
      real*8, intent(out) :: z0hq
      real*8 r0q,beta,fac_smooth_ScPr,fac_rough_ScPr

C**** functional dependence on Sc,Pr for smooth, rough surfaces
      fac_smooth_ScPr = 30.*exp(-13.6d0*kappa*ScPr**twoby3)

      if (ustar.le.0.20d0) then ! smooth regime
        z0hq=nu*fac_smooth_ScPr/ustar + z0min
      else                      ! rough regime
        fac_rough_ScPr = -7.3d0*kappa*sqrt(ScPr)
        r0q=sqrt(sqrt(ustar*z0m/nu))
        z0hq=7.4d0*z0m*exp(fac_rough_ScPr*r0q)
        if (ustar.lt.0.2d0) then ! intermediate regime (lin. interp.)
          beta=(ustar-0.02d0)/0.18d0
          z0hq=(1.-beta)*(nu*fac_smooth_ScPr/ustar+z0min)+beta*z0hq
        endif
      endif

      return
      end subroutine getzhq

      subroutine getcm(z,z0,lmonin,dm,dpsim,cm)
!@sum calculates the drag coefficient (cm) for momentum flux
!@+  (Hartke and Rind, 1997; Zeng et al., 1998; updated by Y. Cheng, 2014).
!@+   It is called from within subroutine dflux.
      implicit none
      real*8, intent(in) :: z,z0,lmonin
      real*8, intent(out) :: dm,dpsim,cm

      real*8 zet,zet0,x,x0,xm

      zet=z/lmonin
      zet0=z0/lmonin

      call find_dpsim(zet,zet0,dpsim)
      dm=max(log(z/z0)-dpsim,1.d-3)
      cm=kappa*kappa/(dm*dm)
      if (cm.gt.cmax) cm=cmax
      if (cm.lt.cmin) cm=cmin
      return
      end subroutine getcm

      subroutine getchq(z,z0,lmonin,dm,dpsih,ch)
!@sum calculates the Stanton number or Dalton number (ch or cq)
!@+  for heat and latent heat fluxes
!@+  (Hartke and Rind, 1997; Zeng et al., 1998; updated by Y. Cheng, 2014).
!@+   It is called from within subroutine dflux.
      implicit none
      real*8, intent(in) :: z,z0,lmonin,dm
      real*8, intent(out) :: dpsih,ch

      real*8 zet,zet0,x,x0,xh,dh

      zet=z/lmonin
      zet0=z0/lmonin

      call find_dpsih(zet,zet0,z,z0,dpsih)
      dh=max(log(z/z0)-dpsih,1.d-3)
      ch=kappa*kappa/(dm*dh)
      if (ch.gt.cmax) ch=cmax
      if (ch.lt.cmin) ch=cmin
      return
      end subroutine getchq

      subroutine simil(z,z0m,z0h,z0q,lmonin,ustar,tstar,qstar,tg,qg
     2                 ,u,t,q,dpsim,dpsih,dpsiq)
!@sum   simil calculates the similarity solutions for wind speed,
!@+     virtual potential temperature, and moisture mixing ratio
!@+     at height z.
!@auth  Ye Cheng/G. Hartke
!@var     z       height above ground at which solution is computed (m)
!@var     z0m     momentum roughness height (m)
!@var     z0h     temperature roughness height (m)
!@var     z0q     moisture roughness height (m)
!@var     lmonin  Monin-Obukhov length scale (m)
!@var     ustar   friction speed (m/sec)
!@var     tg      ground temperature (K)
!@var     qg      ground moisture mixing ratio
!@var     u       computed similarity solution for wind speed (m/sec)
!@var     t       computed similarity solution for virtual potential
!@+               temperature (K)
!@var     q       computed similarity solution for moisture mixing ratio
      implicit none

      real*8,  intent(in) :: z,z0m,z0h,z0q,lmonin,ustar,tstar,qstar
     &                      ,tg,qg
      real*8,  intent(out) :: u,t,q,dpsim,dpsih,dpsiq
      real*8 dm,cm,ch,cq
      call  getcm(z,z0m,lmonin,dm,dpsim,cm)
      call getchq(z,z0h,lmonin,dm,dpsih,ch)
      call getchq(z,z0q,lmonin,dm,dpsiq,cq)
      u=   (ustar/kappa)*(log(z/z0m)-dpsim)
      t=tg+(tstar/kappa)*(log(z/z0h)-dpsih)
      q=qg+(qstar/kappa)*(log(z/z0q)-dpsiq)
      return
      end subroutine simil

      subroutine griddr(z,zhat,xi,xihat,dz,dzh,z1,zn,bgrid,n,ierr)
!@sum griddr computes altitudes on vertical grid. The xi coordinates are
!@+  uniformly spaced and are mapped in a log-linear fashion onto the
!@+  z grid. (The z's are the physical coords.) Also computes the
!@+  altitudes on the secondary grid, zhat, and the derivatives
!@+  dxi/dz evaluated at both all z and zhat. z and zhat are staggered:
!@+  mean quantitied are calculated at z, turbulent kinetic enery
!@+  and fluxes are calculated at zhat.
!@auth  Ye Cheng/G. Hartke

c     Grids:
c
c                n   - - - - - - - - - - - - -
c                    -------------------------  n-1
c                n-1 - - - - - - - - - - - - -
c                    -------------------------  j+1
c     (main,z)   j+1 - - - - - - - - - - - - -
c     (z)            -------------------------  j     (secondary, zhat)
c                j   - - - - - - - - - - - - -
c                    -------------------------  j-1
c                j-1 - - - - - - - - - - - - -
c                    -------------------------    2
c                2   - - - - - - - - - - - - -
c                    -------------------------    1
c                1   - - - - - - - - - - - - -
c
c     dz(j)==zhat(j)-zhat(j-1), dzh(j)==z(j+1)-z(j)
!@var bgrid determines how strongly non-linear the
!@+   mapping is. BGRID=0 gives linear mapping. Increasing BGRID
!@+   packs more points into the bottom of the layer.
!@+   Bgrid is calculated in the beginning of this subroutine every
!@+   time this subroutine is called. zs is the first grid height
!@+   and dzs is the grid separation near the surface. The values of
!@+   parameters byzs(=1/zs) and bydzs(=1/dzs) are obtained by
!@+   calling an old version of this subroutine with bgrid=0.2927
!@+   and has been tested as appropriate for pe(2)=934mb.
!@+   Now we impose that for other pe(2)s, zs and dzs be the same
!@+   as when pe(2)=934mb, so as to maintain the balance between
!@+   the accuracy and stability.
!@+   The new value of bgrid is then calculated below which
!@+   also depends on ztop (i.e., depends on pe(2)).
!@var  z       height of main grids (meter)
!@var  zhat    height of secondary grids (meter)
!@var  xi      an uniformly spaced coordinate mapped to z
!@var  xihat   an uniformly spaced coordinate mapped to zhat
!@var  dz   dxi/(dxi/dz)
!@var  dzh  dxi/(dxi/dzh)
!@var  dxi  (ztop - zbottom)/(n-1)
!@var  ierr Error reporting flag
      use RootFinding_mod, only: NewtonMethod
      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(out) :: z,xi,dz
      real*8, dimension(n-1), intent(out) :: zhat,xihat,dzh
      integer, intent(out) :: ierr
      real*8, intent(in) :: z1,zn
      real*8, intent(out) :: bgrid

      real*8, parameter ::  tolz=1d-3
     &  ,byzs=1.d0/10.d0,bydzs=1.d0/4.7914d0
      real*8 z1pass,znpass,b,xipass,lznbyz1
      common /grids_99/z1pass,znpass,b,xipass,lznbyz1
      external fgrid2
      integer i  !@var i,iter loop variable
      real*8 dxi,zmin,zmax,dxidz,dxidzh

      z1pass=z1
      znpass=zn
      dxi=(zn-z1)/float(n-1)
      bgrid=max((dxi*bydzs-1.)/((zn-z1)*byzs-log(zn/z1)),0.d0)

      b=bgrid
      zmin=z1
      zmax=zn

      do i=1,n-1
        xi(i)=z1+(zn-z1)*float(i-1)/float(n-1)
        xihat(i)=z1+(zn-z1)*(float(i)-0.5)/float(n-1)
      end do
      xi(n)=zn
      z(1)=z1

      lznbyz1 = log(zn/z1)
      dxidz=1.+bgrid*((zn-z1)/z1-lznbyz1)
      dz(1)=dxi/dxidz
      xipass=xihat(1)
      zhat(1)=NewtonMethod(fgrid2,zmin,zmax,tolz,ierr)
      if (ierr.gt.0) return
      dxidzh=1.+bgrid*((zn-z1)/zhat(1)-lznbyz1)
      dzh(1)=dxi/dxidzh

      do i=2,n-1
        xipass=xi(i)
        z(i)=NewtonMethod(fgrid2,zmin,zmax,tolz,ierr)
        if (ierr.gt.0) return
        xipass=xihat(i)
        zhat(i)=NewtonMethod(fgrid2,zmin,zmax,tolz,ierr)
        if (ierr.gt.0) return
        dxidz=1.+bgrid*((zn-z1)/z(i)-lznbyz1)
        dxidzh=1.+bgrid*((zn-z1)/zhat(i)-lznbyz1)
        dz(i)=dxi/dxidz
        dzh(i)=dxi/dxidzh
      end do
      z(n)=zn
      dxidz=1.+bgrid*((zn-z1)/zn-lznbyz1)
      dz(n)=dxi/dxidz

      return
      end subroutine griddr

      subroutine tfix(t,z,ttop,tgrnd,lmonin,tstar,ustar,khs,n)
!@sum tfix linearly interpolates between the ground temperature tgrnd and the
!@+  virtual potential temperature at the middle of the first GCM layer 
!@+  to reset the T(z) profile.
!@+  It is called when the T(z) profile becomes irregular.
!@auth  Ye Cheng/G. Hartke

      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n),intent(in) :: z
      real*8, dimension(n),intent(inout) :: t
      real*8, intent(in) :: ttop,tgrnd,ustar,khs
      real*8, intent(inout) :: lmonin,tstar

      real*8 dtdz
      integer i  !@var i loop variable

      do i=1,n-1
        t(i)=tgrnd+(ttop-tgrnd)*z(i)/z(n)
      end do

      dtdz   = (t(2)-t(1))/(z(2)-z(1))
      tstar  = khs*dtdz/ustar
      if (abs(tstar).gt.smax*abs(t(1)-tgrnd)) tstar=smax*(t(1)-tgrnd)
      if (abs(tstar).lt.smin*abs(t(1)-tgrnd)) tstar=smin*(t(1)-tgrnd)

      lmonin = ustar*ustar*tgrnd/(kappa*grav*tstar)
      if(abs(lmonin).lt.lmonin_min) lmonin=sign(lmonin_min,lmonin)
      if(abs(lmonin).gt.lmonin_max) lmonin=sign(lmonin_max,lmonin)


      return
      end subroutine tfix

      subroutine ccoeff0
!@sum ccoeff0 sets/calculates model coefficients for the
!@+  GISS 2002 turbulence model (Cheng et al., 2002).
!@auth  Ye Cheng
      implicit none

      ! temperary variable
      real*8 :: del,aa,bb,cc,tmp

      prt=     0.82d0
      b1=     19.3d0
      b2=     15.8d0
      b123=b1**(2./3.)
      g1=       .1070d0
      g2=       .0032d0
      g3=       .0864d0
      g4=       .1000d0
      g5=     11.04d0
      g6=       .786d0
      g7=       .643d0
      g8=       .547d0
c
      d1=(7.*g4/3+g8)/g5
      d2=(g3**2-g2**2/3.)-1./(4.*g5**2)*(g6**2-g7**2)
      d3=g4/(3.*g5**2)*(4.*g4+3.*g8)
      d4=g4/(3.*g5**2)*(g2*g6-3.*g3*g7-g5*(g2**2-g3**2))
     &   +g8/g5*(g3**2-g2**2/3.)
      d5=-1./(4.*g5**2)*(g3**2-g2**2/3)*(g6**2-g7**2)
      s0=g1/2.
      s1=-g4/(3.*g5**2)*(g6+g7)+2.*g4/(3.*g5)*(g1-g2/3.-g3)
     &   +g1/(2.*g5)*g8
      s2=-g1/(8.*g5**2)*(g6**2-g7**2)
      s4=2./(3.*g5)
      s5=2.*g4/(3.*g5**2)
      s6=2./(3.*g5)*(g3**2-g2**2/3)-g1/(2.*g5)*(g3-g2/3.)
     &   +g1/(4*g5**2)*(g6-g7)

      s7=3.*g3-g2   !@ useful to calculate w2
      s8=4.*g4      !@ useful to calculate w2

c     find rimax:

      c1=s5+2*d3
      c2=s1-s6-2*d4
      c3=-s2+2.*d5
      c4=s4+2.*d1
      c5=-s0+2.*d2

      rimax=(c2+sqrt(c2**2-4.*c1*c3))/(2.*c1)
      rimax=int(rimax*1000.)/1000.

c     find gm_at_rimax

      aa=c1*rimax*rimax-c2*rimax+c3
      bb=c4*rimax+c5
      cc=2.d0
      if(abs(aa).lt.1d-8) then
         gm_at_rimax= -cc/bb
      else
         tmp=bb*bb-4.*aa*cc
         gm_at_rimax=(-bb-sqrt(tmp))/(2.*aa)
      endif
      ! rimax=.96,  gm_at_rimax=.1366285d6

c     find ghmin,ghmax,gmmax0:

      del=(s4+2.*d1)**2-8.*(s5+2.*d3)
      ghmin=(-s4-2.*d1+sqrt(del))/(2.*(s5+2.*d3))
      ghmin=int(ghmin*10000.)/10000.
      ghmax=(b1*0.53d0)**2
      gmmax0=(b1*0.34d0)**2  ! not in use yet

      ! for level 3 model only:
      g0=2./3.
      d1_3=(7.*g4/3.)/g5
      d2_3=d2
      d3_3=g4/(3.*g5**2)*(4.*g4)
      d4_3=g4/(3.*g5**2)*(g2*g6-3*g3*g7-g5*(g2**2-g3**2))
      d5_3=d5
      s0_3=s0
      s1_3=-g4/(3.*g5**2)*(g6+g7)+2.*g4/(3.*g5)*(g1-g2/3.-g3)
      s2_3=s2
      s3_3=g0*g4/g5*(g3+g2/3.+1./(2.*g5)*(g6+g7))
      s4_3=s4
      s5_3=s5
      s6_3=s6

      return
      end subroutine ccoeff0

      subroutine get_tv(t,q,tv,n)
!@sum get_tv converts temperature T to virtual temperature Tv
      USE CONSTANT, only : deltx
      implicit none
      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: t,q
      real*8, dimension(n), intent(out) :: tv
      integer :: i
      do i=1,n
        tv(i) = t(i)*(1.+deltx*q(i))
      enddo
      return
      end subroutine get_tv

      subroutine getk(km,kh,kq,ke,gma,gha,u,v,t,e,lscale,dzh,n)
!@sum getk computes the turbulent diffusivities for momentum, Km,
!@+  for heat, Kh, for moisture, Kq and for kinetic energy, Ke,
!@+  at the secondary grids,
!@+  using the GISS second order closure model (Cheng et al., 2000).
!@+  u,v,t,q,ke are calculated at the primary grid z, while
!@+  e,lscale,km,kh,gm,gh are calculated at the secondary grid zhat.
!@auth  Ye Cheng/G. Hartke
c     Grids:
c
c                n   - - - - - - - - - - - - -
c                    -------------------------  n-1
c                n-1 - - - - - - - - - - - - -
c                    -------------------------  j+1
c     (main,z)   j+1 - - - - - - - - - - - - -
c     (z)            -------------------------  j     (secondary, zhat)
c                j   - - - - - - - - - - - - -
c                    -------------------------  j-1
c                j-1 - - - - - - - - - - - - -
c                    -------------------------    2
c                2   - - - - - - - - - - - - -
c                    -------------------------    1
c                1   - - - - - - - - - - - - -
c
c     dz(j)==zhat(j)-zhat(j-1), dzh(j)==z(j+1)-z(j)
c     at main: u,v,t,q,ke
c     at edge: e,lscale,km,kh,gm,gh
!@var u,v,t,e,lscale,t_real z-profiles
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var km turbulent viscosity for u and v equations
!@var kh turbulent conductivity for t and q equations
!@var ke turbulent diffusivity for e equation
!@var gma normalized velocity gradient, tau**2*as2
!@var gha normalized temperature gradient, tau**2*an2
!@var tau B1*lscale/sqrt(2*e)
!@var as2 shear squared, (dudz)**2+(dvdz)**2
!@var an2 Brunt-Vaisala frequency, grav/T*dTdz
!@var se stability constant for e, adjustable

      USE CONSTANT, only : teeny

      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t
      real*8, dimension(n-1), intent(in) :: e,lscale,dzh
      real*8, dimension(n-1), intent(out) :: km,kh,kq,ke,gma,gha

      real*8 :: an2,dudz,dvdz,as2,ell,den,qturb,tau,gh,gm,gmmax,sm,sh
     &  ,sq,taue
      integer :: i !@var i loop variable

      do i=1,n-1
        an2=2.*grav*(t(i+1)-t(i))/((t(i+1)+t(i))*dzh(i))
        dudz=(u(i+1)-u(i))/dzh(i)
        dvdz=(v(i+1)-v(i))/dzh(i)
        as2=dudz*dudz+dvdz*dvdz
        ell=lscale(i)
        qturb=sqrt(2.*e(i))
        tau=B1*ell/max(qturb,teeny)
        gh=tau*tau*an2
        gm=tau*tau*as2
        if(gh.lt.ghmin) gh=ghmin
        if(gh.gt.ghmax) gh=ghmax
        gmmax=(1+d1*gh+d3*gh*gh)/(d2+d4*gh)
        if(gm.gt.gmmax) gm=gmmax
        den=1.+d1*gh+d2*gm+d3*gh*gh+d4*gh*gm+d5*gm*gm
        sm=(s0+s1*gh+s2*gm)/den
        sh=(s4+s5*gh+s6*gm)/den
        sq=sh
        taue=tau*e(i)
        km(i)=min(max(taue*sm,kmmin),k_max)
        kh(i)=min(max(taue*sh,khmin),k_max)
        kq(i)=min(max(taue*sq,kqmin),k_max)
        ke(i)=min(max(taue*se,kemin),k_max)
        gma(i)=gm
        gha(i)=gh
      end do
      return
      end subroutine getk

      subroutine e_eqn(esave,e,u,v,t,km,kh,ke,lscale,
     &                     dz,dzh,ustar,dtime,n)
!@sum e_eqn integrates differential eqns for
!@+  the turbulent kinetic energy, e, using tridiagonal method over 
!@+  npbl-1(=7) sublayer edges (i.e., the secondary grids).
!@+  The boundary condition near the bottom is:
!@+  e(1)=(1/2)*B1**(2/3)*ustar**2.
!@+  At the top secondary grid, nearest to the middle of the
!@+  first GCM layer, e is prescribed.
!@auth Ye Cheng/G. Hartke
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var km z-profile of turbulent viscosity
!@var kh z-profile of turbulent conductivity
!@var ke z-profile of turbulent diffusion in eqn for e
!@var lscale z-profile of the turbulent dissipation length scale
!@var z vertical grids (main, meter)
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var dtime time step
!@var ustar friction velocity at the surface
!@var n number of vertical subgrid main layers
      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, intent(in) :: ustar, dtime
      real*8, dimension(n), intent(in) :: u,v,t,dz
      real*8, dimension(n-1), intent(in) :: esave,km,kh,ke,lscale,dzh
      real*8, dimension(n-1), intent(inout) :: e

      real*8 :: an2,dudz,dvdz,as2,qturb,ri,aa,bb,cc,gm,tmp
      integer :: j !@var j loop variable
c
c     sub(j)*e_jm1_kp1+dia(j)*e_j_kp1+sup(j)*e_jp1_kp1 = rhs(j)
c     from kirk:/u/acyxc/papers/2ndOrder/maple/phik.1,
c       sub(j)=-dtime*aj*P1a_jm1_k/dxi^2
c       sup(j)=-dtime*aj*P1a_jp1_k/dxi^2
c       dia(j)=1-(sub(j)+dia(j))+dtime*P3_j_k
c       rhs(j)=PHI_j_k + dtime*P4_j_k
c       where P1a == P1*a
c       aj=dxi/dz(j) if j refers to the primary grid
c       aj=dxi/dzh(j) if j refers to the secondary grid
c       now for e_eqn j refers to the secondary grid, therefore:
c       aj = dxi/dzh(j)
c       a(j-1/2) = dxi/dzh(j-1/2)=dxi/dz(j)
c       a(j+1/2) = dxi/dzh(j+1/2)=dxi/dz(j+1)
c       p1(j-1/2)=0.5*(ke(j)+ke(j-1))
c       p1(j+1/2)=0.5*(ke(j)+ke(j+1))
c       ke(j)=sq*lscale(j)*qturb, qturb=sqrt(2.*e(j))
c
      do j=2,n-2
          qturb=sqrt(2.*e(j))
          sub(j)=-dtime*0.5d0*(ke(j)+ke(j-1))/(dzh(j)*dz(j))
          sup(j)=-dtime*0.5d0*(ke(j)+ke(j+1))/(dzh(j)*dz(j+1))
          dia(j)=1.-(sub(j)+sup(j))+dtime*2*qturb/(b1*lscale(j))
          an2=2*grav*(t(j+1)-t(j))/((t(j+1)+t(j))*dzh(j))
          dudz=(u(j+1)-u(j))/dzh(j)
          dvdz=(v(j+1)-v(j))/dzh(j)
          as2=dudz*dudz+dvdz*dvdz
          rhs(j)=esave(j)+dtime*(km(j)*as2-kh(j)*an2)
       end do

      dia(1)=1.
      sup(1)=0.
      rhs(1)=0.5d0*b123*ustar*ustar

c      j=n-1
c      an2=2.*grav*(t(j+1)-t(j))/((t(j+1)+t(j))*dzh(j))
c      dudz=(u(j+1)-u(j))/dzh(j)
c      dvdz=(v(j+1)-v(j))/dzh(j)
c      as2=max(dudz*dudz+dvdz*dvdz,teeny)
c      ri=an2/as2
c      if(ri.gt.rimax) ri=rimax
c      aa=c1*ri*ri-c2*ri+c3
c      bb=c4*ri+c5
c      cc=2.d0
c      if(abs(aa).lt.1d-8) then
c        gm= -cc/bb
c      else
c        tmp=bb*bb-4.*aa*cc
c        gm=(-bb-sqrt(tmp))/(2.*aa)
c      endif
c      sub(n-1)=0.
c      dia(n-1)=1.
c      rhs(n-1)=max(0.5d0*(B1*lscale(j))**2*as2/max(gm,teeny),teeny)

      sub(n-1)=-1.
      dia(n-1)=1.
      rhs(n-1)=0.

      call TRIDIAG(sub,dia,sup,rhs,e,n-1)

      do j=1,n-1
         e(j)=min(max(e(j),teeny),emax)
      end do

      Return
      end subroutine e_eqn

      subroutine e_les(tstar,ustar,wstar3,dbl,lmonin,zhat,lscale,e,n)
!@sum e_les finds e according to the parameterization of les data
!@Ref Moeng and Sullivan 1994, J. Atmos. Sci., 51, 999-1022.
!@Ref Cheng et al. 2002, J. Atmos. Sci., 59, 1550-1565.
!@auth  Ye Cheng
!@ver   1.0
!@var (see subroutine k_gcm)
      USE CONSTANT, only : by3

      implicit none

      integer, intent(in) :: n   !@var n  array dimension
      real*8, intent(in) :: tstar,ustar,wstar3,dbl,lmonin
      real*8, dimension(n), intent(in) :: zhat,lscale
      real*8, dimension(n), intent(inout) :: e
      integer :: j !@var j loop variable
      real*8 :: tvflx,ustar3,zj,kz,zet,phim,eps,ej

      tvflx=ustar*tstar
      ustar3=ustar*ustar*ustar
      do j=1,n-1
        zj=zhat(j)
        kz=kappa*zj
        if(zj.le.dbl) then
          zet=zj/lmonin
          call find_phim(zet,phim)
          eps=.4d0*wstar3/dbl+ustar3*(1.-zj/dbl)*phim/kz
          ej=.5d0*(19.3d0*lscale(j)*eps)**(2.*by3)
          ej=min(max(ej,emin),emax)
        else
          ej=emin
        endif
        e(j)=max(e(j),ej) ! e(j) on the rhs is an input
c       e(j)=ej
      end do
      return
      end subroutine e_les

      subroutine t_eqn(u,v,t0,t,q,z,kh,kq,dz,dzh,ch,usurf,tgrnd
     &     ,ttop,qtop,dtdt_gcm,dtime,n
     &     ,dpdxr,dpdyr,dpdxr0,dpdyr0,usurf0,gusti,tprime,tdns
     &     ,qdns,ddml_eq_1)
!@sum t_eqn integrates differential eqns for 
!@+  the virtual potential temperature, T, using tridiagonal method 
!@+  over npbl(=8) sublayers between the surface (sublayer 1)
!@+  and the middle of the first GCM layer (sublayer npbl).
!@+  The boundary condition at the bottom is:
!@+   kh * dt/dz = ch * ( usurf*(t1 - tgrnd)
!@+               +(1+xdelt*qtop)*gusti*tprime )
!@+  which includes the effects on the surface flux
!@+  due to the moist convection wind gustiness and the
!@+  downdraft temperature perturbation
!@+  (Redelsperger et al. 2000; Emanuel and Zivkovic 1999),
!@+  where 
!@+  tprime=tdns-ttop/(1+xdelt*qtop),
!@+  t1, q1 are the T and Q at the surface, 
!@+  and tdns is the downdraft temperature in K at (i,j), which is
!@+  calculated in subroutines CONDSE (in CLOUDS2_DRV.f) and 
!@+  PBL (in PBL_DRV.f).
!@+  At the top, i.e., the middle of the first GCM layer, 
!@+  T is prescribed.
!@+  ### the following term was removed from BC at the bottom
!@+  ###   + xdelt * t1/(1+xdelt*q1) * kq * dqdz
!@auth Ye Cheng/G. Hartke
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature ref. to the surface
!@+          (if xdelt=0, t is the actual temperature)
!@var q z-profle of specific humidity
!@var t0 z-profle of t at previous time step
!@var kh z-profile of heat conductivity
!@var kq z-profile of moisture diffusivity
!@var z vertical grids (main, meter)
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var ch  dimensionless heat flux at surface (stanton number)
!@var usurf effective surface velocity
!@var tgrnd virtual potential temperature at the ground
!@+          (if xdelt=0, tgrnd is the actual temperature)
!@var ttop virtual potential temperature at the first GCM layer
!@+          (if xdelt=0, ttop is the actual temperature)
!@var dtdt_gcm ttop tendency from processes other than turbulence
!@var dtime time step
!@var n number of vertical subgrid main layers

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: u,v,t0,q,z,dz
      real*8, dimension(n-1), intent(in) :: dzh,kh,kq
      real*8, dimension(n), intent(inout) :: t
      real*8, intent(in) :: ch,tgrnd
      real*8, intent(in) :: ttop,qtop,dtdt_gcm,dtime,usurf
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0
      real*8, intent(in) ::  usurf0,gusti,tprime,tdns,qdns
      logical, intent(in) :: ddml_eq_1

      real*8 :: facth0,facth,factx,facty
      integer :: i  !@var i loop variable

      do i=2,n-1
         sub(i)=-dtime/(dz(i)*dzh(i-1))*kh(i-1)
         sup(i)=-dtime/(dz(i)*dzh(i))*kh(i)
         dia(i)=1.-(sub(i)+sup(i))
      end do

cccc
c changes for using aturb with SCM
#ifdef SCM
c     factx=(dpdxr-dpdxr0)/(z(n)-z(1))
c     facty=(dpdyr-dpdyr0)/(z(n)-z(1))
      do i=2,n-1
c       rhs(i)=t0(i)-dtime*t(i)*bygrav*(v(i)*facty+u(i)*factx)
        rhs(i) = t0(i)
      end do
#else
      factx=(dpdxr-dpdxr0)/(z(n)-z(1))
      facty=(dpdyr-dpdyr0)/(z(n)-z(1))
      do i=2,n-1
#ifdef PBL_USES_GCM_TENDENCIES
        rhs(i)=t0(i)+dtime*dtdt_gcm ! maybe scale by a function of height
#else
        rhs(i)=t0(i)-dtime*t(i)*bygrav*(v(i)*facty+u(i)*factx)
#endif
      end do
#endif

      facth0  = ch*dzh(1)/kh(1)
      facth  = usurf*facth0
      sup(1) = -1.

      if(ddml_eq_1) then
c        rat = usurf0/(usurf+teeny)
         dia(1) = 1+facth
c        dia(1) = 1+facth*rat
!     &             +xdelt*kq(1)*(q(2)-q(1))/(kh(1)*(1.+xdelt*q(1)))
c        rhs(1) = facth*(tgrnd-(1.d0-rat)*(1.+xdelt*qdns)*tdns)
c        rhs(1) = facth*(tgrnd-gusti/(usurf+teeny)
c    &            *((1.+xdelt*qtop)*tdns-ttop))
         rhs(1) = facth0*(usurf*tgrnd
     &            -gusti*((1.+xdelt*qtop)*tdns-ttop))
      else
         dia(1) = 1+facth
!     &             +xdelt*kq(1)*(q(2)-q(1))/(kh(1)*(1.+xdelt*q(1)))
         rhs(1) = facth*tgrnd
      endif

      dia(n) = 1.
      sub(n) = 0.
      rhs(n) = ttop

      call TRIDIAG(sub,dia,sup,rhs,t,n)

      return
      end subroutine t_eqn

      subroutine q_eqn(q0,q,kq,dz,dzh,cq,usurf,qgrnd,qtop,dtime,n
     &     ,flux_max,fr_sat,usurf0,gusti,qprime,qdns,ddml_eq_1)
!@sum q_eqn integrates differential eqns for
!@+  the specific humidity, Q, using tridiagonal method over npbl(=8)
!@+  sublayers between the surface (sublayer 1)
!@+  and the middle of the first GCM layer (sublayer npbl).
!@+  The boundary condition at the bottom is:
!@+    kq * dq/dz = min ( cq * usurf * (q1 - qgrnd)
!@+                     + cq * gusti * qprime ,
!@+            fr_sat * ( cq * usurf * (q1 - qgrnd)
!@+                     + cq * gusti * qprime )
!@+        - ( 1 - fr_sat ) * flux_max ),
!@+  which includes the effects on the surface flux
!@+  due to the moist convection wind gustiness and the
!@+  downdraft specific humidity  perturbation
!@+  (Redelsperger et al. 2000; Emanuel and Zivkovic 1999),
!@+  where qprime=qdns-qtop, qtop is Q at the first layer
!@+  and qdns is the downdraft humidity in kg/kg, (i,j), which is
!@+  calculated in subroutines CONDSE (in CLOUDS2_DRV.f) 
!@+  and PBL (in PBL_DRV.f).
!@+  At the top, i.e., the middle of the first GCM layer,
!@+  Q is prescribed.
!@auth Ye Cheng/G. Hartke
!@var q z-profle of specific humidity
!@var q0 z-profle of q at previous time step
!@var kq z-profile of moisture diffusivity
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var cq  dimensionless moisture flux at surface (dalton number)
!@var usurf effective surface velocity
!@var qgrnd specific humidity at the ground
!@var qtop specific humidity at the first GCM layer
!@var dtime time step
!@var n number of vertical subgrid main layers
!@var flux_max maximal flux from the unsaturated soil
!@var fr_sat fraction of the saturated soil

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: q0,dz
      real*8, dimension(n-1), intent(in) :: dzh,kq
      real*8, dimension(n), intent(out) :: q
      real*8, intent(in) :: cq,qgrnd,qtop,dtime,usurf
      real*8, intent(in) :: flux_max,fr_sat
      real*8, intent(in) :: usurf0,gusti,qprime,qdns
      logical, intent(in) :: ddml_eq_1

      real*8 :: factq0,factq
      integer :: i  !@var i loop variable

      do i=2,n-1
         sub(i)=-dtime/(dz(i)*dzh(i-1))*kq(i-1)
         sup(i)=-dtime/(dz(i)*dzh(i))*kq(i)
         dia(i)=1.-(sub(i)+sup(i))
      end do

      do i=2,n-1
        rhs(i)=q0(i)
      end do

      factq0  = cq*dzh(1)/kq(1)
      factq  = usurf*factq0
      sup(1) = -1.

      if(ddml_eq_1) then
c        rat = usurf0/(usurf+teeny)
c        dia(1) = 1.+factq*rat
c        rhs(1)= factq*(qgrnd-(1.-rat)*qdns)
         dia(1) = 1.+factq
c        rhs(1)= factq*(qgrnd-gusti/(usurf+teeny)*(qdns-qtop))
         rhs(1)= factq0*(usurf*qgrnd-gusti*(qdns-qtop))
      else
         dia(1) = 1.+factq
         rhs(1) = factq*qgrnd
      endif

      dia(n) = 1.
      sub(n) = 0.
      rhs(n) = qtop

      call TRIDIAG(sub,dia,sup,rhs,q,n)

c**** Now let us check if the computed flux doesnt exceed the maximum
c**** for unsaturated fraction

      if ( fr_sat .ge. 1. ) return   ! all soil is saturated
      if ( cq * usurf * (qgrnd - q(1)) - cq * gusti * qprime
     &   .le. flux_max ) return

c**** Flux is too high, have to recompute with the following boundary
c**** conditions at the bottom:
c**** kq * dq/dz = fr_sat * ( cq * usurf * (q(1) - qgrnd)
c****                       + cq * gusti * qprime )
c****             - ( 1 - fr_sat ) * flux_max

      if(ddml_eq_1) then
c        rat = usurf0/(usurf+teeny)
c        rat = -gusti/(usurf+teeny)
c        dia(1) = 1. + fr_sat*factq*rat
c        rhs(1)= fr_sat*factq*(qgrnd-(1.-rat)*qdns)
c    &            + (1.-fr_sat)*flux_max*dzh(1)/kq(1)
         dia(1) = 1. + fr_sat*factq
c        rhs(1)= fr_sat*factq*(qgrnd-gusti/(usurf+teeny)*(qdns-qtop))
c    &            + (1.-fr_sat)*flux_max*dzh(1)/kq(1)
         rhs(1)= fr_sat*factq0*(usurf*qgrnd-gusti*(qdns-qtop))
     &            + (1.-fr_sat)*flux_max*dzh(1)/kq(1)
      else
         dia(1) = 1. + fr_sat*factq
         rhs(1)= fr_sat*factq*qgrnd
     &            + (1.-fr_sat)*flux_max*dzh(1)/kq(1)
      endif

      call TRIDIAG(sub,dia,sup,rhs,q,n)

      return
      end subroutine q_eqn

      subroutine tr_eqn(tr0,tr,kq,dz,dzh,sfac,constflx,trtop,
#ifdef TRACERS_WATER
     *     tr_evap_max,fr_sat,
#endif
     *     dtime,n)
!@sum tr_eqn integrates differential eqns for
!@+  the tracers, TR, using tridiagonal method over npbl(=8)
!@+  sublayers between the surface (sublayer 1)
!@+  and the middle of the first GCM layer (sublayer npbl).
!@+  The boundary condition at the bottom is:
!@+  kq * dtr/dz = sfac * trs - constflx,
!@+  i.e. for moisture, sfac=cq*usurf, constflx=cq*usurf*qg,
!@+  to get:  kq * dq/dz = cq * usurf * (qs - qg);
!@+  for new moisture (including downdraft effects),
!@+  sfac=cq*(usurf-dusurf), constflx=cq*(usurf*qg + dusurf*qdns),
!@+  or sfac=cq*usurf0, constflx=cq*(usurf*(qg+qdns)-usurf0*qdns),
!@+  to get:  kq * dq/dz = cq*(usurf*(qs-qg) + dusurf*(qdns-qs)).
!@+  This should be flexible enough to deal with most situations.
!@+  At the top, i.e., the middle of the first GCM layer,
!@+  TR is prescribed.
!@auth Ye Cheng/G. Hartke
!@var tr z-profle of tracer concentration
!@var tr0 z-profle of tr at previous time step
!@var kq z-profile of moisture diffusivity
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var sfac factor multiplying surface tracer conc. in b.c.
!@var constflx the constant component of the surface tracer flux
!@var trtop tracer concentration at the first GCM layer
!@var dtime time step
!@var n number of vertical subgrid main layers
#ifdef TRACERS_WATER
!@var tr_evap_max  max amount of possible tracer evaporation
!@var fr_sat fraction of the saturated soil
#endif
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: tr0,dz
      real*8, dimension(n-1), intent(in) :: dzh,kq
      real*8, dimension(n), intent(out) :: tr
      real*8, intent(in) :: sfac,constflx,trtop,dtime
#ifdef TRACERS_WATER
      real*8, intent(in) :: tr_evap_max, fr_sat
#endif
      real*8 :: facttr
      integer :: i  !@var i loop variable

      do i=2,n-1
        sub(i)=-dtime/(dz(i)*dzh(i-1))*kq(i-1)
        sup(i)=-dtime/(dz(i)*dzh(i))*kq(i)
        dia(i)=1.-(sub(i)+sup(i))
      end do

      do i=2,n-1
        rhs(i)=tr0(i)
      end do

      facttr  = dzh(1)/kq(1)

      dia(1) = 1.+sfac*facttr
      sup(1) = -1.
      rhs(1)= facttr*constflx

      dia(n)  = 1.
      sub(n)  = 0.
      rhs(n) = trtop

#ifdef WATER_PROPORTIONAL
      if( .not. force_limit ) then
#endif
      call TRIDIAG(sub,dia,sup,rhs,tr,n)
#ifdef WATER_PROPORTIONAL
      endif
#endif

#ifdef TRACERS_WATER
c**** Check as in q_eqn if flux is limited, and if it is, recalculate
c**** profile

      if ( fr_sat .ge. 1. ) return   ! all soil is saturated

#ifdef WATER_PROPORTIONAL
      if (.not.force_limit) then
#endif

      if ( constflx - sfac * tr(1) .le. tr_evap_max ) return

#ifdef WATER_PROPORTIONAL
      endif
      force_limit = .true.
#endif

c**** Flux is too high, have to recompute with the following boundary
c**** conditions at the bottom:
c**** kq * dq/dz = fr_sat * (sfac * trs - constflx)
c****              + ( 1 - fr_sat ) * tr_evap_max

      dia(1) = 1. + fr_sat*sfac*facttr
      sup(1) = -1.
      rhs(1)= fr_sat*facttr*constflx +
     *     (1.-fr_sat)*tr_evap_max*dzh(1)/kq(1)

      call TRIDIAG(sub,dia,sup,rhs,tr,n)
#endif

      return
      end subroutine tr_eqn

      subroutine uv_eqn(u0,v0,u,v,z,km,dz,dzh,
     2                  ustar,cm,z0m,utop,vtop,utop_old,vtop_old,
     &                  dtime,coriol,
     3                  ug,vg,uocean,vocean,n
     &                  ,dpdxr,dpdyr,dpdxr0,dpdyr0)
!@sum uv_eqn integrates differential eqns for
!@+  mean velocity u and v using tridiagonal method over npbl(=8)
!@+  sublayers between the surface (sublayer 1)
!@+  and the middle of the first GCM layer (sublayer npbl).
!@+  The boundary condition at the bottom is:
!@+  km * du/dz = cm * usurf * u and
!@+  km * dv/dz = cm * usurf * v.
!@+  At the top, i.e., the middle of the first GCM layer,
!@+  u, v are prescribed.
!@auth Ye Cheng/G. Hartke
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var u0 z-profle of u at previous time step
!@var v0 z-profle of v at previous time step
!@var km z-profile of turbulent viscosity
!@var kh z-profile of turbulent conductivity
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var ustar friction velocity at the surface
!@var cm  dimensionless  momentum flux at surface (drag coeff.)
!@var ch  dimensionless heat flux at surface (stanton number)
!@var cq  dimensionless moisture flux at surface (dalton number)
!@var z0m  roughness height for momentum (if itype=1 or 2)
!@var utop due east component of the wind at the first GCM layer
!@var vtop due north component of the wind at the first GCM layer
!@var dtime time step
!@var coriol the Coriolis parameter
!@var ug due east component of the geostrophic wind
!@var vg due north component of the geostrophic wind
!@var utop_old utop after turb. diffusion during previous timestep
!@var vtop_old vtop after turb. diffusion during previous timestep
!@var n number of vertical subgrid main layers
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs,rhs1

      real*8, dimension(n), intent(in) :: u0,v0,z,dz
      real*8, dimension(n), intent(inout) :: u,v
      real*8, dimension(n-1), intent(in) :: km,dzh
      real*8, intent(in) :: ustar,cm,z0m,utop,vtop,utop_old,vtop_old
     &     ,dtime,coriol,ug,vg
     &     ,uocean,vocean
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0

      real*8 :: factx,facty,dpdx,dpdy,usurf,factor
      integer :: i,j,iter  !@var i,j,iter loop variable

      do i=2,n-1
         sub(i)=-dtime/(dz(i)*dzh(i-1))*km(i-1)
         sup(i)=-dtime/(dz(i)*dzh(i))*km(i)
         dia(i)=1.-(sub(i)+sup(i))
      end do

c#ifdef PBL_USES_GCM_TENDENCIES
cc Estimate the sum of non-coriolis non-PBL terms in the
cc momentum equations, termed here "dpdx" and "dpdy":
cc (utop-utop_old)/dtime = -dpdx + coriol*(vtop+vtop_old)/2
cc (vtop-vtop_old)/dtime = -dpdy - coriol*(utop+utop_old)/2
c      dpdx = +coriol*.5*(vtop+vtop_old) - (utop-utop_old)/dtime
c      dpdy = -coriol*.5*(utop+utop_old) - (vtop-vtop_old)/dtime
c#endif

#ifndef SCM
      factx=(dpdxr-dpdxr0)/(z(n)-z(1))
      facty=(dpdyr-dpdyr0)/(z(n)-z(1))
#endif
      do i=2,n-1
ccc for SCM use ug and vg and not dpdx,dpdy
#ifdef SCM
        rhs(i)=u0(i)+dtime*coriol*(v(i)-vg)
        rhs1(i)=v0(i)-dtime*coriol*(u(i)-ug)
#else
c#ifdef PBL_USES_GCM_TENDENCIES
cc        rhs(i)=u0(i)+dtime*(coriol*v(i)-dpdx)
cc        rhs1(i)=v0(i)-dtime*(coriol*u(i)+dpdy)
c        rhs(i) =u0(i)+(utop-utop_old)
c        rhs1(i)=v0(i)+(vtop-vtop_old)
c#else
        dpdx=factx*(z(i)-z(1))+dpdxr0
        dpdy=facty*(z(i)-z(1))+dpdyr0
        rhs(i)=u0(i)+dtime*(coriol*v(i)-dpdx)
        rhs1(i)=v0(i)-dtime*(coriol*u(i)+dpdy)
c#endif /* PBL_USES_GCM_TENDENCIES */
#endif
      end do

      usurf  = sqrt((u(1)-uocean)**2+(v(1)-vocean)**2)
      factor = cm*usurf*dzh(1)/km(1)

      dia(1) = 1.+factor
      sup(1) = -1.
      rhs(1)  = factor*uocean
      rhs1(1) = factor*vocean

      dia(n) = 1.
      sub(n) = 0.
      rhs(n)  = utop
      rhs1(n)  = vtop

      call TRIDIAG(sub,dia,sup,rhs,u,n)
      call TRIDIAG(sub,dia,sup,rhs1,v,n)

      return
      end subroutine uv_eqn

      subroutine t_eqn_sta(t,q,kh,kq,dz,dzh,ch,usurf,tgrnd,ttop,n)
!@sum t_eqn_sta computes the static solutions of
!@+ the virtual potential temperature, T,
!@+ between the surface and the first GCM layer.
!@+ The boundary condition at the bottom is:
!@+ kh * dt/dz = ch * usurf * (t - tg).
!@+ At the top, T is prescribed.
!@+ It is called only at the initialization
!@+ (from within subroutine inits).
!@auth Ye Cheng/G. Hartke
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var q z-profle of specific humidity
!@var kh z-profile of heat conductivity
!@var kq z-profile of specific humidity
!@var z vertical grids (main, meter)
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var ch  dimensionless heat flux at surface (stanton number)
!@var usurf effective surface velocity
!@var tgrnd virtual potential temperature at the ground
!@var ttop virtual potential temperature at the first GCM layer
!@var n number of vertical main subgrid layers
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: q,dz
      real*8, dimension(n), intent(inout) :: t
      real*8, dimension(n-1), intent(in) :: kh,kq,dzh
      real*8, intent(in) :: ch,tgrnd,ttop,usurf

      real*8 :: facth
      integer :: i !@var i loop variable

      do i=2,n-1
         sub(i)=-1./(dz(i)*dzh(i-1))*kh(i-1)
         sup(i)=-1./(dz(i)*dzh(i))*kh(i)
         dia(i)=-(sub(i)+sup(i))
      end do

      do i=2,n-1
        rhs(i)=0.
      end do

      facth  = ch*usurf*dzh(1)/kh(1)

      dia(1) = 1.+facth
     &        +xdelt*kq(1)/kh(1)*(q(2)-q(1))/(1.+xdelt*q(1))
      sup(1) = -1.
      rhs(1) = facth*tgrnd

      dia(n)  = 1.
      sub(n)  = 0.
      rhs(n) = ttop

      call TRIDIAG(sub,dia,sup,rhs,t,n)

      return
      end subroutine t_eqn_sta

      subroutine q_eqn_sta(q,kq,dz,dzh,cq,usurf,qgrnd,qtop,n)
!@sum q_eqn_sta computes the static solutions of
!@+ the specific humidity, Q,
!@+ between the surface and the first GCM layer.
!@+ The boundary condition at the bottom is:
!@+ kq * dq/dz = cq * usurf * (q - qg).
!@+ At the top, Q is prescribed.
!@+ It is called only at the initialization
!@+ (from within subroutine inits).
!@auth Ye Cheng/G. Hartke
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var q z-profle of specific humidity
!@var kq z-profile of moisture diffusivity
!@var z vertical grids (main, meter)
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var cq  dimensionless moisture flux at surface (dalton number)
!@var usurf effective surface velocity
!@var qgrnd specific humidity at the ground
!@var qtop specific humidity at the first GCM layer
!@var n number of vertical main subgrid layers
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: dz
      real*8, dimension(n), intent(inout) :: q
      real*8, dimension(n-1), intent(in) :: kq,dzh
      real*8, intent(in) :: cq,qgrnd,qtop,usurf

      real*8 :: factq
      integer :: i  !@var i loop variable

      do i=2,n-1
         sub(i)=-1./(dz(i)*dzh(i-1))*kq(i-1)
         sup(i)=-1./(dz(i)*dzh(i))*kq(i)
         dia(i)=-(sub(i)+sup(i))
      end do

      do i=2,n-1
        rhs(i)=0.
      end do

      factq  = cq*usurf*dzh(1)/kq(1)

      dia(1) = 1.+factq
      sup(1) = -1.
      rhs(1)= factq*qgrnd

      dia(n)  = 1.
      sub(n)  = 0.
      rhs(n) = qtop

      call TRIDIAG(sub,dia,sup,rhs,q,n)

      return
      end subroutine q_eqn_sta

      subroutine uv_eqn_sta(u,v,z,km,dz,dzh,
     2            ustar,cm,utop,vtop,coriol,uocean,vocean,n
     &            ,dpdxr,dpdyr,dpdxr0,dpdyr0,ug,vg)
!@sum uv_eqn_sta computes the static solutions of the
!@+  wind components, u and v,
!@+  between the surface and the first GCM layer.
!@+  The boundary conditions at the bottom are:
!@+  km * du/dz = cm * usurf * u,
!@+  km * dv/dz = cm * usurf * v.
!@+  At the top, u and v are prescribed.
!@+  It is called only at the initialization
!@+  (from within subroutine inits).
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var z vertical height at the main grids (meter)
!@var zhat vertical height at the secondary grids (meter)
!@var km z-profile of turbulent viscosity
!@var kh z-profile of turbulent conductivity
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var ustar friction velocity at the surface
!@var cm  dimensionless  momentum flux at surface (drag coeff.)
!@var ch  dimensionless heat flux at surface (stanton number)
!@var cq  dimensionless moisture flux at surface (dalton number)
!@var utop due east component of the wind at the first GCM layer
!@var vtop due north component of the wind at the first GCM layer
!@var coriol the Coriolis parameter
!@var n number of vertical subgrid main layers
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs,rhs1

      real*8, dimension(n), intent(in) :: z,dz
      real*8, dimension(n), intent(inout) :: u,v
      real*8, dimension(n-1), intent(in) :: km,dzh
      real*8, intent(in) :: ustar,cm,utop,vtop,coriol,uocean,vocean
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0

      real*8 :: factx,facty,dpdx,dpdy,usurf,factor
      integer :: i,j,iter  !@var i,j,iter loop variable
c**** passed for SCM
      real*8 ug,vg

      do i=2,n-1
          sub(i)=-1./(dz(i)*dzh(i-1))*km(i-1)
          sup(i)=-1./(dz(i)*dzh(i))*km(i)
          dia(i)=-(sub(i)+sup(i))
       end do

      factx=(dpdxr-dpdxr0)/(z(n)-z(1))
      facty=(dpdyr-dpdyr0)/(z(n)-z(1))
c      write(99,*) factx,facty
      do i=2,n-1
ccc if running SCM then use ug and vg instead of dpdx,dpdy
#ifdef SCM
            rhs(i)=coriol*(v(i)-vg)
            rhs1(i)=-coriol*(u(i)-ug)
#else
            dpdx=factx*(z(i)-z(1))+dpdxr0
            dpdy=facty*(z(i)-z(1))+dpdyr0
            rhs(i)=(coriol*v(i)-dpdx)
            rhs1(i)=-(coriol*u(i)+dpdy)
#endif
      end do

      usurf  = sqrt((u(1)-uocean)**2+(v(1)-vocean)**2)
      factor = cm*usurf*dzh(1)/km(1)

      dia(1) = 1.+factor
      sup(1) = -1.
      rhs(1) = factor*uocean
      rhs1(1)= factor*vocean

      dia(n) = 1.
      sub(n) = 0.
      rhs(n)  = utop
      rhs1(n)  = vtop

      call TRIDIAG(sub,dia,sup,rhs,u,n)
      call TRIDIAG(sub,dia,sup,rhs1,v,n)

      return
      end subroutine uv_eqn_sta

      subroutine level2(e,u,v,t,lscale,dzh,n)
!@sum level2 computes the turbulent kinetic energy e (Cheng et al., 2002).
!@auth  Ye Cheng/G. Hartke
!@var e z-profle of turbulent kinetic energy
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var lscale z-profile of the turbulent dissipation length scale
!@var z vertical grids (main, meter)
!@var dzh(j)  z(j+1)-z(j)
!@var n number of vertical subgrid main layers

      USE CONSTANT, only : teeny

      implicit none

      integer, intent(in) :: n     !@var n array dimension
      real*8, dimension(n), intent(in)   :: u,v,t
      real*8, dimension(n-1), intent(in) :: lscale,dzh
      real*8, dimension(n-1), intent(out) :: e
      real*8 :: dudz,dvdz,as2,an2,ri,gm
      real*8 :: aa,bb,cc,tmp
      integer :: i    !@var i loop variable

      do i=1,n-1
        an2=2.*grav*(t(i+1)-t(i))/((t(i+1)+t(i))*dzh(i))
        dudz=(u(i+1)-u(i))/dzh(i)
        dvdz=(v(i+1)-v(i))/dzh(i)
        as2=max(dudz*dudz+dvdz*dvdz,teeny)
        ri=an2/as2
        if(ri.gt.rimax) ri=rimax
        aa=c1*ri*ri-c2*ri+c3
        bb=c4*ri+c5
        cc=2.d0
        if(abs(aa).lt.1d-8) then
          gm= -cc/bb
        else
          tmp=bb*bb-4.*aa*cc
          gm=(-bb-sqrt(tmp))/(2.*aa)
        endif
        e(i)=0.5d0*(B1*lscale(i))**2*as2/max(gm,teeny)
        e(i)=min(max(e(i),teeny),emax)
      end do

      return
      end subroutine level2

      subroutine inits(tgrnd,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2          ttop,qtop,coriol,cm,ch,cq,ustar,lmonin,uocean,vocean,
     3          ilong,jlat,itype
     &          ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &          ,u,v,t,q,e,ug,vg)
!@sum inits initializes the winds, virtual potential temperature,
!@+  and humidity by solving their differential equations for the
!@+  static solutions, using tridiagonal method over npbl(=8)
!@+  sublayers between the surface (sublayer 1)
!@+  and the middle of the first GCM layer (sublayer npbl).
!@+  (Cheng et a., 2002).
!@+  It is called by subroutine init_pbl (in PBL_DRV.f),
!@+  and the latter (init_pbl) is called 
!@+  by subroutine INPUT (in MODELE.f).
!@auth  Ye Cheng/G. Hartke
!@var  n number of sub-grid levels for the PBL
!@var  tgrnd virtual potential temperature of ground,at roughness height
!@var  qgrnd  moisture at the ground, at the roughness height
!@var  zgrnd
!@var  zgs  height of the surface layer (nominally 10 m)
!@var  ztop height of the first model layer, approx 200 m if lm=9
!@var  utop  x component of wind at the top of the layer
!@var  vtop  y component of wind at the top of the layer
!@var  ttop  virtual potential temperature at the top of the layer
!@var  qtop  moisture at the top of the layer
!@var  coriol the Coriolis parameter
!@var  cm  dimensionless  momentum flux at surface (drag coeff.)
!@var  ch  dimensionless heat flux at surface (stanton number)
!@var  cq  dimensionless moisture flux at surface (dalton number)
!@var  bgrid log-linear gridding parameter
!@var  ustar  friction speed
!@var  ilong  longitude identifier
!@var  jlat  latitude identifier
!@var  itype  surface type
!@var  ug,vg  passed for SCM

!@var  iprint longitude for diagnostics
!@var  jprint latitude for diagnostics
      implicit none

      real*8, intent(in) :: tgrnd,qgrnd,zgrnd,zgs,ztop,utop,vtop,ttop
     *     ,qtop,coriol,uocean,vocean
      real*8, intent(out) :: cm,ch,cq,ustar,lmonin
      integer, intent(in) :: ilong,jlat,itype
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0

      real*8, dimension(n-1) :: km,kh,kq,ke,gm,gh
      real*8, dimension(n) :: z,dz,xi,usave,vsave,tsave,qsave
      real*8, dimension(n-1) :: zhat,xihat,dzh,lscale,esave
      real*8 :: lmonin_dry,bgrid,z0m,z0h,z0q,hemi,psi1,psi0,psi
     *     ,usurf,tstar,qstar,ustar0,dtime,test
     *     ,wstar3,wstar2h,usurfq,usurfh,ts

      integer, parameter ::  itmax=100
      integer, parameter ::  iprint=0,jprint=41 ! set iprint>0 to debug
      real*8, parameter ::  w=0.50,tol=1d-3
      integer :: i,iter,ierr  !@var i,iter loop variable
#ifdef TRACERS_SPECIAL_O18
      real*8 :: fac_cq_tr(ntm)   ! not used here
#endif

      real*8 dbl ! I hope it is really a local variable (was global before) I.A

      real*8, dimension(n), intent(out) :: u,v,t,q
      real*8, dimension(n-1), intent(out) :: e
c****  passed for SCM
      real*8  ug,vg
      real*8  dm

      dbl=1000.d0 !initial guess of dbl
      z0m=zgrnd
      z0h=z0m
      z0q=z0m

      call griddr(z,zhat,xi,xihat,dz,dzh,zgs,ztop,bgrid,n,ierr)
      if (ierr.gt.0) then
        print*,"In inits: i,j,itype =",ilong,jlat,itype,tgrnd,qgrnd
        call stop_model("PBL error in inits",255)
      end if

c Initialization for iteration:
      if (coriol.le.0.) then
        hemi=-1.
        else
        hemi= 1.
      endif
      if (tgrnd.gt.ttop) then
        psi1=hemi*5.*radian
        else
        psi1=hemi*15.*radian
      endif
      psi0=atan2(vtop,utop+teeny)
      if(psi0.lt.0.) psi0=psi0+2.*pi
      psi=psi0+psi1
      usurf=0.4d0*sqrt(utop*utop+vtop*vtop)
      u(1)=usurf*cos(psi)
      v(1)=usurf*sin(psi)
      t(1)=tgrnd-0.5d0*(tgrnd-ttop)
      q(1)=qgrnd-0.5d0*(qgrnd-qtop)
      e(1)=1.d-3
        usave(1)=u(1)
        vsave(1)=v(1)
        tsave(1)=t(1)
        qsave(1)=q(1)
        esave(1)=e(1)
      u(n)=utop
      v(n)=vtop
      t(n)=ttop
      q(n)=qtop
      e(n-1)=1.d-3

      do i=2,n-1
        u(i)=u(1)+(u(n)  -u(1))*((z(i)-z(1))/(z(n)  -z(1)))
        v(i)=v(1)+(v(n)  -v(1))*((z(i)-z(1))/(z(n)  -z(1)))
        t(i)=t(1)+(t(n)  -t(1))*((z(i)-z(1))/(z(n)  -z(1)))
        q(i)=q(1)+(q(n)  -q(1))*((z(i)-z(1))/(z(n)  -z(1)))
        e(i)=e(1)+(e(n-1)-e(1))*((z(i)-z(1))/(z(n-1)-z(1)))
        usave(i)=u(i)
        vsave(i)=v(i)
        tsave(i)=t(i)
        qsave(i)=q(i)
        esave(i)=e(i)
      end do

      ustar0=0.
      do iter=1,itmax

        call getl1(e,zhat,dzh,lscale,n)
        call getk(km,kh,kq,ke,gm,gh,u,v,t,e,lscale,dzh,n)
        ts=t(1)/(1+q(1)*xdelt)
        call stars(ustar,tstar,qstar,lmonin,lmonin_dry,tgrnd,qgrnd,ts,
     2             u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,dm,
#ifdef TRACERS_SPECIAL_O18
     *             fac_cq_tr,
#endif
     3             km,kh,kq,dzh,itype,n)
        lmonin = lmonin_dry ! inits always sees dry temps
        !! dbl=.375d0*sqrt(ustar*abs(lmonin)/omega)
        !@+ M.J.Miller et al. 1992, J. Climate, 5(5), 418-434, Eqs(6-7),
        !@+ for heat and mositure

        if(t(2).lt.t(1)) then
          wstar3=-dbl*grav*2.*(t(2)-t(1))*kh(1)/((t(2)+t(1))*dzh(1))
          wstar2h = wstar3**twoby3
        else
          wstar3=0.
          wstar2h = 0.
        endif

        usurfh  = sqrt((u(1)-uocean)**2+(v(1)-vocean)**2+wstar2h)
        usurfq  = usurfh

        call t_eqn_sta(t,q,kh,kq,dz,dzh,ch,usurfh,tgrnd,ttop,n)

        call q_eqn_sta(q,kq,dz,dzh,cq,usurfq,qgrnd,qtop,n)

        call uv_eqn_sta(u,v,z,km,dz,dzh,ustar,cm,utop,vtop,coriol,
     &                  uocean,vocean,n
     &                  ,dpdxr,dpdyr,dpdxr0,dpdyr0,ug,vg)

        call tcheck(t,tgrnd,n)
        call tcheck(q,qgrnd,n)
        call ucheck(u,v,z,ustar,lmonin,z0m,hemi,psi0,psi1,n)

        test=abs(2.*(ustar-ustar0)/(ustar+ustar0))
        if (test.lt.tol) exit
        call level2(e,u,v,t,lscale,dzh,n)

        do i=1,n-1
          u(i)=w*usave(i)+(1.-w)*u(i)
          v(i)=w*vsave(i)+(1.-w)*v(i)
          t(i)=w*tsave(i)+(1.-w)*t(i)
          q(i)=w*qsave(i)+(1.-w)*q(i)
          e(i)=w*esave(i)+(1.-w)*e(i)
          usave(i)=u(i)
          vsave(i)=v(i)
          tsave(i)=t(i)
          qsave(i)=q(i)
          esave(i)=e(i)
        end do
        ustar0=ustar
      end do

c     iter_count=iter_count+min(iter,itmax)
c     write(96,*) "iter_count in inits =",iter_count, iter,dbl

c     call check1(ustar,1,ilong,jlat,1)

      if ((ilong.eq.iprint).and.(jlat.eq.jprint)) then
        call output(u,v,t,q,e,lscale,z,zhat,dzh,
     2              km,kh,kq,ke,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3              ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4              utop,vtop,ttop,qtop,
     5              dtime,bgrid,ilong,jlat,iter,itype,n)
      endif


      return
      end subroutine inits

      subroutine tcheck(t,tgrnd,n)
!@sum tcheck checks for reasonable temperatures
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
c ----------------------------------------------------------------------
c This routine makes sure that the temperature remains within
c  reasonable bounds during the initialization process. (Sometimes the
c  the computed temperature iterated out in left field someplace,
c  *way* outside any reasonable range.) This routine keeps the temp
c  between the maximum and minimum of the boundary temperatures.
c ----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: n    !@var n array dimension
      real*8, dimension(n),intent(inout) ::  t !@var t temperature array
      real*8, intent(in) :: tgrnd  !@var tgrnd ground temperature

      integer, dimension(20) :: imax,imin
      integer nmin,nmax
      real*8 tmin,tmax,dt
      integer i                 !@var i  loop and dummy variables

      if (tgrnd.lt.t(n)) then
        tmin=tgrnd
        tmax=t(n)
      else
        tmin=t(n)
        tmax=tgrnd
      endif
      nmin=0
      nmax=0

      do i=1,n-1
        if (t(i).lt.tmin) then
          nmin=nmin+1
          imin(nmin)=i
        endif
        if (t(i).gt.tmax) then
          nmax=nmax+1
          imax(nmax)=i
        endif
      end do

      dt=0.01*(tmax-tmin)
      if (nmin.gt.0) then
        do i=1,nmin
          t(imin(i))=tmin+float(i)*dt
       end do
      endif

      if (nmax.gt.0) then
        do i=1,nmax
          t(imax(i))=tmax-float(i)*dt
       end do
      endif

      return
      end subroutine tcheck

      subroutine ucheck(u,v,z,ustar,lmonin,z0m,hemi,psi0,psi1,n)
!@sum ucheck makes sure that the winds remain within reasonable
!@+    bounds during the initialization process. (Sometimes the computed
!@+    wind speed iterated out in left field someplace, *way* outside
!@+    any reasonable range.) Tests and corrects both direction and
!@+    magnitude of the wind rotation with altitude. Tests the total
!@+    wind speed via comparison to similarity theory. Note that it
!@+    works from the top down so that it can assume that at level (i),
!@+    level (i+1) displays reasonable behavior.
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
      implicit none
      real*8, parameter :: psistb=15.*radian, psiuns=5.*radian

      integer, intent(in) :: n  !@var n array dimension
      real*8, dimension(n),intent(in) :: z
      real*8, dimension(n),intent(inout) :: u,v

      real*8, intent(in) :: lmonin,z0m,hemi,psi0,psi1,ustar

      real*8 psilim,psirot,angle,utotal,x,x0,dpsim,psiu,zet,zet0,utest
      real*8 xm
      integer i                 !@var i  loop and dummy variables

      if (lmonin.ge.0.) then
        psilim=psistb
        else
        psilim=psiuns
      endif

c First, check the rotation of the wind vector:

      do i=n-1,1,-1
        psiu=atan2(v(i),u(i)+teeny)
        if(psiu.lt.0.) psiu=psiu+2.*pi
        psirot=psiu-psi0
c --------------------------------------------------------------------
c Next, check the magnitude of the wind. If the gradient is incorrect,
c  set the wind magnitude to that given by similarity theory:

        utotal=sqrt(u(i)*u(i)+v(i)*v(i))
        utest =sqrt(u(i+1)*u(i+1)+v(i+1)*v(i+1))
        if (utotal.gt.utest) then
          zet=z(i)/lmonin
          zet0=z0m/lmonin
          call find_dpsim(zet,zet0,dpsim)
          utotal=(ustar/kappa)*(log(z(i)/z0m)-dpsim)
          if (utotal.gt.utest) utotal=0.95d0*utest
          u(i)=utotal*cos(psiu)
          v(i)=utotal*sin(psiu)
        endif
        if (hemi*psirot.lt.0.) then
          angle=psi0+psi1*float(n-i)/float(n-1)
          u(i)=utotal*cos(angle)
          v(i)=utotal*sin(angle)
        elseif (hemi*psirot.gt.psilim) then
          angle=psi0+hemi*psilim*float(n-i)/float(n-1)
          u(i)=utotal*cos(angle)
          v(i)=utotal*sin(angle)
        endif

      end do

      return
      end subroutine ucheck


      subroutine check1(a,n,ilong,jlat,id)
!@sum check1 checks for NaN'S and INF'S in real 1-D arrays.
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
      implicit none

      integer, intent(in) :: n      !@var n  array dimension
      integer, intent(in) :: ilong  !@var ilong  longitude identifier
      integer, intent(in) :: jlat   !@var jlat  latitude identifier
      integer, intent(in) :: id     !@var n  integer id
      real*8, dimension(n), intent(in) ::  a !@var a  real array

      integer i,k !@var i,k loop and dummy variables
      character*16 str  !@var str  output string

      do i=1,n
        write(str,'(e16.8)')a(i)
        k=index(str,'N')+index(str,'n')
        if (k.ne.0) then
          write (99,1000) ilong,jlat,i,id,a(i)
          if (id.lt.100) call stop_model('check1',255)
        endif
      end do

      return
 1000 format (1x,'check1:  ilong = ',6x,i3,/,1x,
     2            '         jlat  = ',6x,i3,/,1x,
     3            '         i     = ',6x,i3,/,1x,
     4            '         id    = ',6x,i3,/,1x,
     5            '         value = ',1pe11.4,/)
      end subroutine check1

      subroutine output(u,v,t,q,e,lscale,z,zhat,dzh,
     2                  km,kh,kq,ke,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3                  ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4                  utop,vtop,ttop,qtop,
     5                  dtime,bgrid,ilong,jlat,iter,itype,n)
!@sum output produces output for diagnostic purposes.
!@auth  Ye Cheng/G. Hartke
!@calls simil
      implicit none
      real*8, parameter :: degree=1./radian

      integer, intent(in) :: n,itype,iter,jlat,ilong
      real*8,  intent(in) :: lscale(n-1),lmonin

      real*8, dimension(n), intent(in) :: u,v,t,q,z
      real*8, dimension(n-1), intent(in) :: e,zhat,dzh,km,kh,kq,ke,gm,gh

      real*8, intent(in) :: cm,ch,cq,z0m,z0h,z0q,ustar,tstar,qstar,tgrnd
     *     ,qgrnd,utop,vtop,ttop,qtop,dtime,bgrid

      integer i  !@var i loop variable

      real*8 :: psitop,psi,utotal,utot1,utot2,sign,shear,tgrad,qgrad
     *     ,uflux,hflux,qflux,bvfrq2,shear2,rich,dqdz,dtdz,dudz
     *     ,phim,phih,dudzs,dtdzs,dqdzs,uratio,tratio,qratio,dbydzh
     *     ,tgradl,prod,utest,ttest,qtest,dpsim,dpsih,dpsiq,zet

      write (99,5000) ilong,jlat,itype,dtime,iter,ustar,tstar,qstar,
     2                lmonin,tgrnd,qgrnd,cm,ch,cq,z0m,z0h,z0q,
     3                utop,vtop,ttop,qtop,
     4                bgrid
      write (99,1000)
      psitop=atan2(v(n),u(n)+teeny)*degree
      if (psitop.lt.0.) psitop=psitop+360.
      do i=1,n
        psi=atan2(v(i),u(i)+teeny)*degree
        if (psi.lt.0.) psi=psi+360.
        psi=psi-psitop
        utotal=sqrt(u(i)*u(i)+v(i)*v(i))
        call simil(z(i),z0m,z0h,z0q,lmonin,ustar,tstar,qstar
     2       ,tgrnd,qgrnd,utest,ttest,qtest,dpsim,dpsih,dpsiq)
        write (99,3000) i,z(i),u(i),v(i),psi,utotal,utest,
     2                    t(i),ttest,q(i),qtest
      end do
      write (99,9000)
      write (99,2000)
      do i=1,n-1
        utot1=u(i+1)*u(i+1)+v(i+1)*v(i+1)
        utot2=u(i)  *u(i)  +v(i)  *v(i)
        if (utot1.ge.utot2) then
          sign=1.
          else
          sign=-1.
        endif
        shear =sqrt((u(i+1)-u(i))*(u(i+1)-u(i))
     2             +(v(i+1)-v(i))*(v(i+1)-v(i)))/dzh(i)
        tgrad =(t(i+1)-t(i))/dzh(i)
        qgrad =(q(i+1)-q(i))/dzh(i)
        uflux=km(i)*shear*sign
        hflux=kh(i)*tgrad
        qflux=kh(i)*qgrad
        bvfrq2=grav*log(t(i+1)/t(i))/dzh(i)+1d-12
        shear2=((u(i+1)-u(i))/dzh(i))**2+
     2         ((v(i+1)-v(i))/dzh(i))**2+1d-8
        rich=bvfrq2/shear2
        write (99,3000) i,zhat(i),e(i),lscale(i),km(i),kh(i),
     2                    gm(i),gh(i),uflux,hflux,qflux,rich
      end do
      write (99,9000)
      write (99,7000)

      do i=1,n-1
        utot1=sqrt(u(i+1)*u(i+1)+v(i+1)*v(i+1))
        utot2=sqrt(  u(i)*  u(i)+  v(i)*  v(i))
        dudz=(utot1-utot2)/dzh(i)
        dtdz=(t(i+1)-t(i))/dzh(i)
        dqdz=(q(i+1)-q(i))/dzh(i)
        zet=zhat(i)/lmonin
        call find_phim(zet,phim)
        call find_phih(zet,phih)
        dudzs=ustar*phim/(kappa*zhat(i))
        dtdzs=tstar*phih/(kappa*zhat(i))
        dqdzs=qstar*phih/(kappa*zhat(i))
        uratio=dudzs/(dudz+teeny)
        tratio=dtdzs/(dtdz+teeny)
        qratio=dqdzs/(dqdz+teeny)

        dbydzh=1./dzh(i)
        shear2=dbydzh*dbydzh*((u(i+1)-u(i))**2+
     2                        (v(i+1)-v(i))**2)
        tgradl=dbydzh*log(t(i+1)/t(i))
        prod=km(i)*shear2-grav*kh(i)*tgradl

        write (99,3000) i,zhat(i),dudz,dudzs,dtdz,dtdzs,dqdz,dqdzs,
     2                            uratio,tratio,qratio,prod
      end do
      write (99,9000)
      write (99,9000)
c ----------------------------------------------------------------------
      return
1000  format (1x,'   i','      z     ','      U     ','      V     ',
     2                   '     psi    ','    Utotal  ','    Utest   ',
     2                   '      T     ','    Ttest   ','      Q     ',
     2                   '    Qtest   ',/)
2000  format (1x,'   i','     zhat   ','      E     ','      L     ',
     2                   '      Km    ','      Kh    ',
     3                   '      Gm    ','      Gh    ','    uflux   ',
     4                   '    hflux   ','    qflux   ','    rich #  ',/)
3000  format (1x,i4,11(1x,1pe11.4))
5000  format (1x,'i      = ',9x,i2,/,1x,
     2            'j      = ',9x,i2,/,1x,
     2            'itype  = ',9x,i2,/,1x,
     2            'dtime  = ',1pe11.4,/,1x,
     3            'iter   = ',7x,i4,/,1x,
     4            'ustar  = ',1pe11.4,/,1x,
     5            'tstar  = ',1pe11.4,/,1x,
     5            'qstar  = ',1pe11.4,/,1x,
     5            'lmonin = ',1pe11.4,/,1x,
     5            'tgrnd  = ',1pe11.4,/,1x,
     5            'qgrnd  = ',1pe11.4,/,1x,
     6            'cm     = ',1pe11.4,/,1x,
     6            'ch     = ',1pe11.4,/,1x,
     6            'cq     = ',1pe11.4,/,1x,
     6            'z0m    = ',1pe11.4,/,1x,
     7            'z0h    = ',1pe11.4,/,1x,
     8            'z0q    = ',1pe11.4,/,1x,
     6            'utop   = ',1pe11.4,/,1x,
     6            'vtop   = ',1pe11.4,/,1x,
     7            'ttop   = ',1pe11.4,/,1x,
     8            'qtop   = ',1pe11.4,/,1x,
     8            'bgrid  = ',1pe11.4,/)
7000  format (1x,'   i','     zhat   ','     dudz   ','   dudz sim ',
     2                   '     dtdz   ','   dtdz sim ','     dqdz   ',
     3                   '   dqdz sim ','    uratio  ','    tratio  ',
     4                   '    qratio  ','  production',/)
9000  format (1x)
      end subroutine output

      subroutine find_phim0(zet,phim)
      implicit none
      ! in:
      real*8 zet
      ! out:
      real*8 phim
      if(zet.ge.0.d0) then ! stable or neutral
        if(zet.le.zet1) then
          phim=1.+gamams*zet
        else
          phim=1.+gamams*zet1+slope1*(zet-zet1)
        endif
      else                ! unstable
c       if(zet.ge.zetm) then
          phim=(1.-gamamu*zet)**(-.25d0)
c       else
c         phim=0.7d0*kappa**(2*by3)*(-zet)**by3
c       endif
      endif
      return
      end subroutine find_phim0

      subroutine find_phim(zet,phim)
      implicit none
      ! in:
      real*8 zet
      ! out:
      real*8 phim
      if(zet.ge.0.d0) then ! stable or neutral
        if(zet.le.zet1) then
          phim=1.+gamams*zet
        else
          phim=1.+gamams*zet1+slope1*(zet-zet1)
        endif
      else                ! unstable
        if(zet.ge.zetm) then
          phim=(1.-gamamu*zet)**(-.25d0)
        else
          phim=0.7d0*kappa**(2*by3)*(-zet)**by3
        endif
      endif
      return
      end subroutine find_phim

      subroutine find_phih(zet,phih)
      implicit none
      ! in:
      real*8 zet
      ! out:
      real*8 phih
      if(zet.ge.0.d0) then ! stable or neutral
        if(zet.le.zet1) then
          phih=sigma*(1.+gamahs*zet)
        else
          phih=sigma*(1.+gamahs*zet1+slope1*(zet-zet1))
        endif
      else                ! unstable
        if(zet.ge.zeth) then
          phih=sigma*(1-gamahu*zet)**(-.5d0)
        else
          phih=.9d0*kappa**(4.d0/3.d0)*(-zet)**(-by3)
        endif
      endif
      return
      end subroutine find_phih

      subroutine find_dpsim(zet,zet0,dpsim)
      implicit none
      ! in:
      real*8 zet,zet0
      ! out:
      real*8 dpsim
      real*8 x,x0,xm
      ! dpsim:
      if(zet.ge.0.d0) then ! stable
        if(zet.le.zet1) then
          dpsim=-gamams*(zet-zet0)
        else
          dpsim=-gamams*(zet1-zet0)
     &          +zet1*(slope1-gamams)*log(zet/zet1)
     &          -slope1*(zet-zet1)
        endif
      else                 ! unstable
        x= (1.-gamamu*zet)**0.25d0
        x0=(1.-gamamu*zet0)**0.25d0
        xm= (1.-gamamu*zetm)**0.25d0
        if(zet.gt.zetm) then
          dpsim=log((1+x)*(1+x)*(1+x*x)/
     1          ((1+x0)*(1+x0)*(1+x0*x0)))-
     2          2.*(atan(x)-atan(x0))
        else
          dpsim=log((1+xm)*(1+xm)*(1+xm*xm)/
     1          ((1+x0)*(1+x0)*(1+x0*x0)))-
     2          2.*(atan(xm)-atan(x0))
     3          +log(zet/zetm)
     4          -1.140125d0*((-zet)**by3-(-zetm)**by3)
                ! 0.7*kappa**(2/3)*3 = 1.140125
        endif
      endif
      return
      end subroutine find_dpsim

      subroutine find_dpsih(zet,zet0,z,z0,dpsih)
      implicit none
      ! in:
      real*8 zet,zet0,z,z0
      ! out:
      real*8 dpsih
      real*8 x,x0,xh,rat
      ! dpsih
      if(zet.ge.0.d0) then ! stable
        if(zet.le.zet1) then
          dpsih=sigma1*log(z/z0)-sigma*gamahs*(zet-zet0)
        else
          dpsih=sigma1*log(zet1/zet0)-sigma*gamahs*(zet1-zet0)
     &         +(1+sigma*(zet1*(slope1-gamahs)-1))*log(zet/zet1)
     &         -sigma*slope1*(zet-zet1)
        endif
      else                 ! unstable
        x= (1.-gamahu*zet)**0.5d0
        x0=(1.-gamahu*zet0)**0.5d0
        xh= (1.-gamahu*zeth)**0.5d0
        if(zet.gt.zeth) then
          if(-gamahu*zet.lt.1d-5) then
            rat=z0/z*(1-.5d0*gamahu*(zet-zet0))
          else
            rat=(1+x)*(1-x0)/((1-x)*(1+x0))
          endif
          dpsih=log(z/z0)+sigma*log(rat)
        else
          dpsih=log(z/z0)
     2          +sigma*log((1+xh)*(1-x0)/((1-xh)*(1+x0)))
     3          -0.7957508d0*((-zeth)**(-by3)-(-zet)**(-by3))
                ! 0.9*kappa**(4/3)*3 = 0.7957508
        endif
      endif
      return
      end subroutine find_dpsih

      subroutine alloc_pbl_args(pbl_args)
      USE CONSTANT, only : IUNDEF_VAL, UNDEF_VAL
#ifdef TRACERS_ON
      use tracer_com, only: gasex_index
#endif
      implicit none
      type (t_pbl_args), intent(inout) :: pbl_args

#ifdef TRACERS_ON
      allocate(pbl_args%trtop(maxNTM))
      allocate(pbl_args%trs(maxNTM))
      allocate(pbl_args%trsfac(maxNTM))
      allocate(pbl_args%trconstflx(maxNTM))
      allocate(pbl_args%trdn1(maxNTM))
      allocate(pbl_args%trprime(maxNTM))
      allocate(pbl_args%trgrnd2(maxNTM))
      allocate(pbl_args%ntix(maxNTM))
      pbl_args%ntix = IUNDEF_VAL
      pbl_args%trtop = UNDEF_VAL
      pbl_args%trs = UNDEF_VAL
      pbl_args%trsfac = UNDEF_VAL
      pbl_args%trconstflx = UNDEF_VAL
      pbl_args%trdn1 = UNDEF_VAL
      pbl_args%trprime = UNDEF_VAL
      pbl_args%trgrnd2 = UNDEF_VAL
#ifdef TRACERS_SPECIAL_O18
      allocate(pbl_args%frack(maxNTM))
      pbl_args%frack = UNDEF_VAL
#endif
      
#ifdef TRACERS_DRYDEP
      allocate(pbl_args%dep_vel(maxNTM))
      allocate(pbl_args%gs_vel(maxNTM))
      pbl_args%dep_vel = 0.d0
      pbl_args%gs_vel = 0.d0 
#endif
#ifdef TRACERS_WATER
      allocate(pbl_args%tr_evap_max(maxNTM))
      pbl_args%tr_evap_max = UNDEF_VAL
#endif
      allocate(pbl_args%kw_gas(gasex_index%getsize()))
      allocate(pbl_args%alpha_gas(gasex_index%getsize()))
      allocate(pbl_args%beta_gas(gasex_index%getsize()))
#endif
      end subroutine alloc_pbl_args

      subroutine dealloc_pbl_args(pbl_args)
      type (t_pbl_args), intent(inout) :: pbl_args

#ifdef TRACERS_ON
      deallocate(pbl_args%trtop)
      deallocate(pbl_args%trs)
      deallocate(pbl_args%trsfac)
      deallocate(pbl_args%trconstflx)
      deallocate(pbl_args%trdn1)
      deallocate(pbl_args%trprime)
      deallocate(pbl_args%trgrnd2)
      deallocate(pbl_args%ntix)
#ifdef TRACERS_SPECIAL_O18
      deallocate(pbl_args%frack)
#endif
      
#ifdef TRACERS_DRYDEP
      deallocate(pbl_args%dep_vel)
      deallocate(pbl_args%gs_vel)
#endif
#ifdef TRACERS_WATER
      deallocate(pbl_args%tr_evap_max)
#endif
#endif
      end subroutine dealloc_pbl_args
      
      END MODULE SOCPBL

      subroutine fgrid2(z,f,df)
!@sum fgrid2 computes functional relationship of z and xi;
!@+  it is used in function NewtonMethod,
!@+  the latter is called by subroutine griddr.
!@auth Ye Cheng/G. Hartke
      implicit none
      real*8, intent(in) :: z
      real*8, intent(out) :: f,df
      real*8 z1,zn,bgrid,xi,lznbyz1
      common /grids_99/z1,zn,bgrid,xi,lznbyz1

      f=z+bgrid*((zn-z1)*log(z/z1)-(z-z1)*lznbyz1)-xi
      df=1.+bgrid*((zn-z1)/z-lznbyz1)

      return
      end subroutine fgrid2

      real*8 function deltaSST(Qnet,Qsol,ustar_oc)
!@sum deltaSST calculate skin-bulk SST difference (deg C)
!@var Qnet Net heat flux (not including solar, +ve dwn) (W/m2)
!@var Qsol Solar heat flux (W/m2)
!@var ustar_oc ocean u* (m/s)
      USE CONSTANT, only : rhows,visc_wtr_kin
      IMPLICIT NONE
      real*8, intent(in) :: Qnet,Qsol,ustar_oc
C**** k thermal cond. 0.596 W/mK (35 psu, 20 deg, 0 press)
      real*8, parameter :: byk= 1.677d0  ! 1/thermal conductivity (K/(W/m))
!@var lam coefficient (non-dim)
      real*8, parameter :: lam=2.4d0  ! Ward (2007)
!@var del micro-layer thickness (m)
      real*8 :: del
!@var fc fraction of solar absorbed in micro-layer
      real*8 :: fc

C**** calculate micro-layer thickness (m)
C**** Cap ustar so that it does not get too small in low wind conditions
C**** ustar_oc > 0.00098 corresponding to tau > 0.001 N/m2
      del = lam*visc_wtr_kin/max(ustar_oc,0.00098d0)

C**** fraction of solar absorbed (Fairall et al, 1996)
      if (del .lt. 1d-8) then
        fc = 0.0545d0 + 11.*del
      else
        fc = 0.137d0 + 11.*del - 6.6d-5*(1.-exp(-del*1250.))/del
      end if
      fc = min(0.24d0,fc)    ! assume max feasible del=1cm?

C**** delta SST = skin - bulk (+ve for flux going down)
C**** includes occurences of warm skin temperatures (generally < 0.01 C)
      deltaSST=(Qnet + fc*Qsol)*del*byk

      end function deltaSST

      real*8 function deltaSnowT(Qnet,Qsol,snow,tg,dQnetdtg,deltatg0)
!@sum deltaSnowT calculate skin-bulk snow T difference (deg C)
      USE SEAICE, only : alams, alami, rhos, ace1i, byrhoi, xsi, rhoi
      IMPLICIT NONE
!@var Qnet Net heat flux (not including solar, +ve dwn) (W/m2)
!@var Qsol Solar heat flux (W/m2)
!@var snow Snow mass (kg/m2)
!@var dQnetdtg derivative of Qnet w.r.t. tg 
!@var deltatg0 starting value of deltaSnowT
      real*8, intent(in) :: Qnet,Qsol,snow,tg,dQnetdtg,deltatg0
!@var ksext solar radiation extinction coeffficent (1/m)
!@var ssi1 salinity in upper layer ice (estimated)
!@var dz1 depth of first thermal layer (m)
!@var hsnow depth of snow (m)
      real*8 :: ksext,dz1,hsnow,delbyk
      real*8 :: ssi1 = 5d0 ! estimate (should be passed?)

!@var fc fraction of solar passing though in micro-layer (LeComte et al, 2013)
      real*8, parameter :: fc = 0.82d0

      hsnow = snow/rhos
      dz1=hsnow+(xsi(1)*ace1i-snow*xsi(2))/rhoi

C**** Effective dz/K
      if (xsi(2)*snow > xsi(1)*ace1i) then ! only snow in first thermal layer
         delbyk = 0.5*hsnow/alams
      else ! snow and ice
         delbyk = hsnow/alams + (0.5*dz1-hsnow)/alami(tg,ssi1)
      end if

C**** snow extinction coefficients. This really depends on fraction of vis and nir
C**** and on snow coniditons (wet or dry) (see solar_ice_frac), but just assume a standard value here
      ksext = 20d0

C**** delta Snow T = skin - bulk (+ve for flux going down)
c      print*,"in1",hsnow,dz1,Qnet,Qsol,delbyk,
c     *     (Qnet+fc*Qsol)*delbyk
c      print*,"in2",(exp(-ksext*hsnow)-1.0+ksext*hsnow)/
c     *     (ksext*ksext*hsnow*alams),(1-fc)*Qsol*
c     *     (exp(-ksext*hsnow)-1.0+ksext*hsnow)/
c     *     (ksext*ksext*hsnow*alams)
c**** explciit is too noisy
c     deltaSnowT= (Qnet + fc*Qsol)*delbyk +
c     *     (1-fc)*Qsol*(exp(-ksext*hsnow)-1.0+ksext*hsnow)/
c     *     (ksext*ksext*hsnow*alams)
c**** implicit 
      deltaSnowT= (Qnet + dQnetdtg*deltatg0 + fc*Qsol
     *    + (1-fc)*Qsol*(exp(-ksext*hsnow)-1.0+ksext*hsnow)/
     *     (ksext*ksext*hsnow*alams))*delbyk
     *     /(1d0 - dQnetdtg*delbyk) 
      
      end function deltaSnowT
