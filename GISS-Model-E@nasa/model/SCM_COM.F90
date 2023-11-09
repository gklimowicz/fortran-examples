! SCM_COM.F90
!@sum Declare SCM global variables and read SCM setup options
!@auth Audrey Wolf (modifications by Fridlind)
!
!--------------------------------------------------------------------------------

  module HorizontalRes
!@sum Trivial horizontal resolution definition for SCM
  implicit none
!@var IM,JM = longitudinal and latitudinal number of grid cells
  integer, parameter :: IM=1,JM=1
  end module HorizontalRes

!--------------------------------------------------------------------------------

  module SCM_COM
  use VerticalRes, only : LM
  implicit none
  save

!@var nstepSCM current SCM time step
  integer :: nstepSCM=0

!@type SCMoptions type for SCM setup options
  type SCMoptions
    logical :: sflx,Tskin,Qskin,Ps,z0m,ustar,alb
    logical :: wind,geo,temp,theta,wvmr,rh
    logical :: ozone,omega,w,VadvHwind
    logical :: ls_v,ls_h,ls_h_UV,Qrad
    logical :: nudge,Fnudge,sfcQrad
    logical :: BeersLaw,PlumeDiag
    logical :: TopHat,allowMC,allowCTEI
    real*8 :: lat,lon,area,tau
    integer :: sfc
  end type SCMoptions
  type(SCMoptions) SCMopt
!@var SCMopt%sflx = T:use prescribed sensible and latent heat fluxes
!@var SCMopt%Tskin = T:use prescribed skin T for radiation and/or ocean sensible heat flux
!@var SCMopt%Qskin = T:use prescribed skin Q (requires ocean surface setting)
!@var SCMopt%Ps = T:use prescribed surface pressure
!@var SCMopt%wind = T:specify winds
!@var SCMopt%geo = T:use geostrophic winds for Coriolis forcing
!@var SCMopt%temp = T:specify absolute temperature
!@var SCMopt%theta = T:specify potential temperature with 1000-mb ref
!@var SCMopt%wvmr = T:specify water vapor mixing ratio
!@var SCMopt%rh = T:specify relative rather than specific humidity
!@var SCMopt%z0m = T:specify surface roughness height
!@var SCMopt%ustar = T:specify surface friction velocity
!@var SCMopt%alb = T:specify surface mid-visible albedo
!@var SCMopt%omega = T:specify omega for qv, theta vertical forcings
!@var SCMopt%ozone = T:specify ozone profile
!@var SCMopt%w = T:specify large-scale vertical wind
!@var SCMopt%VadvHwind = T:vertical forcing of horizontal wind (using Omega or W)
!@var SCMopt%ls_v = T:specify qv and dry static energy / Cp vert adv flux divergence
!@var SCMopt%ls_h = T:specify qv and dry static energy / Cp horiz adv flux divergence
!@var SCMopt%ls_h_UV = T:horizontal forcing of horizontal wind (specified)
!@var SCMopt%Qrad = T:specify fixed radiative heating profile
!@var SCMopt%nudge = T:nudge qv and T with timescale tau
!@var SCMopt%Fnudge = T:apply scale factor profile to qv and T nudging
!@var SCMopt%sfcQrad = T:allow longwave atmospheric heating associated with the surface
!@var SCMopt%lat,SCMopt%lon = SCM latitude and longitude
!@var SCMopt%area = SCM nominal area (m2)
!@var SCMopt%tau = nudging time constant (s) for qv and T
!@var SCMopt%sfc = -1:GCM land only, 0(default):GCM   1:land, 2:ocean 
!@var SCMopt%TopHat = T:assume top-hat distributions for input profiles (else piecewise linear)
!@var SCMopt%allowMC = T:allow moist convection
!@var SCMopt%allowCTEI = T:allow cloud-top entrainment instability
!@var SCMopts%BeersLaw = T:use Beer's Law treatment (only) for radiative heating
!@var SCMopts%PlumeDiag = T:report moist convection plume diagnostics

!@var SCMinputs type for SCM inputs at each time step
  type SCMinputs
    real*8 U(LM),V(LM),Ug(LM),Vg(LM)
    real*8 T(LM),TH(LM),Q(LM),Omega(LM),W(LM),O3(LM)
    real*8 SadvV(LM),QadvV(LM),TadvH(LM),QadvH(LM)
    real*8 UadvH(LM),VadvH(LM)
    real*8 Qrad(LM),Fnudge(LM)
    real*8 time,lhf,shf,Qskin,Tskin,Ps,z0m,ustar,alb
    real*8 BeersLaw_f0,BeersLaw_f1,BeersLaw_kappa
  end type SCMinputs
  type(SCMinputs) SCMin
!@var SCMin%U SCM input zonal wind at GCM sigma levels (m/s)
!@var SCMin%V SCM input meridional wind at GCM sigma levels (m/s)
!@var SCMin%Ug SCM input geostrophic zonal wind at GCM sigma levels (m/s)
!@var SCMin%Vg SCM input geostrophic meridional wind at GCM sigma levels (m/s)
!@var SCMin%T SCM input absolute T at GCM sigma levels (K)
!@var SCMin%Q SCM input water vapor mixing ratio at GCM sigma levels (kg/kg)
!@var SCMin%Omega SCM input pressure tendency at GCM sigmal levels (mb/s)
!@var SCMin%W SCM input vertical wind at GCM sigmal levels (m/s)
!@var SCMin%O3 SCM input ozone at GCM sigmal levels (molec/atm-cm)
!@var SCMin%SadvV input vertical flux divergence of dry static energy / Cp (K/s)
!@var SCMin%QadvV input water vapor mixing ratio vertical flux div at GCM sigma levels (kg/kg/s)
!@var SCMin%TadvH input absolute T horizontal flux div at GCM sigma levels (K/s)
!@var SCMin%QadvH input water vapor mixing ratio horizontal flux div at GCM sigma levels (kg/kg/s)
!@var SCMin%UadvH input zonal wind div at GCM sigma levels (m/s/s)
!@var SCMin%VadvH input meridional wind div at GCM sigma levels (m/s/s)
!@var SCMin%Qrad input radiative heating rate profile at GCM sigma levels (W/m2)
!@var SCMin%Fnudge input nudging scale factor profile for qv and T (-)
!@var SCMin%time SCM input time (d)
!@var SCMin%lhf SCM input surface turbulent latent heat flux (W/m2)
!@var SCMin%shf SCM input surface turbulent sensible heat flux (W/m2)
!@var SCMin%Qskin SCM input surface skin water vapor mixing ratio (kg/kg)
!@var SCMin%Tskin SCM input surface skin temperature (K)
!@var SCMin%Ps SCM input surface pressure (mb)
!@var SCMin%z0m SCM input surface roughness height (m)
!@var SCMin%ustar SCM input surface friction velocity (m/s)
!@var SCMin%alb SCM input surface albedo (-)
!@var SCMin%BeersLaw_f0 SCM input cloud-top longwave cooling asymptote (W/m2)
!@var SCMin%BeersLaw_f1 SCM input cloud-base longwave heating asymptote (W/m2)
!@var SCMin%BeersLaw_kapp SCM input longwave absorption coefficient (m2/kg)

  end module SCM_COM

!-------------------------------------------------------------------------------

  subroutine alloc_SCM_COM()

  use filemanager, only : file_exists
  use Dictionary_mod, only : is_set_param,get_param
  use SCM_COM, only : SCMopt,SCMin
  implicit none
  real*8 dum_array(3)
  integer idum

  ! required and optional inputs

  SCMopt%Ps = file_exists('SCM_PS')

  SCMopt%wvmr = file_exists('SCM_WVMR')

  SCMopt%temp = file_exists('SCM_TEMP')
  SCMopt%theta = file_exists('SCM_THETA')
  if( SCMopt%temp .and. SCMopt%theta ) &
    call stop_model('alloc_SCM_COM: either T or theta required',255)

  SCMopt%Tskin = file_exists('SCM_TSKIN')
  SCMopt%Qskin = file_exists('SCM_QSKIN')
  SCMopt%ustar = file_exists('SCM_USTAR')
  SCMopt%sflx = file_exists('SCM_SFLUX')
  SCMopt%wind = file_exists('SCM_WIND')
  SCMopt%geo = file_exists('SCM_GEO')

  SCMopt%ozone = file_exists('SCM_OZONE')

  SCMopt%omega = file_exists('SCM_OMEGA')
  SCMopt%w = file_exists('SCM_W')
  SCMopt%ls_v = file_exists('SCM_LS_V')
  if( SCMopt%omega .and. SCMopt%w ) &
    call stop_model('alloc_SCM_COM: at most one of omega and w',255)
  if( SCMopt%ls_v .and. .not. SCMopt%omega ) &
    call stop_model( 'alloc_SCM_COM: omega needed for convergence',255)

  SCMopt%ls_h = file_exists('SCM_LS_H')
  SCMopt%ls_h_UV = file_exists('SCM_LS_H_UV')
  SCMopt%Qrad = file_exists('SCM_QRAD')
  SCMopt%Fnudge = file_exists('SCM_FNUDGE')

  ! F(default): assume piecewise linear relation between input layers when reading input profiles
  ! T(optional): assume top-hat distribution within input layers when reading input profiles
  !              (which requires uniform pressure grid)
  call get_param('SCM_TopHat',SCMopt%TopHat,default=.false.)

  ! T(default): standard operation
  ! F(optional): skip moist convection
  call get_param('SCM_allowMC',SCMopt%allowMC,default=.true.)

  ! T(default): standard operation
  ! F(optional): skip cloud-top entrainment instability
  call get_param('SCM_allowCTEI',SCMopt%allowCTEI,default=.true.)

  ! optional Beer's Law radiative heating (which takes 3 input parameters)
  SCMopt%BeersLaw = is_set_param('SCM_BeersLaw')
  if( SCMopt%BeersLaw )then
    call get_param('SCM_BeersLaw',dum_array,3)
    SCMin%BeersLaw_f0    = dum_array(1)
    SCMin%BeersLaw_f1    = dum_array(2)
    SCMin%BeersLaw_kappa = dum_array(3)
  endif
     
  if( SCMopt%Qrad .and. SCMopt%BeersLaw ) &
    call stop_model( 'alloc_SCM_COM: at most one of Qrad or BeersLaw',255)

  ! ignore longwave atmospheric heating associated with the surface
  ! when using parameterized or specified profile of radiative heating
  SCMopt%sfcQrad = .not. ( SCMopt%Qrad .or. SCMopt%BeersLaw )

  ! 0(default): surface is set by GCM input files in run deck
  call get_param('SCM_sfc',SCMopt%sfc,default=0)

  ! optional nudging
  SCMopt%nudge = is_set_param('SCM_tau')
  if( SCMopt%nudge ) call get_param('SCM_tau',SCMopt%tau)

  ! optional (gray) surface albedo
  SCMopt%alb = is_set_param('SCM_alb')
  if( SCMopt%alb ) call get_param('SCM_alb',SCMin%alb)

  ! optional surface roughness length for momentum
  SCMopt%z0m = is_set_param('SCM_z0m')
  if( SCMopt%z0m ) call get_param('SCM_z0m',SCMin%z0m)

  ! if not time-varying surface friction velocity above,
  ! optional fixed surface friction velocity
  if( .not. SCMopt%ustar )then
    SCMopt%ustar = is_set_param('SCM_ustar')
    if( SCMopt%ustar ) call get_param('SCM_ustar',SCMin%ustar)
  else
    if( is_set_param('SCM_ustar') ) &
      call stop_model('alloc_SCM_COM: redundant ustar values',255)
  endif

  if( SCMopt%z0m .and. SCMopt%ustar ) &
    call stop_model('alloc_SCM_COM: at most one of z0m or ustar',255)

  ! if not time-varying surface turbulent heat fluxes above,
  ! optional fixed sensible and latent heat fluxes
  if( .not. SCMopt%sflx )then
    SCMopt%sflx = is_set_param('SCM_shf')
    if( SCMopt%sflx )then 
      call get_param('SCM_shf',SCMin%shf)
      call get_param('SCM_lhf',SCMin%lhf)
    endif
  endif

  ! land surface currently requires specified heat fluxes, and ocean
  ! setting must be used for setting skin water vapor mixing ratio
  if( SCMopt%sfc.eq.1 )then
    if( .not.SCMopt%sflx ) &
      call stop_model('alloc_SCM_COM: land surface fluxes required',255)
    if( SCMopt%Qskin ) &
      call stop_model('alloc_SCM_COM: setting Qskin requires ocean setting',255)
  endif

  ! ocean surface requires specified heat fluxes or skin temperature
  if( SCMopt%sfc.eq.2 )then
    if( .not.SCMopt%sflx .and. .not.SCMopt%Tskin ) &
      call stop_model('alloc_SCM_COM: ocean surface fluxes required',255)
  endif

  ! optional vertical forcing of horizontal winds
  SCMopt%VadvHwind = is_set_param('SCM_VadvHwind')

  if( SCMopt%VadvHwind .and. .not. ( SCMopt%geo .and. ( SCMopt%omega .or. SCMopt%w ))) &
    call stop_model('alloc_SCM_COM: SCM_VadvHwind makes no sense')
        
  if( SCMopt%ls_h_UV .and. .not. SCMopt%geo ) &
    call stop_model('alloc_SCM_COM: SCMopt%ls_h_UV requires geostrophic wind')
        
  call get_param('SCM_PlumeDiag',SCMopt%PlumeDiag,default=.false.)

  end subroutine alloc_SCM_COM

!--------------------------------------------------------------------------------

  module GEOM
!@sum  GEOM contains geometric variables and arrays
!@auth M. Kelley

  implicit none
  save

  real*8, dimension(1,1) :: &
!@var  lat2d latitude of mid point of primary grid box (radians)
!@var  lat2d_dg latitude of mid point of primary grid box (degrees)
    lat2d,lat2d_dg &
!@var  lon2d longitude of mid point of primary grid box (radians)
!@var  lon2d_dg longitude of mid point of primary grid box (degrees)
   ,lon2d,lon2d_dg &
   ,lon_dg,lat_dg &
   ,sinlat2d, coslat2d

!@var  axyp,byaxyp area of grid box (+inverse) (m^2)
!@+    used for consistent conversion to and from extensive
!@+    units and in scale-aware cloud parameterization elements
!@+    (only the latter can influence SCM results).
  real*8, dimension(1,1) :: axyp, byaxyp

  integer, dimension(1) :: imaxj

  contains

  subroutine GEOM_ATM
  use CONSTANT, only : PI,TWOPI,radian
  use Dictionary_mod, only : get_param
  use SCM_COM, only : SCMopt
  implicit none

  ! mandatory rundeck parameters: lon and lat of target point
  call get_param('SCM_lon',SCMopt%lon)
  call get_param('SCM_lat',SCMopt%lat)

  if(abs(SCMopt%lon).gt.180d0 .or. abs(SCMopt%lat).gt.90d0) &
    call stop_model('geom_atm: invalid SCM_lon,SCM_lat in rundeck',255)

  lon2d_dg(1,1) = SCMopt%lon
  lat2d_dg(1,1) = SCMopt%lat

  call get_param('SCM_area',SCMopt%area)
  axyp(1,1) = SCMopt%area
  byaxyp(1,1) = 1d0/axyp(1,1)

  lon2d(1,1) = lon2d_dg(1,1)*radian
  lat2d(1,1) = lat2d_dg(1,1)*radian

  sinlat2d(1,1) = sin(lat2d(1,1))
  coslat2d(1,1) = cos(lat2d(1,1))
  lon2d(1,1) = lon2d(1,1) + pi ! IDL has a value of zero
  if(lon2d(1,1) .lt. 0.) lon2d(1,1)= lon2d(1,1) + twopi

  imaxj = 1

  lon_dg = lon2d_dg
  lat_dg = lat2d_dg

  end subroutine GEOM_ATM

  subroutine lonlat_to_ij(ll,ij)
  implicit none
  real*8, intent(in) :: ll(2)
  integer, intent(out) :: ij(2)
  ij = (/ 1, 1 /)
  end subroutine lonlat_to_ij

  end module GEOM
