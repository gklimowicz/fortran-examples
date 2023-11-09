! SCM.F90
!@sum Read SCM input files, initialize SCM, and update SCM inputs
!@auth Wolf (modifications by Fridlind)
!
!-------------------------------------------------------------------------------

  module SCM_mod
!@sum  SCM_mod contains variable declarations used by subroutines in SCM.f only
!@auth Fridlind
  implicit none

!@type SCMin_tProfile structured type for SCM input profiles on native
!vertical grid
  type SCMin_tProfile
!@var SCMin_tP%ntime number of simulation time steps
!@var SCMin_tP%nlev number of pressure levels
    integer ntime,nlev
!@var SCMin_tP%L input layer center pressure (mb) or altitude (m)
    real*8, dimension(:), allocatable :: L
!@var SCMin_tP%value input variable value, same units as SCMin
    real*8, dimension(:,:), allocatable :: value
  end type SCMin_tProfile
  type(SCMin_tProfile) SCMin_tU,SCMin_tV,SCMin_tUg,SCMin_tVg
  type(SCMin_tProfile) SCMin_tT,SCMin_tTH,SCMin_tQ
  type(SCMin_tProfile) SCMin_tOzone
  type(SCMin_tProfile) SCMin_tOmega,SCMin_tW
  type(SCMin_tProfile) SCMin_tSadvV,SCMin_tQadvV
  type(SCMin_tProfile) SCMin_tTadvH,SCMin_tQadvH
  type(SCMin_tProfile) SCMin_tUadvH,SCMin_tVadvH
  type(SCMin_tProfile) SCMin_tQrad,SCMin_tFnudge

!@type SCMin_tSscalar structured type for SCM input scalars
  type SCMin_tScalar
!@var SCMin_tS%ntime number of simulation time steps
    integer ntime
!@var SCMin_tS%value input scalar value in correct units
!@var SCMin_tS%f_scale scale factor to achieve correct units
!@var SCMin_tS%f_offset offset factor to achieve correct units
    real*8 :: f_scale,f_offset
    real*8, dimension(:), allocatable :: value
  end type SCMin_tScalar
  type(SCMin_tScalar) SCMin_tLHF,SCMin_tSHF
  type(SCMin_tScalar) SCMin_tTskin,SCMin_tQskin,SCMin_tPs
  type(SCMin_tScalar) SCMin_tUstar

!@type SCMreadXscalar structured type for reading SCM input scalars
  type SCMreadXscalar
!@var SCMreadX%ntime number of input time steps
    integer ntime
!@var SCMreadXscalar%year input calendar year
!@var SCMreadXscalar%month input calendar month
!@var SCMreadXscalar%day input calendar day
!@var SCMreadXscalar%hour input calendar hour
!@var SCMreadXscalar%value input scalar in provided units
    integer, dimension(:), allocatable :: year,month,day,hour
    real*8, dimension(:), allocatable :: value
  end type SCMreadXscalar
  type(SCMreadXscalar) SCMreadX

!@type SCMreadZprofile structured type for reading SCM input profiles
  type SCMreadZprofile
!@var SCMreadZ%ntime number of input time steps
!@var SCMreadZ%nlev number of vertical levels
    integer ntime,nlev
!@var SCMin_tP%L input layer center pressure (mb) or altitude (m)
    real*8, dimension(:), allocatable :: L
!@var SCMreadZ%year input calendar year
!@var SCMreadZ%month input calendar month
!@var SCMreadZ%day input calendar day
!@var SCMreadZ%hour input calendar hour
!@var SCMreadZ%value input profiles in provided units (see SCM.f)
    integer, dimension(:), allocatable :: year,month,day,hour
    real*8, dimension(:,:), allocatable :: value
  end type SCMreadZprofile
  type(SCMreadZprofile) SCMreadZ

!@type SCMnml_type structured type for reading SCM input variable namelist
  type SCMnml_type
!@var SCMnml_type%nvar number of input variables in namelist
    integer nvar
!@var SCMnml_type%model_varname reference variable name throughout SCM
!@var SCMnml_type%file_varname variable name in netCDF input file
    character(len=80) :: model_varname,file_varname
!@var SCMnml_type%f_scale applied to input data
!@var SCMnml_type%f_offset applied to input data
    real*8 f_scale,f_offset
  end type SCMnml_type
  type(SCMnml_type), dimension(:), allocatable :: SCMreadNML

! constant parameters used for SCM

!@param SCMp_zero   set to zero above topmost input level
!@param SCMp_one    set to one above topmost input level
!@param SCMp_const  extend upward using the topmost input level 
!@param SCMp_append append McClatchey et al. (1972) above topmost input level
  integer, parameter :: SCMp_zero   = 0, &
                        SCMp_one    = 1, &
                        SCMp_const  = 2, &
                        SCMp_append = 3

  end module SCM_mod

!-------------------------------------------------------------------------------

  subroutine init_SCM
  ! initialize SCM

  use SCM_com, only : SCMopt,SCMin
  use resolution, only : LM
  use atm_com, only : T,PK,Q,QCL,U,V,PMID
  use fluxes, only : atmocn,atmlnd
  use fluxes, only : FLAND,FOCEAN,FLICE,FLAKE0,FEARTH0
  use ghy_com, only : FEARTH
  use lakes_com, only : FLAKE
  use radpar, only : KEEPAL
  use constant, only : TF,LHE,SHA
  use CLOUDS_COM, only : CLDSAV 
  implicit none
  integer L
  real*8 dqsum,fcond

  ! report initialization
  write(6,*) 'SCM initializing ... '

  ! read input data, check for required fields, do all time interpolation
  write(6,*) ' ... SCM reading inputs ...'
  call read_SCM_inputs

  ! fill SCMin fields for the first time step, on the GCM vertical grid
  write(6,*) ' ... SCM initializing input variables ...'
  call update_SCM_inputs

  ! initialize atmospheric state
  write(6,*) ' ... SCM initializing atmospheric state ...'
  do L = 1,LM
    if( SCMopt%theta .or. SCMopt%temp )then
      T(1,1,L) = SCMin%T(L)/PK(L,1,1) ! temperature -> potential temperature
    endif
    if( SCMopt%wvmr )then
      Q(1,1,L) = SCMin%Q(L)
      ! saturation adjustment (relieve any supersaturation w/r/t liquid in
      ! initial state)
      call get_dq_cond( T(1,1,L), Q(1,1,L), PK(L,1,1), 1d0, lhe, pmid(L,1,1), dqsum, fcond )
      if( dqsum > 0d0 ) CLDSAV(L,1,1) = 1d0
      QCL(1,1,L) = dqsum
      Q(1,1,L) = Q(1,1,L)-dqsum
      T(1,1,L) = T(1,1,L)+dqsum*LHE/SHA/PK(L,1,1)
    endif
    ! horizontal winds may be specified throughout or else geostrophic
    ! throughout; if geostrophic, may be initialized using SCMopt%wind
    if( SCMopt%wind )then
      U(1,1,L) = SCMin%U(L)
      V(1,1,L) = SCMin%V(L)
    else if( SCMopt%geo )then
      U(1,1,L) = SCMin%Ug(L)
      V(1,1,L) = SCMin%Vg(L)
    endif
  enddo 

  ! initialize surface from run deck if requested
  !      Note: for SCM we are overwriting surface coverage fractions for
  !            land,ocean,landice that are read from GCM input files.
  !            If your box could have land or ocean ice you can consider
  !            how you want to treat that. For Ocean(or Lake) ice see
  !            si_atm%rsi or si_ocn%rsi   where rsi is the ratio of ice 
  !            coverage to water coverage.   
  write(6,*) ' ... SCM initializing surface state ...'
  FLICE(1,1)=0.0
  if( SCMopt%sfc > 0 )then
    if( SCMopt%sfc==1 )then
      FLAND(1,1) = 1.
      FOCEAN(1,1) = 0.
      FLAKE(1,1) = 0.
    else if( SCMopt%sfc==2 )then
      FLAND(1,1) = 0.
      FOCEAN(1,1) = 1.
      FLAKE(1,1) = 0.
    endif
  else if (SCMopt%sfc == -1) then
    FLAND(1,1) = 1. 
    FOCEAN(1,1) = 0.  
    FLAKE(1,1) = 0.
  endif
  FLAKE0(1,1) = FLAKE(1,1)
  FEARTH(1,1) = FLAND(1,1)-FLICE(1,1) 
  FEARTH0(1,1) = FEARTH(1,1)

  if( SCMopt%Tskin )then
    atmocn%GTEMP(1,1)  = SCMin%Tskin - TF
    atmocn%GTEMPR(1,1) = SCMin%Tskin
    atmlnd%GTEMP(1,1)  = SCMin%Tskin - TF
    atmlnd%GTEMPR(1,1) = SCMin%Tskin
  endif
  if( SCMopt%alb ) KEEPAL = 1 ! use surface albedo in run deck

  ! report initialization
  write(6,*) '... SCM initialization complete.'

  end subroutine init_SCM

!-------------------------------------------------------------------------------

  subroutine read_SCM_inputs
  ! read each input variable, convert units as specified in
  ! SCM input variable namelist, and interpolate to SCM simulation 
  ! times, extrapolating where required

  ! by design this routine relies on SCMopt switches making sense, 
  ! so any potential inconsistencies in those switches must be
  ! handled in subroutine alloc_SCM_COM

  use SCM_com, only : SCMopt
  use SCM_mod
  use filemanager, only : file_exists
  implicit none

  ! read variable namelist, which provides input file variable names
  ! if they differ from model names, and any unit conversions required
  ! (all input variables must be included in the namelist; 
  ! included variables that are not used will be ignored)
  call read_SCM_namelist('SCM_NML')

  ! specified sensible and latent heat fluxes (W/m2),
  if( SCMopt%sflx .and. file_exists('SCM_SFLUX') )then
    call read_SCM_scalar('SCM_SFLUX','LHF',SCMin_tLHF)
    call read_SCM_scalar('SCM_SFLUX','SHF',SCMin_tSHF)
  endif

  ! specified skin temperature (K)
  if( SCMopt%Tskin )then
    call read_SCM_scalar('SCM_TSKIN','Tskin',SCMin_tTskin)
  endif

  ! specified skin water vapor mixing ratio (kg/kg),
  ! currently overwriting values of ocean surface skin
  ! in order to represent dry idealized boundary layers
  if( SCMopt%Qskin )then
    call read_SCM_scalar('SCM_QSKIN','Qskin',SCMin_tQskin)
  endif

  ! specified surface friction velocity (m/s)
  if( SCMopt%ustar .and. file_exists('SCM_USTAR') )then
    call read_SCM_scalar('SCM_USTAR','Ustar',SCMin_tUstar)
  endif

  ! specified surface pressure (mb)
  if( SCMopt%Ps )then
    call read_SCM_scalar('SCM_PS','Ps',SCMin_tPs)
  endif

  ! specified atmospheric temperature or potential temperature (K)
  if( SCMopt%temp )then
    call read_SCM_profile('SCM_TEMP','T',SCMin_tT)
  endif
  if( SCMopt%theta )then
    call read_SCM_profile('SCM_THETA','TH',SCMin_tTH)
  endif

  ! specified atmospheric water vapor (kg/kg)
  if( SCMopt%wvmr )then
    call read_SCM_profile('SCM_WVMR','Q',SCMin_tQ)
  endif

  ! initial horizontal and geostrophic winds (m/s)
  if( SCMopt%wind )then
    call read_SCM_profile('SCM_WIND','U',SCMin_tU)
    call read_SCM_profile('SCM_WIND','V',SCMin_tV)
  endif
  if( SCMopt%geo )then
    call read_SCM_profile('SCM_GEO','Ug',SCMin_tUg)
    call read_SCM_profile('SCM_GEO','Vg',SCMin_tVg)
  endif

  ! specified ozone profile (kg/kg), to be used only where non-zero
  if( SCMopt%ozone )then
    call read_SCM_profile('SCM_OZONE','O3',SCMin_tOzone)
  endif

  ! specified large-scale vertical wind (mb/s) or forcing terms
  if( SCMopt%omega )then
    call read_SCM_profile('SCM_OMEGA','Omega',SCMin_tOmega)
  endif
  if( SCMopt%w )then
    call read_SCM_profile('SCM_W','W',SCMin_tW)
  endif
  if( SCMopt%ls_v )then
    call read_SCM_profile('SCM_LS_V','SadvV',SCMin_tSadvV)
    call read_SCM_profile('SCM_LS_V','QadvV',SCMin_tQadvV)
  endif

  ! specified large-scale horizontal tendences of temperature (K/s)
  ! and water vapor (kg/kg/s)
  if( SCMopt%ls_h )then
    call read_SCM_profile('SCM_LS_H','TadvH',SCMin_tTadvH)
    call read_SCM_profile('SCM_LS_H','QadvH',SCMin_tQadvH)
  endif

  ! specified large-scale horizontal tendences of winds (m/s/s)
  if( SCMopt%ls_h_UV )then
    call read_SCM_profile('SCM_LS_H_UV','UadvH',SCMin_tUadvH)
    call read_SCM_profile('SCM_LS_H_UV','VadvH',SCMin_tVadvH)
  endif

  ! specified radiative heating rate profile (K/s)
  if( SCMopt%Qrad )then
    call read_SCM_profile('SCM_QRAD','Qrad',SCMin_tQrad)
  endif

  ! specified thermodynamic nudging scale factor profile (-)
  if( SCMopt%Fnudge )then
    call read_SCM_profile('SCM_FNUDGE','Fnudge',SCMin_tFnudge)
  endif

  end subroutine read_SCM_inputs

!-------------------------------------------------------------------------------

  subroutine read_SCM_namelist(file_name)
  ! read all variable names and any conversion factors to expected units

  use SCM_mod
  use filemanager, only : openunit,closeunit
  implicit none
  character(len=*), intent(in) :: file_name
  character(len=80) :: model_varname,file_varname
  real*8 scale_factor,offset
  namelist/scm_input_var/model_varname,file_varname,scale_factor,offset
  integer iu_nml,nvar_nml,ivar,ios

  ! report reading
  write(6,*) ' ... SCM reading ',file_name

  ! read namelist once just to check and count variables
  call openunit('SCM_NML',iu_nml,.false.,.true.)
  nvar_nml = 0
  do
    model_varname = 'notset'
    read(iu_nml,nml=scm_input_var,iostat=ios)
    if( ios.ne.0 ) exit ! end of file
    if( trim(model_varname).eq.'notset' ) &
      call stop_model('SCM: namelist requires model_varname',255)
    nvar_nml = nvar_nml + 1
  enddo
  call closeunit(iu_nml)

  ! allocate and read all elements
  allocate(SCMreadNML(nvar_nml))
  call openunit('SCM_NML',iu_nml,.false.,.true.)
  do ivar = 1,nvar_nml
    ! set defaults
    offset = 0.
    scale_factor = 1.
    model_varname = 'notset'
    file_varname = 'notset'
    ! read values provided
    read(iu_nml,nml=scm_input_var,iostat=ios)
    SCMreadNML(ivar)%model_varname = trim(model_varname)
    if(trim(file_varname).eq.'notset')then
      SCMreadNML(ivar)%file_varname = SCMreadNML(ivar)%model_varname
    else
      SCMreadNML(ivar)%file_varname = trim(file_varname)
    endif
    SCMreadNML(ivar)%f_scale = scale_factor
    SCMreadNML(ivar)%f_offset = offset
    SCMreadNML(ivar)%nvar = nvar_nml ! redundant
    write(6,*) ' ... namelist variable names: ',trim(model_varname)
  enddo
  call closeunit(iu_nml)

  end subroutine read_SCM_namelist

!-------------------------------------------------------------------------------

  subroutine read_SCM_scalar(file_name,model_vname,SCMin_tS)
  ! read a single scalar and interpolate to SCM simulation times

  use SCM_mod
  implicit none
  character(len=*), intent(in) :: file_name,model_vname
  type(SCMin_tScalar), intent(inout) :: SCMin_tS
  character(len=80) file_vname
  real*8 factor_scale,factor_offset
  integer, dimension(:), allocatable :: read_intarray
  real*8, dimension(:), allocatable :: read_fltarray
  integer ncid,dimid,varid,stat,ivar
  include 'netcdf.inc'

  ! report reading
  write(6,*) ' ... SCM reading ',file_name,': ',model_vname

  ! obtain file variable name, unit offset, conversion from namelist
  file_vname='notset'
  do ivar =1,SCMreadNML(1)%nvar
    if( SCMreadNML(ivar)%model_varname.eq.model_vname )then
      file_vname = SCMreadNML(ivar)%file_varname
      factor_scale = SCMreadNML(ivar)%f_scale
      factor_offset = SCMreadNML(ivar)%f_offset
    endif
  enddo
  if( file_vname=='notset' ) call stop_model('SCM: scalar not found in namelist',255)

  ! read time dimension
  stat = nf_open(file_name,0,ncid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'nf_open')
  stat = nf_inq_dimid(ncid,'time',dimid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_time')
  stat = nf_inq_dimlen(ncid,dimid,SCMreadX%ntime)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'time')

  ! allocate temporary storage variables
  allocate(read_intarray(SCMreadX%ntime))
  allocate(read_fltarray(SCMreadX%ntime))
  allocate(SCMreadX%year(SCMreadX%ntime))
  allocate(SCMreadX%month(SCMreadX%ntime))
  allocate(SCMreadX%day(SCMreadX%ntime))
  allocate(SCMreadX%hour(SCMreadX%ntime))
  allocate(SCMreadX%value(SCMreadX%ntime))

  ! read time variables (integer)

  stat = nf_inq_varid(ncid,'year',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_year')
  stat = nf_get_var_int(ncid,varid,read_intarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'year')
  SCMreadX%year = read_intarray

  stat = nf_inq_varid(ncid,'month',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_mon')
  stat = nf_get_var_int(ncid,varid,read_intarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'month')
  SCMreadX%month = read_intarray

  stat = nf_inq_varid(ncid,'day',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_day')
  stat = nf_get_var_int(ncid,varid,read_intarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'day')
  SCMreadX%day = read_intarray

  stat = nf_inq_varid(ncid,'hour',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_hour')
  stat = nf_get_var_int(ncid,varid,read_intarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'hour')
  SCMreadX%hour = read_intarray

  ! read scalar time series (double precision)

  stat = nf_inq_varid(ncid,file_vname,varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,file_vname)
  stat = nf_get_var_double(ncid,varid,read_fltarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,file_vname)
  read_fltarray = read_fltarray * factor_scale
  read_fltarray = read_fltarray + factor_offset
  SCMreadX%value = read_fltarray

  stat = nf_close(ncid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'close')

  ! interpolate scalar time series to SCM time steps
  call interp_t_SCM_scalar(SCMreadX,SCMin_tS)

  deallocate(read_intarray)
  deallocate(read_fltarray)
  deallocate(SCMreadX%year)
  deallocate(SCMreadX%month)
  deallocate(SCMreadX%day)
  deallocate(SCMreadX%hour)
  deallocate(SCMreadX%value)

  end subroutine read_SCM_scalar

!-------------------------------------------------------------------------------

  subroutine interp_t_SCM_scalar(SCMreadXt,SCMin_tS)
  ! read a single scalar and interpolate to SCM simulation times,
  ! using fractional year as common time unit for interpolation

  use SCM_mod
  use model_com, only : ItimeI,ItimeE,Iyear1,hourI,DTsrc,modelEclock
  use JulianCalendar_mod, only : JDendOfM
  use TimeConstants_mod, only : SECONDS_PER_YEAR
  implicit none
  type(SCMreadXscalar), intent(inout) :: SCMreadXt
  type(SCMin_tScalar), intent(inout) :: SCMin_tS
  real*8, allocatable :: SCMin_time(:),SCMread_time(:)
  real*8 ft1,ft2,dt
  integer day,it_SCM,it_read

  ! allocate scalar size to total number of simulation time steps,
  ! where the first time step is time zero by convention
  SCMin_tS%ntime = ItimeE - ItimeI + 1
  allocate(SCMin_tS%value(0:SCMin_tS%ntime-1))

  ! SCM times to interpolate to, using fractional years as common unit
  allocate(SCMin_time(SCMin_tS%ntime))
  call modelEclock%get(dayOfYear=day)
  SCMin_time(1) = Iyear1 + (day+hourI*1d0/24.)/365.
  dt = DTsrc*1d0/3600./24./365.
  do it_SCM = 2,SCMin_tS%ntime
    SCMin_time(it_SCM) = SCMin_time(1) + dt*(it_SCM-1)
  enddo

  ! input times to interpolate from, fractional years as common unit
  allocate(SCMread_time(SCMreadXt%ntime))
  do it_read = 1,SCMreadXt%ntime
    if( SCMreadXt%month(it_read) == 2 .and. SCMreadXt%day(it_read) > 28 ) &
      call stop_model('SCM: no leap years in ModelE',255)
    SCMread_time(it_read) = SCMreadXt%year(it_read)*1d0 + &
                          ( JDendOfM(SCMreadXt%month(it_read)-1) + &
                          SCMreadXt%day(it_read) + &
                          SCMreadXt%hour(it_read)*1d0/24. )/365.
  enddo

  ! interpolate from read times to input times, extrapolate if needed
  do it_SCM = 1,SCMin_tS%ntime
    if( SCMin_time(it_SCM) <= SCMread_time(1) )then
      SCMin_tS%value(it_SCM-1) = SCMreadXt%value(1)
    else if( SCMin_time(it_SCM) > SCMread_time(SCMreadXt%ntime) )then
      SCMin_tS%value(it_SCM-1) = SCMreadXt%value(SCMreadXt%ntime)
    else
      do it_read = 1,SCMreadXt%ntime-1
        if( SCMin_time(it_SCM) > SCMread_time(it_read) .and. &
            SCMin_time(it_SCM) <= SCMread_time(it_read+1) )then
          ft2 = (SCMin_time(it_SCM)-SCMread_time(it_read))/ &
                (SCMread_time(it_read+1)-SCMread_time(it_read))
          ft1 = 1. - ft2
          SCMin_tS%value(it_SCM-1) = ft1*SCMreadXt%value(it_read) + &
                                     ft2*SCMreadXt%value(it_read+1)
          exit
        endif
      enddo
    endif
  enddo

  deallocate(SCMin_time)
  deallocate(SCMread_time)

  end subroutine interp_t_SCM_scalar

!-------------------------------------------------------------------------------

  subroutine read_SCM_profile(file_name,model_vname,SCMin_tP)
  ! read a single profile and interpolate to SCM simulation times
  ! (must interpolate to GCM pressure grid later, at each SCM time step)

  use SCM_mod
  implicit none
  character(len=*), intent(in) :: file_name,model_vname
  character(len=80) file_vname
  type(SCMin_tProfile), intent(inout) :: SCMin_tP
  integer, dimension(:), allocatable :: read_intarray
  real*8, dimension(:), allocatable :: read_fltarray,read_levarray
  real*8, dimension(:,:), allocatable :: read_fltprofs
  integer ncid,dimid,varid,stat,Ldata,ivar
  real*8 factor_scale,factor_offset
  include 'netcdf.inc'

  ! report reading
  write(6,*) ' ... SCM reading ',file_name,': ',model_vname

  ! obtain file variable name, unit offset, conversion from namelist
  file_vname = 'notset'
  do ivar =1,SCMreadNML(1)%nvar
    if( SCMreadNML(ivar)%model_varname.eq.model_vname )then
      file_vname = SCMreadNML(ivar)%file_varname
      factor_scale = SCMreadNML(ivar)%f_scale
      factor_offset = SCMreadNML(ivar)%f_offset
    endif
  enddo
  if( file_vname=='notset' ) call stop_model('SCM: profile not found in namelist',255)

  ! read time and height dimensions
  stat = nf_open(file_name,0,ncid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'nf_open')
  stat = nf_inq_dimid(ncid,'time',dimid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_time')
  stat = nf_inq_dimlen(ncid,dimid,SCMreadZ%ntime)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'time')
  ! vertical grid must be pressure coordinate
  stat = nf_inq_dimid(ncid,'lev',dimid)
  stat = nf_inq_dimlen(ncid,dimid,SCMreadZ%nlev)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'lev')

  ! allocate temporary storage variables
  allocate(read_intarray(SCMreadZ%ntime))
  allocate(read_levarray(SCMreadZ%nlev))
  allocate(read_fltarray(SCMreadZ%ntime))
  allocate(read_fltprofs(SCMreadZ%nlev,SCMreadZ%ntime))
  allocate(SCMreadZ%year(SCMreadZ%ntime))
  allocate(SCMreadZ%month(SCMreadZ%ntime))
  allocate(SCMreadZ%day(SCMreadZ%ntime))
  allocate(SCMreadZ%hour(SCMreadZ%ntime))
  allocate(SCMreadZ%value(SCMreadZ%ntime,SCMreadZ%nlev))
  allocate(SCMreadZ%L(SCMreadZ%nlev))

  ! read time variables and profile time series
  ! as a function of pressure (mb) or altitude (m)
  ! coordinates that do not change in time

  ! read time variables (integer)

  stat = nf_inq_varid(ncid,'year',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_year')
  stat = nf_get_var_int(ncid,varid,read_intarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'year')
  SCMreadZ%year = read_intarray

  stat = nf_inq_varid(ncid,'month',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_mon')
  stat = nf_get_var_int(ncid,varid,read_intarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'month')
  SCMreadZ%month = read_intarray

  stat = nf_inq_varid(ncid,'day',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_day')
  stat = nf_get_var_int(ncid,varid,read_intarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'day')
  SCMreadZ%day = read_intarray

  stat = nf_inq_varid(ncid,'hour',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_hour')
  stat = nf_get_var_int(ncid,varid,read_intarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'hour')
  SCMreadZ%hour = read_intarray

  ! read pressure profile (double precision)

  stat = nf_inq_varid(ncid,'lev',varid)
  if(stat.ne.NF_NOERR) stat = nf_inq_varid(ncid,'P',varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'id_P')
  stat = nf_get_var_double(ncid,varid,read_levarray)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'P')
  SCMreadZ%L = read_levarray
  ! check for progression from surface upwards
  do Ldata = 2,SCMreadZ%nlev
    if( SCMreadZ%L(Ldata) > SCMreadZ%L(Ldata-1) ) &
      call stop_model('SCM: now requires decreasing P grid',255)
  enddo

  ! read profile time series (double precision)

  stat = nf_inq_varid(ncid,file_vname,varid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,file_vname)
  stat = nf_get_var_double(ncid,varid,read_fltprofs)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,file_vname)
  read_fltprofs = read_fltprofs * factor_scale
  read_fltprofs = read_fltprofs + factor_offset
  ! reorder time and level indices
  do Ldata = 1,SCMreadZ%nlev
    SCMreadZ%value(1:SCMreadZ%ntime,Ldata) = read_fltprofs(Ldata,1:SCMreadZ%ntime)
  enddo

  stat = nf_close(ncid)
  if(stat.ne.NF_NOERR) call handle_err(stat,file_name,'close')

  ! interpolate profile time series to SCM time steps
  call interp_t_SCM_profile(SCMreadZ,SCMin_tP)

  deallocate(read_intarray)
  deallocate(read_fltarray)
  deallocate(read_fltprofs)
  deallocate(SCMreadZ%year)
  deallocate(SCMreadZ%month)
  deallocate(SCMreadZ%day)
  deallocate(SCMreadZ%hour)
  deallocate(SCMreadZ%L)
  deallocate(SCMreadZ%value)

  end subroutine read_SCM_profile

!-------------------------------------------------------------------------------

  subroutine interp_t_SCM_profile(SCMreadZt,SCMin_tP)
  ! read a profile time series and interpolate to SCM simulation times,
  ! using fractional year as common time unit for interpolation;
  ! assumes input vertical grid is constant in time

  use SCM_mod
  use model_com, only : ItimeI,ItimeE,Iyear1,hourI,DTsrc,modelEclock
  use JulianCalendar_mod, only : JDendOfM
  use TimeConstants_mod, only : SECONDS_PER_YEAR
  implicit none
  type(SCMreadZprofile), intent(inout) :: SCMreadZt
  type(SCMin_tProfile), intent(inout) :: SCMin_tP
  real*8, allocatable :: SCMin_time(:),SCMread_time(:)
  real*8 ft1,ft2,dt
  integer day,it_SCM,it_read

  ! number of levels and times for output
  SCMin_tP%nlev = SCMreadZt%nlev
  SCMin_tP%ntime = ItimeE - ItimeI + 1
  if( allocated(SCMin_tP%value) )then
    deallocate(SCMin_tP%value)
    deallocate(SCMin_tP%L)
  endif
  allocate(SCMin_tP%value(0:SCMin_tP%ntime-1,SCMin_tP%nlev))

  ! output vertical grid
  allocate(SCMin_tP%L(SCMin_tP%nlev))
  SCMin_tP%L = SCMreadZt%L

  ! SCM times to interpolate to, using fractional years as common unit
  allocate(SCMin_time(SCMin_tP%ntime))
  call modelEclock%get(dayOfYear=day)
  SCMin_time(1) = Iyear1 + (day+hourI*1d0/24.)/365.
  dt = DTsrc*1d0/3600./24./365.
  do it_SCM = 2,SCMin_tP%ntime
    SCMin_time(it_SCM) = SCMin_time(1) + dt*(it_SCM-1)
  enddo

  ! input times to interpolate from, fractional years as common unit
  allocate(SCMread_time(SCMreadZt%ntime))
  do it_read = 1,SCMreadZt%ntime
    SCMread_time(it_read) = SCMreadZt%year(it_read)*1d0 + &
      ( JDendOfM(SCMreadZt%month(it_read)-1) + &
                 SCMreadZt%day(it_read) + &
                 SCMreadZt%hour(it_read)*1d0/24. )/365.
  enddo

  ! interpolate from read times to input times, extrapolate if needed
  do it_SCM = 1,SCMin_tP%ntime
    if( SCMin_time(it_SCM) <= SCMread_time(1) )then
      SCMin_tP%value(it_SCM-1,:) = SCMreadZt%value(1,:)
    else if( SCMin_time(it_SCM) > SCMread_time(SCMreadZt%ntime) )then
      SCMin_tP%value(it_SCM-1,:) = SCMreadZt%value(SCMreadZt%ntime,:)
    else
      do it_read = 1,SCMreadZt%ntime-1
        if( SCMin_time(it_SCM) > SCMread_time(it_read) .and. &
            SCMin_time(it_SCM) <= SCMread_time(it_read+1) )then
          ft2 = (SCMin_time(it_SCM)-SCMread_time(it_read))/ &
                (SCMread_time(it_read+1)-SCMread_time(it_read))
              ft1 = 1. - ft2
          SCMin_tP%value(it_SCM-1,:) = &
               ft1*SCMreadZt%value(it_read,:) + &
               ft2*SCMreadZt%value(it_read+1,:)
          exit
        endif
      enddo
    endif
  enddo

  deallocate(SCMin_time)
  deallocate(SCMread_time)

  end subroutine interp_t_SCM_profile

!-------------------------------------------------------------------------------

  subroutine interp_p_SCM_profile(SCMin_tP,SCMinP,model_vname,i_above)
  ! interpolate from input profile grid to GCM vertical grid
  ! (input profiles already interpolated to SCM time steps);
  ! treatment above input grid depends on i_above;
  ! force potential temperature to modestly increase with height 
  ! in the appended layers

  use SCM_com, only : nstepSCM, SCMin, SCMopt
  use SCM_mod
  use constant, only : LHE
  use resolution, only : LM
  use atm_com, only : PMID,PEDN,PK,PDSIG
  implicit none
  character(len=*), intent(in) :: model_vname
  type(SCMin_tProfile), intent(in) :: SCMin_tP
  real*8, intent(inout) :: SCMinP(LM)
  integer, intent(in) :: i_above
  real*8 fp,dPsum,dPadd,dPtot
  real*8 rh_below,th_min
  real*8, parameter :: stab_min = 0.02 ! minimum appended dth/dp (K/mb)
  integer L,Lgcm,Lbot,Ltop
  real*8, external :: QSAT
  real*8 p_input(SCMin_tP%nlev),psi_input(SCMin_tP%nlev)
  real*8 pe_input(SCMin_tP%nlev+1)
  integer n_input,method
  real*8 Pgcm,Pgcm_up,Pgcm_dn,p_up,p_dn,dp,psi_up,psi_dn,total,weight,slope

  ! interpolate from input on pressure grid, using pressure-weighted averages
  ! to conserve vertical integrals, with one of two methods:
  ! (a) piecewise linear between input points 
  ! (b) top-hat distribution within each input layer (requires input dP=constant)

  ! copy some input variables for ease of notation
  p_input(:)   = SCMin_tP%L(:)
  psi_input(:) = SCMin_tP%value(nstepSCM,:)
  n_input      = SCMin_tP%nlev

  if( SCMopt%TopHat )then
  ! create pressure edges for top-hat method
    dp = SCMin_tP%L(1)-SCMin_tP%L(2)
    do L = 2, n_input-1
      if( abs(dp-(SCMin_tP%L(L)-SCMin_tP%L(L+1))) > epsilon(1d0) )then
        print*,'top-hat method requires constant dP in input'
        call stop_model('bad pressure grid in interp_p_SCM_profile',255)
      endif
    enddo
    pe_input(1) = SCMin_tP%L(1)+dp/2.
    do L = 1, n_input
      pe_input(L+1) = SCMin_tP%L(L)-dp/2.
    enddo
  endif

  ! loop over GCM layers
  do Lgcm = 1, LM

    if( model_vname=='Omega' .or. model_vname=='W' )then
    ! vertical winds located at GCM lower edges
      if( Lgcm==1 )then ! zero vertical wind at surface
        SCMinP(Lgcm) = 0.
        cycle
      endif
      Pgcm    = Pedn(Lgcm,  1,1)
      Pgcm_up = Pmid(Lgcm,  1,1)
      Pgcm_dn = Pmid(Lgcm-1,1,1)
    else
    ! other variables located at GCM layer centers
      Pgcm    = Pmid(Lgcm,  1,1)
      Pgcm_up = Pedn(Lgcm+1,1,1)
      Pgcm_dn = Pedn(Lgcm,  1,1)
    endif

    if(( .not. SCMopt%TopHat .and. Pgcm_dn > p_input(1)  ) .or. &
       ( SCMopt%TopHat       .and. Pgcm_dn > pe_input(1) ))then
    ! use lowermost specfied value when relevant GCM pressure below input grid

      SCMinP(Lgcm) = psi_input(1)

    else if(( .not. SCMopt%TopHat .and. Pgcm_up >= p_input(n_input)    ) .or. &
            (  SCMopt%TopHat      .and. Pgcm_up >= pe_input(n_input+1) ))then
    ! relevant GCM pressures within input grid

      if( .not. SCMopt%TopHat )then
      ! integrate using piecewise linear assumption

        ! locate upper and lower GCM pressures on input grid
        do L = 1, n_input-1
          if( Pgcm_dn >= p_input(L+1) ) exit
        enddo
        Lbot=L
        do L = Lbot, n_input-1
          if( Pgcm_up >= p_input(L+1) ) exit
        enddo
        Ltop=L

        if( Lbot==Ltop )then
        ! GCM upper and lower pressures fall between same adjacent points on input grid

          psi_dn = psi_input(Lbot)
          psi_up = psi_input(Lbot+1)
          p_dn   = p_input(Lbot)
          p_up   = p_input(Lbot+1)
          slope  = (psi_dn-psi_up)/(p_dn-p_up)
          SCMinP(Lgcm) = psi_up + slope*(Pgcm-p_up)

        else
        ! GCM pressures span at least one point on input grid

          total  = 0.
          weight = 0.
          do L = Lbot, Ltop
            psi_dn = psi_input(L)
            psi_up = psi_input(L+1)
            p_dn   = p_input(L)
            p_up   = p_input(L+1)
            slope  = (psi_dn-psi_up)/(p_dn-p_up)
            if( L==Ltop )then
              dp = p_dn-Pgcm_up
              total = total + 0.5*dp*(psi_dn+psi_up+slope*(Pgcm_up-p_up))
            else if( L==Lbot )then
              dp = Pgcm_dn-p_up
              total = total + 0.5*dp*(psi_up+psi_up+slope*dp)
            else
              dp = p_dn-p_up
              total = total + 0.5*dp*(psi_up+psi_dn)
            endif
            weight = weight + dp
          enddo
          SCMinP(Lgcm) = total/weight

        endif ! Lbot==Ltop or not

      else ! SCMopt%TopHat
      ! integrate using top-hat assumption

        ! locate upper and lower GCM pressures on input grid
        do L = 1, n_input
          if( Pgcm_dn >= pe_input(L+1) ) exit
        enddo
        Lbot=L
        do L = Lbot, n_input
          if( Pgcm_up >= pe_input(L+1) ) exit
        enddo
        Ltop=L

        total  = 0.
        weight = 0.
        p_dn   = pe_input(L)
        p_up   = pe_input(L+1)
        do L = Lbot, Ltop
          if( L==Ltop )then
            dp = p_dn-Pgcm_up
          else if( L==Lbot )then
            dp = Pgcm_dn-p_up
          else
            dp = p_dn-p_up
          endif
          total = total + dp*psi_input(L)
          weight = weight + dp
        enddo
        SCMinP(Lgcm) = total/weight

      endif ! SCMopt%TopHat or not

    else if( .not. SCMopt%TopHat .and. Pgcm_dn >= p_input(n_input) )then
    ! use uppermost specfied value when bottom of GCM layer within input grid,
    ! only when not assuming top-hat distribution (for backwards compatibility)

      SCMinP(Lgcm) = psi_input(n_input)

    else
    ! relevant GCM pressures above input grid:
    ! fill overlying data according to i_above

      select case (i_above)
      case (SCMp_zero) 
        SCMinP(Lgcm) = 0.
      case (SCMp_one) 
        SCMinP(Lgcm) = 1.
      case (SCMp_const)
        SCMinP(Lgcm) = psi_input(n_input)
      case (SCMp_append)
        call append_SCM_profile(Pgcm,SCMinP(Lgcm),model_vname)
      case default
        print*,'stop in interp_p_SCM_profile: no such i_above=',i_above
        call stop_model('bad i_above in interp_p_SCM_profile')
      end select

      ! force potential temperature to modestly increase with height
      ! to avoid triggering moist convection

      if( model_vname=='T' )then
        th_min = SCMinP(Lgcm-1)/PK(Lgcm-1,1,1) + stab_min*(Pmid(Lgcm-1,1,1)-Pmid(Lgcm,1,1))
        SCMinP(Lgcm) = PK(Lgcm,1,1)*max(SCMinP(Lgcm)/PK(Lgcm,1,1),th_min)
      endif
      if( model_vname=='TH' )then
       th_min = SCMinP(Lgcm-1) + stab_min*(Pmid(Lgcm-1,1,1)-Pmid(Lgcm,1,1))
        SCMinP(Lgcm) = max(SCMinP(Lgcm),th_min)
      endif

      ! prevent relative humidity from increasing with height
      ! (requires this routine to be called for T before Q,
      ! for which there is no check)

      if( model_vname=='Q' )then
        rh_below = SCMinP(Lgcm-1)/QSAT(SCMin%T(Lgcm-1),LHE,Pmid(Lgcm-1,1,1))
        SCMinP(Lgcm) = min( SCMinP(Lgcm),rh_below*QSAT(SCMin%T(Lgcm),LHE,Pmid(Lgcm,1,1)) )
      endif

    endif ! location of GCM pressure levels
  enddo   ! loop over GCM pressure layers

  end subroutine interp_p_SCM_profile

!-------------------------------------------------------------------------------

  subroutine append_SCM_profile(P_layer,X_value,model_vname)
  ! add climatological T and qv profiles above input data layers:
  ! McClatchey data (1972) for standard profile using PHATMO subroutine
  ! (see phatmo subroutine header for full documentation; only
  ! inputs and outputs used here are defined here)

  use SCM_com, only : SCMopt
  use SCM_mod
  use model_com, only : modelEclock
  use constant, only : KAPA
  implicit none
  character(len=*), intent(in) :: model_vname
  real*8, intent(in) :: P_layer    ! layer pressure (mb)
  real*8, intent(inout) :: X_value ! q (kg/kg) or T (K) value
  real*8 H,D,T,O,Q,S,OCM,WCM       ! only T and Q used here
  integer :: NPHD=1                ! input P, return T and Q option
  integer NATM                     ! latitude, season option
  integer INDATM(12,3)             ! NATM indexed by month, latitude
  integer month                    ! current month
  integer INDLAT                   ! latitude index
  data  INDATM/ 1,1,1,1,1,1,1,1,1,1,1,1, &
                3,3,3,2,2,2,2,2,2,3,3,3, &
                5,5,5,4,4,4,4,4,4,5,5,5/

  ! find closest match from season and latitude options
  call modelEclock%get(month=month)
  if( abs(SCMopt%lat) <= 25. ) then
    INDLAT = 1
  else if( abs(SCMopt%lat) <= 50. ) then
    INDLAT = 2
  else
    INDLAT = 3
  endif
  NATM = INDATM(month,INDLAT)

  call PHATMO(P_layer,H,D,T,O,Q,S,OCM,WCM,NPHD,NATM)

  select case (model_vname)
  case ('Q') 
    X_value = Q
  case ('T') 
    X_value = T
  case ('TH') 
    X_value = T/(P_layer/1000.)**KAPA
  case default
    print*,'variable not available in append_SCM_profile:',model_vname
    call stop_model('no such variable in append_SCM_profile',255)
  end select

  end subroutine append_SCM_profile

!-------------------------------------------------------------------------------

  subroutine update_SCM_inputs
  ! fill SCMin variables for the current time step

  ! by design this routine relies on SCMopt switches making sense, 
  ! so any potential inconsistencies in those switches must be
  ! handled in subroutine alloc_SCM_COM

  use SCM_com, only : SCMopt,SCMin,nstepSCM
  use SCM_mod
  use resolution, only : LM
  use atm_com, only : P,PMID,PEDN,PK,T,PDSIG
  use fluxes, only : atmocn,atmlnd
  use constant, only : TF,KAPA,GRAV,RGAS
  use filemanager, only : file_exists
  implicit none
  integer L
  real*8 DZ(LM)

  ! specified surface pressure
  if( SCMopt%Ps )then
    SCMin%Ps = SCMin_tPs%value(nstepSCM)
    PEDN(1,1,1) = SCMin%Ps
    call CALC_AMPK(LM)
  endif

  ! specified surface fluxes
  if( SCMopt%sflx .and. file_exists('SCM_SFLUX') )then
    SCMin%lhf = SCMin_tLHF%value(nstepSCM)
    SCMin%shf = SCMin_tSHF%value(nstepSCM)
  endif

  ! specified skin temperature (K)
  if( SCMopt%Tskin )then
    SCMin%Tskin = SCMin_tTskin%value(nstepSCM)
    atmocn%GTEMP(1,1)  = SCMin%Tskin - TF
    atmocn%GTEMPR(1,1) = SCMin%Tskin
    atmlnd%GTEMP(1,1)  = SCMin%Tskin - TF
    atmlnd%GTEMPR(1,1) = SCMin%Tskin
  endif

  ! specified skin water vapor mixing ratio (kg/kg)
  if( SCMopt%Qskin ) SCMin%Qskin = SCMin_tQskin%value(nstepSCM)

  if( SCMopt%ustar .and. file_exists('SCM_USTAR') ) &
    SCMin%ustar = SCMin_tUstar%value(nstepSCM)

  ! specified temperature or potential temperature with 1000-mb ref
  if( SCMopt%temp )then
    call interp_p_SCM_profile(SCMin_tT,SCMin%T,'T',SCMp_append)
  endif
  if( SCMopt%theta )then
    call interp_p_SCM_profile(SCMin_tTH,SCMin%TH,'TH',SCMp_append)
    do L = 1,LM
      SCMin%T(L) = SCMin%TH(L)*(PMID(L,1,1)/1000.)**KAPA
    enddo
  endif

  ! specified water vapor mixing ratio
  ! (this call needs SCMin%T to be current)
  if( SCMopt%wvmr )then
    call interp_p_SCM_profile(SCMin_tQ,SCMin%Q,'Q',SCMp_append)
  endif

  ! initial horizontal and geostrophic winds
  if( SCMopt%wind )then
    call interp_p_SCM_profile(SCMin_tU,SCMin%U,'U',SCMp_const)
    call interp_p_SCM_profile(SCMin_tV,SCMin%V,'V',SCMp_const)
  endif
  if( SCMopt%geo )then
    call interp_p_SCM_profile(SCMin_tUg,SCMin%Ug,'Ug',SCMp_const)
    call interp_p_SCM_profile(SCMin_tVg,SCMin%Vg,'Vg',SCMp_const)
  endif

  ! specified ozone profile (molec/atm-cm), to be used only where non-zero
  ! (convert from input units of kg/kg)
  if( SCMopt%ozone )then
    call interp_p_SCM_profile(SCMin_tOzone,SCMin%O3,'O3',SCMp_zero)
    do L = 1,LM
      DZ(L) = PDSIG(L,1,1)/PMID(L,1,1) &
            * (RGAS/GRAV)*T(1,1,L)*PK(L,1,1)
      SCMin%O3(L) = SCMin%O3(L)*DZ(L)*PMID(L,1,1)*100. &
                  / (RGAS*T(1,1,L)*2.14E-2)
    enddo
  endif

  ! large-scale forcing terms
  ! vertical velocity always applied through omega, in pressure units (mb/s)
  if( SCMopt%omega )then
    call interp_p_SCM_profile(SCMin_tOmega,SCMin%Omega,'Omega',SCMp_zero)
  endif
  if( SCMopt%w )then
    call interp_p_SCM_profile(SCMin_tW,SCMin%W,'W',SCMp_zero)
    do L = 1,LM
      SCMin%Omega(L) = -SCMin%W(L)*GRAV*PMID(L,1,1)/(RGAS*T(1,1,L)*PK(L,1,1))
    enddo
  endif
  if( SCMopt%ls_v )then
    call interp_p_SCM_profile(SCMin_tSadvV,SCMin%SadvV,'SadvV',SCMp_zero)
    call interp_p_SCM_profile(SCMin_tQadvV,SCMin%QadvV,'QadvV',SCMp_zero)
  endif
  if( SCMopt%ls_h )then
    call interp_p_SCM_profile(SCMin_tTadvH,SCMin%TadvH,'TadvH',SCMp_zero)
    call interp_p_SCM_profile(SCMin_tQadvH,SCMin%QadvH,'QadvH',SCMp_zero)
  endif
  if( SCMopt%ls_h_UV )then
    call interp_p_SCM_profile(SCMin_tUadvH,SCMin%UadvH,'UadvH',SCMp_zero)
    call interp_p_SCM_profile(SCMin_tVadvH,SCMin%VadvH,'VadvH',SCMp_zero)
  endif

  ! specified radiative heating profile
  if( SCMopt%Qrad )then
    call interp_p_SCM_profile(SCMin_tQrad,SCMin%Qrad,'Qrad',SCMp_zero)
  endif

  ! specified thermodynamic nudging scale factor profile
  if( SCMopt%Fnudge )then
    call interp_p_SCM_profile(SCMin_tFnudge,SCMin%Fnudge,'Fnudge',SCMp_one)
  endif

  end subroutine update_SCM_inputs

!-------------------------------------------------------------------------------

  subroutine handle_err(errcode,file_name,file_string)
  implicit none
  character(len=*), intent(in) :: file_name,file_string
  include 'netcdf.inc'
  integer errcode

  print*,'Error reading file: ',file_name
  print*,'Error string: ',file_string
  print*,'Error code: ', nf_strerror(errcode)
  call stop_model('SCM handle_err: netcdf stop here',255)

  end subroutine handle_err
