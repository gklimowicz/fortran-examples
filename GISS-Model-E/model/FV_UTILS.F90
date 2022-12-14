#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort_core(__LINE__,rc)

module FV_UTILS

  USE ESMF
  USE CONSTANT, only: KAPA

  implicit none
  private

  public :: FV_CORE
  public :: abort_core
  public :: ReverseLevels
  public :: ConvertPotTemp_GISS2FV
  public :: ConvertPotTemp_FV2GISS
  public :: ConvertPressure_GISS2FV
  public :: ConvertPressure_FV2GISS
  public :: EdgePressure_GISS
  public :: DryTemp_GISS
  public :: DeltPressure_DryTemp_GISS
  public :: DeltPressure_GISS
  public :: compute_phi
  public :: init_app_clock
  public :: MoveFile
  public :: AllocateFvExport3D
  public :: DumpState
  public :: load_configuration
  public :: SetupForESMF
  public :: HumidityInit
  public :: Copy_modelE_to_FV_import
  public :: Tendency
  public :: ClearTendencies
  public :: SaveTendencies
  public :: clear_accumulated_mass_fluxes

  Interface ReverseLevels
     Module Procedure reverse_3d_r8
     Module Procedure reverse_3d_r4
  End Interface

  Interface ConvertPotTemp_GISS2FV
     Module Procedure CnvPotTemp_GISS2FV_r8
     Module Procedure CnvPotTemp_GISS2FV_r4
  End Interface

  Interface ConvertPotTemp_FV2GISS
     Module Procedure CnvPotTemp_FV2GISS_r4
  End Interface

  Interface ConvertPressure_GISS2FV
     Module Procedure ConvertPressure_GISS2FV_r4
     Module Procedure ConvertPressure_GISS2FV_r8
  End Interface

  ! Public parameters
  character(len=*), parameter, public :: FVCORE_INTERNAL_RESTART = 'dyncore_internal_restart'
  character(len=*), parameter, public :: FVCORE_IMPORT_RESTART   = 'dyncore_import_restart'
  character(len=*), parameter, public :: FVCORE_LAYOUT           = 'fvcore_layout.rc'

  character(len=*), parameter, public :: TENDENCIES_FILE = 'tendencies_checkpoint'

  Real*8, parameter, public :: PRESSURE_UNIT_GISS  =  100 ! 1 mb
  Real*8, parameter, public :: PRESSURE_UNIT_FV    =    1 ! 1 pa
  Real*8, parameter, public :: PRESSURE_UNIT_RATIO = PRESSURE_UNIT_GISS/PRESSURE_UNIT_FV

  integer, parameter, public :: INITIAL_START = 2 ! ISTART=2: cold start, no FV restart files
  integer, parameter, public :: EXTEND_RUN = 3    ! ISTART>2: FV restart files are read

  ! This data structure is o convenient entity for storing persistent data between
  ! calls to this module.  In addition to
  Type FV_CORE

     type (esmf_gridcomp) :: gc   ! This is the handle for the fv dynamical core

     type (esmf_grid)     :: grid ! Although modelE is not an ESMF component, it does have an ESMF_Grid
     type (esmf_vm)       :: vm   ! Should be eliminated ... Only used to obtain NPES for config file.

     ! Import and Export states for FV dycore
     type(esmf_state) :: import   ! Allocated within FV component
     type(esmf_state) :: export   ! Allocated within FV component

     ! The following pointers can be re-extracted from the import state at each iteration,
     ! but it is convenient to have a simpler means of access.
     real*4, pointer, dimension(:,:,:) :: dudt, dvdt, dtdt, dpedt  ! Tendencies
     real*4, pointer, dimension(:,:,:) :: Q  ! Humidity
     real*4, pointer, dimension(:,:,:) :: Qtr  ! other tracers
     real*4, pointer, dimension(:,:)   :: phis

     ! modelE does not work directly with tendencies.  Instead, tendencies are derived
     ! by differencing before and after physics.  Therefore, the final dynamical state
     ! must be preserved for the following
     ! of modelE fields
     real*8, pointer, dimension(:,:,:) :: U_old, V_old, dPT_old, PE_old, dT_old

  END Type FV_CORE

  ! private data of convenience
  Integer :: rc ! return code from ESMF

  ! The following parameters address the fact that FV and modelE use
  ! different units and different reference pressures for potential temperature.
  ! Superficially there is redundancy between the two sets, but in some sense
  ! this is merely coincidental.
  Real*8, parameter :: REF_PRESSURE_GISS = 100 ! 1 mb = 100 pa
  Real*8, parameter :: REF_PRESSURE_FV   =   1 ! pa
  Real*8, parameter :: REF_RATIO = REF_PRESSURE_FV / REF_PRESSURE_GISS

contains

  !----------------------------------------------------------------
  ! The following routine interfaces the VERIFY_ macro with the
  ! GISS termination routine.   The line number and return code are
  ! written into a buffer which is passed to stop_model().
  !----------------------------------------------------------------
  Subroutine abort_core(line,rc)
    Implicit None
    Integer, Intent(In) :: line
    Integer, Intent(In) :: rc

    Character(len=100) :: buf

    Write(buf,*)'FV core failure at line',line,'\n error code:', rc
    Call stop_model(Trim(buf),99)

  End Subroutine abort_core

  ! Reverse the order of vertical levels.
  Function reverse_3d_r8(A) Result(B)
    Real*8, Intent(In) :: A(:,:,:)
    Real*8             :: B(Size(A,1),Size(A,2),size(A,3))

    Integer, parameter :: K_IDX = 3
    Integer :: k, n


    n = Size(A, K_IDX)

    Do k = 1, n
       B(:,:,k) = A(:,:,1+n-k)
    End Do

  End Function reverse_3d_r8

  ! Single precision variant of Reverse()
  Function reverse_3d_r4(A) Result(B)
    Real*4, Intent(In) :: A(:,:,:)
    Real*4             :: B(Size(A,1),Size(A,2),size(A,3))

    Integer, parameter :: K_IDX = 3
    Integer :: k, n

    n = Size(A, K_IDX)

    Do k = 1, n
       B(:,:,k) = A(:,:,1+n-k)
    End Do

  End Function reverse_3d_r4

  function EdgePressure_GISS() Result(PE)
    USE RESOLUTION, only: IM, LM, LS1, Ptop, PSFMPT
    Use DYNAMICS, only : SIG, SIGE
    use atm_com, only : P
    use domain_decomp_atm, only: grid, getDomainBounds

    REAL*8 :: PE(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM+1)

    INTEGER :: L, i_0, i_1, j_0, j_1

    call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

    Do L = 1, LM+1

       If (L < LS1) THEN
          PE(:,:,L) = SIGE(L)*P(I_0:I_1,J_0:J_1) + Ptop
       Else
          PE(:,:,L)  = SIGE(L)*PSFMPT + Ptop
       End IF

    End Do

  end function EdgePressure_GISS

  ! Compute Delta-pressure for GISS model
  function DeltPressure_GISS() Result(dP)
    USE RESOLUTION, only: IM, LM, LS1
    Use ATM_COM, only: T
    USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds

    REAL*8 :: dP(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM)
    REAL*8 :: PE(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM+1)

    INTEGER :: k

    PE = EdgePressure_GISS()
    do k = 1, LM
       dP(:,:,k) = (PE(:,:,k)-PE(:,:,k+1))
    end do

  end function DeltPressure_GISS

  ! Convert Potential Temperature into (dry) Temperature
  function DeltPressure_DryTemp_GISS() Result(dPT)
    USE RESOLUTION, only: IM, LM, LS1
    Use ATM_COM, only: T
    USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds

    REAL*8 :: dPT(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM)
    REAL*8 :: PE(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM+1)
    REAL*8 :: T_dry(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM)

    INTEGER :: J_0,J_1,k
    Call GetDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

    T_dry = DryTemp_GISS()
    PE = EdgePressure_GISS()
    do k = 1, LM
       dPT(:,:,k) = (PE(:,:,k)-PE(:,:,k+1)) * T_dry(:,:,k)
    end do

  end function DeltPressure_DryTemp_GISS

  ! Convert Potential Temperature into (dry) Temperature
  function DryTemp_GISS() Result(T_dry)
    USE RESOLUTION, only: IM, LM, LS1
    Use ATM_COM, only: T
    USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds

    REAL*8 :: T_dry(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM)
    REAL*8 :: PKZ(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM)

    INTEGER :: I_0, I_1, J_0,J_1
    Call GetDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

    PKZ = PKZ_GISS()
    T_dry = PKZ * T(I_0:I_1,J_0:J_1,:)

  end function DryTemp_GISS

  function PKZ_GISS() Result(PKZ)
    USE RESOLUTION, only: IM, LM, LS1, Ptop, PSFMPT
    Use DYNAMICS, only : SIG
    use atm_com, only : P
    use domain_decomp_atm, only: grid, getDomainBounds

    REAL*8 :: PKZ(grid % I_STRT:grid % I_STOP,grid % J_STRT:grid % J_STOP,LM)

    INTEGER :: L
    INTEGER :: I_0, I_1, J_0,J_1
    Call GetDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

    Do L = 1, LM

       If (L < LS1) THEN
          PKZ(:,:,L) = (SIG(L)*P(I_0:I_1,J_0:J_1) + Ptop) ** KAPA
       Else
          PKZ(:,:,L) = (SIG(L)*PSFMPT + Ptop) ** KAPA
       End IF

    End Do
  end function PKZ_GISS

  subroutine ConvertPressure_GISS2FV_r4(P_giss, P_fv)
    Real*8, intent(in) :: P_giss(:,:,:)
    Real*4, intent(out) :: P_fv(:,:,:)   ! no halo in this case

    P_fv = ReverseLevels(P_giss * PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_GISS2FV_r4

  subroutine ConvertPressure_GISS2FV_r8(P_giss, P_fv)
    Real*8, intent(in) :: P_giss(:,:,:)
    Real*8, intent(out) :: P_fv(:,:,:)   ! no halo in this case

    P_fv = ReverseLevels(P_giss * PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_GISS2FV_r8

  ! Convert pressure from GISS representation to FV representation.
  ! Both the order of levels and the units must be adjusted.
  subroutine ConvertPressure_FV2GISS(P_fv, P_giss)
    Real*4, intent(in) :: P_fv(:,:,:)
    Real*8, intent(out) :: P_giss(:,:,:)

    P_giss = ReverseLevels(P_fv / PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_FV2GISS

  ! Convert potential temperature between the two representations.
  subroutine CnvPotTemp_GISS2FV_r8(PT_giss, PT_fv)
    Real*8, intent(in) :: PT_giss(:,:,:)
    Real*8, intent(out) :: PT_fv(:,:,:)

    PT_fv = ReverseLevels(PT_giss * REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_GISS2FV_r8

  ! Convert potential temperature between the two representations.
  subroutine CnvPotTemp_GISS2FV_r4(PT_giss, PT_fv)
    Real*8, intent(in) :: PT_giss(:,:,:)
    Real*4, intent(out) :: PT_fv(:,:,:)

    PT_fv = ReverseLevels(PT_giss * REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_GISS2FV_r4

  !------------------------------------------------------------------------
  ! Convert potential temperature as exported by FV to the corresponding
  ! potential temperauture within modelE.  Note that FV export uses
  ! a reference pressure of 10^5 pa.
  !------------------------------------------------------------------------
  subroutine CnvPotTemp_FV2GISS_r4(PT_fv, PT_giss)
    Real*4, intent(in) :: PT_fv(:,:,:)
    Real*8, intent(out) :: PT_giss(:,:,:)

    ! As an export, FV provides Pot Temp with a reference
    ! pressure of 10^5 Pa, whereas GISS uses 100 Pa (1 mb)
    REAL*8, PARAMETER :: REF_RATIO = 100000 / 100
    integer :: m,n

    ! PT_GISS has a halo, while PT_fv does not.
    m = size(PT_GISS,1) ! exclude halo
    n = size(PT_GISS,2) ! exclude halo
    PT_giss(2:m-1,2:n-1,:) = ReverseLevels(PT_fv / REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_FV2GISS_r4

  function compute_phi(P, T, SZ, zatmo) result(phi)
    USE CONSTANT, only : rgas,bykapa,bykapap1,bykapap2
    USE DYNAMICS, only: DSIG, SIG, SIGE
    USE RESOLUTION, only: IM, JM, LM, LS1, LM, PTOP, PSFMPT
    USE DOMAIN_DECOMP_ATM, Only: grid, getDomainBounds
    implicit none

    real*8, intent(in), dimension(grid%i_strt_halo:grid%i_stop_halo, &
                                  grid%j_strt_halo:grid%j_stop_halo) :: p,zatmo
    real*8, intent(in), dimension(grid%i_strt_halo:grid%i_stop_halo, &
                                  grid%j_strt_halo:grid%j_stop_halo,lm) :: t,sz
    real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
                      grid%j_strt_halo:grid%j_stop_halo,lm) :: phi

    REAL*8 :: PKE(LS1:LM+1)
    REAL*8 :: PIJ, PHIDN
    REAL*8 :: DP, BYDP, P0, TZBYDP, X
    REAL*8 :: PUP, PKUP, PKPUP, PKPPUP
    REAL*8 :: PDN, PKDN, PKPDN, PKPPDN
    INTEGER :: I, J, L

    INTEGER :: I_0, I_1, J_0, J_1
    LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

    call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1, &
         & HAVE_NORTH_POLE = HAVE_NORTH_POLE, &
         & HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)

    DO L=LS1,LM+1
       PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
    END DO

!$OMP  PARALLEL DO PRIVATE(I,J,L,DP,P0,PIJ,PHIDN,TZBYDP,X,
!$OMP*             BYDP,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP)
    DO J=J_0,J_1

       DO I=I_0,I_1

          PIJ=P(I,J)

          PDN=PIJ+PTOP
          PKDN=PDN**KAPA
          PHIDN=ZATMO(I,J)

          !**** LOOP OVER THE LAYERS
          DO L=1,LM
             PKPDN=PKDN*PDN
             PKPPDN=PKPDN*PDN
             IF(L.GE.LS1) THEN
                DP=DSIG(L)*PSFMPT
                BYDP=1./DP
                P0=SIG(L)*PSFMPT+PTOP
                TZBYDP=2.*SZ(I,J,L)*BYDP
                X=T(I,J,L)+TZBYDP*P0
                PUP=SIGE(L+1)*PSFMPT+PTOP
                PKUP=PKE(L+1)
                PKPUP=PKUP*PUP
                PKPPUP=PKPUP*PUP
             ELSE
                DP=DSIG(L)*PIJ
                BYDP=1./DP
                P0=SIG(L)*PIJ+PTOP
                TZBYDP=2.*SZ(I,J,L)*BYDP
                X=T(I,J,L)+TZBYDP*P0
                PUP=SIGE(L+1)*PIJ+PTOP
                PKUP=PUP**KAPA
                PKPUP=PKUP*PUP
                PKPPUP=PKPUP*PUP
             END IF
             !**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
             PHI(I,J,L)=PHIDN+RGAS*(X*PKDN*BYKAPA-TZBYDP*PKPDN*BYKAPAP1 &
                  &      -(X*(PKPDN-PKPUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2) &
                  &      *BYDP*BYKAPAP1)
             !**** CALULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
             PHIDN=PHIDN+RGAS*(X*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPDN-PKPUP) &
                  &     *BYKAPAP1)
             PDN=PUP
             PKDN=PKUP
          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO

    !**** SET POLAR VALUES FROM THOSE AT I=1
    IF (HAVE_SOUTH_POLE) THEN
       DO L=1,LM
          PHI(2:IM,1,L)=PHI(1,1,L)
       END DO
    END IF

    IF (HAVE_NORTH_POLE) THEN
       DO L=1,LM
          PHI(2:IM,JM,L)=PHI(1,JM,L)
       END DO
    END IF

  end function compute_phi

  function init_app_clock(start_time, end_time, interval) Result(clock)
    integer :: start_time(6)
    integer :: end_time(6)
    integer :: interval
    type (esmf_clock)              :: clock

    type (esmf_time) :: startTime
    type (esmf_time) :: stopTime
    type (esmf_timeinterval) :: timeStep
    type (esmf_calendar) :: gregorianCalendar

    ! initialize calendar to be Gregorian type
    gregorianCalendar = esmf_calendarcreate(ESMF_CALKIND_GREGORIAN, name="ApplicationCalendar", rc=rc)
    VERIFY_(rc)

    call ESMF_CalendarSetDefault(ESMF_CALKIND_GREGORIAN, rc=rc)
    VERIFY_(rc)   

    ! initialize start time
    write(*,*)'Time Set Start: ',START_TIME
    call esmf_timeset(startTime, YY=START_TIME(1), MM= START_TIME(2), &
         & DD=START_TIME(3), H=START_TIME(4),                         &
         & M=START_TIME(5),  S=START_TIME(6),                         &
         & calendar=gregorianCalendar, rc=rc)
    VERIFY_(rc)

    ! initialize stop time
    write(*,*)'Time Set End: ',END_TIME
    call esmf_timeset(stopTime, YY=END_TIME(1), MM= END_TIME(2), &
         & DD=END_TIME(3), H=END_TIME(4),                        &
         & M=END_TIME(5),  S=END_TIME(6),                        &
         & calendar=gregorianCalendar, rc=rc)
    VERIFY_(rc)

    ! initialize time interval
    call esmf_timeintervalset(timeStep, &
         & S=INT(interval), rc=rc)
    VERIFY_(rc)

    ! initialize the clock with the above values
    clock = esmf_clockcreate(timeStep, startTime, stoptime=stopTime, name="ApplClock",rc=rc)
    VERIFY_(rc)

    call ESMF_ClockSet ( clock, CurrTime=startTime, rc=rc)
    VERIFY_(rc)

  end function init_app_clock

  subroutine MoveFile(fv_fname, fv_dfname, suffix)
    character(len=*), intent(in) :: fv_fname, fv_dfname, suffix

    call system('mv ' // FVCORE_INTERNAL_RESTART // suffix // ' ' // trim(fv_fname) )
    call system('mv ' // FVCORE_IMPORT_RESTART // ' ' // trim(fv_dfname) )
  end subroutine MoveFile

  subroutine AllocateFvExport3D ( state, name )
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    type(ESMF_State),  intent(INOUT) :: state
    character(len=*),  intent(IN   ) :: name

    real, pointer :: ptr(:,:,:)
    integer       :: status
    logical       :: alloc

    alloc = .true.
    call ESMFL_StateGetPointerToData ( state, ptr , name , alloc , rc=status )
    VERIFY_(status)

  end subroutine AllocateFvExport3D
#ifdef CUBED_SPHERE
#define USE_MAPL
#endif

#ifdef USE_MAPL
  subroutine DumpState(fv, clock, fv_fname, fv_dfname, suffix, isFinalize)
    use MAPL_mod, only: MAPL_MetaComp
    use MAPL_mod, only: MAPL_Get
    use MAPL_mod, only: MAPL_GetResource
    use MAPL_mod, only: MAPL_GetObjectFromGC
    use MAPL_mod, only: MAPL_GetResource
    use domain_decomp_atm, only: AM_I_ROOT

    Type (FV_CORE),    intent(inout) :: fv
    type (esmf_clock) :: clock
    character(len=*), intent(in) :: fv_fname, fv_dfname, suffix
    logical, intent(in) :: isFinalize

    type (MAPL_MetaComp), pointer :: internalState
    type (ESMF_State) :: esmfInternalState
    integer :: hdr,rc

    call SaveTendencies(fv, FVCORE_IMPORT_RESTART)

    if(isFinalize) then
       call ESMF_GridCompFinalize( fv % gc, importstate=fv%import, exportstate=fv%export, clock=clock, rc=rc)
    else
       call MAPL_GetObjectFromGC( fv%gc, internalSTate)
       call MAPL_Get(internalState, internal_ESMF_state=ESMFInternalState, rc=rc)
       VERIFY_(rc)
       call MAPL_GetResource( internalState   , hdr,         &
                               default=0, &
                               LABEL="INTERNAL_HEADER:", &
                               RC=rc)
       VERIFY_(rc)
       call StateWriteToFile(ESMFInternalState, clock, &
            & FVCORE_INTERNAL_RESTART // suffix,       &
            & 'binary',internalState, hdr==1, rc=rc)
       VERIFY_(rc)
    endif

    ! Now move the file into a more useful name
    if (AM_I_ROOT()) then
      call MoveFile(fv_fname, fv_dfname, suffix)
    end if

  end subroutine DumpState

!-------------------------------------------------------------------------------
  subroutine StateWriteToFile(STATE,CLOCK,FILENAME,FILETYPE,MPL,HDR,RC)
!-------------------------------------------------------------------------------
    use MAPL_mod, only: MAPL_MetaComp
    Use MAPL_IOMod, only: GETFILE, Free_file, MAPL_VarWrite, Write_parallel
    USE RESOLUTION, only: IM, JM, LM
    type(ESMF_State),                 intent(INOUT) :: STATE
    type(ESMF_Clock),                 intent(IN   ) :: CLOCK
    character(len=*),                 intent(IN   ) :: FILENAME
    character(LEN=*),                 intent(IN   ) :: FILETYPE
    type(MAPL_MetaComp),              intent(INOUT) :: MPL
    logical,                          intent(IN   ) :: HDR
    integer, optional,                intent(  OUT) :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="StateWriteToFile"

    type (ESMF_StateItem_Flag), pointer   :: ITEMTYPES(:)
    character(len=ESMF_MAXSTR ), pointer  :: ITEMNAMES(:)
    integer                               :: ITEMCOUNT
    integer                               :: UNIT
    integer                               :: I, J, K, L, N
    integer                               :: YYYY, MM, DD, H, M, S
    type(ESMF_Time)                       :: currentTime
    integer                               :: HEADER(6)

! Get information from state
!---------------------------

    call ESMF_StateGet(STATE,ITEMCOUNT=ITEMCOUNT,RC=rc)
    VERIFY_(rc)

    allocate(ITEMNAMES(ITEMCOUNT),STAT=rc)
    VERIFY_(rc)
    allocate(ITEMTYPES(ITEMCOUNT),STAT=rc)
    VERIFY_(rc)

    call ESMF_StateGet(STATE,ITEMNAMELIST=ITEMNAMES, &
         & ITEMTYPELIST=ITEMTYPES,RC=rc)
    VERIFY_(rc)

! Open file
!----------

    if (filetype == 'binary' .or. filetype == 'BINARY') then
       UNIT = GETFILE(FILENAME, form="unformatted", rc=rc)
       VERIFY_(rc)
    elseif(filetype=="formatted".or.filetype=="FORMATTED") then
       UNIT = GETFILE(FILENAME, form="formatted", rc=rc)
       VERIFY_(rc)
    else
       UNIT=0
    end if

! Write data
!-----------

    call ESMF_ClockGet (clock, currTime=currentTime, rc=rc)
    VERIFY_(rc)
    call ESMF_TimeGet(CurrentTime, &
      &   YY=YYYY, MM=MM, DD=DD,   &
      &   H=H, M=M, S=S, rc=rc)
    VERIFY_(rc)

    HEADER(1) = YYYY
    HEADER(2) = MM
    HEADER(3) = DD
    HEADER(4) = H
    HEADER(5) = M
    HEADER(6) = S
    
    call Write_Parallel(HEADER, UNIT, RC=rc)
    VERIFY_(rc)
    
    HEADER(1) = IM
    HEADER(2) = JM
    HEADER(3) = LM
    HEADER(4) = 0
    HEADER(5) = 0
    
    call Write_Parallel(HEADER(1:5), UNIT, RC=rc)
    VERIFY_(rc)

    if(UNIT/=0) then
       do J = 1, ITEMCOUNT
          call MAPL_VarWrite(UNIT=UNIT, STATE=STATE, NAME=ITEMNAMES(J), rc=rc)
          VERIFY_(rc)
       end do
       call FREE_FILE(UNIT)
    else
       rc = -1  ! not yet
       VERIFY_(rc)
    endif
    
    deallocate(ITEMNAMES) 
    deallocate(ITEMTYPES)

  end subroutine StateWriteToFile

#else

  subroutine DumpState(fv, clock, fv_fname, fv_dfname, suffix, isFinalize)
    use GEOS_mod, only: RecordPhase=>GEOS_RecordPhae
    use domain_decomp_atm, only: AM_I_ROOT

    Type (FV_CORE),    intent(inout) :: fv
    type (esmf_clock), intent(in) :: clock
    character(len=*), intent(in) :: fv_fname, fv_dfname, suffix
    logical, intent(in) :: isFinalize

    call SaveTendencies(fv, FVCORE_IMPORT_RESTART)

    if(isFinalize) then
       call ESMF_GridCompFinalize( fv % gc, fv%import, fv%export, clock, rc=rc)
    else
       call ESMF_GridCompFinalize( fv % gc, fv % import, fv % export, clock, &
            &  phase=RecordPhase, rc=rc)
    endif

    ! Now move the file into a more useful name
    if (AM_I_ROOT()) then
      call MoveFile(fv_fname, fv_dfname, suffix)
    end if

  end subroutine DumpState
#endif

  function load_configuration(config_file) result( config )
    use ESMF
    use FILEMANAGER
    Use MODEL_COM,  only: DT=>DTsrc
    character(len=*), parameter :: Iam="FV_INTERFACE::loadconfiguration"
    character(len=*), intent(in) :: config_file
    type (esmf_config)           :: config

    integer :: iunit
    type (ESMF_VM) :: vm

    config = esmf_configcreate(rc=rc)
    VERIFY_(rc)

    call openunit(config_file, iunit, qbin=.false., qold=.false.)
    write(iunit,*)'FVCORE_INTERNAL_CHECKPOINT_FILE:  ', FVCORE_INTERNAL_RESTART
    write(iunit,*)'FVCORE_INTERNAL_RESTART_FILE:     ', FVCORE_INTERNAL_RESTART
!!$ write(iunit,*)'FVCORE_IMPORT_CHECKPOINT_FILE:    ', TENDENCIES_FILE
!!$ write(iunit,*)'FVCORE_IMPORT_RESTART_FILE:       ', TENDENCIES_FILE
    write(iunit,*)'FVCORE_LAYOUT:                    ', FVCORE_LAYOUT
    write(iunit,*)'RUN_DT:                           ', DT
    call closeUnit(iunit)

    Call ESMF_VMGetGlobal(vm, rc=rc)
    call esmF_VMbarrier(vm, rc=rc)
    call esmf_configloadfile(config, config_file, rc=rc)
    VERIFY_(rc)

  end function load_configuration

  subroutine SetupForESMF(fv, vm, grid, cf, config_file)
    use DOMAIN_DECOMP_ATM, only : AM_I_ROOT
#ifdef CUBED_SPHERE
    use FVCubed_dycore_GridCompMod, only: SetServices
#else
    use fvdycore_gridcompmod, only: SetServices
#endif

    Type(FV_CORE), intent(inout) :: fv
    type (esmf_vm),   intent(in) :: vm
    type (esmf_grid) :: grid
    type (esmf_config), intent(inout) :: cf
    character(len=*),  intent(in) :: config_file ! filename for resource file

#ifdef CUBED_SPHERE
    character(len=*), parameter :: gridCompName = 'FVCORE'
#else
    character(len=*), parameter :: gridCompName = 'FV dynamics'
#endif

    fv % vm  =vm
    fv % grid=grid

    ! Load configuration information from resource file
    !  ------------------------------------------------
    cf = load_configuration(config_file)

    !  Create the dynamics component, using same layout as application
    !  ------------------------------------------------------------------
    fv % gc = ESMF_GridCompCreate ( name=gridCompName, &
         & grid=grid, &
         & config=cf, rc=rc)
    VERIFY_(rc)

    ! Create couplings
    fv % import = ESMF_StateCreate ( name='fv dycore imports', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc )
    VERIFY_(rc)
    fv % export = ESMF_StateCreate ( name='fv dycore exports', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc )
    VERIFY_(rc)

    !  Register services for components
    !  --------------------------------
    !   print*,'calling set services'
    call ESMF_GridCompSetServices ( fv % gc, SetServices, rc=rc )
    VERIFY_(rc)

    ! 1) Create layout resource file - independent of actual restart file.
    If (AM_I_ROOT()) Then
       Call Write_Layout(FVCORE_LAYOUT, fv)
    End If

    contains

      Subroutine write_layout(fname, fv)
#ifdef CUBED_SPHERE
        Use MAPL_IOMod, only: GETFILE, Free_file
#else
        Use GEOS_IOMod, only: GETFILE, Free_file
#endif

      USE RESOLUTION, only: IM, JM, LM, LS1
      Use MODEL_COM,  only: DT=>DTsrc
      USE DOMAIN_DECOMP_ATM, only : grid
      character(len=*), intent(in) :: fname
      type (fv_core), intent(in)   :: fv

      Integer :: unit
      Integer :: npes
      Integer :: mppnx, mppny


      npes = grid%nproc
      unit = GetFile(fname, form="formatted", rc=rc)

      mppnx = 0
      mppny = 0
      if(mod(NPES,6) == 0) then
         mppnx = int(floor(sqrt(real(NPES/6))))
         mppny = (NPES / mppnx) / 6
      endif

      ! Lat-Lon parameters
      write(unit,*)' # empty line #'
      Write(unit,*)'xy_yz_decomp:',1,npes,npes,1
      Write(unit,*)'    im: ',IM
      Write(unit,*)'    jm: ',JM
      Write(unit,*)'    km: ',LM
      Write(unit,*)'    dt: ',DT
      Write(unit,*)' ntotq: ',1
      Write(unit,*)'  imxy: ',IM
      Write(unit,'(a,100(1x,I3))')'  jmxy: ',1+grid%jer(:)-grid%jsr(:) ! the original code is wrong.
      Write(unit,'(a,100(1x,I3))')'  jmyz: ',1+grid%jer(:)-grid%jsr(:) ! why are we even writing this.
      Write(unit,*)'   kmyz: ',LM

      ! Common parameters
      Write(unit,*)'nsplit: ',0
      Write(unit,*)'    nq: ',1

      ! Cubed-Sphere parameters
      Write(unit,*)'      npx: ',IM
      Write(unit,*)'      npy: ',JM
      Write(unit,*)'      npz: ',LM
      Write(unit,*)'       dt: ',DT
      Write(unit,*)'   npes_x: ',mppnx
      Write(unit,*)'   npes_y: ',mppny
      Write(unit,*)' consv_te: ',0.
      Write(unit,*)' dgrid_imports: ',.true.

      Call Free_File(unit)

    End Subroutine write_layout

  end subroutine SetupForESMF

  subroutine HumidityInit(fv)
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData

    Type(FV_CORE), intent(inout) :: fv

    type (ESMF_FieldBundle)        :: bundle
    type (ESMF_Field)              :: Qfield

    ! Specific Humidity - need to reserve space in FV import state
    call ESMFL_StateGetPointerToData (fv%export, fv%q, 'Q', rc=rc)
    VERIFY_(rc)
    fv%q = 0.0

    Qfield = ESMF_FieldCreate( fv%grid, fv%q, name='Q', rc=rc)
    VERIFY_(rc)

    !call stop_model('HumidityInit: use TRADV string?',255)

    ! First obtain a reference field from the state (avoids tedious details)
    call ESMF_StateGet(fv%import, 'TRADV'  , bundle, rc=rc)
    VERIFY_(rc)

    Call ESMF_FieldBundleAdd(bundle, (/Qfield/), rc=rc)
    VERIFY_(rc)

  end subroutine HumidityInit

  Subroutine Copy_modelE_to_FV_import(fv)
    USE ATM_COM, only:  Q     ! Secific Humidity
    USE ATM_COM, only:  ZATMO ! Geopotential Height?
    USE ATM_COM, Only : U, V, T
    Use Domain_decomp_atm, Only: grid, getDomainBounds
    Type (FV_CORE) :: fv

    Integer :: nq

    Integer :: I_0, I_1, j_0, j_1

    Call GetDomainBounds(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1)

    ! 1) Link/copy modelE data to import state
    fv % phis=ZATMO(I_0:I_1,j_0:j_1)
#ifdef NO_FORCING
    fv % phis = 0
#endif

    ! Moisture
#ifndef ADIABATIC
    fv % Q = ReverseLevels(Q(I_0:I_1,j_0:j_1,:))
#ifdef NO_FORCING
    fv % Q = 0
#endif
#endif

  End Subroutine Copy_modelE_to_FV_import

  !-----------
  ! d/dt A
  ! Note that FV needs tendency with the order of vertical levels reversed.
  !-----------
  Function Tendency(A, A_old) Result(tend)
    USE MODEL_COM, Only: DTsrc
    REAL*8, INTENT(IN) :: A(:,:,:)
    REAL*8, INTENT(IN) :: A_old(:,:,:)
    REAL*8 :: tend(size(a,1),size(a,2),size(a,3))

    tend = (A - A_old)/DTsrc

  End Function Tendency

  subroutine ClearTendencies(fv)
    Type (FV_CORE) :: fv

    fv % dudt = 0
    fv % dvdt = 0
    fv % dtdt = 0
    fv % dpedt = 0

  end subroutine ClearTendencies

  subroutine SaveTendencies(fv, fv_dfname)
    use domain_decomp_atm, only: AM_I_ROOT
    use FILEMANAGER

    type (fv_core),    intent(inout) :: fv
    integer :: iunit
    character(len=*), intent(in) :: fv_dfname

    call OpenUnit(trim(fv_dfname) , iunit, qbin=.true.)

    call SaveArr(iunit, fv % U_old)
    call SaveArr(iunit, fv % V_old)
    call SaveArr(iunit, fv % dPT_old)
    call SaveArr(iunit, fv % PE_old)
    call SaveArr(iunit, fv % dT_old)

    call closeunit(iunit)

  contains

    subroutine SaveArr(iunit, arr)
      use domain_decomp_atm, only: grid, dwrite8_parallel, getDomainBounds
      integer, intent(in) :: iunit
      real*8, intent(in) :: arr(:,:,:)
      real*8, allocatable :: padArr(:,:,:)
      integer :: I_0,  I_1,  J_0,  J_1
      integer :: I_0H, I_1H, J_0H, J_1H

      Call GetDomainBounds(grid, i_strt=I_0, i_stop=I_1, j_strt=J_0, j_stop=J_1, &
           & i_strt_halo = I_0H, i_stop_halo = I_1H, &
           & j_strt_halo = J_0H, j_stop_halo = J_1H)
      allocate(padArr(I_0H:I_1H,J_0H:J_1H, size(arr,3)))

      padArr(I_0:I_1,J_0:J_1,:) = arr(:,:,:)
      call dwrite8_parallel(grid, iunit, nameunit(iunit), padArr)

      deallocate(padArr)

    end subroutine SaveArr

  end subroutine SaveTendencies

  subroutine clear_accumulated_mass_fluxes()
    USE ATM_COM, ONLY: MUS,MVs,MWs

    MUs(:,:,:) = 0.0
    MVs(:,:,:) = 0.0
    MWs(:,:,:) = 0.0

  end subroutine clear_accumulated_mass_fluxes

end module FV_UTILS
