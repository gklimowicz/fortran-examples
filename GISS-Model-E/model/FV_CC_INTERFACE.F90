#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort_core(__LINE__,rc)
#define RETURN_(status) If (Present(rc)) rc=status; return

!--------------------------------------------------------------------------------
! The following module provides an interface for modelE to use the ESMF-wrapped 
! version of the FV dynamical core from GEOS-5.  Many elements closely resemble
! an ESMF coupler component, which is to be expected given the intended role for
! this layer.
!--------------------------------------------------------------------------------

module FV_INTERFACE_MOD

  use ESMF_mod
  use FV_UTILS

#ifdef CUBED_SPHERE
  use FV_CS_Mod
#else
  use FV_LatLon_Mod
#endif

  implicit none
  save
  private

  ! except for

  public :: Initialize          ! Uses modelE data to initialize the FV gridded component
  public :: Run                 ! Execute the FV method to integrate the core forward in time
  public :: Checkpoint      
  public :: Finalize            ! Uses modelE data to finalize the FV gridded component

  public :: FV_CORE_WRAPPER     ! Derived type to encapsulate FV + modelE data

  public :: Compute_Tendencies  ! Temporarily a separate method, but will move within Run() soon.

  Interface Initialize
    module procedure Initialize_fv
  End Interface

  Interface Run
    module procedure run_fv
  End Interface

  Interface Compute_Tendencies
    module procedure compute_tendencies_external
    module procedure compute_tendencies_internal
  End Interface

  ! private data of convenience
  Integer :: rc ! return code from ESMF

  Type FV_CORE_WRAPPER
    PRIVATE
    type(FV_CORE) :: fv
  END Type FV_CORE_WRAPPER

  Type (FV_CORE_WRAPPER), public :: fvstate

  Character(Len=*), Parameter :: config_file = 'fv_config.rc'

contains

!-------------------------------------------------------------------------------
  subroutine Initialize_fv(fv_wrapper, istart, kdisk)
!-------------------------------------------------------------------------------
    use atm_com, only : clock=>atmclock
    Use Domain_decomp_atm, only: grid, AM_I_ROOT
    type (fv_core_wrapper),    intent(inout) :: fv_wrapper
    integer,           intent(in) :: istart,kdisk
    type (ESMF_config)            :: cf
    character(len=1) :: suffix

! if warm start, copy FV restart files to expected names
    if(AM_I_ROOT() .and. istart.ge.8) then
      if(istart.lt.10) then
        call system('cp AICfv  dyncore_internal_restart')
        call system('cp AICdfv tendencies_checkpoint')
      else
        write(suffix,'(i1)') KDISK
        call system('cp  fv.'// suffix // ' dyncore_internal_restart')
        call system('cp dfv.'// suffix // ' tendencies_checkpoint')
      endif
    endif

    call gridCompSetup(fv_wrapper%fv, grid%esmf_grid, cf, config_file)

    ! The FV components requires its own restart file for managing
    ! its internal state.  We check to see if the file already exists, and if not
    ! create one based upon modelE's internal state.
    call createInternalRestart(fv_wrapper%fv, istart, cf, clock)

    call gridCompInit(fv_wrapper%fv, clock)

    Call allocateTendencyStorage(fv_wrapper%fv, istart)

    call addQfieldToFVImport(fv_wrapper%fv)

  end subroutine Initialize_fv

!-------------------------------------------------------------------------------
  subroutine run_fv(fv_wrapper)
!-------------------------------------------------------------------------------
!@vers 2013/03/25
    use atm_com, only : clock=>atmclock
    USE RESOLUTION, Only : IM, JM, LM
    USE ATM_COM, Only : U, V, T, P, ZATMO
    USE MODEL_COM, only : DT => DTSRC
    USE SOMTQ_COM, only: QMOM, TMOM, MZ
#ifdef FVCUBED_SKIPPED_THIS
    USE ATMDYN, only:  COMPUTE_MASS_FLUX_DIAGS
#endif
    USE DYNAMICS, ONLY: PU, PV, SD, CONV
    USE ATM_COM, ONLY: MMA, PHI, GZ, MUs, MVs, MWs

    Type (FV_CORE_WRAPPER)    :: fv_wrapper
    type (ESMF_TimeInterval) :: timeInterval

    integer :: L,istep, NS, NIdyn_fv
    integer :: addIncsPhase

! Phase number used to invoke the core's RunAddIncs routine changes
! between latlon and cubedsphere cores
#ifdef CUBED_SPHERE
    addIncsPhase = 1
#else
    addIncsPhase = 91
#endif

!@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
    Call CALC_AMP (P, MMA)

    Call Copy_modelE_to_FV_import(fv_wrapper%fv)

    call clear_accumulated_mass_fluxes()

    ! Run dycore
    NIdyn_fv = 1 !DTsrc / (DT)
    do istep = 1, NIdyn_fv

      call ESMF_GridCompRun(fv_wrapper%fv%gc, fv_wrapper%fv%import, &
        fv_wrapper%fv%export, clock, phase=addIncsPhase, rc=rc )

      call clearTendencies(fv_wrapper%fv)
      call ESMF_GridCompRun(fv_wrapper%fv%gc, fv_wrapper%fv%import, &
        fv_wrapper%fv%export, clock, rc=rc )

      call ESMF_TimeIntervalSet(timeInterval, s = nint(DT), rc=rc)
      call ESMF_ClockAdvance(clock, timeInterval, rc=rc)

      call accumulate_mass_fluxes(fv_wrapper%fv)
      call Copy_FV_export_to_modelE(fv_wrapper%fv) ! inside loop to accumulate MUs,MVs,MWs

      call reset_tmom
#if defined(USE_FV_Q)
      call reset_qmom
#endif

      phi = compute_phi(P, T, TMOM(MZ,:,:,:), ZATMO)
#ifdef FVCUBED_SKIPPED_THIS
      call compute_mass_flux_diags(phi, pu, pv, dt)
#endif

    end do

    call compute_cp_vvel(MUs,MVs,MWs,p)

    gz  = phi
    
  contains

  !---------------------------------------------------------------------------
    subroutine reset_tmom
  !---------------------------------------------------------------------------

      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : t,pmid,pedn
      USE DOMAIN_DECOMP_ATM, ONLY: grid
      USE SOMTQ_COM, only : mz,tmom,qmom
      implicit none

      integer :: i,j,l
      REAL*8 :: rdsig

      INTEGER :: I_0, I_1, J_1, J_0

!**** Extract useful local domain parameters from "grid"
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP

      tmom = 0 ! except for vertical slopes:
      DO J=J_0,J_1
        DO I=I_0,I_1
          RDSIG=(PMID(1,I,J)-PEDN(2,I,J))/(PMID(1,I,J)-PMID(2,I,J))
          TMOM(MZ,I,J,1)=(T(I,J,2)-T(I,J,1))*RDSIG
          DO L=2,LM-1
            RDSIG=(PMID(L,I,J)-PEDN(L+1,I,J))/(PMID(L-1,I,J)-PMID(L+1,I,J))
            TMOM(MZ,I,J,L)=(T(I,J,L+1)-T(I,J,L-1))*RDSIG
          END DO
          RDSIG=(PMID(LM,I,J)-PEDN(LM+1,I,J))/(PMID(LM-1,I,J)-PMID(LM,I,J))
          TMOM(MZ,I,J,LM)=(T(I,J,LM)-T(I,J,LM-1))*RDSIG
        END DO
      END DO

      return

    end subroutine reset_tmom

  !---------------------------------------------------------------------------
    subroutine reset_qmom
  !---------------------------------------------------------------------------

      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : q,pmid,pedn
      USE DOMAIN_DECOMP_ATM, ONLY: grid
      USE SOMTQ_COM, only : mz,tmom,qmom
      implicit none

      integer :: i,j,l
      REAL*8 :: rdsig

      INTEGER :: I_0, I_1, J_1, J_0

!**** Extract useful local domain parameters from "grid"
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP

      qmom = 0 ! except for vertical slopes:
      DO J=J_0,J_1
        DO I=I_0,I_1
          RDSIG=(PMID(1,I,J)-PEDN(2,I,J))/(PMID(1,I,J)-PMID(2,I,J))
          QMOM(MZ,I,J,1)=(Q(I,J,2)-Q(I,J,1))*RDSIG
          IF(Q(I,J,1)+QMOM(MZ,I,J,1).LT.0.) QMOM(MZ,I,J,1)=-Q(I,J,1)
          DO L=2,LM-1
            RDSIG=(PMID(L,I,J)-PEDN(L+1,I,J))/(PMID(L-1,I,J)-PMID(L+1,I,J))
            QMOM(MZ,I,J,L)=(Q(I,J,L+1)-Q(I,J,L-1))*RDSIG
            IF(Q(I,J,L)+QMOM(MZ,I,J,L).LT.0.) QMOM(MZ,I,J,L)=-Q(I,J,L)
          END DO
          RDSIG=(PMID(LM,I,J)-PEDN(LM+1,I,J))/(PMID(LM-1,I,J)-PMID(LM,I,J))
          QMOM(MZ,I,J,LM)=(Q(I,J,LM)-Q(I,J,LM-1))*RDSIG
          IF(Q(I,J,LM)+QMOM(MZ,I,J,LM).LT.0.) QMOM(MZ,I,J,LM)=-Q(I,J,LM)
        END DO
      END DO

      return

     end subroutine reset_qmom

  end subroutine run_fv

!-------------------------------------------------------------------------------
  subroutine Checkpoint(fv_wrapper, modelEfilename)
!-------------------------------------------------------------------------------
    use model_com, only : rsf_file_name,xlabel,lrunid
    use atm_com, only : clock=>atmclock
    Type (FV_Core_Wrapper) :: fv_wrapper
    character(len=*), intent(in) :: modelEfilename

    character(len=*), parameter :: SUFFIX_TEMPLATE = '.YYYYMMDD_HHMMz.bin'
    character(len=len(SUFFIX_TEMPLATE)) :: suffix
    Type (ESMF_Time)  :: currentTime
    integer :: year, month, day, hour, minute, second
    logical :: isFinalize

    character(len=80) :: fv_fname, fv_dfname
    character(len=1)  :: c1
    integer :: m

    if( any(rsf_file_name(:)==trim(modelEfilename)) ) then
      ! standard checkpoint file
      do m=1,size(rsf_file_name)
        if(rsf_file_name(m) == trim(modelEfilename)) exit
      enddo
      write(c1,'(i1)') m
      fv_fname='fv.'//c1
      fv_dfname='dfv.'//c1
    elseif( index(modelEfilename,'.rsf')==9 ) then
      ! end-of-month file of the form 1MONYEAR.rsfRUNNAME
      fv_fname  = modelEfilename(1:8)//'.fv'//XLABEL(1:LRUNID)
      fv_dfname = modelEfilename(1:8)//'.dfv'//XLABEL(1:LRUNID)
    else
      call stop_model('checkpoint_fv: unrecognized modelEfilename',255)
    endif

    ! dyncore_internal_restart.YYYYMMDD_HHMMz.bin
    call ESMF_ClockGet(clock, currTime=currentTime, rc=rc)
    call ESMF_TimeGet(currentTime, YY=year, MM= month, &
         DD=day, H=hour, M=minute, &
         S=second, rc=rc)
    write(suffix,'(".",i4.4,i2.2,i2.2,"_",i2.2,i2.2,"z.bin")'), &
         & year, month, day, hour, minute

    isFinalize = .false.
    call DumpState(fv_wrapper%fv, clock, fv_fname, fv_dfname, suffix, isFinalize)

  end subroutine Checkpoint

!-------------------------------------------------------------------------------
  subroutine Finalize(fv_wrapper, kdisk)
!-------------------------------------------------------------------------------
    use atm_com, only : clock=>atmclock
    Type (FV_Core_Wrapper) :: fv_wrapper
    integer, intent(in) :: kdisk

    character(len=0) :: suffix = ''
    logical :: isFinalize
    character(len=28) :: fv_fname, fv_dfname
    character(len=1)  :: c1

    write(c1,'(i1)') kdisk
    fv_fname='fv.'//c1
    fv_dfname='dfv.'//c1

    isFinalize = .true.
    call DumpState(fv_wrapper%fv, clock, fv_fname, fv_dfname, suffix, isFinalize)

    Deallocate(fv_wrapper%fv%U_old, fv_wrapper%fv%V_old, &
          &    fv_wrapper%fv%dPT_old, fv_wrapper%fv%dT_old, &
          &    fv_wrapper%fv%PE_old)

  end subroutine Finalize

!-------------------------------------------------------------------------------
  Subroutine allocateTendencyStorage(fv, istart)
!-------------------------------------------------------------------------------
    Use Domain_decomp_atm, only: GRID, getDomainBounds, AM_I_ROOT
    USE RESOLUTION, only: IM, LM, LS1
    USE MODEL_COM, only: DTsrc
    USE ATM_COM, only: U, V, T
    use FILEMANAGER
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
#ifdef CUBED_SPHERE
    use FV_StateMod, only: INTERP_AGRID_TO_DGRID
#endif

    Type (FV_Core) :: fv
    integer, intent(in) :: istart

    Integer :: I_0, I_1, J_0, J_1, J_0H, J_1H
    integer :: iunit
!local
      real*8,allocatable :: U_temp(:,:,:)
      real*8,allocatable :: V_temp(:,:,:)
      real*8,allocatable :: U_d(:,:,:)
      real*8,allocatable :: V_d(:,:,:)

    call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, &
         & J_STRT=J_0, J_STOP=J_1, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

    ! 1) Link/copy modelE data to import state
    call ESMFL_StateGetPointerToData ( fv%import,fv%dudt,'DUDT',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv%import,fv%dvdt,'DVDT',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv%import,fv%dtdt,'DTDT',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv%import,fv%dpedt,'DPEDT',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv%import,fv%phis, 'PHIS', rc=rc)
    VERIFY_(rc)

    ! 2) Allocate space for storing old values from which to compute tendencies
    Allocate(fv%U_old(I_0:I_1,J_0:J_1,LM), &
         &   fv%V_old(I_0:I_1,J_0:J_1,LM), &
         &   fv%dPT_old(I_0:I_1,J_0:J_1,LM), &
         &   fv%dT_old(I_0:I_1,J_0:J_1,LM), &
         &   fv%PE_old(I_0:I_1,J_0:J_1,LM+1))

    select case (istart)
    case (:initial_start)
      ! Do a cold start.  Set Old = Current.
#ifdef CUBED_SPHERE
      Allocate(U_d(I_0:I_1,J_0:J_1+1,LM), &
        V_d(I_0:I_1+1,J_0:J_1,LM))
      Allocate(U_temp(I_0:I_1,J_0:J_1,LM), &
        V_temp(I_0:I_1,J_0:J_1,LM))

      U_temp(I_0:I_1,J_0:J_1,:) = U(I_0:I_1,J_0:J_1,:)
      V_temp(I_0:I_1,J_0:J_1,:) = V(I_0:I_1,J_0:J_1,:)
      Call INTERP_AGRID_TO_DGRID(U_temp, V_temp, U_d, V_d)
      U(I_0:I_1,J_0:J_1,:) = U_d(I_0:I_1,J_0:J_1,:)
      V(I_0:I_1,J_0:J_1,:) = V_d(I_0:I_1,J_0:J_1,:)

      Deallocate(U_d, V_d)
      Deallocate(U_temp, V_temp)

      ! the tendency must be scaled by DT for use by the core's ADD_INCS routine
      fv%dudt=ReverseLevels(U(I_0:I_1,J_0:J_1,:))/DTsrc
      fv%dvdt=ReverseLevels(V(I_0:I_1,J_0:J_1,:))/DTsrc
#else
      fv%dudt=0
      fv%dvdt=0
#endif
      fv%dtdt=0
      fv%dpedt=0

      fv%U_old = U(I_0:I_1,J_0:J_1,:)
      fv%V_old = V(I_0:I_1,J_0:J_1,:)
      fv%dPT_old = DeltPressure_DryTemp_GISS()
      fv%dT_old = DryTemp_GISS()
      fv%PE_old   = EdgePressure_GISS()
    case (extend_run:)
      ! input couplings are in a file somewhere and already read in
      call openunit(TENDENCIES_FILE, iunit, qbin=.true.,qold=.true.)
      call readArr(iunit, fv%U_old)
      call readArr(iunit, fv%V_old)
      call readArr(iunit, fv%dPT_old)
      call readArr(iunit, fv%PE_old)
      call readArr(iunit, fv%dT_old)
      call closeunit(iunit)

      if (AM_I_ROOT()) then
        call system('rm ' // TENDENCIES_FILE ) ! clean up
      end if

      call compute_tendencies(fv)
      !!  case default
      !!     call stop_model('ISTART option not supported',istart)
    end select

  contains

    subroutine readArr(iunit, arr)
      use domain_decomp_atm, only: grid, dread8_parallel, getDomainBounds
      integer, intent(in) :: iunit
      real*8, intent(out) :: arr(:,:,:)
      real*8, allocatable :: padArr(:,:,:)
      integer :: I_0, I_1, j_0, j_1
      integer :: i_0h, i_1h
      integer :: j_0h, j_1h

      call getDomainBounds(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1, &
           & i_strt_halo = i_0h, i_stop_halo = i_1h, &
           & j_strt_halo = j_0h, j_stop_halo = j_1h)

      allocate(padArr(i_0h:i_1h,j_0h:j_1h, size(arr,3)))

      call dread8_parallel(grid, iunit, nameunit(iunit), padArr)
      arr(:,:,:) = padArr(I_0:I_1,j_0:j_1,:)

      deallocate(padArr)

    end subroutine readArr

  End Subroutine allocateTendencyStorage

!-------------------------------------------------------------------------------
  subroutine compute_tendencies_internal(fv)
!-------------------------------------------------------------------------------
    USE RESOLUTION, only: IM, LM
    USE ATM_COM, Only: U, V, T
    USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
    USE DYNAMICS, only: DUT, DVT

    Implicit None
    Type (FV_CORE) :: fv

    Integer :: I_0, I_1, J_0, J_1
    Integer :: J_0h, J_1h
    call getDomainBounds(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1)

    ! U, V
    DUT(I_0:I_1,J_0:J_1,:) = Tendency(U(I_0:I_1,J_0:J_1,:), &
      fv%U_old(I_0:I_1,J_0:J_1,:))
    DVT(I_0:I_1,J_0:J_1,:) = Tendency(V(I_0:I_1,J_0:J_1,:), &
      fv%V_old(I_0:I_1,J_0:J_1,:))
    call ConvertUV_GISS2FV(ReverseLevels(DUT(I_0:I_1,J_0:J_1,:)), &
      ReverseLevels(DVT(I_0:I_1,J_0:J_1,:)), &
      fv%dudt, fv%dvdt)

    ! delta pressure weighted Temperature
    fv % dtdt = ReverseLevels(DeltPressure_GISS() * &
      Tendency(DryTemp_GISS(), fv%dT_old)) * &
      (PRESSURE_UNIT_RATIO)

    ! Edge Pressure
    Call ConvertPressure_GISS2FV( Tendency(EdgePressure_GISS(), fv%PE_old), &
      fv%dpedt)

#ifdef NO_FORCING
    call clearTendencies(fv)
#endif

  end subroutine compute_tendencies_internal

  ! This routine exposes compute_tendencies to the main program, using the fv_wrapper
!-------------------------------------------------------------------------------
  subroutine compute_tendencies_external(fv_wrapper)
!-------------------------------------------------------------------------------
    Type (FV_CORE_WRAPPER) :: fv_wrapper

    call compute_tendencies(fv_wrapper%fv)
  end subroutine compute_tendencies_external

end module FV_INTERFACE_MOD
