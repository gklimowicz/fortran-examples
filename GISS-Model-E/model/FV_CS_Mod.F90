#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort_core(__LINE__,rc)

module FV_CS_Mod

  use ESMF
  use FV_UTILS

  implicit none
  private

  public :: Create_Restart_File
  public :: GridSpecificInit
  public :: copy_fv_export_to_modele
  public :: ConvertUV_GISS2FV
  public :: accumulate_mass_fluxes

  ! private data of convenience
  Integer :: rc ! return code from ESMF

  contains

  Subroutine Create_Restart_File(fv, istart, cf, clock)
    USE DOMAIN_DECOMP_ATM, ONLY: GRID, getDomainBounds, AM_I_ROOT, ESMF_GRID_ATM
    Use MAPL_IOMod, only: GETFILE, Free_file, GEOS_VarWrite=>MAPL_VarWrite, Write_parallel
    USE RESOLUTION, only: IM, JM, LM, LS1, PMTOP, PTOP, PSFMPT
    Use DYNAMICS, only: sige,sig
    USE MODEL_COM, only: DT=>DTsrc
    Use ATM_COM, only: U, V, T, P, Q
    Use Constant, only: omega, radius, grav, rgas, kapa, deltx

    Type (FV_Core), Intent(InOut) :: fv
    integer, intent(in) :: istart
    Type (ESMF_Config), Intent(InOut) :: cf
    Type (ESMF_Clock),  Intent(In) :: clock

    Character(Len=ESMF_MAXSTR) :: rst_file
    Integer :: unit

    Integer, Parameter :: N_TRACERS = 0
    Integer :: I_0, I_1, j_0, j_1, j_0h, j_1h, L
    Logical :: exist

    Real*8 :: ak(size(sige)), bk(size(sige))

    real*8, allocatable, dimension(:,:,:) :: U_d
    real*8, allocatable, dimension(:,:,:) :: V_d
    real*8, allocatable, dimension(:,:,:) :: PE, PKZ, PT

    Call ESMF_ConfigGetAttribute(cf, value=rst_file, label='FVCORE_INTERNAL_RESTART_FILE:', &
         & default=FVCORE_INTERNAL_RESTART,rc=rc)

    if(istart .ge. extend_run) then
    ! Check to see if restart file exists
       inquire(file=FVCORE_INTERNAL_RESTART, EXIST=exist)
       if (exist) then
          if (AM_I_ROOT()) then
             print*,'Using checkpoint file: ',FVCORE_INTERNAL_RESTART
          end if
          return
       end if
       ! Uh oh
       call stop_model('fv part of restart file not found',255)
       return
    end if

    ! If we got to here, then this means then we'll have to create a restart file
    ! from scratch.
    unit = GetFile(rst_file, form="unformatted", rc=rc)
    VERIFY_(rc)

    ! 1) Start date
    Call write_start_date(clock, unit)

    ! 2) Grid size
    Call WRITE_PARALLEL( (/ IM, JM, LM, LM+1-LS1, N_TRACERS /), unit )

    ! 3) Pressure coordinates
    ! Keep in mind that L is reversed between these two models

    Call Compute_ak_bk(ak, bk, sige, Ptop, PSFMPT, unit)
    Call WRITE_PARALLEL( ak, unit)
    Call WRITE_PARALLEL( bk, unit)

    Call GetDomainBounds(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1, j_strt_halo=j_0h, j_stop_halo=j_1h)

    ! 4) 3D fields velocities
    Allocate(U_d(I_0:I_1, J_0:J_1, LM))
    Allocate(V_d(I_0:I_1, J_0:J_1, LM))

    U_d = 0.d0
    V_d = 0.d0

!!$      call set_zonal_flow(U_d, V_d, j_0, j_1)

    ! Compute potential temperature from modelE (1 mb -> 1 pa ref)
    Allocate(PT(I_0:I_1, J_0:J_1, LM))
    Call ConvertPotTemp_GISS2FV(VirtualTemp(T(I_0:I_1,J_0:J_1,:), Q(I_0:I_1,J_0:J_1,:)), PT)

    ! Compute PE, PKZ from modelE
    Allocate(PKZ(I_0:I_1, J_0:J_1, LM))
    Allocate(PE(I_0:I_1, J_0:J_1, LM+1))
    Call ComputePressureLevels(unit, grid, VirtualTemp(T, Q), P, SIG, SIGE, Ptop, KAPA, PE, PKZ )

    ! ESMF5 seems to register variables in alphabetical order, so we write them that way.
    ! note AK and BK have already been written.
    Call GEOS_VarWrite(unit, ESMF_GRID_ATM, PE)
    Call GEOS_VarWrite(unit, ESMF_GRID_ATM, PKZ)
    Call GEOS_VarWrite(unit, ESMF_GRID_ATM, PT)
    Call GEOS_VarWrite(unit, ESMF_GRID_ATM, U_d(I_0:I_1,J_0:J_1,:))
    Call GEOS_VarWrite(unit, ESMF_GRID_ATM, V_d(I_0:I_1,J_0:J_1,:))

    Deallocate(PE)
    Deallocate(PKZ)
    Deallocate(V_d)
    Deallocate(U_d)
    Deallocate(PT)


    Call Free_File(unit)

  CONTAINS

    ! Computes virtual pot. temp. from pot. temp. and specific humidity
    !------------------------------------------------------------------
    Function VirtualTemp(T, Q) Result(T_virt)
      Use Constant, only: DELTX
      Real*8, Intent(In) :: T(:,:,:)
      Real*8, Intent(In) :: Q(:,:,:)
      Real*8             :: T_virt(size(T,1), size(T,2), size(T,3))

      T_virt = T * (1 + deltx * Q)

    End Function VirtualTemp

    Subroutine ComputePressureLevels(unit, grid, T_virt, P, sig, sige, ptop, kapa, PE, PKZ)
      USE DOMAIN_DECOMP_ATM, only: dist_grid, getDomainBounds
      USE RESOLUTION, only: IM, LM
      Integer, intent(in) :: unit
      type (dist_grid) :: grid
      real*8, dimension(grid % i_strt_halo:,grid % j_strt_halo:,:) :: T_virt
      real*8, dimension(grid % i_strt_halo:,grid % j_strt_halo:) :: P
      real*8, dimension(grid % i_strt:grid % i_stop,grid % j_strt:grid % j_stop,LM) :: PKZ
      real*8, dimension(grid % i_strt:grid % i_stop,grid % j_strt:grid % j_stop,LM+1) :: PE
      real*8 :: sig(:), sige(:)
      real*8 :: ptop, kapa

      Integer :: I_0, I_1, i_0h, i_1h, j_0, j_1, j_0h, j_1h
      Integer :: I,J,L,L_fv
      Real*8, Allocatable :: PK(:,:,:), PELN(:,:,:), PE_trans(:,:,:)

      !    Request local bounds from modelE grid.
      Call GetDomainBounds(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1, &
            & i_strt_halo=i_0h, i_stop_halo=i_1h, j_strt_halo=j_0h, j_stop_halo=j_1h)

      PE = -99999
      PKZ = -99999

      Call ConvertPressure_GISS2FV(EdgePressure_GISS(), PE(i_0:i_1,j_0:j_1,:))

      Allocate(pk(i_0h:i_1h,j_0h:j_1h,LM+1))
      PK(I_0:I_1,J_0:J_1,:) = PE**KAPA
      do l=1,LM
        do j=j_0,j_1
          do i=I_0,I_1
            if (PE(i,j,l+1)-PE(i,j,l) /= 0.0) then
               PKZ(i,j,l) = ( PK(i,j,l+1)-PK(i,j,l) ) / &
                            ( KAPA*log( PE(i,j,l+1)/PE(i,j,l) ) )
            endif
          enddo
        enddo
      enddo
      deallocate(pk)

    End Subroutine ComputePressureLevels

    Subroutine Compute_ak_bk(ak, bk, sige, Ptop, PSFMPT, unit)
      USE RESOLUTION, only: LM, LS1
      Real*8 :: sige(:)
      Real*8 :: Ptop, PSFMPT
      Integer :: unit

      Real*8, intent(out) :: ak(size(sige)), bk(size(sige))
      Integer :: L, L_fv

      Do L = 1, LM+1
         L_fv = LM+2-L
         Select Case(L)
         Case (:LS1-1)
            ak(L_fv) = Ptop*(1-sige(L)) * PRESSURE_UNIT_RATIO ! convert from hPa
            bk(L_fv) = sige(L)
         Case (LS1:)
            ak(L_fv)   = (sige(L)*PSFMPT + Ptop) * PRESSURE_UNIT_RATIO! convert from hPa
            bk(L_fv)   = 0
         End Select
      End Do

    End Subroutine Compute_ak_bk


    Subroutine write_start_date(clock, unit)
      Type (ESMF_Clock), Intent(In) :: clock
      Integer          , Intent(In) :: unit

      Type (ESMF_Time)  :: currentTime
      Integer :: int_pack(6), YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

      call ESMF_ClockGet(clock, currTime=currentTime, rc=rc)

      call ESMF_TimeGet(currentTime, YY=year, MM= month, &
           DD=day, H=hour, M=minute, &
           S=second, rc=rc)

      INT_PACK(1) = YEAR
      INT_PACK(2) = MONTH
      INT_PACK(3) = DAY
      INT_PACK(4) = HOUR
      INT_PACK(5) = MINUTE
      INT_PACK(6) = SECOND
      Call WRITE_PARALLEL(INT_PACK(1:6), unit=UNIT)
    End Subroutine write_start_date

  End Subroutine Create_Restart_File

  subroutine GridSpecificInit(fv, clock)
    use FV_StateMod, only:FV_RESET_CONSTANTS
    use FV_Control_Mod, only : z_tracer
    use CONSTANT, only: pi, omega, sha, radius, rvap, grav, lhe, rgas, kapa

    Type(FV_CORE), intent(inout) :: fv
    type (esmf_clock) :: clock

    call ESMF_GridCompInitialize ( fv % gc, importState=fv % import, exportState=fv % export, clock=clock, &
!         & phase=ESMF_SINGLEPHASE, rc=rc )  ! todo: figure out what happened to ESMF_SINGLEPHASE
         & phase=1, rc=rc )
    VERIFY_(rc)

! z_tracer=.false. is necessary to get correct mass flux exports.
! Eventually, fix the FV code that is used when z_tracer=.true. 
    z_tracer = .false.

    call  FV_RESET_CONSTANTS( FV_PI=pi, &
                            & FV_OMEGA=omega ,&
                            & FV_CP=rgas/kapa ,&
                            & FV_RADIUS=radius ,&
                            & FV_RGAS=rgas ,&
                            & FV_RVAP=rvap ,&
                            & FV_KAPPA=kapa ,&
                            & FV_GRAV=grav ,&
                            & FV_HLV=lhe ,&
                            & FV_ZVIR=rvap/rgas-1  )

    ! Initialize component (and import/export states)
    !  ----------------------------------------------
    call allocateFvExport3D ( fv % export,'U_DGRID' )
    call allocateFvExport3D ( fv % export,'V_DGRID' )
    call allocateFvExport3D ( fv % export,'TH' )
    call allocateFvExport3D ( fv % export,'PLE' )
    call allocateFvExport3D ( fv % export,'Q' )
    call allocateFvExport3D ( fv % export,'MFX' )
    call allocateFvExport3D ( fv % export,'MFY' )
    call allocateFvExport3D ( fv % export,'MFZ' )
  end subroutine GridSpecificInit

  Subroutine Copy_FV_export_to_modelE(fv)
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    Use Resolution, only: IM,JM,LM,LS1,Ptop
    Use ATM_COM, only: U, V, T, P, Q
    USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
    USE GEOM

    Type (FV_CORE) :: fv
    real*4, Dimension(:,:,:), Pointer :: T_fv, PLE, U_d, V_d

    Integer :: unit

    Integer :: i,j,k
    Integer :: i_0, i_1, j_0, j_1

    Call GetDomainBounds(grid, i_strt=i_0, i_stop=i_1, j_strt=j_0, j_stop=j_1)

    ! First compute updated values for modelE.  Then capture the
    ! new state in fv % *_old for computing tendencies.
    ! ----------------------------------------------------------
    call ESMFL_StateGetPointerToData ( fv % export,U_d,'U_DGRID',rc=rc)
      VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % export,V_d,'V_DGRID',rc=rc)
      VERIFY_(rc)
    U(I_0:I_1,J_0:J_1,:) = ReverseLevels(U_d)
    V(I_0:I_1,J_0:J_1,:) = ReverseLevels(V_d)
    fv % U_old = U(I_0:I_1,J_0:J_1,:)
    fv % V_old = V(I_0:I_1,J_0:J_1,:)

    ! Potential temperature (save dry Temperature for computing tendencies)
    !----------------------------------------------------------------------
    call ESMFL_StateGetPointerToData ( fv % export,T_fv,'TH',rc=rc)
    VERIFY_(rc)
    call ConvertPotTemp_FV2GISS(T_fv, T)

    ! Use edge pressure export to update surface pressure
    !----------------------------------------------------------------------
    call ESMFL_StateGetPointerToData ( fv % export,PLE,'PLE',rc=rc)
    VERIFY_(rc)

    call ConvertPressure_FV2GISS(PLE, fv % PE_old)
    ! Just need surface pressure - Ptop
    P(I_0:I_1,J_0:J_1) = fv % PE_old(:,:,1) - Ptop
    CALL CALC_AMPK(LS1-1)

    ! Preserve state information for later computation of tendencies.
    fv % dT_old = DryTemp_GISS()
    fv % dPT_old = DeltPressure_DryTemp_GISS()

#if defined(USE_FV_Q)
    Q(I_0:I_1,j_0:j_1,:) = ReverseLevels(fv % Q)
#endif

  End Subroutine Copy_FV_export_to_modelE


  subroutine ConvertUV_GISS2FV(U_orig, V_orig, U_d, V_d)
! Cubed-sphere modelE works with native-grid winds, so this routine is trivial.
    Use Resolution, only : LM
    Use Domain_decomp_atm, only : grid, getDomainBounds
    Real*8, intent(in), Dimension(grid % I_STRT:,grid % J_STRT:,:) :: U_orig, V_orig
    Real*4, intent(out), Dimension(grid % I_STRT:,grid % J_STRT:,:) :: U_d, V_d

    Integer :: i,j,k
    integer :: I_0, I_1, j_0, j_1

    Call GetDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

    Do k = 1, LM
       Do j = j_0,j_1
          Do i = i_0,i_1
             U_d(i,j,k) = U_orig(i,j,k)
             V_d(i,j,k) = V_orig(i,j,k)
          End do
       end do
    end do
  end subroutine ConvertUV_GISS2FV

  subroutine accumulate_mass_fluxes(fv)
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    Use Resolution, only: LM
    USE GEOM, ONLY: AXYP
    USE ATM_COM, ONLY: MUs,MVs,MWs
    USE DYNAMICS, ONLY: PU,PV,CONV,SD,PIT
    Use MODEL_COM, only: DT=>DTsrc
    USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
    Use Constant, only: grav
    implicit none
    type (FV_core) :: fv
    real*4, Dimension(:,:,:), Pointer :: mfx_X, mfx_Y, mfx_Z
    integer :: I_0, I_1, J_0, J_1
    integer :: i,j,l,k
    integer :: rc
    real*8 :: DTfac

    DTfac = DT

    Call GetDomainBounds(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1)

    ! Horizontal and Vertical mass fluxes
    !---------------
    !     \item {\tt MFX}:       Mass-Weighted U-Wind on C-Grid (Pa m^2/s)
    !     \item {\tt MFY}:       Mass-Weighted V-wind on C-Grid (Pa m^2/s)
    !     \item {\tt MFZ}:       Vertical mass flux (kg/(m^2*s))

    call ESMFL_StateGetPointerToData ( fv % export,mfx_X,'MFX',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % export,mfx_Y,'MFY',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % export,mfx_Z,'MFZ',rc=rc)
    VERIFY_(rc)

#ifdef NO_MASS_FLUX
    mfx_X = 0
    mfx_Y = 0
    mfx_Z = 0
#endif

    PU(I_0:I_1,J_0:J_1,1:LM) = ReverseLevels(mfx_X)
    PV(I_0:I_1,J_0:J_1,1:LM) = ReverseLevels(mfx_Y)
    mfx_Z = ReverseLevels(mfx_Z)
    SD(I_0:I_1,J_0:J_1,1:LM-1) = mfx_Z(:,:,1:LM-1)

!
! Add missing factor of grav to FV exports, change units
!
    pu = pu*grav/PRESSURE_UNIT_RATIO
    pv = pv*grav/PRESSURE_UNIT_RATIO


! Change Units of vertical mass fluxes to mb m^2/s
    do l=1,lm-1
       do j=j_0,j_1
          do i=I_0,I_1
             sd(i,j,l) = grav*grav*AXYP(i,j)*sd(i,j,l)/PRESSURE_UNIT_RATIO
          enddo
       enddo
    enddo

    ! Surface Pressure tendency - vert integral of horizontal convergence
    PIT(I_0:I_1,J_0:J_1) = mfx_Z(:,:,0) + sum(SD(I_0:I_1,J_0:J_1,1:LM-1),3)

    ! Recopy into CONV to support prior usage
    CONV(I_0:I_1,J_0:J_1,1) = PIT(I_0:I_1,J_0:J_1)
    CONV(I_0:I_1,J_0:J_1,2:LM) = SD(I_0:I_1,J_0:J_1,1:LM-1)

    MUs(I_0:I_1,J_0:J_1,:) = MUs(I_0:I_1,J_0:J_1,:) + PU(I_0:I_1,J_0:J_1,:)*DTfac
    MVs(I_0:I_1,J_0:J_1,:) = MVs(I_0:I_1,J_0:J_1,:) + PV(I_0:I_1,J_0:J_1,:)*DTfac
    MWs(I_0:I_1,J_0:J_1,1:LM-1) = MWs(I_0:I_1,J_0:J_1,1:LM-1) + SD(I_0:I_1,J_0:J_1,1:LM-1)*DTfac

  end subroutine accumulate_mass_fluxes

end module FV_CS_Mod
