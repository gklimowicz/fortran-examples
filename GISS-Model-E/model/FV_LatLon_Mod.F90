#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort_core(__LINE__,rc)

module FV_LatLon_Mod

  use ESMF_MOD
  use FV_UTILS

  implicit none
  private

  public :: createInternalRestart
  public :: gridCompInit
  public :: copy_fv_export_to_modele
  public :: ConvertUV_GISS2FV
  public :: accumulate_mass_fluxes

  ! private data of convenience
  Integer :: rc ! return code from ESMF

  contains

  Subroutine createInternalRestart(fv, istart, cf, clock)
    USE DOMAIN_DECOMP_1D, ONLY: GRID, GET, AM_I_ROOT
    Use GEOS_IOMod, only: GETFILE, Free_file, GEOS_VarWrite, Write_parallel
    USE RESOLUTION, only: IM, JM, LM, LS1
    Use MODEL_COM, only: sige, sig, Ptop, DT, PMTOP
    Use MODEL_COM, only: U, V, T, P, PSFMPT, Q
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

    call getDomainBounds(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1, j_strt_halo=j_0h, j_stop_halo=j_1h)

    ! 4) 3D fields velocities
    Allocate(U_d(I_0:I_1, J_0:J_1, LM))
    Allocate(V_d(I_0:I_1, J_0:J_1, LM))

    write(*,*)'Calling ComputeRestartVelocities()'
    Call ComputeRestartVelocities(unit, grid, U, V, U_d, V_d)

!!$      call set_zonal_flow(U_d, V_d, j_0, j_1)

    Call GEOS_VarWrite(unit, grid % ESMF_GRID, U_d(I_0:I_1,J_0:J_1,:))
    Call GEOS_VarWrite(unit, grid % ESMF_GRID, V_d(I_0:I_1,J_0:J_1,:))

    Deallocate(V_d)
    Deallocate(U_d)

    ! Compute potential temperature from modelE (1 mb -> 1 pa ref)
    Allocate(PT(I_0:I_1, J_0:J_1, LM))
    Call ConvertPotTemp_GISS2FV(VirtualTemp(T(I_0:I_1,J_0:J_1,:), Q(I_0:I_1,J_0:J_1,:)), PT)
    Call GEOS_VarWrite(unit, grid % ESMF_GRID, PT)
    Deallocate(PT)

    ! Compute PE, PKZ from modelE
    Allocate(PKZ(I_0:I_1, J_0:J_1, LM))
    Allocate(PE(I_0:I_1, J_0:J_1, LM+1))
    Call ComputePressureLevels(unit, grid, VirtualTemp(T, Q), P, SIG, SIGE, Ptop, KAPA, PE, PKZ )

    Call GEOS_VarWrite(unit, grid % ESMF_GRID, PE)
    Call GEOS_VarWrite(unit, grid % ESMF_GRID, PKZ)

    Deallocate(PE)
    Deallocate(PKZ)

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
      USE DOMAIN_DECOMP_1D, only: dist_grid, get
      USE RESOLUTION, only: IM, LM, LS1
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
      call getDomainBounds(grid, &
    &      i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1, &
    &      i_strt_halo=i_0h, i_stop_halo=i_1h, &
    &      j_strt_halo=j_0h, j_stop_halo=j_1h)

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

    Subroutine ComputeRestartVelocities(unit, grid, U_b, V_b, U_d, V_d)
      use domain_decomp_1d, only: DIST_GRID, NORTH, HALO_UPDATE

      Integer, intent(in) :: unit
      Type (Dist_Grid) :: grid
      real*8, dimension(:,grid % j_strt_halo:,:) :: U_b, V_b

      real*8, intent(out) :: U_d(:,:,:), V_d(:,:,:)

      Call HALO_UPDATE(grid, U_b, FROM=NORTH)
      Call HALO_UPDATE(grid, V_b, FROM=NORTH)
      Call Regrid_B_to_D(ReverseLevels(U_b), ReverseLevels(V_b), U_d, V_d)

    End Subroutine ComputeRestartVelocities

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

  End Subroutine createInternalRestart

  subroutine gridCompInit(fv, clock)
    Type(FV_CORE), intent(inout) :: fv
    Type(esmf_clock), intent(in) :: clock

    call ESMF_GridCompInitialize ( fv%gc, importState=fv%import, &
         exportState=fv%export, clock=clock, &
         & phase=ESMF_SINGLEPHASE, rc=rc )
    VERIFY_(rc)

    ! Initialize component (and import/export states)
    !  ----------------------------------------------
    call allocateFvExport3D ( fv % export,'U' )
    call allocateFvExport3D ( fv % export,'V' )
    call allocateFvExport3D ( fv % export,'TH' )
    call allocateFvExport3D ( fv % export,'PLE' )
    call allocateFvExport3D ( fv % export,'Q' )
    call allocateFvExport3D ( fv % export,'MFX' )
    call allocateFvExport3D ( fv % export,'MFY' )
    call allocateFvExport3D ( fv % export,'MFZ' )
  end subroutine gridCompInit

  Subroutine Copy_FV_export_to_modelE(fv)
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    Use Resolution, only: IM,JM,LM,LS1
    Use DYNAMICS,   Only: PU, PV, SD
    Use MODEL_COM,  Only: U, V, T, Q, P, PTOP,PSFMPT
    USE DOMAIN_DECOMP_1D, only: grid, GET
    USE GEOM

    Type (FV_CORE) :: fv
    real*4, Dimension(:,:,:), Pointer :: U_a, V_a, T_fv, PLE

    Integer :: unit

    Integer :: i,j,k
    Integer :: i_0, i_1, j_0, j_1

    call getDomainBounds(grid, i_strt=i_0, i_stop=i_1, j_strt=j_0, j_stop=j_1)

    ! First compute updated values for modelE.  Then capture the
    ! new state in fv % *_old for computing tendencies.
    ! ----------------------------------------------------------

    ! Velocity field
    !---------------
    call ESMFL_StateGetPointerToData ( fv % export,U_a,'U',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % export,V_a,'V',rc=rc)
    VERIFY_(rc)
    call Regrid_A_to_B(ReverseLevels(U_a), ReverseLevels(V_a), U(I_0:I_1,J_0:J_1,:), V(I_0:I_1,J_0:J_1,:))
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


  !------------------------------------------------------------------------
  ! This following routine interpolates U, V from the Arakawa A grid
  ! to the Arakawa B grid.  B-grid velocities correspond to GISS modelE,
  ! while A-grid velocities are used for import/export of the FV component.
  !------------------------------------------------------------------------
  subroutine Regrid_A_to_B(U_a, V_a, U_b, V_b)
    Use Resolution, only : IM, LM
    Use Domain_decomp_1d, only : grid, get, SOUTH, HALO_UPDATE
    Real*4, intent(in), Dimension(:,grid % J_STRT:,:) :: U_a, V_a
    Real*8, intent(out),Dimension(:,grid % J_STRT:,:) :: U_b, V_b

    Real*8, allocatable, dimension(:,:,:) :: Ua_halo, Va_halo
    Integer :: i,j,k,ip1
    integer :: I_0, I_1, j_0stgr, j_1stgr,J_0h,J_1h,J_0,J_1

    call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, J_STRT_STGR=j_0stgr, J_STOP_STGR=j_1stgr, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H, &
         & J_STRT=J_0, J_STOP=J_1)
    Allocate(Ua_halo(IM,J_0h:J_1h,LM), Va_halo(IM,J_0h:J_1h,LM))

    Ua_halo(I_0:I_1,J_0:J_1,:) = U_a
    Va_halo(I_0:I_1,J_0:J_1,:) = V_a

    Call Halo_Update(grid, Ua_halo,FROM=SOUTH)
    Call Halo_Update(grid, Va_halo,FROM=SOUTH)

    Do k = 1, LM
       Do j = j_0STGR, j_1STGR
          i = IM
          Do ip1 = 1, IM

             u_b(i,j,k) = (Ua_halo(ip1,j-1,k) + Ua_halo(i,j-1,k) + Ua_halo(ip1,j,k) + Ua_halo(i,j,k))/4
             v_b(i,j,k) = (Va_halo(ip1,j-1,k) + Va_halo(i,j-1,k) + Va_halo(ip1,j,k) + Va_halo(i,j,k))/4
             i = ip1

          End do
       end do
    end do

    Deallocate(Ua_halo, Va_halo)

  end subroutine Regrid_A_to_B


  !--------------------------------------------------------------------
  ! This following routine interpolates U, V from the Arakawa B grid
  ! to the Arakawa D grid.  B-grid velocities correspond to GISS modelE,
  ! while D-grid velocities are used for the initial conditions of FV
  !--------------------------------------------------------------------
  Subroutine regrid_B_to_D(u_b, v_b, u_d, v_d)
    USE RESOLUTION, only: IM, JM, LM
    Use Domain_decomp_1d, only: grid, get
    Implicit None
    Real*8, intent(in) :: U_b(:,grid % j_strt_halo:,:)
    Real*8, intent(in) :: V_b(:,grid % j_strt_halo:,:)
    Real*8, intent(out) :: U_d(:,grid % j_strt:,:)
    Real*8, intent(out) :: V_d(:,grid % j_strt:,:)

    Integer :: i, j, k, im1
    Integer :: j_0, j_1, j_0s, j_1s
    Logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

    call getDomainBounds(grid, j_strt=j_0, j_stop=j_1, j_strt_skp=j_0s, j_stop_skp=j_1s, &
         & HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE)

    If (HAVE_SOUTH_POLE) Then
       ! 1st interpolate to 'A' grid
       Call FixPole(U_b(:,2,:), V_b(:,2,:), U_d(:,1,:))
     ! V_d(:,1,:) = (V_b(:,2,:) + CSHIFT(V_b(:,2,:),1,1))/2
       U_d(:,1,:)=0 ! not used (but needs legal value)
    End If

    If (HAVE_NORTH_POLE) Then
       ! 1st interpolate to 'A' grid
       Call FixPole(U_b(:,JM,:), V_b(:,JM,:), U_d(:,JM,:))
     ! V_d(:,JM,:) = (V_b(:,JM,:) + CSHIFT(V_b(:,JM,:),1,1))/2
    End If

! V-grid
    Do k = 1, LM
       Do j = j_0s, j_1s
          Do i = 1, IM
             V_d(i,j,k) = (V_b(i,j,k) + V_b(i,j+1,k)) /2
          End Do
       !  V_d(:,j,:) = (V_b(:,j,:) + CSHIFT(V_b(:,j,:),1,1))/2
       End Do
    End Do
! U-grid
    Do k = 1, LM
       Do j = j_0s, j_1
          im1 = IM
          Do i = 1, IM
             U_d(i,j,k) = (U_b(i,j,k) + U_b(im1,j,k)) /2
             im1 = i
          End Do
       !  U_d(:,j,:) = (U_b(:,j,:) + CSHIFT(U_b(:,j,:),1,1))/2
       End Do
    End Do

  Contains

    Subroutine FixPole(Ub, Vb, Ud)
      USE GEOM, only : SINIP, COSIP ! trig of lon and stgr lon
      Real*8, intent(in) :: Ub(IM, LM), Vb(IM,LM)
      Real*8, intent(out) :: Ud(IM, LM)

      Real*8 :: US, VS
      Integer :: k

      do k = 1, LM

         us = Sum( -ub(:,k) * SINIP - vb(:,k) * COSIP) / IM
         vs = Sum(  ub(:,k) * COSIP - vb(:,k) * SINIP) / IM

         ud(:,k) = (ub(:,k) + (-us * SINIP + vs * COSIP))/2

      end do

    End Subroutine FixPole
  End Subroutine regrid_B_to_D

  subroutine set_zonal_flow(U_d, V_d, j0, j1)
    use GEOM, only: cosp
    use domain_decomp_1d, only: grid
    Real*8, intent(out) :: U_d(:,grid % j_strt:,:)
    Real*8, intent(out) :: V_d(:,grid % j_strt:,:)
    integer :: j0, j1
    integer :: j

    do j = j0, j1
       U_d(:,j,:) = cosp(j)
       V_d(:,j,:) = 0
    end do
  end subroutine set_zonal_flow

  !------------------------------------------------------------------------
  ! This following routine interpolates U, V from the Arakawa B grid
  ! to the Arakawa A grid.  B-grid velocities correspond to GISS modelE,
  ! while A-grid velocities are used for import/export of the FV component.
  !------------------------------------------------------------------------
  subroutine ConvertUV_GISS2FV(U_b, V_b, U_a, V_a)
    Use Resolution, only : IM, JM, LM
    Use Domain_decomp_1d, only : grid, get, NORTH, HALO_UPDATE
    Real*8, intent(in), Dimension(:,grid % J_STRT:,:) :: U_b, V_b
    Real*4, intent(out), Dimension(:,grid % J_STRT:,:) :: U_a, V_a

    Real*8, allocatable, dimension(:,:,:) :: Ub_halo, Vb_halo

    Integer :: i,j,k,im1
    integer :: I_0, I_1, j_0, j_1, j_0s, j_1s, j_0h, j_1h
    logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

    call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1, J_STRT_SKP=J_0s, J_STOP_SKP=J_1S, &
         & HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE, &
         & J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

    Allocate(Ub_halo(IM,J_0h:J_1h,LM), Vb_halo(IM,J_0h:J_1h,LM))

    Ub_halo(I_0:I_1,J_0:J_1,:) = U_b
    Vb_halo(I_0:I_1,J_0:J_1,:) = V_b

    Call Halo_Update(grid, Ub_halo,FROM=NORTH)
    Call Halo_Update(grid, Vb_halo,FROM=NORTH)

    Do k = 1, LM
       Do j = j_0s,j_1s
          im1 = IM
          Do i = 1, IM

             u_a(i,j,k) = (ub_halo(im1,j,k) + ub_halo(i,j,k) + ub_halo(im1,j+1,k) + ub_halo(i,j+1,k))/4
             v_a(i,j,k) = (vb_halo(im1,j,k) + vb_halo(i,j,k) + vb_halo(im1,j+1,k) + vb_halo(i,j+1,k))/4
             im1 = i

          End do
       end do
    end do

    deallocate(ub_halo, Vb_halo)

    ! Polar conditions are a bit more complicated
    ! First determine an "absolute" U, V at the pole, then use sin/cos to
    ! map to the individual longitudes.
    If (HAVE_SOUTH_POLE) then
       Call FixPole(U_b(:,2,:), V_b(:,2,:), U_a(:,1,:), V_a(:,1,:))
    End If

    If (HAVE_NORTH_POLE) then
       Call FixPole(U_b(:,JM,:), V_b(:,JM,:), U_a(:,JM,:), V_a(:,JM,:))
    End If
  Contains
    ! This routine works for both poles.
    ! The sign is correct for the NP, but (-1) appears as a factor twice
    ! in the derivation for the south pole.
    !----------------------------------------
    Subroutine FixPole(Ub, Vb, Ua, Va)
      USE GEOM, only : SINIV, COSIV, SINIP, COSIP ! trig of lon and stgr lon
      Real*8, intent(in) :: Ub(IM, LM), Vb(IM,LM)
      Real*4, intent(out) :: Ua(IM, LM), Va(IM,LM)

      Real*8 :: us, vs
      Integer :: k

      do k = 1, LM

         us = Sum( -ub(:,k) * SINIP - vb(:,k) * COSIP) / IM
         vs = Sum(  ub(:,k) * COSIP - vb(:,k) * SINIP) / IM

         ua(:,k) = -us * SINIV + vs * COSIV
         va(:,k) = -us * COSIV - vs * SINIV

      end do

    End Subroutine FixPole
  end subroutine ConvertUV_GISS2FV

  subroutine accumulate_mass_fluxes(fv)
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    Use Resolution, only: IM,JM,LM,LS1
    USE GEOM, ONLY: AXYP
    USE GEOM, ONLY: DYP, DXV
    USE DYNAMICS, ONLY: MUs,MVs,MWs
    USE DYNAMICS, ONLY: PU,PV,CONV,SD,PIT
    USE MODEL_COM, only: DTsrc,DT,DSIG
    USE DOMAIN_DECOMP_1D, only: get, grid

!   Use Constant, only: radius,pi,grav
    implicit none
    type (FV_core) :: fv
    real*4, Dimension(:,:,:), Pointer :: PLE, mfx_X, mfx_Y, mfx_Z
    integer :: I_0, I_1, J_0, J_1
    integer :: J_0S, J_1S
    integer :: J_0H, J_1H
    integer :: i,im1,j,l,k
    integer :: rc
    logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
    real*8 :: DTLF, DTfac
    real*8 :: area,dlon,dlat,acap,rcap
    real*8 :: rX,rY,rZ
    real*8 :: sum1, sum2, mySum
    real*8, allocatable :: sine(:),cosp(:),cose(:)

    REAL*8 PVS,PVN
    REAL*8 PVSA(LM),PVNA(LM)
    REAL*8 TMPim1(IM)

    real*8, parameter :: grav=9.80
    real*8, parameter :: pi=3.14159265358979323846
    real*8, parameter :: radius= 6376000.00000000

    DTfac = DT

    call getDomainBounds(grid, &
    &    i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1, &
    &    J_STRT_SKP=J_0S, J_STOP_SKP=J_1S, &
    &    J_STRT_HALO=J_0H, J_STOP_HALO = J_1H, &
    &    HAVE_NORTH_POLE = HAVE_NORTH_POLE, &
    &    HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)

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
    call ESMFL_StateGetPointerToData ( fv % export,PLE,'PLE',rc=rc)
    VERIFY_(rc)


    PU(I_0:I_1,J_0:J_1,1:LM) = mfx_X
    PV(I_0:I_1,J_0:J_1,1:LM) = mfx_Y

#ifdef NO_MASS_FLUX
    mfx_X = 0
    mfx_Y = 0
    mfx_Z = 0
#endif
    PU = ReverseLevels(PU)/PRESSURE_UNIT_RATIO
    PV = ReverseLevels(PV)/PRESSURE_UNIT_RATIO
    mfx_Z = ReverseLevels(mfx_Z)/PRESSURE_UNIT_RATIO

    allocate ( sine(jm) )
    allocate ( cosp(jm) )
    allocate ( cose(jm) )
    dlon = 2.0*pi/im
    dlat = pi/(jm-1)
    do j=2,jm
         sine(j) = sin(-0.5*pi + ((j-1)-0.5)*(pi/(jm-1)))
    enddo
    cosp( 1) =  0.0
    cosp(jm) =  0.0
    do j=2,jm-1
         cosp(j) = (sine(j+1)-sine(j)) / dlat
    enddo
    do j=2,jm
       cose(j) = 0.5 * (cosp(j-1) + cosp(j))
    enddo
       cose(1) = cose(2)
! Adjust area scale factors between FV and GISS
    do l=1,lm
       do j=j_0,j_1
          do i=I_0,I_1
             PU(i,j,l) = PU(i,j,l)*DYP(j)/(radius*dlat)
          enddo
          do i=I_0,I_1
             PV(i,j,l) = PV(i,j,l)*DXV(j)/(cose(j)*radius*dlon)
          enddo
       enddo
    enddo
! Shift C-grid PU to Eastward orientation for GISS
    do l=1,lm
       do j=j_0,j_1
          im1 = I_0
          do i=I_1,I_0,-1
             TMPim1(i) = PU(im1,j,l)
             im1 = i
          enddo
          PU(:,j,l) = TMPim1
       enddo
    enddo
    deallocate ( sine )
    deallocate ( cosp )
    deallocate ( cose )

! Change Units of vertical mass fluxes
    do l=0,lm
       do j=j_0,j_1
          do i=I_0,I_1
             area = AXYP(i,j)
             mfx_Z(i-i_0+1,j-j_0+1,l) = grav*area*mfx_Z(i-i_0+1,j-j_0+1,l) ! convert to (mb m^2/s)
          enddo
       enddo
    enddo
    SD(I_0:I_1,J_0:J_1,1:LM-1) = (mfx_Z(:,:,1:LM-1)) ! SD only goes up to LM-1
    ! Surface Pressure tendency - vert integral of horizontal convergence
    PIT(I_0:I_1,J_0:J_1) = mfx_Z(:,:,0) + sum(SD(I_0:I_1,J_0:J_1,1:LM-1),3)

    ! Recopy into CONV to support prior usage
    CONV(I_0:I_1,J_0:J_1,1) = PIT(I_0:I_1,J_0:J_1)
    CONV(I_0:I_1,J_0:J_1,2:LM) = SD(I_0:I_1,J_0:J_1,1:LM-1)

    MUs(I_0:I_1,J_0:J_1,:) = MUs(I_0:I_1,J_0:J_1,:) + PU(I_0:I_1,J_0:J_1,:)*DTfac
    MVs(I_0:I_1,J_0:J_1,:) = MVs(I_0:I_1,J_0:J_1,:) + PV(I_0:I_1,J_0:J_1,:)*DTfac
    MWs(I_0:I_1,J_0:J_1,1:LM-1) = MWs(I_0:I_1,J_0:J_1,1:LM-1) + SD(I_0:I_1,J_0:J_1,1:LM-1)*DTfac
  end subroutine accumulate_mass_fluxes

end module FV_LatLon_Mod
