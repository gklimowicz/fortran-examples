module GlobalSum_mod
  use MpiSupport_mod, only: ROOT_PROCESS
  use dist_grid_mod
  use GatherScatter_mod, only: gatherReal8
  implicit none
  private

!@var GLOBALSUM output a bit-reproducible global-hemisphere-zonal sum for an array
  public :: globalsum

  interface globalsum
    module procedure globalsum_int_reduce
    module procedure globalsum_j
    module procedure globalsum_ij
    module procedure globalsum_ijk
    module procedure globalsum_ijk_ik
    module procedure globalsum_other_ijk
    module procedure globalsum_other_ijk_ik ! not used?
    module procedure globalsum_jk
    module procedure globalsum_xxxj_xxx
    module procedure globalsum_xxxij_xxx
  end interface

!@var checksum output a bit-reproducible checksum for an array
  public :: checksum ! Communicate overlapping portions of subdomains
  public :: checksumJ! Communicate overlapping portions of subdomains
  public :: checksum_column ! (k, i, j)

  interface checksum
    module procedure checksum_1d
    module procedure checksum_2d
    module procedure checksum_3d
  end interface

  interface checksumj
    module procedure checksumj_2d
    module procedure checksumj_3d
  end interface

  interface checksum_column
    module procedure checksum_column_2d
    module procedure checksum_column_3d
    module procedure int_checksum_column_3d
    module procedure checksum_column_4d
  end interface

#ifdef USE_MPI
  include 'mpif.h'
#endif

contains

  !---------------------------
  subroutine globalsum_int_reduce(grd_dum, ivar, isum, all)
    type(dist_grid), intent(in) :: grd_dum
    integer, intent(in)  :: ivar
    integer, intent(out) :: isum
    logical, optional, intent(in) :: all

    logical :: all_
    integer :: ierr

    all_ = .false.
    if (present(all)) all_ = all

#ifdef USE_MPI
    if (all_) then
      call MPI_Allreduce(ivar, isum, 1, MPI_integer, MPI_SUM, &
           &        getMpiCommunicator(grd_dum), ierr)
    else
      call MPI_Reduce(ivar, isum, 1, MPI_integer, MPI_SUM, ROOT_PROCESS, &
           &        getMpiCommunicator(grd_dum), ierr)
    end if
#else
    isum = ivar
#endif

  end subroutine globalsum_int_reduce

  subroutine globalSum_J(grd_dum, arr, gsum, hsum, istag, iskip, all, jband)
    type (dist_grid),  intent(IN) :: grd_dum
    real*8,            intent(IN) :: arr(grd_dum%j_strt_halo:)
    real*8,            intent(OUT):: gsum
    real*8, optional,  intent(OUT):: hsum(2)
    integer,optional,  intent(IN) :: istag
    integer,optional,  intent(IN) :: iskip
    logical,optional,  intent(IN) :: all
    integer, optional, intent(IN) :: jband(2)
    
    integer :: i_0, i_1, j_0, j_1, IM, JM, ierr
    real*8  :: garr(grd_dum%jm_world)
    logical :: istag_, iskip_
    
    integer :: JSTRT, JSTOP
    ! now local
    
    call getDomainBounds(grd_dum, J_STRT=j_0, J_STOP=j_1)
    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD
    
    istag_ = .false.
    if (present(istag)) then
      if (istag == 1) istag_ = .true.
    end if
    
    iskip_ = .false.
    if (present(iskip)) then
      if (iskip == 1) iskip_ = .true.
    end if
    
    JSTRT=1
    JSTOP=JM
    
    if (present(jband)) then
      JSTRT=jband(1)
      JSTOP=jband(2)
    end if
    
#ifdef DEBUG_DECOMP
    if (size(arr) /= grd_dum%j_stop_halo-grd_dum%j_strt_halo+1) then
#ifdef USE_MPI
      call MPI_ABORT(ierr)
#else
      stop
#endif
    end if
#endif
    
    call gatherReal8(grd_dum, arr, garr, shape(arr), 1, .false.)
    
    if (AM_I_ROOT()) then
      if (istag_) then
        gsum = sum(garr(2:JM),1)
      elseif (iskip_) then
        gsum = sum(garr(2:JM-1),1)
      else
        gsum = sum(garr(JSTRT:JSTOP),1)
      endif
      if (present(hsum)) then
        if (istag_) then
          hsum(1)   = sum( garr(2     :JM/2),1   )
          hsum(2)   = sum( garr(2+JM/2:JM  ),1   )
          hsum(1)   = hsum(1) + 0.5*garr(1+JM/2)
          hsum(2)   = hsum(2) + 0.5*garr(1+JM/2)
        else
          hsum(1)   = Sum( garr(1     :JM/2),1   )
          hsum(2)   = sum( garr(1+JM/2:JM  ),1   )
        endif
      endif
    endif

#ifdef USE_MPI
    if (present(all)) then
      if (all) then
        call MPI_BCAST(gsum,1,MPI_DOUBLE_PRECISION,ROOT_PROCESS, &
             &           getMpiCommunicator(grd_dum), ierr)
        if (present(hsum)) call MPI_BCAST(hsum,2, MPI_DOUBLE_PRECISION, &
             & ROOT_PROCESS, getMpiCommunicator(grd_dum), ierr)
      end if
    end if
#endif
  end subroutine globalSum_J
  
  subroutine globalSum_IJ(grd_dum, arr, gsum, hsum, zsum, istag, all, iskip)
    type (dist_grid),  intent(IN) :: grd_dum
    real*8,            intent(IN) :: arr(grd_dum%i_strt_halo:, &
         &                                     grd_dum%j_strt_halo:)
    real*8,            intent(OUT):: gsum
    real*8, OPTIONAL,  INTENT(OUT):: hsum(2)
    real*8, OPTIONAL,  INTENT(OUT):: zsum(grd_dum%j_strt: &
         &                                      grd_dum%j_stop)
    integer,OPTIONAL,  INTENT(IN) :: istag
    integer,OPTIONAL,  INTENT(IN) :: iskip

    integer :: i_0, i_1, j_0, j_1, IM, JM, J_0STG, ierr, J_0S,J_1S
    real*8, allocatable  :: zon(:)
    real*8  :: garr(grd_dum%jm_world)
    logical :: istag_,iskip_
    logical,OPTIONAL,  INTENT(IN) :: all

    ! now local

    call getDomainBounds(grd_dum, j_strt=j_0, j_stop=j_1)
    allocate(zon(j_0:j_1))
    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    j_0STG  = grd_dum%j_strt_stgr
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD
    j_0S = grd_dum%j_strt_skp
    j_1S  = grd_dum%j_stop_skp

    istag_ = .false.
    if (present(istag)) Then
      if (istag == 1) istag_ = .true.
    end if

    iskip_ = .false.
    if (present(iskip)) Then
      if (iskip == 1) iskip_ = .true.
    end if

    if (istag_) Then
      zon(J_0STG:J_1)  = sum(arr(:,j_0STG:j_1),1)
    Else if (iskip_) Then
      zon(J_0S:J_1S)  = sum(arr(:,j_0S:j_1S),1)
    Else
      zon  = sum(arr(:,j_0:j_1),1)
    end if

    Call gatherReal8(grd_dum, zon, garr, shape(zon), 1, .false.)

    if (AM_I_ROOT()) then
      if (istag_) then
        gsum = sum(garr(2:JM),1)
      Else if(iskip_) then
        gsum = sum(garr(2:JM-1),1)
      Else
        gsum = sum(garr(1:JM),1)
      endif
      if (present(hsum)) then
        if (istag_) then
          hsum(1)   = Sum( garr(2     :JM/2),1   )
          hsum(2)   = Sum( garr(2+JM/2:JM  ),1   )
          hsum(1)   = hsum(1) + 0.5*garr(1+JM/2)
          hsum(2)   = hsum(2) + 0.5*garr(1+JM/2)
        Else
          hsum(1)   = Sum( garr(1     :JM/2),1   )
          hsum(2)   = Sum( garr(1+JM/2:JM  ),1   )
        endif
      endif
    endif
    if (present(zsum)) zsum = zon

#ifdef USE_MPI
    if (present(all)) Then
      if (all) Then
        Call MPI_BCAST(gsum,1,MPI_DOUBLE_PRECISION,ROOT_PROCESS, &
             &           getMpiCommunicator(grd_dum), ierr)
        if (present(hsum))     Call MPI_BCAST(hsum,2, &
             &           MPI_DOUBLE_PRECISION,ROOT_PROCESS, getMpiCommunicator(grd_dum), ierr)
      end if
    end if
#endif

  end subroutine globalSum_IJ

  subroutine globalSum_IJK(grd_dum, arr, gsum, hsum, zsum, istag, iskip)
    type (dist_grid),   INTENT(IN) :: grd_dum
    real*8,             INTENT(IN) :: arr(grd_dum%i_strt_halo:, &
         &                                      grd_dum%j_strt_halo:,:)
    real*8,             INTENT(OUT):: gsum(size(arr,3))
    real*8, OPTIONAL,   INTENT(OUT):: hsum(2,size(arr,3))
    real*8, OPTIONAL,   INTENT(OUT):: zsum(grd_dum%j_strt: &
         &                                       grd_dum%j_stop, &
         &                                       size(arr,3))
    integer,OPTIONAL,   INTENT(IN) :: istag,iskip

    integer :: i_0, i_1, j_0, j_1, j_0S, j_1S, IM, JM, j_0STG
    real*8  :: zon(grd_dum%j_strt:grd_dum%j_stop,size(arr,3))
    real*8  :: garr(grd_dum%jm_world,size(arr,3))
    logical :: istag_,iskip_
    ! now local

    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    j_0  = grd_dum%j_strt
    j_0S = grd_dum%j_strt_skp
    j_0STG  = grd_dum%j_strt_stgr
    j_1  = grd_dum%j_stop
    j_1S = grd_dum%j_stop_skp
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD

    istag_ = .false.
    if (present(istag)) Then
      if (istag == 1) istag_ = .true.
    end if

    iskip_ = .false.
    if (present(iskip)) Then
      if (iskip == 1) iskip_ = .true.
    end if

    if (istag_) Then
      zon(j_0STG:j_1,:)  = sum(arr(:,j_0STG:j_1,:),1)
    Else if (iskip_) Then
      zon(j_0S:j_1S,:)  = sum(arr(:,j_0S:j_1S,:),1)
    Else
      zon  = sum(arr(:,j_0:j_1,:),1)
    end if

    Call gatherReal8(grd_dum, zon, garr, shape(zon), 1, .false.)

    if (AM_I_ROOT()) then
      if (istag_) then
        gsum = sum(garr(2:JM,:),1)
      Else if(iskip_) then
        gsum = sum(garr(2:JM-1,:),1)
      Else
        gsum = sum(garr(1:JM,:),1)
      endif
      if (present(hsum)) then
        if (istag_) then
          hsum(1,:)   = Sum( garr(2     :JM/2,:),1   )
          hsum(2,:)   = Sum( garr(2+JM/2:JM  ,:),1   )
          hsum(1,:)   = hsum(1,:) + 0.5*garr(1+JM/2,:)
          hsum(2,:)   = hsum(2,:) + 0.5*garr(1+JM/2,:)
        Else
          hsum(1,:)   = Sum( garr(1     :JM/2,:),1   )
          hsum(2,:)   = Sum( garr(1+JM/2:JM  ,:),1   )
        endif
      endif
    endif
    if (present(zsum)) zsum = zon

  end subroutine globalSum_IJK

  subroutine globalSum_OTHER_IJK(grd_dum, arr, gsum, jband, all)
    type (dist_grid),   intent(IN) :: grd_dum
    real*8,             intent(IN) :: arr(:,grd_dum%j_strt_halo:)
    real*8,             intent(OUT):: gsum(size(arr,1))
    integer,           intent(IN) :: jband(2)
    logical,optional,   intent(IN) :: all
    
    integer :: i_0, i_1, j_0, j_1, IM, JM, jb1, jb2
    logical :: all_
    integer :: ierr
#ifdef USE_MPI
    real*8  :: garr(size(arr,1),grd_dum%jm_world)
#endif
    ! now local
    
    all_ = .false.
    if (present(all)) all_ = all
    
    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    j_0  = grd_dum%j_strt
    j_1  = grd_dum%j_stop
    IM = SIZE(arr,1)
    JM   = grd_dum%JM_WORLD

    jb1 = jband(1)
    jb2 = jband(2)

#ifdef USE_MPI
    call gatherReal8(grd_dum, arr, garr, shape(arr), 2, .false.)
    if (AM_I_ROOT()) gsum = sum(garr(:,jb1:jb2),2)
    if (all_) then
      call MPI_BCAST(gsum, size(gsum), MPI_DOUBLE_PRECISION, ROOT_PROCESS, &
           &        getMpiCommunicator(grd_dum), ierr)
    end if
#else
    gsum = sum(arr(:, max(jb1,j_0):min(jb2,j_1) ),2)
#endif
  end subroutine globalSum_OTHER_IJK
  
  subroutine globalSum_OTHER_IJK_IK(grd_dum, arr, gsum, jband, all)
    type (dist_grid),   intent(IN) :: grd_dum
    real*8,             intent(IN) :: arr(:,grd_dum%j_strt_halo:,:)
    real*8,             intent(OUT):: gsum(size(arr,1), size(arr,3))
    integer,            intent(IN) :: jband(2)
    logical,optional,   intent(IN) :: all
    
    integer :: i_0, i_1, j_0, j_1, IM, JM, jb1, jb2
    logical :: all_
    integer :: ierr
#ifdef USE_MPI
    real*8  :: garr(size(arr,1),grd_dum%jm_world,size(arr,3))
#endif
    ! now local

    all_ = .false.
    if (present(all)) all_ = all

    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    j_0  = grd_dum%j_strt
    j_1  = grd_dum%j_stop
    IM = SIZE(arr,1)
    JM   = grd_dum%JM_WORLD

    jb1 = jband(1)
    jb2 = jband(2)

#ifdef USE_MPI
    call gatherReal8(grd_dum, arr, garr, shape(arr), 2, .false.)
    if (AM_I_ROOT()) gsum = sum(garr(:,jb1:jb2,:),2)
    if (all_) then
      call MPI_BCAST(gsum, size(gsum), MPI_DOUBLE_PRECISION, ROOT_PROCESS, &
           &        getMpiCommunicator(grd_dum), ierr)
    end if
#else
    gsum = sum(arr(:, max(jb1,j_0):min(jb2,j_1) ,:),2)
#endif
  end subroutine globalSum_OTHER_IJK_IK
  
  subroutine globalSum_IJK_IK(grd_dum, arr, gsum, all)
    type (dist_grid),   intent(IN) :: grd_dum
    real*8,             intent(IN) :: arr(:,grd_dum%j_strt_halo:,:)
    real*8,             intent(OUT):: gsum(size(arr,1), size(arr,3))
    logical,optional,   intent(IN) :: all
    
    integer :: i_0, i_1, j_0, j_1, IM, JM, LM
    logical :: all_
    integer :: ierr
#ifdef USE_MPI
    real*8  :: garr(size(arr,1),grd_dum%jm_world,size(arr,3))
#endif
    ! now local
#ifdef USE_MPI
    !     type (ESMF_Grid)                           :: GRID
    Integer :: scnts(0:npes_world-1), sdspl(0:npes_world-1)
    integer :: rcnts(0:npes_world-1), rdspl(0:npes_world-1)
    integer ::  dik_map(0:npes_world-1), dik, dik_sum
    integer :: npes, nik, i,k,j,p, ijk, iremain
    real*8, allocatable :: tsum(:)
    real*8, allocatable :: send_buf(:)
    real*8, allocatable :: recv_buf(:,:)
#endif
    
    all_ = .false.
    if (present(all)) all_ = all
    
    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    j_0  = grd_dum%j_strt
    j_1  = grd_dum%j_stop
    IM = SIZE(arr,1)
    JM   = grd_dum%JM_WORLD
    LM =  size(arr,3)
    
#ifdef USE_MPI
    npes = getNumProcesses(grd_dum)
    ! Number of sums each processor computes is dik
    do p = 0, npes-1
      dik_map(p)=(IM*LM)/npes
      iremain=mod(IM*LM, npes)
      if (iremain > 0 .and. iremain > p)  dik_map(p)= dik_map(p)+1
    end do
    
    dik=dik_map(rank)
    
    allocate(tsum(dik))
    allocate(send_buf(maxval(dik_map) *(grd_dum%dj)*npes))
    allocate(recv_buf(dik, JM))
    
    ! ugly packing for transpose
    ijk = 0
    do j = 1, (grd_dum%dj)
      p = 0
      dik_sum=0
      nik = 0
      ijk = dik_map(p)*(j-1)
      do k = 1, lm
        do i = 1, im
          nik = nik+1
          ijk=ijk+1
          send_buf(ijk+dik_sum) = arr(i,j+grd_dum%j_strt-1,k)
          if (nik == dik_map(p)) then
            dik_sum=dik_sum+dik_map(p)*(grd_dum%dj)
            p = p+1
            if (p == npes) exit
            nik= 0
            ijk = dik_map(p)*(j-1)
          end if
        end do
      end do
    end do
    
    scnts=dik_map*(grd_dum%dj)
    sdspl(0)=0
    do p = 1, npes-1
      sdspl(p)=sdspl(p-1)+(grd_dum%dj)*dik_map(p-1)
    end do
    
    rcnts=dik*(grd_dum%dj_map)
    rdspl(0)=0
    do p = 1, npes-1
      rdspl(p)=rdspl(p-1)+rcnts(p-1)
    end do
    
    call MPI_ALLTOALLV(send_buf, scnts, sdspl, mpi_double_precision, &
         &             recv_buf, rcnts, rdspl, mpi_double_precision, &
         &             getMpiCommunicator(grd_dum), ierr)
    
    tsum=sum(recv_buf,2)
    
    rcnts=dik_map
    rdspl(0)=0
    do p = 1, npes-1
      rdspl(p)=rdspl(p-1)+rcnts(p-1)
    end do
    
    call MPI_GatherV(tsum, dik, mpi_double_precision, &
         & gsum, dik_map, rdspl, mpi_double_precision, &
         & ROOT_PROCESS, getMpiCommunicator(grd_dum), ierr)
    
    deallocate(recv_buf)
    deallocate(send_buf)
    deallocate(tsum)
    
    if (all_) then
      call MPI_BCAST(gsum, size(gsum), MPI_DOUBLE_PRECISION, ROOT_PROCESS, &
           &        getMpiCommunicator(grd_dum), ierr)
    end if
#else
    gsum = sum(arr(:,j_0:j_1,:),2)
#endif
  end subroutine globalSum_IJK_IK
  
  subroutine globalSum_JK(grd_dum, arr, gsum, hsum, istag, all)
    type (dist_grid),  intent(IN) :: grd_dum
    real*8,            intent(IN) :: arr(grd_dum%j_strt_halo:,:)
    real*8,            intent(OUT):: gsum(size(arr,2))
    real*8, optional,  intent(OUT):: hsum(2,size(arr,2))
    integer,optional,  intent(IN) :: istag
    logical,optional,  intent(IN) :: all

    integer :: ierr
    integer :: i_0, i_1, j_0, j_1, IM, JM
    real*8  :: garr(grd_dum%jm_world,size(arr,2))
    logical :: istag_
    
    ! now local
    
    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    j_0  = grd_dum%j_strt
    j_1  = grd_dum%j_stop
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD
    istag_ = .false.
    if (present(istag)) then
      if (istag == 1) istag_ = .true.
    end if
    
    call gatherReal8(grd_dum, arr, garr, shape(arr), 1, .false.)

    if (AM_I_ROOT()) then
      if (istag_) then
        gsum = sum(garr(2:JM,:),1)
      else
        gsum = sum(garr(1:JM,:),1)
      endif
      if (present(hsum)) then
        if (istag_) then
          hsum(1,:)   = sum( garr(2     :JM/2,:),1   )
          hsum(2,:)   = sum( garr(2+JM/2:JM  ,:),1   )
          hsum(1,:)   = hsum(1,:) + 0.5*garr(1+JM/2,:)
          hsum(2,:)   = hsum(2,:) + 0.5*garr(1+JM/2,:)
        else
          hsum(1,:)   = sum( garr(1     :JM/2,:),1   )
          hsum(2,:)   = sum( garr(1+JM/2:JM  ,:),1   )
        endif
      endif
    endif
    
#ifdef USE_MPI
    if (present(all)) then
      if (all) then
        call MPI_BCAST(gsum,size(gsum),MPI_DOUBLE_PRECISION,ROOT_PROCESS, &
             &           getMpiCommunicator(grd_dum), ierr)
        if (present(hsum)) call MPI_BCAST(hsum,size(hsum), &
             &           MPI_DOUBLE_PRECISION, ROOT_PROCESS, getMpiCommunicator(grd_dum), ierr)
      end if
    end if
#endif
  end subroutine globalSum_JK
  
  subroutine globalSum_XXXJ_XXX(grd_dum, arr, gsum, all)
    type (dist_grid), intent(IN) :: grd_dum
    real*8, intent(IN) :: arr(:,:,:,grd_dum%j_strt_halo:)
    real*8, intent(Out) :: gsum(:,:,:)
    logical, optional, intent(IN) :: all
    
    integer :: i_0, i_1, j_0, j_1, IM, JM
    real*8  :: garr(size(arr,1),size(arr,2),size(arr,3), &
         &     grd_dum%jm_world)
    logical :: all_
    
    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    j_0  = grd_dum%j_strt
    j_1  = grd_dum%j_stop
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD

#ifdef USE_MPI
    Call gatherReal8(grd_dum, arr, garr, shape(arr), 4, all=all)
    all_=.false.
    if (present(all)) all_=all
    if (AM_I_ROOT() .or. all_)  gsum = SUM(garr,4)
#else
    gsum = SUM(arr(:,:,:,J_0:J_1),4)
#endif

  end subroutine globalSum_XXXJ_XXX

  subroutine globalSum_XXXIJ_XXX(grd_dum, arr, gsum, all)
    type (dist_grid), INTENT(IN) :: grd_dum
    Real*8, INTENT(IN) :: &
         &     arr(:,:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
    Real*8, INTENT(Out) :: gsum(:,:,:)
    Logical, Optional, INTENT(IN) :: all

    integer :: i_0, i_1, j_0, j_1, IM, JM
    real*8  :: garr(size(arr,1),size(arr,2),size(arr,3), &
         &     grd_dum%jm_world)
    real*8  :: larr(size(arr,1),size(arr,2),size(arr,3), &
         &     grd_dum%j_strt_halo:grd_dum%j_stop_halo)
    logical :: all_

    i_0  = grd_dum%i_strt
    i_1  = grd_dum%i_stop
    j_0  = grd_dum%j_strt
    j_1  = grd_dum%j_stop
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD

#ifdef USE_MPI
    larr = sum(arr,4)
    Call gatherReal8(grd_dum, larr, garr, shape(larr), 4, all=all)
    all_=.false.
    if (present(all)) all_=all
    if (AM_I_ROOT() .or. all_)  gsum = SUM(garr,4)
#else
    gsum = SUM(SUM(arr(:,:,:,I_0:I_1,J_0:J_1),4),4)
#endif

  end subroutine globalSum_XXXIJ_XXX

  subroutine checksum_1D(grd_dum, arr, line, file, unit, STGR, SKIP)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    real*8,            INTENT(IN) :: &
         &                arr(grd_dum%j_strt_halo:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit
    logical, OPTIONAL, INTENT(IN) :: stgr
    logical, OPTIONAL, INTENT(IN) :: skip

#ifdef DEBUG_DECOMP
    integer :: unit_
    real*8  :: asum, L1norm
    real*8  :: t_arr(grd_dum%j_strt_halo:grd_dum%j_stop_halo)
    integer :: J_0, J_1
    integer :: stgr_, skip_

    J_0 = grd_dum%J_STRT
    J_1 = grd_dum%J_STOP

    unit_ = checksum_UNIT ! default
    if (present(unit)) unit_ = unit

    stgr_ = 0
    if (present(stgr)) THEN
      if (stgr) stgr_=1
    end if

    skip_ = 0
    if (present(skip)) THEN
      if (skip) skip_=1
    end if

    t_arr = arr
    Call globalSum(grd_dum, t_arr, asum,istag=stgr_,iskip=skip_)
    t_arr(J_0:J_1) = ABS(t_arr(J_0:J_1))
    Call globalSum(grd_dum, t_arr, L1norm,istag=stgr_,iskip=skip_)

    if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') &
         &     file,line, asum, L1norm

#endif

  end subroutine checksum_1D

  subroutine checksum_2D(grd_dum, arr, line, file, unit, stgr, &
       &     skip)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    real*8,            INTENT(IN) :: &
         &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit
    logical, OPTIONAL, INTENT(IN) :: stgr
    logical, OPTIONAL, INTENT(IN) :: skip

    integer :: unit_
    real*8  :: asum, L1norm
    real*8 :: &
         &  t_arr(size(arr,1),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
    integer :: J_0, J_1
    integer :: stgr_,skip_

    J_0 = grd_dum%J_STRT
    J_1 = grd_dum%J_STOP

    unit_ = checksum_UNIT ! default
    if (present(unit)) unit_ = unit

    stgr_ = 0
    if (present(stgr)) THEN
      if (stgr) stgr_=1
    end if

    skip_ = 0
    if (present(skip)) THEN
      if (skip) skip_=1
    end if

    t_arr(:,J_0:J_1) = arr(:,J_0:J_1)
    Call globalSum(grd_dum,      t_arr,asum, istag=stgr_,iskip=skip_)
    t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
    Call globalSum(grd_dum,      t_arr,L1norm,istag=stgr_,iskip=skip_)

    if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') &
         &     file,line, asum, L1norm


  end subroutine checksum_2D

  subroutine checksum_3D(grd_dum, arr, line, file, unit, stgr, skip)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    real*8,            INTENT(IN) :: &
         &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit
    logical, OPTIONAL, INTENT(IN) :: stgr
    logical, OPTIONAL, INTENT(IN) :: skip


    integer :: unit_
    integer :: k
    real*8, DIMENSION(Size(arr,3))  :: asum, L1norm

    real*8 :: &
         &  t_arr(size(arr,1),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
    integer :: J_0, J_1
    Integer :: stgr_,skip_


    J_0 = grd_dum%J_STRT
    J_1 = grd_dum%J_STOP

    unit_ = checksum_UNIT ! default
    if (present(unit)) unit_ = unit

    stgr_ = 0
    if (present(stgr)) THEN
      if (stgr) stgr_=1
    end if

    skip_ = 0
    if (present(skip)) THEN
      if (skip) skip_=1
    end if

    Do k = 1, Size(arr, 3)
      t_arr(:,J_0:J_1) = arr(:,J_0:J_1,k)
      Call globalSum(grd_dum, t_arr, asum(k), istag=stgr_,iskip=skip_)
      t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
      Call globalSum(grd_dum, t_arr,L1norm(k),istag=stgr_,iskip=skip_)
    end Do
    if (AM_I_ROOT()) Then
      Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') &
           &       file,line, Sum(asum), Sum(L1norm)
    end if

  end subroutine checksum_3D

  subroutine checksumj_2D(grd_dum, arr, line, file, unit, stgr, &
       &     skip)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    real*8,            INTENT(IN) :: &
         &                arr(grd_dum%j_strt_halo:,:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit
    logical, OPTIONAL, INTENT(IN) :: stgr
    logical, OPTIONAL, INTENT(IN) :: skip

    integer :: unit_
    real*8  :: asum, L1norm
    real*8 :: &
         &  t_arr(size(arr,2),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
    integer :: J_0, J_1
    integer :: stgr_,skip_

    J_0 = grd_dum%J_STRT
    J_1 = grd_dum%J_STOP

    unit_ = checksum_UNIT ! default
    if (present(unit)) unit_ = unit

    stgr_ = 0
    if (present(stgr)) THEN
      if (stgr) stgr_=1
    end if

    skip_ = 0
    if (present(skip)) THEN
      if (skip) skip_=1
    end if

    t_arr(:,J_0:J_1) = Transpose(arr(J_0:J_1,:))
    Call globalSum(grd_dum,      t_arr,asum, istag=stgr_,iskip=skip_)
    t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
    Call globalSum(grd_dum,      t_arr,L1norm,istag=stgr_,iskip=skip_)

    if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') &
         &     file,line, asum, L1norm


  end subroutine checksumj_2D

  subroutine checksumj_3D(grd_dum, arr, line, file, unit, stgr)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    real*8,            INTENT(IN) :: &
         &                arr(grd_dum%j_strt_halo:,:,:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit
    logical, OPTIONAL, INTENT(IN) :: stgr

    real*8, DIMENSION(Size(arr,1))  :: asum, L1norm
    real*8 :: &
         &  t_arr(size(arr,2),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
    integer :: unit_
    integer :: k
    integer :: J_0, J_1
    Integer :: stgr_


    J_0 = grd_dum%J_STRT
    J_1 = grd_dum%J_STOP

    unit_ = checksum_UNIT ! default
    if (present(unit)) unit_ = unit

    stgr_ = 0
    if (present(stgr)) THEN
      if (stgr) stgr_=1
    end if

    Do k = 1, Size(arr, 3)
      t_arr(:,J_0:J_1) = Transpose(arr(J_0:J_1,:,k))
      Call globalSum(grd_dum, t_arr, asum(k), istag=stgr_)
      t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
      Call globalSum(grd_dum, t_arr, L1norm(k), istag=stgr_)
    end Do
    if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') &
         &     file,line, Sum(asum), Sum(L1norm)


  end subroutine checksumj_3D

  subroutine checksum_COLUMN_2D(grd_dum,arr,line,file,unit,stgr)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    real*8,            INTENT(IN) :: &
         &                arr(:,grd_dum%j_strt_halo:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit
    logical, OPTIONAL, INTENT(IN) :: stgr

    integer :: unit_
    integer :: k
    integer :: stgr_
    real*8, DIMENSION(Size(arr,1))  :: asum, L1norm

    unit_ = checksum_UNIT ! default
    if (present(unit)) unit_ = unit
    stgr_ = 0
    if (present(stgr)) THEN
      if (stgr) stgr_=1
    end if

    Do k = 1, Size(arr, 1)
      Call globalSum(grd_dum,      arr(k,:), asum(k), istag=stgr_)
      Call globalSum(grd_dum, Abs(arr(k,:)), L1norm(k), istag=stgr_)
    end DO

    if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') &
         &     file,line, Sum(asum), Sum(L1norm)

  end subroutine checksum_COLUMN_2D

  subroutine checksum_COLUMN_3D(grd_dum,arr,line,file,unit,stgr)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    real*8,            INTENT(IN) :: &
         &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit
    logical, OPTIONAL, INTENT(IN) :: stgr


    integer :: unit_
    integer :: k
    real*8, DIMENSION(Size(arr,1))  :: asum, L1norm
    integer :: stgr_

    unit_ = checksum_UNIT ! default
    if (present(unit)) unit_ = unit
    stgr_ = 0
    if (present(stgr)) THEN
      if (stgr) stgr_=1
    end if

    Do k = 1, Size(arr, 1)
      Call globalSum(grd_dum,      arr(k,:,:), asum(k), istag=stgr_)
      Call globalSum(grd_dum, Abs(arr(k,:,:)), L1norm(k), istag=stgr_)
    end DO

    if (AM_I_ROOT()) Then
      Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') &
           &       file,line, Sum(asum), Sum(L1norm)
      CALL SYS_FLUSH(unit_)
    end if

  end subroutine checksum_COLUMN_3D

  subroutine INT_checksum_COLUMN_3D(grd_dum, arr, line, file, unit)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    integer,           INTENT(IN) :: &
         &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit

    Call checksum_COLUMN(grd_dum, Real(arr,KIND=KIND(1.0D+0)), line, &
         &     file, unit)

  end subroutine INT_checksum_COLUMN_3D

  subroutine checksum_COLUMN_4D(grd_dum,arr,line,file,unit,stgr, &
       &                              skip)
    TYPE (dist_grid),   INTENT(IN) :: grd_dum
    real*8,            INTENT(IN) :: &
         &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
    integer,           INTENT(IN) :: line
    character(LEN=*),  INTENT(IN) :: file
    integer, OPTIONAL, INTENT(IN) :: unit
    logical, OPTIONAL, INTENT(IN) :: stgr
    logical, OPTIONAL, INTENT(IN) :: skip


    integer :: unit_
    integer :: i_0, i_1, j_0, j_1, k
    real*8, DIMENSION(Size(arr,1),Size(arr,4)) :: asum, L1norm
    Integer :: stgr_,skip_

    unit_ = checksum_UNIT ! default
    if (present(unit)) unit_ = unit

    i_0 = grd_dum%i_strt
    i_1 = grd_dum%i_stop
    j_0 = grd_dum%j_strt
    j_1 = grd_dum%j_stop

    stgr_ = 0
    if (present(stgr)) THEN
      if (stgr) stgr_=1
    end if

    skip_ = 0
    if (present(skip)) THEN
      if (skip) skip_=1
    end if

    Do k = 1, Size(arr,1)
      CALL globalSum(grd_dum,    arr(k,:,:,:),   asum(k,:), &
           &        istag=stgr_, iskip=skip_)
      CALL globalSum(grd_dum,Abs(arr(k,:,:,:)),L1norm(k,:), &
           &       istag=stgr_, iskip=skip_)
    end Do

    if (AM_I_ROOT()) Then
      Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') &
           &       file,line, SUM(asum), SUM(L1norm)
    end if

  end subroutine checksum_COLUMN_4D

end module GlobalSum_mod
