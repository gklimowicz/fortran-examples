#include "rundeck_opts.h"
#ifdef MPI_DEFS_HACK
#include "mpi_defs.h"
#endif

#define FILL(N) IAND(USABLE_FROM,N)==N
  
module Halo_mod
  use dist_grid_mod
#ifdef USE_MPI
  use MpiSupport_mod, only: getNeighbors
  use MpiSupport_mod, only: createDecompMpiType
  use MpiSupport_mod, only: freeMpiType
#endif
  implicit none
  private

!@var halo_update Update data in halo for local domain using data from
!@+   neighbouring processes
  public :: halo_update_south_pole
  public :: halo_update_north_pole
  public :: halo_update ! Communicate overlapping portions of subdomains
  public :: halo_updateJ ! jx
  public :: halo_update_column ! K, I, J
  public :: halo_update_block ! K, L, I, J
  public :: halo_update_mask ! K, L, I, J

  integer, parameter :: ALL = NORTH + SOUTH ! no east/west for now

  interface halo_update_south_pole
    module procedure halo_update_south_pole_1d
    module procedure halo_update_south_pole_2d
    module procedure halo_update_south_pole_3d
  end interface
  interface halo_update_north_pole
    module procedure halo_update_north_pole_1d
    module procedure halo_update_north_pole_2d
    module procedure halo_update_north_pole_3d
  end interface

!@var halo_update generic wrapper for 2D and 3D routines
  interface halo_update
    module procedure halo_update_1d  ! j
    module procedure halo_update_2d  ! i,j
    module procedure halo_update_3d  ! i,j,k
    module procedure halo_update_2dint  ! i,j,k
  end interface

  interface halo_updatej
    module procedure halo_updatej_2d
  end interface

  interface halo_update_column
    module procedure halo_update_column_2d  ! m,j
    module procedure halo_update_column_3d  ! m,i,j
    module procedure int_halo_update_column_3d  ! m,i,j
    module procedure halo_update_column_4d  ! m,i,j,k
    module procedure halo_update_column_7d  ! m1,m2,m3,m4,i,j
  end interface

  interface halo_update_block
    module procedure halo_update_block_4d  ! k,l,i,j
  end interface

#ifdef USE_MPI
#include "mpif.h"
#endif

contains

  subroutine halo_update_south_pole_1d(grid, arr)
    type (dist_grid),   intent(inout)   :: grid
    real*8,   intent(inout) :: arr(grid%j_strt_halo:)
#ifdef USE_MPI
    integer :: status(MPI_STATUS_SIZE), ierr

    call incrementMpiTag(grid)
    if (haveLatitude(grid, 1) .and. (.not. haveLatitude(grid, 2))) then
      call mpi_send( arr(1), 1, MPI_DOUBLE_PRECISION, 1, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(2), 1, MPI_DOUBLE_PRECISION, 1, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if 
    if (haveLatitude(grid, 2) .and. (.not. haveLatitude(grid, 1))) then
      call mpi_send( arr(2), 1, MPI_DOUBLE_PRECISION, 0, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(1), 1, MPI_DOUBLE_PRECISION, 0, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if
#endif

  end subroutine halo_update_south_pole_1d

  subroutine halo_update_north_pole_1d(grid, arr)
    type (dist_grid),   intent(inout)   :: grid
    real*8,   intent(inout) :: arr(grid%j_strt_halo:)
#ifdef USE_MPI
    integer :: status(MPI_STATUS_SIZE), ierr, jm, npes

    jm = grid%jm_world
    npes = grid%npes_used

    call incrementMpiTag(grid)
    if (haveLatitude(grid, jm) .and. (.not. haveLatitude(grid, jm-1))) then
      call mpi_send( arr(jm), 1, MPI_DOUBLE_PRECISION, npes-2, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(jm-1), 1, MPI_DOUBLE_PRECISION, npes-2, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if 
    if (haveLatitude(grid, jm-1) .and. (.not. haveLatitude(grid, jm))) then
      call mpi_send( arr(jm-1), 1, MPI_DOUBLE_PRECISION, npes-1, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(jm), 1, MPI_DOUBLE_PRECISION, npes-1, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if
#endif

  end subroutine halo_update_north_pole_1d

  subroutine halo_update_south_pole_2d(grid, arr)
    type (dist_grid),   intent(inout)   :: grid
    real*8,   intent(inout) :: arr(grid%i_strt_halo:,grid%j_strt_halo:)
#ifdef USE_MPI
    integer :: status(MPI_STATUS_SIZE), ierr, m

    m  = size(arr,1)

    call incrementMpiTag(grid)
    if (haveLatitude(grid, 1) .and. (.not. haveLatitude(grid, 2))) then
      call mpi_send( arr(1,1), m, MPI_DOUBLE_PRECISION, 1, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(1,2), m, MPI_DOUBLE_PRECISION, 1, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if 
    if (haveLatitude(grid, 2) .and. (.not. haveLatitude(grid, 1))) then
      call mpi_send( arr(1,2), m, MPI_DOUBLE_PRECISION, 0, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(1,1), m, MPI_DOUBLE_PRECISION, 0, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
   end if
#endif

  end subroutine halo_update_south_pole_2d

  subroutine halo_update_north_pole_2d(grid, arr)
    type (dist_grid),   intent(inout)   :: grid
    real*8,   intent(inout) :: arr(grid%i_strt_halo:,grid%j_strt_halo:)
#ifdef USE_MPI
    integer :: status(MPI_STATUS_SIZE), ierr, jm, npes, m

    m  = size(arr,1)
    jm = grid%jm_world
    npes = grid%npes_used

    call incrementMpiTag(grid)
    if (haveLatitude(grid, jm) .and. (.not. haveLatitude(grid, jm-1))) then
      call mpi_send( arr(:,jm), m, MPI_DOUBLE_PRECISION, npes-2, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(:,jm-1), m, MPI_DOUBLE_PRECISION, npes-2, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if 
    if (haveLatitude(grid, jm-1) .and. (.not. haveLatitude(grid, jm))) then
      call mpi_send( arr(:,jm-1), m, MPI_DOUBLE_PRECISION, npes-1, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(:,jm), m, MPI_DOUBLE_PRECISION, npes-1, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if
#endif

  end subroutine halo_update_north_pole_2d

  subroutine halo_update_south_pole_3d(grid, arr)
    type (dist_grid),   intent(inout)   :: grid
    real*8,   intent(inout) :: arr(grid%i_strt_halo:,grid%j_strt_halo:,:)
#ifdef USE_MPI
    integer :: status(MPI_STATUS_SIZE), ierr
    integer :: new_type

    ! create a new mpi type for use in communication
    !-------------------------------
    new_type = createDecompMpiType(MPI_DOUBLE_PRECISION, shape(arr), 2)

    call incrementMpiTag(grid)
    if (haveLatitude(grid, 1) .and. (.not. haveLatitude(grid, 2))) then
      call mpi_send( arr(1,1,1), 1, new_type, 1, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(1,2,1), 1, new_type, 1, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if 
    if (haveLatitude(grid, 2) .and. (.not. haveLatitude(grid, 1))) then
      call mpi_send( arr(1,2,1), 1, new_type, 0, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(1,1,1), 1, new_type, 0, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if

    call freeMpiType(new_type)
#endif

  end subroutine halo_update_south_pole_3d

  subroutine halo_update_north_pole_3d(grid, arr)
    type (dist_grid),   intent(inout)   :: grid
    real*8,   intent(inout) :: arr(grid%i_strt_halo:,grid%j_strt_halo:,:)
#ifdef USE_MPI
    integer :: status(MPI_STATUS_SIZE), ierr, jm, npes
    integer :: new_type

    ! create a new mpi type for use in communication
    !-------------------------------
    new_type = createDecompMpiType(MPI_DOUBLE_PRECISION, shape(arr), 2)

    jm = grid%jm_world
    npes = grid%npes_used

    call incrementMpiTag(grid)
    if (haveLatitude(grid, jm) .and. (.not. haveLatitude(grid, jm-1))) then
      call mpi_send( arr(1,jm,1), 1, new_type, npes-2, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(1,jm-1,1), 1, new_type, npes-2, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if 
    if (haveLatitude(grid, jm-1) .and. (.not. haveLatitude(grid, jm))) then
      call mpi_send( arr(1,jm-1,1), 1, new_type, npes-1, &
        getMpiTag(grid), getMpiCommunicator(grid), ierr )
      call mpi_recv(arr(1,jm,1), 1, new_type, npes-1, &
        getMpiTag(grid), getMpiCommunicator(grid), status, ierr )
    end if

    call freeMpiType(new_type)
#endif

  end subroutine halo_update_north_pole_3d

  subroutine halo_update_1d(grid, arr, from)
    type (dist_grid),   intent(in)   :: grid
    real*8,   intent(inout) :: arr(grid%j_strt_halo:)
    integer, optional, intent(in)    :: from

    call sendrecv(grid, arr, shape(arr), 1, from &
         &     ,hasperiodicbc(grid))

  end subroutine halo_update_1d

  subroutine halo_update_2d(grid, arr, from)
    type (dist_grid),   intent(in)   :: grid
    real*8,   intent(inout) :: arr(grid%i_strt_halo:,grid%j_strt_halo:)
    integer, optional, intent(in)    :: from

    call sendrecv(grid, arr, shape(arr), 2, from, hasperiodicbc(grid))
  end subroutine halo_update_2d

  subroutine halo_update_2dint(grid, arr, from)
    type (dist_grid),   intent(in)    :: grid
    integer,   intent(inout) :: arr(grid%i_strt_halo:,grid%j_strt_halo:)
    integer, optional, intent(in)    :: from

    call sendrecv_int(grid, arr, shape(arr), 2, from, hasperiodicbc(grid))
  end subroutine halo_update_2dint

  subroutine halo_updatej_2d(grid, arr, from)
    type (dist_grid),   intent(in)   :: grid
    real*8,   intent(inout) :: arr(grid%j_strt_halo:,:)
    integer, optional, intent(in)    :: from

    call sendrecv(grid, arr, shape(arr), 1, from, hasperiodicbc(grid))

  end subroutine halo_updatej_2d

  subroutine halo_update_3d(grid, arr, from, jdim)
    type (dist_grid),   intent(in)    :: grid
    real*8,   intent(inout) :: arr(:,:,:)
    integer, optional, intent(in)    :: from
    integer, optional, intent(in)    :: jdim

    integer :: jd

    if(present(jdim)) then
      jd = jdim
    else
      jd = 2
    endif
    call sendrecv(grid, arr, shape(arr), jd, from, hasperiodicbc(grid))

  end subroutine halo_update_3d

#ifdef USE_MPI
  subroutine halo_update_mask(grid, sbufs, sbufn, rbufs, rbufn)
    type (dist_grid),   intent(in)    :: grid
    real*8,   intent(in) :: sbufn(:)
    real*8,   intent(in) :: sbufs(:)
    real*8,   intent(out) :: rbufn(:)
    real*8,   intent(out) :: rbufs(:)

    integer :: numSendSouth, numSendNorth
    integer :: numRecvSouth, numRecvNorth
    integer :: pe_south, pe_north
    integer, save :: tag = 1
    integer, parameter :: NUM_TAGS = 100
    integer :: status(MPI_STATUS_SIZE)
    integer :: ier

    numSendSouth = size(sBufS)
    numSendNorth = size(sBufN)

    numRecvSouth = size(rBufS)
    numRecvNorth = size(rBufN)

    call getNeighbors(rank, npes_world, pe_south, pe_north, .false.)

    tag = 1 + mod(tag - 1, NUM_TAGS)
    call MPI_sendrecv(sBufS, numSendSouth, MPI_DOUBLE_PRECISION, &
         &     pe_south, tag, &
         &     rBufN, numRecvNorth, MPI_DOUBLE_PRECISION, pe_north, tag, &
         &     getMpiCommunicator(grid), status, ier)

    tag = 1 + mod(tag - 1, NUM_TAGS)

    call MPI_sendrecv(sBufN, numSendNorth, MPI_DOUBLE_PRECISION, &
         &     pe_north, tag, &
         &     rBufS, numRecvSouth, MPI_DOUBLE_PRECISION, pe_south, tag, &
         &     getMpiCommunicator(grid), status, ier)

  end subroutine halo_update_mask
#endif

#ifndef USE_MPI
  subroutine halo_update_mask(grid, sbufs, sbufn, rbufs, rbufn)
    type (dist_grid),   intent(in)    :: grid
    real*8,   intent(in) :: sbufn(:)
    real*8,   intent(in) :: sbufs(:)
    real*8,   intent(out) :: rbufn(:)
    real*8,   intent(out) :: rbufs(:)
    ! NOOP
  end subroutine halo_update_mask
#endif

  subroutine halo_update_column_2d(grid, arr, from)
    type (dist_grid),   intent(in)   :: grid
    real*8,   intent(inout) :: arr(:,grid%j_strt_halo:)
    integer, optional, intent(in)    :: from

    call sendrecv(grid, arr, shape(arr), 2, from, hasperiodicbc(grid))

  end subroutine halo_update_column_2d

  subroutine halo_update_column_3d(grid, arr, from)
    type (dist_grid),   intent(in)    :: grid
    real*8,   intent(inout) :: arr(:,grid%i_strt_halo:,grid%j_strt_halo:)
    integer, optional, intent(in)    :: from

    call sendrecv(grid, arr, shape(arr), 3, from, hasperiodicbc(grid))

  end subroutine halo_update_column_3d

  subroutine int_halo_update_column_3d(grid, arr, from)
    type (dist_grid),   intent(in)   :: grid
    integer,  intent(inout) :: arr(:,grid%i_strt_halo:,grid%j_strt_halo:)
    integer, optional, intent(in)    :: from

    real*8 :: foo(size(arr,1),size(arr,2),size(arr,3))

    ! cheat
    foo = arr
    call halo_update_column(grid, foo, from)
    arr = nint(foo)

  end subroutine int_halo_update_column_3d

  subroutine halo_update_column_4d(grid, arr, from)
    type (dist_grid),   intent(in)    :: grid
    real*8,   intent(inout) :: arr(:,grid%i_strt_halo:,grid%j_strt_halo:,:)
    integer, optional, intent(in)    :: from

    call sendrecv(grid, arr, shape(arr), 3, from)

  end subroutine halo_update_column_4d

  subroutine halo_update_block_4d(grid, arr, from)
    type (dist_grid),   intent(in)   :: grid
    real*8,   intent(inout) :: arr(:,:,grid%i_strt_halo:,grid%j_strt_halo:)
    integer, optional, intent(in)    :: from

    call sendrecv(grid, arr, shape(arr), 4, from)

  end subroutine halo_update_block_4d

  subroutine halo_update_column_7d(grid, arr, from)
    type (dist_grid),   intent(in)   :: grid
    real*8, intent(inout) :: arr(:,:,:,:,grid%i_strt_halo:,grid%j_strt_halo:)
    integer, optional, intent(in)    :: from

    call sendrecv(grid, arr, shape(arr), 6, from, hasperiodicbc(grid))

  end subroutine halo_update_column_7d

#ifdef USE_MPI
  Subroutine sendrecv(grid, arr, shp, dist_idx, from, bc_periodic_)
    Type (Dist_Grid) :: grid
    Real(Kind=8) :: arr(*)
    Integer :: shp(:)
    Integer :: dist_idx
    Integer, optional :: from
    Logical, optional :: bc_periodic_

    Integer :: new_type
    Integer :: npy, npx, px, py, pe_south, pe_north
    Integer :: off(4)
    Integer :: USABLE_FROM
    Integer :: status(MPI_STATUS_SIZE), ierr
    Integer :: n, sz
    Logical :: bc_periodic

    USABLE_FROM = usableFrom(from)
    bc_periodic = isPeriodic(bc_periodic_)

    ! create a new mpi type for use in communication
    !-------------------------------
    new_type = createDecompMpiType(MPI_DOUBLE_PRECISION, shp,dist_idx)

    ! Determine neigboring processes
    !-------------------------------
    call gridRootPELocation(grid,  px,  py)
    call gridPELayout(grid, npx, npy)
    Call GetNeighbors(py, npy, pe_south, pe_north, bc_periodic)
    
    sz = Product(shp(1:dist_idx-1))
    n  = shp(dist_idx)
    off(1) = 1
    off(2) = 1+sz*1
    off(3) = 1+sz*(n-2)
    off(4) = 1+sz*(n-1)

    IF(FILL(NORTH)) THEN
      call incrementMpiTag(grid)
      Call MPI_SendRecv( &
           &       arr(off(2)), 1, new_type, pe_south, getMpiTag(grid), &
           &       arr(off(4)), 1, new_type, pe_north, getMpiTag(grid), &
           &       getMpiCommunicator(grid), status, ierr)
    End If

    IF(FILL(SOUTH)) THEN
      call incrementMpiTag(grid)
      Call MPI_SendRecv( &
           &       arr(off(3)), 1, new_type, pe_north, getMpiTag(grid), &
           &       arr(off(1)), 1, new_type, pe_south, getMpiTag(grid), &
           &       getMpiCommunicator(grid), status, ierr)
    End If

    call freeMpiType(new_type)

  end subroutine SendRecv

  subroutine sendrecv_int(grid, arr, shp, dist_idx, from, &
       &     bc_periodic_)
    Type (Dist_Grid) :: grid
    Integer :: arr(*)
    Integer :: shp(:)
    Integer :: dist_idx
    Integer, optional :: from
    Logical, optional :: bc_periodic_

    Integer :: new_type
    Integer :: npy, npx, px, py, pe_south, pe_north
    Integer :: off(4)
    Integer :: USABLE_FROM
    Integer :: status(MPI_STATUS_SIZE), ierr
    Integer :: n, sz
    Logical :: bc_periodic

    USABLE_FROM = usableFrom(from)
    bc_periodic = isPeriodic(bc_periodic_)
    ! create a new mpi type for use in communication
    !-------------------------------
    new_type = createDecompMpiType(MPI_INTEGER, shp,dist_idx)

    ! Determine neigboring processes
    !-------------------------------
    call gridRootPELocation(grid,  px,  py)
    call gridPELayout(grid, npx, npy)
    Call GetNeighbors(py, npy, pe_south, pe_north, bc_periodic)

    sz = Product(shp(1:dist_idx-1))
    n  = shp(dist_idx)
    off(1) = 1
    off(2) = 1+sz*1
    off(3) = 1+sz*(n-2)
    off(4) = 1+sz*(n-1)

    IF(FILL(NORTH)) THEN
      call incrementMpiTag(grid)
      Call MPI_SendRecv( &
           &       arr(off(2)), 1, new_type, pe_south, getMpiTag(grid), &
           &       arr(off(4)), 1, new_type, pe_north, getMpiTag(grid), &
           &       getMpiCommunicator(grid), status, ierr)
    End If

    IF(FILL(SOUTH)) THEN
      call incrementMpiTag(grid)
      Call MPI_SendRecv( &
           &       arr(off(3)), 1, new_type, pe_north, getMpiTag(grid), &
           &       arr(off(1)), 1, new_type, pe_south, getMpiTag(grid), &
           &       getMpiCommunicator(grid), status, ierr)
    End If

    call freeMpiType(new_type)

  end subroutine SendRecv_int
#endif

#ifndef USE_MPI
  subroutine sendrecv(grid, arr, shp, dist_idx, from, bc_periodic_)
    Type (Dist_Grid) :: grid
    Real(Kind=8) :: arr(*)
    Integer :: shp(:)
    Integer :: dist_idx
    Integer, optional :: from
    Logical, optional :: bc_periodic_
  end subroutine sendrecv

  subroutine sendrecv_int(grid, arr, shp, dist_idx, from, bc_periodic_)
    Type (Dist_Grid) :: grid
    Integer :: arr(*)
    Integer :: shp(:)
    Integer :: dist_idx
    Integer, optional :: from
    Logical, optional :: bc_periodic_
  end subroutine sendrecv_int
#endif

! Helper function to handle optional arguments related to halo directions
  integer function usableFrom(fromDirection)
    integer, optional, intent(in) :: fromDirection
    usableFrom = ALL
    if (present(fromDirection)) usableFrom = fromDirection
  end function usableFrom

end module Halo_mod
