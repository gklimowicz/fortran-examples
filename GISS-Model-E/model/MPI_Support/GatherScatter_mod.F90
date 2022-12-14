#include "rundeck_opts.h"

module GatherScatter_mod
  use MpiSupport_mod
  use dist_grid_mod
  implicit none
  private

  public :: pack_data
  public :: pack_block
  public :: pack_dataj
  public :: pack_j

  public :: unpack_data
  public :: unpack_block
  public :: unpack_dataj
  public :: unpack_j

  public :: pack_column
  public :: unpack_column

  public :: gatherReal8
  public :: scatterReal8

#ifdef USE_MPI
  public :: integerGather
  public :: integerScatter
#endif

!@var PACK Generic routine to pack  a global array
!@+   with the data from the corresponding distributed array.
  interface pack_data
    module procedure pack_1d       ! (i)
    module procedure ipack_1d      ! (i)
    module procedure pack_2d       ! (i,j)
    module procedure ipack_2d      ! (i,j)
    module procedure lpack_2d      ! (i,j)
    module procedure pack_3d       ! (i,j,l)
    module procedure ipack_3d      ! (i,j,l)
    module procedure pack_4d       ! (i,j,l,m)
    module procedure pack_5d       ! (i,j,l,m,n)
  end interface

  interface pack_dataj
    module procedure packj_2d     ! (j,k)
    module procedure ipackj_2d    ! (j,k)
    module procedure packj_3d     ! (j,k,l)
    module procedure packj_4d     ! (j,k,l,m)
  end interface

!@var UNPACK Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
  interface unpack_data
    module procedure unpack_1d      ! (i)
    module procedure unpack_2d      ! (i,j)
    module procedure iunpack_2d     ! (i,j)
    module procedure lunpack_2d     ! (i,j)
    module procedure unpack_3d      ! (i,j,l)
    module procedure iunpack_3d     ! (i,j,l)
    module procedure unpack_4d      ! (i,j,l,m)
    module procedure unpack_5d      ! (i,j,l,m,n)
  end interface

  interface unpack_dataj
    module procedure unpackj_2d     ! (j,k)
    module procedure unpackj_3d     ! (j,k,l)
    module procedure unpackj_4d     ! (j,k,l,m)
  end interface

!@var PACK_COLUMN Generic routine to pack  a global array
!@+   with the data from the corresponding distributed array.
  interface pack_column
    module procedure pack_column_1d  ! (k,  j  )
    module procedure pack_column_2d  ! (k,i,j  )
    module procedure pack_column_i2d ! (k,i,j  )
    module procedure pack_column_3d  ! (k,i,j,l)
  end interface

!@var UNPACK_COLUMN Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
  interface unpack_column
    module procedure unpack_column_1d  ! (k,  j  )
    module procedure unpack_column_2d  ! (k,i,j  )
    module procedure iunpack_column_2d ! (k,i,j  )
    module procedure unpack_column_3d  ! (k,i,j,l)
  end interface

!@var PACK_BLOCK  Generic routine to pack  a global array
!@+   with the data from the corresponding distributed array.
  interface pack_block
    module procedure ipack_block_2d    ! (k,l,i,j  )
    module procedure  pack_block_2d    ! (k,l,i,j  )
    module procedure  pack_block_3d    ! (k,l,m,i,j)
  end interface

!@var UNPACK_BLOCK  Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
  interface unpack_block
    module procedure iunpack_block_2d    ! (k,l,i,j  )
    module procedure  unpack_block_2d    ! (k,l,i,j  )
    module procedure  unpack_block_3d    ! (k,l,m,i,j)
  end interface

!@var PACK_J Generic routine to pack  a global array
!@+   with the data from the corresponding distributed array.
  interface pack_j
    module procedure unpackj_2d     ! (j,k)
    module procedure unpackj_3d     ! (j,k,l)
    module procedure unpackj_4d     ! (j,k,l,m)
  end interface
!@var UNPACK Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
  interface unpack_j
    module procedure unpackj_2d     ! (j,k)
    module procedure unpackj_3d     ! (j,k,l)
    module procedure unpackj_4d     ! (j,k,l,m)
  end interface

  interface localCopy
    module procedure localCopyReal8
    module procedure localCopyInteger
    module Procedure localCopyLogical
  end interface

#ifdef USE_MPI
  include 'mpif.h'
#endif

contains

#ifdef USE_MPI

  subroutine getCountsAndDisplacements(AI, counts, displacements)
    type (AxisIndex), intent(in) :: AI(0:,:)
    integer, intent(out) :: counts(0:)
    integer, intent(out) :: displacements(0:)

    integer :: i1, in, count
    integer :: p, dim
    integer :: numProcesses, numDims

    numProcesses = size(AI,1)
    numDims = size(AI,2)

    displacements(0) = 0
    do p = 0, numProcesses - 1
      count = 1
      do dim = 1, size(AI,2)
        i1 = AI(p,dim)%min
        in = AI(p,dim)%max
        count = count*(in - i1 + 1)
      end do
      counts(p) = count
      displacements(p+1) = displacements(p) + counts(p)

    enddo

  end subroutine getCountsAndDisplacements

  subroutine integerGather(grid, arr_loc, arr_glob, shp, numDims, dist_idx, baseType, baseSize, all)
    type (dist_grid), intent(in) :: grid
    integer, intent(in) :: arr_loc(*)
    integer, intent(out) :: arr_glob(*)
    integer, intent(in) :: numDims
    integer, intent(in) :: shp(numDims)
    integer, intent(in) :: dist_idx

    integer, intent(in) :: baseType
    integer, intent(in) :: baseSize
    logical, intent(in) :: all

    integer :: globalType, localType
    integer :: ierr, rc
    integer :: scount, offset, npes
    integer, allocatable :: rcounts(:), displs(:)
    type (AxisIndex), pointer :: ai(:,:)

    ! create new mpi types for use in communication
    !----------------------------------------------
    localType = createDecompMpiType(baseType,shp,dist_idx)
    globalType = getGlobalMpiType(baseType, shp, dist_idx, grid%JM_WORLD)

    npes = getNumAllProcesses(grid)

    allocate (rcounts(0:NPES-1), displs(0:NPES), stat=rc)
    allocate(AI(NPES,3))
    call getAxisIndex(grid, AI)
    call getCountsAndDisplacements(AI(:,2:2), rcounts, displs)
    deallocate(AI)

    scount = rcounts(rank)
    offset = baseSize*getOffset(shp, scount, dist_idx)

    if (all) then
      call MPI_AllGatherV(arr_loc(1 + offset), scount, localType, &
           &       arr_glob(1), rcounts, displs, globalType,  &
           &       getMpiCommunicator(grid), ierr)
    else
      call MPI_GatherV(arr_loc(1 + offset), scount, localType,   &
           &       arr_glob(1), rcounts, displs, globalType, &
           &       ROOT_PROCESS, getMpiCommunicator(grid), ierr)
    end if

    deallocate(rcounts,displs)

    call freeMpiType(localType)
    call freeMpiType(globalType)

  end subroutine integerGather

  subroutine integerScatter(grid, arr_glob, arr_loc, shp, numDims, dist_idx, baseType, baseSize)
    type (dist_grid) :: grid
    integer, intent(in) :: arr_glob(*)
    integer, intent(out) :: arr_loc(*)
    integer, intent(in) :: numDims
    integer, intent(in) :: shp(numDims)
    integer, intent(in) :: dist_idx
    integer, intent(in) :: baseType
    integer, intent(in) :: baseSize

    integer :: shp_glob(size(shp))
    integer :: localType, globalType
    integer :: rc, ierr
    integer :: p, rcount, offset, npes
    integer, allocatable :: scounts(:), displs(:)
    type (AxisIndex), pointer :: AI(:,:)

    ! create a new mpi type for use in communication
    !-------------------------------
    localType = createDecompMpiType(baseType,shp, dist_idx)
    globalType = getGlobalMpiType(baseType, shp, dist_idx, grid%JM_WORLD)

    npes = getNumAllProcesses(grid)

    allocate (scounts(0:NPES-1), displs(0:NPES), stat=rc)
    allocate(AI(NPES,3))
    call getAxisIndex(grid, AI)
    call getCountsAndDisplacements(AI(:,2:2), scounts, displs)
    deallocate(AI)

    rcount = scounts(rank)
    offset = baseSize*getOffset(shp, rcount, dist_idx)

    Call MPI_ScatterV(arr_glob(1), scounts, displs, globalType, &
         &     arr_loc(1+offset), rcount, localType, &
         &     ROOT_PROCESS, getMpiCommunicator(grid), ierr)

    call freeMpiType(localType)
    call freeMpiType(globalType)
    
    deallocate(scounts,displs)

  end subroutine integerScatter

  subroutine gatherReal8(grid, arr_loc, arr_glob, shp, dist_idx, all)
    type (dist_grid), intent(in) :: grid
    real*8, intent(in) :: arr_loc(*)
    real*8, intent(out) :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional :: all

    logical :: all_
    ! f77 interface handles integer, real, logical, ...
    external :: genericTypeGather 

    all_=.false.
    If (Present(all)) all_=all

    call genericTypeGather(grid, arr_loc, arr_glob, shp, size(shp), dist_idx, &
           & MPI_DOUBLE_PRECISION, 2, all_)

  end subroutine gatherReal8

  subroutine gatherInteger(grid, arr_loc, arr_glob, shp, dist_idx, all)
    type (dist_grid), intent(in) :: grid
    integer, intent(in) :: arr_loc(*)
    integer, intent(out) :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional :: all

    logical :: all_
    ! f77 interface handles integer, real, logical, ...
    external :: genericTypeGather 

    all_=.false.
    If (Present(all)) all_=all

    call genericTypeGather(grid, arr_loc, arr_glob, shp, size(shp), dist_idx, &
           & MPI_INTEGER, 1, all_)

  end subroutine gatherInteger

  subroutine gatherLogical(grid, arr_loc, arr_glob, shp, dist_idx, all)
    type (dist_grid), intent(in) :: grid
    logical, intent(in) :: arr_loc(*)
    logical, intent(out) :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional :: all

    logical :: all_
    ! f77 interface handles integer, real, logical, ...
    external :: genericTypeGather 

    all_=.false.
    If (Present(all)) all_=all

    call genericTypeGather(grid, arr_loc, arr_glob, shp, size(shp), dist_idx, &
           & MPI_INTEGER, 1, all_)

  end subroutine gatherLogical

  integer function getGlobalMpiType(baseType, shp, dist_idx, jm) result(newType)
    integer, intent(in) :: baseType
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    integer, intent(in) :: jm
    
    integer :: globalShape(size(shp))
    
    globalShape = shp
    globalShape(dist_idx) = JM
    newType = createDecompMpiType(baseType, globalShape, dist_idx)
  end function getGlobalMpiType
  
  integer function getOffset(shp, numInterior, dist_idx) result(offset)
!@sum Computes the offset of the interior of local arrays vs their
!@+ exterior.   Should be just 1, when there is no halo.

!@var shp is the shape of a local array (including halos)
    integer, intent(in) :: shp(:)
!@var numInterior is the size of local interior of distributed axis
    integer, intent(in) :: numInterior
!@var dist_idx specifies which axis is distributed (only supports 1d dist)
    integer, intent(in) :: dist_idx
    
    integer :: haloWidth

    haloWidth = (shp(dist_idx) - numInterior)/2
    offset = product(shp(:dist_idx-1)) * haloWidth
  end function getOffset

  subroutine scatterReal8(grid, arr_glob, arr_loc, shp, dist_idx, local)
    type (dist_grid) :: grid
    real(kind=8) :: arr_glob(*)
    real(kind=8) :: arr_loc(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: local
    logical :: local_

    local_ = .false.
    if (present(local)) local_ = local

    if (local_) then
      call localCopyReal8(grid, arr_glob, arr_loc, getGlobalShape(grid, shp, dist_idx))
    else
      call genericTypeScatter(grid, arr_glob, arr_loc, shp, size(shp), dist_idx, &
           & MPI_DOUBLE_PRECISION, 2)
    end if
  end subroutine scatterReal8

  subroutine scatterInteger(grid, arr_glob, arr_loc, shp, dist_idx, local)
    type (dist_grid) :: grid
    integer :: arr_glob(*)
    integer :: arr_loc(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: local
    logical :: local_

    local_ = .false.
    if (present(local)) local_ = local

    if (local_) then
      call localCopyInteger(grid, arr_glob, arr_loc, getGlobalShape(grid, shp, dist_idx))
    else
      call genericTypeScatter(grid, arr_glob, arr_loc, shp, size(shp), dist_idx, &
           & MPI_INTEGER, 1)
    end if

  end subroutine scatterInteger

  subroutine scatterLogical(grid, arr_glob, arr_loc, shp, dist_idx, local)
    type (dist_grid) :: grid
    logical :: arr_glob(*)
    logical :: arr_loc(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: local
    logical :: local_

    local_ = .false.
   if (present(local)) local_ = local

    if (local_) then
      call localCopyLogical(grid, arr_glob, arr_loc, getGlobalShape(grid, shp, dist_idx))
    else
      call genericTypeScatter(grid, arr_glob, arr_loc, shp, size(shp), dist_idx, &
           & MPI_INTEGER, 1)
    end if
  end subroutine scatterLogical

#endif

#ifndef USE_MPI
  subroutine gatherReal8(grid, arr_loc, arr_glob, shp, dist_idx, all)
    type (dist_grid) :: grid
    real(kind=8) :: arr_loc(*)
    real(kind=8) :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: all

    integer :: nn

    nn = product(shp)
    arr_glob(1:nn) = arr_loc(1:nn)

  end subroutine gatherReal8

  subroutine gatherInteger(grid, arr_loc, arr_glob, shp, dist_idx, all)
    type (dist_grid) :: grid
    integer :: arr_loc(*)
    integer :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: all

    integer :: nn

    nn = product(shp)
    arr_glob(1:nn) = arr_loc(1:nn)

  end subroutine gatherInteger

  subroutine gatherLogical(grid, arr_loc, arr_glob, shp, dist_idx, all)
    type (dist_grid) :: grid
    logical :: arr_loc(*)
    logical :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: all

    integer :: nn

    nn = product(shp)
    arr_glob(1:nn) = arr_loc(1:nn)

  end subroutine gatherLogical
  
  subroutine scatterReal8(grid, arr_glob, arr_loc, shp, dist_idx, local)
    type (dist_grid) :: grid
    real(kind=8) :: arr_loc(*)
    real(kind=8) :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: local ! unused in serial mode

    call localCopyReal8(grid, arr_glob, arr_loc, getGlobalShape(grid, shp, dist_idx))

  end subroutine scatterReal8

  subroutine scatterInteger(grid, arr_glob, arr_loc, shp, dist_idx, local)
    type (dist_grid) :: grid
    integer :: arr_loc(*)
    integer :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: local ! unused in serial mode

    call localCopyInteger(grid, arr_glob, arr_loc, getGlobalShape(grid, shp, dist_idx))

  end subroutine scatterInteger

  subroutine scatterLogical(grid, arr_glob, arr_loc, shp, dist_idx, local)
    type (dist_grid) :: grid
    logical :: arr_loc(*)
    logical :: arr_glob(*)
    integer, intent(in) :: shp(:)
    integer, intent(in) :: dist_idx
    logical, optional, intent(in) :: local ! unused in serial mode

    call localCopyLogical(grid, arr_glob, arr_loc, getGlobalShape(grid, shp, dist_idx))

  end subroutine scatterLogical
#endif

  subroutine localCopyReal8(grid, arr_glob, arr_loc, globalShape)
    type (dist_grid), intent(in) :: grid
    integer, intent(in) :: globalShape(3)
    real(kind=8) :: arr_glob(globalShape(1),globalShape(2),globalShape(3))
    real(kind=8) :: arr_loc(globalShape(1),grid%j_strt_halo:grid%j_stop_halo,globalShape(3))
      
    integer :: j_0, j_1
      
    call getDomainBounds(grid, j_strt=j_0, j_stop=j_1)
    arr_loc(:,j_0:j_1,:) = arr_glob(:,j_0:j_1,:)
    
  end subroutine localCopyReal8

  subroutine localCopyInteger(grid, arr_glob, arr_loc, globalShape)
    type (dist_grid), intent(in) :: grid
    integer, intent(in) :: globalShape(3)
    integer :: arr_glob(globalShape(1),globalShape(2),globalShape(3))
    integer :: arr_loc(globalShape(1),grid%j_strt_halo:grid%j_stop_halo,globalShape(3))
      
    integer :: j_0, j_1
      
    call getDomainBounds(grid, j_strt=j_0, j_stop=j_1)
    arr_loc(:,j_0:j_1,:) = arr_glob(:,j_0:j_1,:)
    
  end subroutine localCopyInteger

  subroutine localCopyLogical(grid, arr_glob, arr_loc, globalShape)
    type (dist_grid), intent(in) :: grid
    integer, intent(in) :: globalShape(3)
    logical :: arr_glob(globalShape(1),globalShape(2),globalShape(3))
    logical :: arr_loc(globalShape(1),grid%j_strt_halo:grid%j_stop_halo,globalShape(3))
      
    integer :: j_0, j_1
      
    call getDomainBounds(grid, j_strt=j_0, j_stop=j_1)
    arr_loc(:,j_0:j_1,:) = arr_glob(:,j_0:j_1,:)
    
  end subroutine localCopyLogical
    
  function getGlobalShape(grid, localShape, dist_idx) result(globalShape)
    type (dist_grid), intent(in) :: grid
    integer, intent(in) :: localShape(:)
    integer, intent(in) :: dist_idx

    integer :: globalShape(3)

    globalShape(1) = product(localShape(1:dist_idx-1))
    globalShape(2) = grid%jm_world
    globalShape(3) = product(localShape(dist_idx+1:))

  end function getGlobalShape

  subroutine pack_1d(grd_dum,ARR,ARR_GLOB)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(grd_dum%j_strt_halo:)
    real*8, intent(out) :: arr_glob(grd_dum%jm_world)

    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)

  end subroutine pack_1d

  subroutine ipack_1d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr(grd_dum%j_strt_halo:)
    integer, intent(out) :: arr_glob(grd_dum%jm_world)

    call gatherInteger(grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)

  end subroutine ipack_1d

  subroutine pack_2d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum

    real*8, intent(in) :: arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
    real*8, intent(inout) :: arr_glob(:,:)

    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)

  end subroutine pack_2d

  subroutine ipack_2d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
    integer, intent(inout) :: arr_glob(:,:)

    call gatherInteger(grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)

  end subroutine ipack_2d

  subroutine lpack_2d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    logical, intent(in) :: arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
    logical, intent(inout) :: arr_glob(:,:)

    call gatherLogical(grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)

  end subroutine lpack_2d

  subroutine pack_3d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
    real*8, intent(out) :: arr_glob(:,:,:)

    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)

  end subroutine pack_3d

  subroutine ipack_3d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
    integer, intent(out) :: arr_glob(:,:,:)

    call gatherInteger(grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)

  end subroutine ipack_3d

  subroutine pack_4d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:,:)
    real*8, intent(out) :: arr_glob(:,:,:,:)

    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)

  end subroutine pack_4d

  subroutine pack_5d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:,:,:)
    real*8, intent(out) :: arr_glob(:,:,:,:,:)

    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)

  end subroutine pack_5d

  subroutine packj_2d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(grd_dum%j_strt_halo:,:)
    real*8, intent(inout) :: arr_glob(:,:)

    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)

  end subroutine packj_2d

  subroutine ipackj_2d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr(grd_dum%j_strt_halo:,:)
    integer, intent(inout) :: arr_glob(:,:)

    call gatherInteger(grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)

  end subroutine ipackj_2d

  subroutine packj_3d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(grd_dum%j_strt_halo:,:,:)
    real*8, intent(inout) :: arr_glob(:,:,:)

    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)

  end subroutine packj_3d

  subroutine packj_4d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(grd_dum%j_strt_halo:,:,:,:)
    real*8, intent(inout) :: arr_glob(:,:,:,:)

    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)

  end subroutine packj_4d

  subroutine unpack_1d(grd_dum,arr_glob,arr)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(grd_dum%jm_world)
    real*8, intent(out) :: arr(grd_dum%j_strt_halo:)

    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 1)

  end subroutine unpack_1d

  subroutine unpack_2d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:)
    real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:)
    logical, optional :: local

    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 2, local=local)

  end subroutine unpack_2d

  subroutine iunpack_2d(grd_dum, arr_glob, arr, local)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr_glob(:,:)
    integer, intent(out) :: arr(:,grd_dum%j_strt_halo:)
    logical, optional :: local

    call scatterInteger(grd_dum, arr_glob, arr, shape(arr), 2, local=local)

  end subroutine iunpack_2d

  subroutine lunpack_2d(grd_dum, arr_glob, arr, local)
    type (dist_grid),  intent(in) :: grd_dum
    logical, intent(in) :: arr_glob(:,:)
    logical, intent(out) :: arr(:,grd_dum%j_strt_halo:)
    logical, optional :: local

    call scatterLogical(grd_dum, arr_glob, arr, shape(arr), 2, local=local)

  end subroutine lunpack_2d

  subroutine unpack_3d(grd_dum, arr_glob, arr, local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:)
    real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:,:)
    logical, optional :: local

    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 2, local=local)

  end subroutine unpack_3d

  subroutine iunpack_3d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr_glob(:,:,:)
    integer, intent(out) :: arr(:,grd_dum%j_strt_halo:,:)
    logical, optional :: local

    call scatterInteger(grd_dum, arr_glob, arr, shape(arr), 2, local=local)

  end subroutine iunpack_3d

  subroutine unpack_4d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:,:)
    real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:,:,:)
    logical, optional :: local
    
    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 2, local=local)

  end subroutine unpack_4d
  
  subroutine unpack_5d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:,:,:)
    real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:,:,:,:)
    logical, optional :: local
    
    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 2, local=local)
    
  end subroutine unpack_5d
  
  subroutine unpackj_2d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:)
    real*8, intent(out) :: arr(grd_dum%j_strt_halo:,:)
    logical, optional :: local
    
    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 1, local=local)
    
  end subroutine unpackj_2d
  
  subroutine unpackj_3d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:)
    real*8, intent(out) :: arr(grd_dum%j_strt_halo:,:,:)
    logical, optional :: local
    
    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 1, local=local)
    
  end subroutine unpackj_3d
  
  subroutine unpackj_4d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:,:)
    real*8, intent(out) :: arr(grd_dum%j_strt_halo:,:,:,:)
    logical, optional :: local

    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 1, local=local)
    
  end subroutine unpackj_4d
  
  subroutine pack_column_1d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(:,grd_dum%j_strt_halo:)
    real*8, intent(out) :: arr_glob(:,:)
    
    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)

  end subroutine pack_column_1d
  
  subroutine pack_column_2d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(:,:,grd_dum%j_strt_halo:)
    real*8, intent(inout) :: arr_glob(:,:,:)
    
    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 3, all=.false.)
    
  end subroutine pack_column_2d
  
  subroutine pack_column_i2d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr(:,:,grd_dum%j_strt_halo:)
    integer, intent(inout) :: arr_glob(:,:,:)
    
    call gatherInteger(grd_dum, arr, arr_glob, shape(arr), 3, all=.false.)
    
  end subroutine pack_column_i2d
  
  subroutine pack_column_3d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(:,:,grd_dum%j_strt_halo:,:)
    real*8, intent(out) :: arr_glob(:,:,:,:)
    
    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 3, all=.false.)
    
  end subroutine pack_column_3d
  
  subroutine pack_block_3d(grd_dum,arr,arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr(:,:,:,:,grd_dum%j_strt_halo:)
    real*8, intent(out) :: arr_glob(:,:,:,:,:)
    
    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 5, all=.false.)
    
  end subroutine pack_block_3d
  
  subroutine unpack_column_1d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:)
    real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:)
    logical, optional :: local

    call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 2, local=local)

  end subroutine unpack_column_1d
  
  subroutine unpack_column_2d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:)
    real*8, intent(out) :: arr(:,:,grd_dum%j_strt_halo:)
    logical, optional, intent(in) :: local
    
    call scatterReal8(grd_dum , arr_glob, arr, shape(arr), 3, local=local)

  end subroutine unpack_column_2d
  
  subroutine iunpack_column_2d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr_glob(:,:,:)
    integer, intent(out) :: arr(:,:,grd_dum%j_strt_halo:)
    logical, optional :: local
    
    call scatterInteger(grd_dum, arr_glob, arr, shape(arr), 3, local=local)

  end subroutine iunpack_column_2d
  
  subroutine unpack_column_3d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:,:)
    real*8, intent(out) :: arr(:,:,grd_dum%j_strt_halo:,:)
    logical, optional, intent(in) :: local
    
    call scatterReal8(grd_dum , arr_glob, arr, shape(arr), 3, local=local)
    
  end subroutine unpack_column_3d

  subroutine unpack_block_3d(grd_dum, arr_glob, arr, local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:,:,:)
    real*8, intent(out) :: arr(:,:,:,:,grd_dum%j_strt_halo:)
    logical, optional, intent(in) :: local
    
    call scatterReal8(grd_dum , arr_glob, arr, shape(arr), 5, local=local)
    
  end subroutine unpack_block_3d
  
  subroutine ipack_block_2d(grd_dum, arr, arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr(:,:,:,grd_dum%j_strt_halo:)
    integer, intent(inout) :: arr_glob(:,:,:,:)

    call gatherInteger(grd_dum, arr, arr_glob, shape(arr), 4)

  end subroutine ipack_block_2d
  
  subroutine pack_block_2d(grd_dum, arr, arr_glob)
    type (dist_grid),  intent(in) :: grd_dum
    real*8 , intent(in) :: arr(:,:,:,grd_dum%j_strt_halo:)
    real*8 , intent(inout) :: arr_glob(:,:,:,:)
    
    call gatherReal8(grd_dum, arr, arr_glob, shape(arr), 4, all=.false.)
    
  end subroutine pack_block_2d
  
  subroutine iunpack_block_2d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    integer, intent(in) :: arr_glob(:,:,:,:)
    integer, intent(out) :: arr(:,:,:,grd_dum%j_strt_halo:)
    logical, optional :: local
    
    call scatterInteger(grd_dum, arr_glob, arr, shape(arr), 4, local=local)

  end subroutine iunpack_block_2d
  
  subroutine unpack_block_2d(grd_dum,arr_glob,arr,local)
    type (dist_grid),  intent(in) :: grd_dum
    real*8, intent(in) :: arr_glob(:,:,:,:)
    real*8, intent(out) :: arr(:,:,:,grd_dum%j_strt_halo:)
    logical, optional, intent(in) :: local

    call scatterReal8(grd_dum , arr_glob, arr, shape(arr), 4, local=local)

  end subroutine unpack_block_2d
  
end module GatherScatter_mod

#ifdef USE_MPI
subroutine genericTypeGather(grid, arr_loc, arr_glob, shp, numDims, dist_idx, &
     & baseType, baseSize, all)
  use dist_grid_mod, only: dist_grid
  use GatherScatter_mod, only: integerGather
  implicit none
  type (dist_grid), intent(in) :: grid
  integer, intent(in) :: arr_loc(*)
  integer, intent(out) :: arr_glob(*)
  integer, intent(in) :: numDims
  integer, intent(in) :: shp(numDims)
  integer, intent(in) :: dist_idx
  integer, intent(in) :: baseType
  integer, intent(in) :: baseSize
  logical, intent(in) :: all

  call integerGather(grid, arr_loc, arr_glob, shp, numDims, dist_idx, baseType, baseSize, all)

end subroutine genericTypeGather

subroutine genericTypeScatter(grid, arr_glob, arr_loc, shp, numDims, dist_idx, &
     & baseType, baseSize)
  use dist_grid_mod, only: dist_grid
  use GatherScatter_mod, only: integerScatter
  implicit none
  type (dist_grid), intent(in) :: grid
  integer, intent(in) :: arr_glob(*)
  integer, intent(out) :: arr_loc(*)
  integer, intent(in) :: numDims
  integer, intent(in) :: shp(numDims)
  integer, intent(in) :: dist_idx
  integer, intent(in) :: baseType
  integer, intent(in) :: baseSize

  call integerScatter(grid, arr_glob, arr_loc, shp, numDims, dist_idx, baseType, baseSize)

end subroutine genericTypeScatter
#endif
