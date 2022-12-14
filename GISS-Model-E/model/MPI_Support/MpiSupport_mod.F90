#include "rundeck_opts.h"
  
module MpiSupport_mod
!@sum Collection of low-level routines that support higher-level
!@+   interfaces for domain decomposition.
!@auth Tom Clune <thomas.l.clune@nasa.gov>
  implicit none
  private

#ifdef USE_MPI
  public :: getNeighbors
  public :: createDecompMpiType
  public :: freeMpiType
#endif

  public :: am_i_root
  interface am_i_root
    module procedure amRootWorld
    module procedure amRootComm
  end interface

#ifdef USE_MPI
  include 'mpif.h'
#endif

  type MpiType_type
    integer :: numBlocks
    integer :: blockLen
    integer :: stride
    integer :: baseType
    integer :: mpiType
  end type MpiType_type
  ! Avoid allocation performance penalty - use safe max size
  integer, parameter :: MAX_MPI_TYPES = 1024
  integer, parameter :: NOT_FOUND = -1

  type (mpiType_type) :: mpiTypes(MAX_MPI_TYPES)
  integer, save :: numMpiTypes = 0

  integer, parameter, public ::  ROOT_PROCESS = 0

contains

#ifdef USE_MPI
  subroutine getNeighbors(rank, npes, pe_south,pe_north,bc_periodic)
    integer, intent(in)  :: rank
    integer, intent(in)  :: npes
    integer, intent(out) :: pe_south
    integer, intent(out) :: pe_north
    logical, intent(in)  :: bc_periodic

    pe_south = mod( npes + rank - 1, npes)
    pe_north = mod(        rank + 1, npes)

    if ( .not. bc_periodic ) then
      if ( rank == 0      ) pe_south = MPI_PROC_NULL
      if ( rank == npes-1 ) pe_north = MPI_PROC_NULL
    endif

  end subroutine getNeighbors

  function createDecompMpiType(baseType, counts, decompIdx) result(newType)
    integer, intent(in) :: baseType
    integer, intent(in) :: counts(:)
    integer, intent(in) :: decompIdx
    integer :: newType

    integer :: stride, numBlocks, blockLen
    integer :: new_len
    integer(kind=MPI_ADDRESS_KIND) disp(2), ext_lb, baseByteLen
    integer :: mpiType1, mpiType2
    integer :: ierr

#ifdef MPITYPE_LOOKUP_HACK
    integer m
#endif

    numBlocks = product(counts(decompIdx+1:))
    blockLen = product(counts(:decompIdx-1))
    stride = counts(decompIdx) * blockLen
    ext_lb = 0

#ifdef MPITYPE_LOOKUP_HACK
    m = getTypeIndex(numBlocks, blockLen, stride, baseType)
    if (m /= NOT_FOUND) then
      newType = mpiTypes(m)%mpiType 
      return
    end if
#endif

    call MPI_Type_vector(numBlocks, blockLen, stride, baseType, &
         &     mpiType1, ierr)
    !Call MPI_Type_extent(baseType, baseByteLen, ierr)
    call MPI_Type_get_extent(baseType, ext_lb, baseByteLen, ierr)
    new_len = baseByteLen * blockLen
    !Call MPI_Type_struct(2, (/ 1, 1 /), (/ 0, new_len /),
    disp(1)=0; disp(2)=new_len;
    call MPI_Type_create_struct(2, (/ 1, 1 /), disp, &
         &     (/ mpiType1, MPI_UB /), mpiType2, ierr)
    call MPI_Type_Free(mpiType1, ierr)

    call MPI_Type_Commit(mpiType2, ierr)

    newType = mpiType2

#ifdef MPITYPE_LOOKUP_HACK
    numMpiTypes = numMpiTypes + 1
    if ( numMpiTypes > MAX_MPI_TYPES ) call stop_model("numMpiTypes > 1024",255)
    mpiTypes(numMpiTypes) = MpiType_type(numBlocks, blockLen, stride, baseType, newType)
#endif

  end function createDecompMpiType

  subroutine freeMpiType(mpiType)
    integer, intent(in) :: mpiType

    integer :: ierr
#ifndef MPITYPE_LOOKUP_HACK
    call MPI_Type_Free(mpiType, ierr)
#endif
  end subroutine freeMpiType
#endif

  integer function getTypeIndex(numBlocks, blockLen, stride, baseType)
    integer, intent(in) :: numBlocks
    integer, intent(in) :: blockLen
    integer, intent(in) :: stride
    integer, intent(in) :: baseType

    integer :: m

    do m = 1, numMpiTypes
      if ( numBlocks == mpiTypes(m)%numBlocks .and. &
           & blockLen == mpiTypes(m)%blockLen .and. &
           & stride == mpiTypes(m)%stride   .and. &
           & baseType == mpiTypes(m)%baseType ) then
        getTypeIndex = m
        return
      endif
    enddo

    getTypeIndex = NOT_FOUND

  end function getTypeIndex

  logical function amRootWorld()
    integer :: rank
    integer :: ierr

#ifdef USE_MPI
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    amRootWorld = (rank == ROOT_PROCESS)
#else
    amRootWorld = .true.
#endif

  end function amRootWorld

  logical function amRootComm(mpiCommunicator)
    integer, intent(in) :: mpiCommunicator
    integer :: rank
    integer :: ierr

#ifdef USE_MPI
    call MPI_Comm_rank(mpiCommunicator, rank, ierr)
    amRootComm = (rank == ROOT_PROCESS)
#else
    amRootComm = .true.
#endif

  end function amRootComm

end module MpiSupport_mod
