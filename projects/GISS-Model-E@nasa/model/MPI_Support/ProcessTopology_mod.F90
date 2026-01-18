module ProcessTopology_mod
  implicit none
  private

  public :: ProcessTopology_type
  public :: newProcessTopology
  public :: getNeighborRank
  public :: getMpiRank
  public :: getRankTopology
  public :: getPartition
  public :: useMpi
  public :: getMpiCommunicator
  public :: flatRank
  public :: amRoot

  integer, parameter :: NUM_DIMENSIONS = 3
  type ProcessTopology_type
    private

    logical :: useMPI
    integer :: communicator
    integer :: rank
    integer :: npes
    integer :: partition(NUM_DIMENSIONS)
    integer :: rankTopology(NUM_DIMENSIONS)

  end type ProcessTopology_type

  interface getNeighborRank
    module procedure getNeighborRank_util
    module procedure getNeighborRank_method
  end interface

  interface newProcessTopology
    module procedure newProcessTopology_serial
#ifdef USE_MPI
    module procedure newProcessTopology_mpi
#endif
  end interface

  interface amRoot
    module Procedure amRoot_serial
#ifdef USE_MPI
    module Procedure amRoot_mpi
#endif
  end interface

contains

  function newProcessTopology_serial() result (this)
    type (ProcessTopology_type) ::this

    this%communicator = -1
    this%useMPI = .false.
    this%rank = 0
    this%npes = 1
    this%partition = (/1,1,1/)
    this%rankTopology = (/0,0,0/)

  end function newProcessTopology_serial

#ifdef USE_MPI
  function newProcessTopology_mpi(communicator, partition) result (this)
    type (ProcessTopology_type) :: this
    integer, intent(in) :: communicator
    integer, intent(in) :: partition(NUM_DIMENSIONS)
    integer :: ier

    this%useMPI = .true.
    this%communicator = communicator

    call MPI_Comm_rank(communicator, this%rank, ier)
    call MPI_Comm_size(communicator, this%npes, ier)

    this%partition = partition

    this%rankTopology(1) = mod(this%rank, partition(1))
    this%rankTopology(2) = mod(this%rank/partition(1), partition(2))
    this%rankTopology(3) = this%rank/(partition(1)*partition(2))

  end function newProcessTopology_mpi
#endif

  logical function amRoot_serial() result(amRoot)
    amRoot = .true.
  end function amRoot_serial

#ifdef USE_MPI
  logical function amRoot_mpi(this) result(amRoot)
    type (ProcessTopology_type), intent(in) :: this
    !    if (this%useMPI) then
    amRoot = (this%rank == 0)
    !    else
    !       amRoot = .true.
    !    end if

  end function amRoot_mpi
#endif

  function getNeighborRank_util(partition, rank, direction) result(neighborRank)
    integer, intent(in) :: partition(NUM_DIMENSIONS)
    integer, intent(in) :: rank(NUM_DIMENSIONS)
    integer, intent(in) :: direction(NUM_DIMENSIONS)
    integer :: neighborRank(NUM_DIMENSIONS)

    neighborRank = mod(rank + direction + partition, partition)

  end function getNeighborRank_util

  function getNeighborRank_method(this, direction) result(neighborRank)
    type (ProcessTopology_type), intent(in) :: this
    integer, intent(in) :: direction(NUM_DIMENSIONS)
    integer :: neighborRank(NUM_DIMENSIONS)

    neighborRank = getNeighborRank(this%partition, this%rankTopology, direction)

  end function getNeighborRank_method

  function getMpiRank(this) result(mpiRank) 
    type (ProcessTopology_type), intent(inout) :: this
    integer mpiRank

    mpiRank = this%rank

  end function getMpiRank

  function getRankTopology(this) result(topology) 
    type (ProcessTopology_type), intent(in) :: this
    integer :: topology(NUM_DIMENSIONS)

    topology = this%rankTopology

  end function getRankTopology

  logical function useMpi(this)
    type (ProcessTopology_type), intent(in) :: this
    useMPI = this%useMPI
  end function useMpi

  integer function getMpiCommunicator(this) result(communicator)
    type (ProcessTopology_type), intent(in) :: this
    communicator = this%communicator
  end function getMpiCommunicator

  integer function flatRank(this, rankTopology) result(rank)
    type (ProcessTopology_type), intent(in) :: this
    integer, intent(in) :: rankTopology(NUM_DIMENSIONS)

    rank = rankTopology(1) + this%partition(1)*(rankTopology(2) + this%partition(2)*rankTopology(3))

  end function flatRank

  function getPartition(this) result(partition)
    type (ProcessTopology_type), intent(in) :: this
    integer :: partition(NUM_DIMENSIONS)

    partition = this%partition

  end function getPartition

end module ProcessTopology_mod
