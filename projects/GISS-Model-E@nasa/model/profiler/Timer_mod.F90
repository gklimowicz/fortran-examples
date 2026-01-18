! This module implements a minimal capability for profiling fortran
! applications.  It is meant to be directly included in other source
! trees and thereby be very easily used.  Products which are far
! superior in their capabilities are available, but generally require
! a larger investment to integrate with a user's application.
! Perhaps eventually, this module could be used as an interface to
! other libraries so that the user need not modify their interfaces.

module Timer_mod
   implicit none
   private

   public :: Timer_type ! derived type

   ! Free functions
   public :: getWTime

   ! Methods
   public :: start
   public :: stop
   public :: reset
   public :: resetGlobal
   public :: getAverageTripTime

   ! Accessor methods
   public :: getNumTrips   
   public :: isActive
   public :: getInclusiveTime
   public :: getExclusiveTime
   public :: getMaximumTime
   public :: getMinimumTime

   ! private
   public :: addTrip ! still public for testing purposes
   private :: addTime ! for testing purposes
   private :: setActive
   private :: setInactive

   public :: r64

   integer, parameter :: r64 = selected_real_kind(14)


   public :: gather
   public :: getMinProcess
   public :: getMaxProcess

   public :: setSynchronous

   type Timer_type
      private
      real (kind=r64) :: inclusiveTime  = 0
      real (kind=r64) :: exclusiveTime  = 0
      real (kind=r64) :: maximumTime  = 0
      real (kind=r64) :: minimumTime  = huge(1._r64)
      integer         :: numTrips     = 0
      real (kind=r64) :: startTime    = 0.
      real (kind=r64) :: startExclusiveTime = 0
      logical         :: isActive     = .false.

      logical :: synchronize = .false.
      integer :: communicator = -1
      integer :: minProcess = 0 ! only used with MPI
      integer :: maxProcess = 0 ! only used with MPI

   end type Timer_type

   interface start
      module procedure start_
      module procedure startAtTime
   end interface

   interface stop
      module procedure stop_
      module procedure stopAtTime
   end interface

   interface reset
      module procedure reset_
   end interface

   interface gather
      module procedure gather_timer
   end interface

   ! private shared variable for computing exclusive
   ! time
   real (kind=r64), save :: globalExclusiveTime = 0

contains

   subroutine setSynchronous(this, flag, communicator)
      type (Timer_type), intent(inOut) :: this
      logical, intent(in) :: flag
      integer, intent(in) :: communicator
      this%synchronize = flag
      this%communicator = communicator
   end subroutine setSynchronous

   integer function getNumTrips(this)
      type (Timer_type), intent(in) :: this
      getNumTrips = this%numTrips
   end function getNumTrips

   real(kind=r64) function getWTime() result(time)
#ifdef USE_MPI
      include 'mpif.h'
      time = mpi_Wtime()
#else      
      integer :: counter, rate
      call system_clock(counter, rate)
      time = real(counter,kind=r64) / rate
#endif
   end function getWTime

   subroutine start_(this)
      type (Timer_type), intent(inout) :: this
#ifdef USE_MPI
      integer :: ier
      include 'mpif.h'
#endif
#ifdef USE_MPI
      if (this%synchronize) call mpi_barrier(this%communicator, ier)
#endif
      call startAtTime(this, getWTime())
   end subroutine start_

   subroutine startAtTime(this, time)
      type (Timer_type), intent(inout) :: this
      real(kind=r64), intent(in) :: time

      this%isActive = .true.
      this%startTime = time
      this%startExclusiveTime = globalExclusiveTime
      call addTrip(this)

   end subroutine startAtTime

   subroutine stop_(this)
      type (Timer_type), intent(inout) :: this

      call stopAtTime(this, getWTime())

   end subroutine stop_

   subroutine stopAtTime(this, time)
      type (Timer_type), intent(inout) :: this
      real(kind=r64), intent(in) :: time

      real(kind=r64) :: dtInclusive, dtExclusive

      dtInclusive = (time - this%startTime)
      dtExclusive = dtInclusive - (globalExclusiveTime - this%startExclusiveTime)
      call addTime_(this, dtInclusive, dtExclusive)
      this%isActive = .false.

   end subroutine stopAtTime

   function clockTick()
      real (kind=r64) :: clockTick
#ifdef USE_MPI
      include 'mpif.h'
      clockTick = mpi_Wtick()
      write(*,*)'MPI: tick = ', clockTick
#else      
      integer :: clockRate
      call system_clock(count_rate = clockRate)
      clockTick = 1 / real(clockRate, kind=r64)
      write(*,*)'serial: tick = ', clockTick
#endif
   end function clockTick

   subroutine setActive(this)
      type (Timer_type), intent(inout) :: this
      this%isActive = .true.
   end subroutine setActive

   subroutine setInactive(this)
      type (Timer_type), intent(inout) :: this
      this%isActive = .false.
   end subroutine setInactive

   subroutine reset_(this)
      type (Timer_type), intent(inout) :: this
      call setInactive(this)
      this%numTrips = 0
      this%inclusiveTime = 0
      this%exclusiveTime = 0
      this%maximumTime = 0
      this%minimumTime = huge(1._r64)
      this%synchronize = .false.
   end subroutine reset_

   ! Needed for testing purposes only - otherwise can contribute roundoff that
   ! plagues later tests.
   subroutine resetGlobal()
     globalExclusiveTime = 0
   end subroutine resetGlobal

   subroutine addTime(this, dtInclusive, dtExclusive)
      type (Timer_type), intent(inout) :: this
      real (kind=r64),   intent(in)    :: dtInclusive
      real (kind=r64), optional, intent(in) :: dtExclusive

      real (kind=r64) :: dtExclusive_

      dtExclusive_ = dtInclusive
      if (present(dtExclusive)) dtExclusive_ = dtExclusive
      call addTime_(this, dtInclusive, dtExclusive_)
      
   end subroutine addTime

   subroutine addTime_(this, dtInclusive, dtExclusive)
      type (Timer_type), intent(inout) :: this
      real (kind=r64),   intent(in)    :: dtInclusive
      real (kind=r64),   intent(in)    :: dtExclusive

      this%inclusiveTime = this%inclusiveTime + dtInclusive
      this%exclusiveTime = this%exclusiveTime + dtExclusive

      this%maximumTime = max(this%maximumTime, dtExclusive)
      this%minimumTime = min(this%minimumTime, dtExclusive)

      globalExclusiveTime = globalExclusiveTime + dtExclusive

   end subroutine addTime_

   function getInclusiveTime(this) result (inclusiveTime)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: inclusiveTime

      inclusiveTime = this%inclusiveTime

   end function getInclusiveTime

   function getExclusiveTime(this) result (exclusiveTime)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: exclusiveTime

      exclusiveTime = this%exclusiveTime

   end function getExclusiveTime

   logical function isActive(this)
      type (Timer_type), intent(in) :: this

      isActive = this%isActive

   end function isActive

   function getMaximumTime(this) result(time)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: time

      time = this%maximumTime

   end function getMaximumTime

   function getMinimumTime(this) result(time)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: time

      if (this%numTrips > 0) then
         time = this%minimumTime
      else
         time = 0
      end if

   end function getMinimumTime

   function getAverageTripTime(this) result(time)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: time

      if (this%numTrips == 0) then
         time = 0
      else
         time = this%exclusiveTime / this%numTrips
      end if

   end function getAverageTripTime

   subroutine addTrip(this)
      type (Timer_type), intent(inout) :: this
      this%numTrips = this%numTrips + 1
   end subroutine addTrip

   function countMinMaxAvg(count, minTime, maxTime, avgTime) result (string)
      use TimeFormatUtilities_mod, only: formatSeconds
      integer, intent(in) :: count
      real (kind=r64), intent(in) :: minTime, maxTime, avgTime

      character(len=100) :: string
      character(len=20) :: avgStr, maxStr, minStr

      avgStr = formatSeconds(avgTime, decimalsAfterPoint=6)
      maxStr = formatSeconds(maxTime, decimalsAfterPoint=6)
      minStr = formatSeconds(minTime, decimalsAfterPoint=6)
      
      write(string,'(i4,3(2x,a12))') count, &
           & avgStr, maxStr, minStr

   end function countMinMaxAvg

   integer function getMinProcess(this)
      type (Timer_type), intent(in) :: this
      getMinProcess = this%minProcess
   end function getMinProcess

   integer function getMaxProcess(this)
      type (Timer_type), intent(in) :: this
      getMaxProcess = this%maxProcess
   end function getMaxProcess

#ifdef USE_MPI
   function gather_timer(this, comm) result(globalTimer)
      type (Timer_type), intent(in) :: this
      integer, intent(in) :: comm
      type (Timer_type) :: globalTimer

      integer :: npes, rank, ier
      integer, parameter :: root = 0

      real (kind=r64) :: localInfo(2)
      real (kind=r64):: globalInfo(2)

      include 'mpif.h'

      call MPI_Comm_size(comm, npes, ier)
      call MPI_Comm_rank(comm, rank, ier)

      call MPI_Allreduce(getNumTrips(this),        globalTimer%numTrips,        1, MPI_INTEGER, &
           & MPI_SUM, comm, ier)

      call MPI_Allreduce(getInclusiveTime(this),   globalTimer%inclusiveTime,   1, MPI_DOUBLE_PRECISION, &
           & MPI_SUM, comm, ier)
      call MPI_Allreduce(getexclusiveTime(this),   globalTimer%exclusiveTime,   1, MPI_DOUBLE_PRECISION, &
           & MPI_SUM, comm, ier)

      localInfo = (/ getInclusiveTime(this), real(rank,kind=r64) /)

      call MPI_Allreduce(localInfo, globalInfo, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm, ier)
      globalTimer%maximumTime = globalInfo(1)
      globalTimer%maxProcess  = nint(globalInfo(2))

      call MPI_Allreduce(localInfo, globalInfo, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, comm, ier)
      globalTimer%minimumTime = globalInfo(1)
      globalTimer%minProcess  = nint(globalInfo(2))

      globalTimer%inclusiveTime = globalTimer%inclusiveTime / npes
      globalTimer%exclusiveTime = globalTimer%exclusiveTime / npes

   end function gather_timer

#else
   function gather_timer(this, comm) result(globalTimer)
      type (Timer_type), intent(in) :: this
      integer, intent(in) :: comm
      type (Timer_type) :: globalTimer

      integer :: npes, rank, ier
      integer, parameter :: root = 0

      real (kind=r64) :: localInfo(2)
      real (kind=r64):: globalInfo(2)

      npes = 1
      rank = 0
      globalTimer = this

   end function gather_timer
#endif

end module Timer_mod
