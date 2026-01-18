! This module implements a singleton intended to manage the set
! of timers associated with an entire application.  A special
! timer "main" is created at initialization. 
module TimerList_mod
   use Timer_mod, only: Timer_type
   implicit none
   private

   public :: TimerList_type
   public :: initialize
   public :: finalize

   public :: getNumTimers
   public :: isTimer
   public :: addTimer
   public :: getTimer
   public :: getName

   public :: start
   public :: stop

   public :: reset
   public :: getDefaultList

   public :: gather

   integer, parameter :: MAX_NAME_LENGTH = 30

   type NamedTimer_type
      private
      type (Timer_type) :: timer
      character(len=MAX_NAME_LENGTH) :: name
   end type NamedTimer_type

   type TimerList_type
      private
      type (NamedTimer_type), pointer :: list(:) => null()
   end type TimerList_type

   type (TimerList_type), target, save :: defaultList ! singleton

   integer, parameter :: MAX_TIMERS = 100
   integer, parameter :: NOT_FOUND =-1

   interface initialize
      module procedure initialize_
      module procedure initializeDefault
   end interface

   interface finalize
      module procedure finalize_
      module procedure finalizeDefault
   end interface

   interface start
      module procedure startByName
      module procedure startByNameAtTime
      module procedure startByNameDefault
      module procedure startByNameAtTimeDefault 
   end interface

   interface stop
      module procedure stopByName
      module procedure stopByNameAtTime
      module procedure stopByNameDefault
      module procedure stopByNameAtTimeDefault
   end interface

   interface getName
      module procedure getNameByIndex
   end interface

   interface getTimer
      module procedure getTimer_
      module procedure getTimerDefault
      module procedure getTimerByIndex
   end interface

   interface isTimer
      module procedure isTimer_
      module procedure isTimerDefault
   end interface

   interface getNumTimers
      module procedure getNumTimers_
      module procedure getNumTimersDefault
   end interface

   interface addTimer
      module procedure addTimer_
      module procedure addTimerDefault
   end interface

   interface reset
      module procedure reset_
      module procedure resetDefault
      module procedure resetTimer
   end interface

   interface gather
      module procedure gather_list
   end interface

   integer, parameter :: r64 = selected_real_kind(14)
   
   logical, save :: test_ = .false.

contains

   function getDefaultList()
      type(TimerList_type), pointer :: getDefaultList
      getDefaultList => defaultList
   end function getDefaultList

   subroutine initialize_(this, mockTime)
      type (TimerList_type), intent(inout) :: this
      real(kind=r64), optional, intent(in) :: mockTime
      real(kind=r64) :: mockTime_

      if (present(mockTime)) mockTime_ = mockTime

      if (associated(this%list)) deallocate(this%list)
      allocate(this%list(0))

      call addTimer(this, 'main')
      if (present(mockTime)) then
         call start(this, 'main', time = mockTime_)
      else
         call start(this, 'main')
      end if
   end subroutine initialize_
      
   ! test enables
   subroutine initializeDefault(mockTime)
      real(kind=r64), optional, intent(in) :: mockTime
      call initialize(defaultList, mockTime)
   end subroutine initializeDefault
   
   subroutine finalize_(this, mockTime)
      type (TimerList_type), target, intent(inout) :: this
      real(kind=r64), optional, intent(in) :: mockTime

      real(kind=r64) :: mockTime_

      if (present(mockTime)) mockTime_ = mockTime

      if (present(mockTime)) then
         call stop(this, 'main', time = mockTime_)
      else
         call stop(this, 'main')
      end if

      call checkTimerConsistenncy()

   contains

      subroutine checkTimerConsistenncy()
         use Timer_mod, only: isActive
         type (NamedTimer_type), pointer :: namedTimer
         integer :: i
         character(len=50) :: message
         do i = 1, getNumTimers(this)
            namedTimer => this%list(i)
            if (isActive(namedTimer%timer)) then
               message = 'Unbalanced start/stop for timer <'//trim(namedTimer%name)//'>.'
               call stop_model(message, 256)
               return
            end if
         end do
      end subroutine checkTimerConsistenncy
       
   end subroutine finalize_

   subroutine finalizeDefault(mockTime)
      real(kind=r64), optional, intent(in) :: mockTime
      call finalize(defaultList, mockTime)
   end subroutine finalizeDefault

   subroutine reset_(this)
      type (TimerList_type), intent(inOut) :: this
      if (associated(this%list)) deallocate(this%list)
      allocate(this%list(0))
   end subroutine reset_

   subroutine resetDefault()
      call reset(defaultList)
   end subroutine resetDefault

   subroutine resetTimer(name)
      use Timer_mod, only: reset
      character(len=*), intent(in) :: name
      type (Timer_type), pointer :: timer
      timer => getTimer(name)
      if (associated(timer)) call reset(timer)
   end subroutine resetTimer

   logical function isTimer_(this, name)
      type (TimerList_type), intent(in) :: this
      character(len=*), intent(in) :: name

      isTimer_ = getIndex(this, name) /= NOT_FOUND

   end function isTimer_
   
   logical function isTimerDefault(name)
      character(len=*), intent(in) :: name
      isTimerDefault = isTimer(defaultList, name)
   end function isTimerDefault
   
   integer function getNumTimers_(this)
      type (TimerList_type), intent(in) :: this
      getNumTimers_ = size(this%list)
   end function getNumTimers_

   integer function getNumTimersDefault()
      getNumTimersDefault = getNumTimers(defaultList)
   end function getNumTimersDefault
   
   subroutine addTimer_(this, name, synchronize, communicator)
      use Timer_mod, only: reset
      use Timer_mod, only: setSynchronous
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name
      logical, optional, intent(in) :: synchronize
      integer, optional, intent(in) :: communicator
      type (NamedTimer_type), allocatable :: tmpList(:)

      integer :: n

      n = size(this%list)

      if (n > 0) then
         allocate(tmpList(n))
         tmpList = this%list
      end if
      if (associated(this%list)) deallocate(this%list)
      allocate(this%list(n+1))
      if (n > 0) then
         this%list(:n) = tmpList
         deallocate(tmpList)
      end if
      this%list(n+1)%name = name
      call reset(this%list(n+1)%timer)
      if (present(synchronize)) then
         call setSynchronous(this%list(n+1)%timer, synchronize, communicator)
      end if

   end subroutine addTimer_

   subroutine addTimerDefault(name, synchronize)
      character(len=*), intent(in) :: name
      logical, optional, intent(in) :: synchronize
      call addTimer(defaultList, name, synchronize)
   end subroutine addTimerDefault

   subroutine startByName(this, name)
      use Timer_mod, only: start
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name
      integer :: index

      index = getIndex(this, name)
      if (index == NOT_FOUND) then
         call addTimer(this, name)
         index = getNumTimers(this)
      end if

      call start(this%list(index)%timer)

   end subroutine startByName

   subroutine startByNameDefault(name)
      character(len=*), intent(in) :: name
      call startByName(defaultList, name)
   end subroutine startByNameDefault

   subroutine startByNameAtTime(this, name, time)
      use Timer_mod, only: start
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name
      real(kind=r64), intent(in) :: time
      call start(this%list(getIndex(this, name))%timer, time)
   end subroutine startByNameAtTime

   subroutine startByNameAtTimeDefault(name, time)
      use Timer_mod, only: start
      character(len=*), intent(in) :: name
      real(kind=r64), intent(in) :: time
      call startByNameAtTime(defaultList, name, time)
   end subroutine startByNameAtTimeDefault

   subroutine stopByName(this, name)
      use Timer_mod, only: stop
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name

      character(len=70) :: message
      integer :: index

      index = getIndex(this, name)
      if (index /= NOT_FOUND) then
         call stop(this%list(index)%timer)
      else
         message = 'Timer <'//trim(name)//'> has not been declared prior to use.'
         call stop_model(message, 256)
         return
      end if

   end subroutine stopByName

   subroutine stopByNameDefault(name)
      character(len=*), intent(in) :: name
      call stopByName(defaultList, name)
   end subroutine stopByNameDefault

   subroutine stopByNameAtTime(this, name, time)
      use Timer_mod, only: stop
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name
      real(kind=r64), intent(in) :: time

      integer :: index
      index = getIndex(this, name)
      call stop(this%list(index)%timer, time)

   end subroutine stopByNameAtTime

   subroutine stopByNameAtTimeDefault(name, time)
      use Timer_mod, only: stop
      character(len=*), intent(in) :: name
      real(kind=r64), intent(in) :: time
      call stopByNameAtTime(defaultList, name, time)
   end subroutine stopByNameAtTimeDefault

   integer function getIndex(this, name)
      type (TimerList_type), intent(in) :: this
      character(len=*), intent(in) :: name
      integer :: i

      do i = 1, getNumTimers(this)
         if (trim(adjustl(name)) == trim(adjustl(this%list(i)%name))) then
            getIndex = i
            return
         end if
      end do

      getIndex = NOT_FOUND

   end function getIndex

   function getTimer_(this, name) result(timer)
      type (TimerList_type), target, intent(in) :: this
      character(len=*), intent(in) :: name
      type (Timer_type), pointer :: timer

      timer => this%list(getIndex(this, name))%timer

   end function getTimer_

   function getTimerDefault(name) result(timer)
      character(len=*), intent(in) :: name
      type (Timer_type), pointer :: timer
      timer => getTimer(defaultList, name)
   end function getTimerDefault

   function getTimerByIndex(this, index) result(timer)
      type (TimerList_type), target, intent(in) :: this
      integer, intent(in) :: index
      type (Timer_type), pointer :: timer

      timer => this%list(index)%timer

   end function getTimerByIndex

   function getNameByIndex(this, index) result(name)
      type (TimerList_type), target, intent(in) :: this
      integer, intent(in) :: index
      character(len=MAX_NAME_LENGTH) :: name

      name = this%list(index)%name

   end function getNameByIndex

   function gather_list(this, communicator) result(globalList)
      use Timer_mod, only: gather
      type (TimerList_type), intent(in) :: this
      integer, intent(in) :: communicator

      type (TimerList_type) :: globalList
      type (Timer_type), pointer :: globalTimer
      type (Timer_type), pointer :: localTimer
      integer :: i

      call initialize(globalList)

      do i = 1, getNumTimers(this)
         call addTimer(globalList, trim(getName(this, i)))
         localTimer => getTimer(this, i)
         globalTimer => getTimer(globalList, i)
         globalTimer = gather(localTimer, communicator)
      end do
   end function gather_list

end module TimerList_mod
