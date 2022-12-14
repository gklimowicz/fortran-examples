module ModelClock_mod
  use BaseTime_mod
  use TimeInterval_mod
  use Time_mod
  implicit none
  private

  public :: ModelClock

  type :: ModelClock
!!$    private
    type (Time) :: currentTime
    type (TimeInterval) :: dt

    ! modelE legacy representation
    integer :: tick

  contains
    procedure :: getTimeAtBeginningOfCurrentDay
    procedure :: getCurrentTime
    procedure :: getTimeInSecondsFromDate
    procedure :: getAbsoluteTimeInSeconds
    procedure :: getTimeTick
    procedure :: getDt
    procedure :: isBeginningOfDay
    procedure :: nextTick
    procedure :: get
    procedure :: getYear
    procedure :: getMonth
    procedure :: getDate
    procedure :: getDayOfYear
    procedure :: getHour
    procedure :: getAbbrev ! month abbreviation
    procedure :: toString => toString_clock
  end type ModelClock

  interface modelClock
     module procedure newModelClock_time
     module procedure newModelClock_string
  end interface modelClock

  ! used to force keyword use in get()
  type UnusedType
  end type UnusedType

contains

  ! constructor
  function newModelClock_time(startTime, dt, startTick) result(clock)
    use AbstractCalendar_mod
    type (ModelClock) :: clock
    type (Time), intent(in) :: startTime
    type (TimeInterval), intent(in) :: dt
    integer, intent(in) :: startTick

    clock%currentTime = startTime
    clock%tick = startTick
    clock%dt = dt

  end function newModelClock_time


  ! constructor
  function newModelClock_string(string, calendar, dt) result(clock)
     use AbstractCalendar_mod
     use BaseTime_mod
     type (ModelClock) :: clock
     character(len=*), intent(in) :: string
     class (AbstractCalendar), intent(in) :: calendar
     type (TimeInterval), intent(in) :: dt

     type (Time) :: t
     integer :: startTick
     character(len=80) :: tmpString

     t = newTime(calendar)

     read(string,'(i8,1x,a)') startTick, tmpString

     call t%setBaseTime(newBaseTime(tmpString))

     clock = ModelClock(t, dt, startTick)
 
  end function newModelClock_string

  subroutine nextTick(this)
    class (ModelClock), intent(inout) :: this

    this%tick = this%tick + 1
    call this%currentTime%setBaseTime(newBaseTime(this%currentTime + this%dt))

  end subroutine nextTick

  integer function getTimeTick(this)
    class (ModelClock), intent(in) :: this
    getTimeTick = this%tick
  end function getTimeTick

  type (TimeInterval) function getDt(this) result(dt)
    class (ModelClock), intent(in) :: this
    dt = this%dt
  end function getDt

  logical function isBeginningOfDay(this)
    class (ModelClock), intent(in) :: this

    type (Time) :: timeAtPreviousStep
    
    timeAtPreviousStep = newTime(this%currentTime%calendar)
    call timeAtPreviousStep%setBaseTime(newBaseTime(this%currentTime - this%dt))

    isBeginningOfDay = (this%getHour() == 0) .and. (timeAtPreviousStep%getHour() /= 0)

  end function isBeginningOfDay

  function getAbsoluteTimeInSeconds(this) result (secs)
    integer*8 :: secs
    class (ModelClock), intent(inout) :: this
    secs = this%currentTime%getWhole()
  end function getAbsoluteTimeInSeconds

  function getTimeAtBeginningOfCurrentDay(this) result(t)
    type (Time) :: t
    class (ModelClock), intent(in) :: this
    integer :: year
    integer :: month
    integer :: date
    t     = this%getCurrentTime()
    year  = t%getYear()
    month = t%getMonth()
    date  = t%getDate()
    call t%setByDate(year, month, date, 0)
  end function getTimeAtBeginningOfCurrentDay

  function getCurrentTime(this) result(t)
    type (Time) :: t
    class (ModelClock), intent(in) :: this
    t = this%currentTime
  end function getCurrentTime

  function getTimeInSecondsFromDate(this, year, month, date, hour) result (seconds)
    use AbstractCalendar_mod, only: AbstractCalendar
    use BaseTime_mod
    type (BaseTime) :: seconds
    class (ModelClock), intent(inout) :: this
    integer, intent(in) :: year, month, date, hour
    type (Time) :: aTime
    class (AbstractCalendar), pointer :: pCalendar

    pCalendar => this%currentTime%calendar
    aTime = newTime(pCalendar)
    call aTime%setByDate(year, month, date, hour)
    seconds = newBaseTime(this%currentTime - aTime)

  end function getTimeInSecondsFromDate


  subroutine get(this, unused, year, month, dayOfYear, date, hour, amn)
!@sum  getDate gets Calendar info from internal timing info
!@auth Gavin Schmidt (updated by Tom CLune)
    use TimeConstants_mod, only: INT_SECONDS_PER_HOUR
    use JulianCalendar_mod, only: JULIAN_MONTHS
    use CalendarMonth_mod, only: LEN_MONTH_ABBREVIATION, CalendarMonth

    class (ModelClock), intent(in) :: this
    type (UnusedType), optional :: unused
    integer, optional, intent(out) :: year
    integer, optional, intent(out) :: month
    integer, optional, intent(out) :: dayOfYear
    integer, optional, intent(out) :: date
    integer, optional, intent(out) :: hour
    character(len=LEN_MONTH_ABBREVIATION), optional, intent(out) :: amn
    integer :: mnth

    if (present(year)) year = this%currentTime%getYear()
    if (present(dayOfYear)) dayOfYear = this%currentTime%getDayOfYear()
    if (present(month)) month = this%currentTime%getMonth()

    if (present(amn)) then
      mnth = this%currentTime%getMonth()
      amn = this%currentTime%getAbbreviation()
    end if

    if (present(date)) date = this%currentTime%getDate()
    if (present(hour)) hour = this%currentTime%getHour()

    return
  end subroutine get

  integer function getYear(this) result(year)
    class (ModelClock), intent(in) :: this
    year = this%currentTime%getYear()
  end function getYear

  integer function getMonth(this) result(month)
    class (ModelClock), intent(in) :: this
    month = this%currentTime%getMonth()
  end function getMonth

  integer function getDayOfYear(this) result(dayOfYear)
    class (ModelClock), intent(in) :: this
    dayOfYear = this%currentTime%getDayOfYear()
  end function getDayOfYear


  integer function getDate(this) result(date)
    class (ModelClock), intent(in) :: this
    date = this%currentTime%getDate()
  end function getDate


  integer function getHour(this) result(hour)
    class (ModelClock), intent(in) :: this

    hour = this%currentTime%getHour()
  end function getHour

  function getAbbrev(this) result(abbrev)
    use CalendarMonth_mod, only: LEN_MONTH_ABBREVIATION
    character(len=LEN_MONTH_ABBREVIATION) abbrev
    class (ModelClock), intent(in) :: this

    abbrev = this%currentTime%getAbbreviation()
  end function getAbbrev

  function toString_clock(this) result(string)
     use StringUtilities_mod, only: toString
     character(len=:), allocatable :: string
     class (ModelClock), intent(in) :: this

     string = toString(this%tick) // ' ' // this%currentTime%toString()

  end function toString_clock

end module ModelClock_mod
