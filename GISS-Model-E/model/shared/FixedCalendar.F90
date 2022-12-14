module FixedCalendar_mod
  use AbstractCalendar_mod, only: AbstractCalendar
  use AbstractCalendar_mod, only: MONTHS_PER_YEAR
  use AbstractCalendar_mod, only: HOURS_PER_DAY
  use CalendarMonth_mod, only: CalendarMonth, LEN_MONTH_ABBREVIATION
  use TimeInterval_mod, only: TimeInterval
  use Rational_mod, only: Rational, modulo, floor, ceiling, nint
  use BaseTime_mod, only: BaseTime, newBaseTime
  implicit none
  private

  public :: FixedCalendar

  type, abstract, extends(AbstractCalendar) :: FixedCalendar
     private
     integer :: daysPerYear
     type (TimeInterval) :: secondsPerDay  ! primary
     type (TimeInterval) :: secondsPerHour ! derived
     type (CalendarMonth) :: months(0:MONTHS_PER_YEAR+1)
   contains

     ! Setters
     procedure :: setDaysPerYear
     procedure :: getDaysPerYear
     procedure :: setSecondsPerDay
     procedure :: setNthCalendarMonth

     procedure :: getYear
     procedure :: getDayOfYear
     procedure :: getMonth
     procedure :: getAbbrev
     procedure :: getCalendarMonth

     procedure :: getDate
     procedure :: getHour
     procedure :: getSeconds

     procedure :: convertToTime

     procedure :: getDaysInYear_year
     procedure :: getDaysInYear_constant
     generic :: getDaysInYear => getDaysInYear_constant
     procedure :: getMaxDaysInYear

     procedure :: getDaysInMonth_monthAndYear
     procedure :: getDaysInMonth_month
     generic :: getDaysInMonth => getDaysInMonth_month

     procedure :: getSecondsInYear_year
     procedure :: getSecondsInYear_constant
     generic :: getSecondsInYear => getSecondsInYear_constant

     procedure :: getSecondsPerHour
     procedure :: getSecondsPerDay

  end type FixedCalendar

  ! Time is expressed as seconds since January 01 0h in BASE_YEAR
  integer, parameter :: BASE_YEAR = 1

contains


  subroutine setDaysPerYear(this, daysPerYear)
    class (FixedCalendar), intent(inout) :: this
    integer, intent(in) :: daysPerYear
    this%daysPerYear = daysPerYear
  end subroutine setDaysPerYear


  function getDaysPerYear(this) result(daysPerYear)
    integer :: daysPerYear
    class (FixedCalendar), intent(in) :: this

    daysPerYear = this%daysPerYear
  end function getDaysPerYear


  subroutine setSecondsPerDay(this, secondsPerDay)
    class (FixedCalendar), intent(inout) :: this
    type (TimeInterval), intent(in) :: secondsPerDay
    this%secondsPerDay = secondsPerDay
    this%secondsPerHour = TimeInterval(secondsPerDay / HOURS_PER_DAY)
  end subroutine setSecondsPerDay


  subroutine setNthCalendarMonth(this, n, month)
    class (FixedCalendar), intent(inout) :: this
    integer, intent(in) :: n
    type (CalendarMonth), intent(in) :: month
    this%months(n) = month
  end subroutine setNthCalendarMonth


  ! Must be careful not to divide by secondsPerYear ... can result in
  ! rational denominator larger than 64 bit.

  integer function getYear(this, t) result(year)
    class (FixedCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    integer :: hour, day
    
    
    hour = floor(t / this%secondsPerHour)
    day = hour / HOURS_PER_DAY
    year = BASE_YEAR + day / this%daysPerYear

  end function getYear


  integer function getDayOfYear(this, t) result(dayOfYear)
    class (FixedCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    type (Rational) :: days
    type (Rational) :: secondsPerYear

    secondsPerYear = this%secondsPerDay * this%daysPerYear
    days = modulo(t, secondsPerYear) / this%secondsPerDay
    dayOfYear = 1 + floor(days) ! convention starts from "1"

  end function getDayOfYear


  integer function getMonth(this, t) result(month)
    class (FixedCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    integer :: dayOfYear
    integer :: m

    dayOfYear = this%getDayOfYear(t)

    m = 1
    do while (dayOfYear > this%months(m)%lastDayInMonth)
       m = m + 1
    end do

    month = m

  end function getMonth

  integer function getDate(this, t) result(date)
    class (FixedCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    integer :: dayOfYear
    integer :: month

    dayOfYear = this%getDayOfYear(t)
    month = this%getMonth(t)
    date = dayOfYear - this%months(month)%firstDayInMonth + 1

  end function getDate

  integer function getHour(this, t) result(hour)
    class (FixedCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    type (Rational) :: hours
    type (Rational) :: secondsPerDay


    hours = modulo(t, this%secondsPerDay) / this%secondsPerHour
    hour = floor(hours)

  end function getHour
  
  function getAbbrev(this, t) result(abbrev)
    character(len=LEN_MONTH_ABBREVIATION) :: abbrev
    class (FixedCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    abbrev = this%months(this%getMonth(t))%abbreviation

  end function getAbbrev

  function getSeconds(this, t) result(seconds)
    type (TimeInterval) :: seconds
    class (FixedCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t


    seconds = TimeInterval(modulo(t, this%secondsPerHour))

  end function getSeconds

  function convertToTime(this, year, month, date, hour) result(t)
    type (BaseTime) :: t
    class (FixedCalendar), intent(in) :: this
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: date
    integer, intent(in) :: hour

    integer :: dayOfYear
    integer :: numDays
    integer :: numYears
    integer :: numHours

    numYears = year - BASE_YEAR
    dayOfYear = this%months(month)%firstDayInMonth + date - 2
    numDays = numYears * this%daysPerYear + dayOfYear
    numHours = numDays*HOURS_PER_DAY + hour

    t = newBaseTime(numHours * this%secondsPerHour)

  end function convertToTime

  function getCalendarMonth(this, month, year) result(m)
    type (CalendarMonth) :: m
    class (FixedCalendar), intent(in) :: this
    integer, intent(in) :: month
    integer, intent(in) :: year ! unused for fixed calendar
    m = this%months(month)
  end function getCalendarMonth


  integer function getDaysInYear_year(this, year) result(numDays)
    class (FixedCalendar), intent(in) :: this
    integer, intent(in) :: year  ! ignore - no leap days in Fixed calendar
    numDays = this%getDaysInYear()
  end function getDaysInYear_year


  integer function getDaysInYear_constant(this) result(numDays)
    class (FixedCalendar), intent(in) :: this
    numDays = this%daysPerYear
  end function getDaysInYear_constant


  integer function getMaxDaysInYear(this) result(numDays)
    class (FixedCalendar), intent(in) :: this
    numDays = this%daysPerYear
  end function getMaxDaysInYear


  integer function getDaysInMonth_month(this, month) result(numDays)
    class (FixedCalendar), intent(in) :: this
    integer, intent(in) :: month
    numDays = this%months(month)%daysInMonth
  end function getDaysInMonth_month

  integer function getDaysInMonth_monthAndYear(this, month, year) result(numDays)
    class (FixedCalendar), intent(in) :: this
    integer, intent(in) :: month
    integer, intent(in) :: year  ! no leap days in fixed calendar
    numDays = this%getDaysInMonth(month)
  end function getDaysInMonth_monthAndYear


  function getSecondsInYear_year(this, year) result(secondsInYear)
    type (TimeInterval) :: secondsInYear
    class (FixedCalendar), intent(in) :: this
    integer, intent(in) :: year ! unused for fixed calendars

    secondsInYear = this%getSecondsInYear()

  end function getSecondsInYear_year

  function getSecondsInYear_constant(this) result(secondsInYear)
    type (TimeInterval) :: secondsInYear
    class (FixedCalendar), intent(in) :: this

    secondsInYear = &
         & TimeInterval((this%secondsPerHour * HOURS_PER_DAY) * this%daysPerYear)

  end function getSecondsInYear_constant


  function getSecondsPerDay(this) result(secondsPerDay)
    type (TimeInterval) :: secondsPerDay
    class (FixedCalendar), intent(in) :: this
    secondsPerDay = this%secondsPerDay
  end function getSecondsPerDay


  function getSecondsPerHour(this) result(secondsPerHour)
    type (TimeInterval) :: secondsPerHour
    class (FixedCalendar), intent(in) :: this
    secondsPerHour = this%secondsPerHour
  end function getSecondsPerHour


end module FixedCalendar_mod
