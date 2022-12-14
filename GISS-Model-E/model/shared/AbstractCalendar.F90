!---------------------------------------------------------------------
! Calendar classes support conversion of raw times to and from 
! human-comprehensible descriptions of time such as year, month
! date, hour, etc.
!
! Although a climate model could function with a relatively simple
! implementation of a Calendar, the requirement to support simulations
! of other planets introduces the need for additional abstraction.
!
! The abstraction defined in this module should be sufficient for a
! wide variety of calendars, but may need some minor modification for
! treatment of non-constant years (e.g. leap-years in Gregorian
! calendar).
!---------------------------------------------------------------------

module AbstractCalendar_mod
  use AbstractTimeStamp_mod
  use CalendarDate_mod
  use AbstractTimeStampHashMap_mod
  use TimeConstants_mod, only: MONTHS_PER_YEAR => INT_MONTHS_PER_YEAR
  use TimeConstants_mod, only: HOURS_PER_DAY => INT_HOURS_PER_DAY
  implicit none
  private

  public :: AbstractCalendar
  public :: MONTHS_PER_YEAR
  public :: HOURS_PER_DAY

  type, abstract :: AbstractCalendar
     private

     type (AbstractTimeStampHashMap) :: transitionDates
     logical :: verbose

   contains

     ! Time dependent quantities: (accessors)
     procedure(get), deferred :: getYear
     procedure(get), deferred :: getDayOfYear
     procedure(get), deferred :: getMonth
     procedure(get), deferred :: getDate
     procedure(get), deferred :: getHour
     procedure(getAbbrev), deferred :: getAbbrev ! abbrev of current month
     procedure :: getSeconds ! in current hour
     procedure :: getCalendarDate
     procedure :: getAnniversaryDate

     ! Convert calendar info into raw time
     procedure(convertToTime), deferred :: convertToTime


     ! delegate to hash map of TimeStamp's
     procedure :: initTransitionDates
     procedure :: addTransitionDate
     procedure :: findTransitionDate
     procedure :: printTransitionDates ! to unit

     procedure :: print
     procedure :: printFancy
     procedure(getCalendarMonth), deferred :: getCalendarMonth

     ! Things that one might think are time independent but are not:
     procedure(getDaysInYear_year), deferred :: getDaysInYear_year
     procedure :: getDaysInYear_time
     procedure(getMaxDaysInYear), deferred :: getMaxDaysInYear
     generic :: getDaysInYear => getDaysInYear_year, getDaysInYear_time

     procedure(getDaysInMonth_monthAndYear), deferred :: getDaysInMonth_monthAndYear
     procedure :: getDaysInMonth_time
     generic :: getDaysInMonth => getDaysInMonth_monthAndYear, getDaysInMonth_time

     procedure(getSecondsInYear_year), deferred :: getSecondsInYear_year
     procedure :: getSecondsInYear_time
     generic :: getSecondsInYear => getSecondsInYear_year, getSecondsInYear_time

     ! Time independent items:
     procedure, nopass :: getMonthsPerYear
     procedure, nopass :: getHoursPerDay
     procedure(getConstInterval), deferred :: getSecondsPerDay
     procedure(getConstInterval), deferred :: getSecondsPerHour

    procedure :: setVerbose
    procedure :: getVerbose

  end type AbstractCalendar

  abstract interface

     ! Common interface for many methods that extract year/month/day/date/hour.
     integer function get(this, t) result (n)
       use BaseTime_mod, only: BaseTime
       import AbstractCalendar
       class (AbstractCalendar), intent(in) :: this
       class (BaseTime), intent(in) :: t
     end function get

     function getAbbrev(this, t) result (abbrev)
       use CalendarMonth_mod, only: LEN_MONTH_ABBREVIATION
       use BaseTime_mod, only: BaseTime
       import AbstractCalendar
       character(len=LEN_MONTH_ABBREVIATION) :: abbrev
       class (AbstractCalendar), intent(in) :: this
       class (BaseTime), intent(in) :: t
     end function getAbbrev

     function getCalendarMonth(this, month, year) result (m)
       use CalendarMonth_mod, only: CalendarMonth
       import AbstractCalendar
       type (CalendarMonth) :: m
       class (AbstractCalendar), intent(in) :: this
       integer, intent(in) :: month
       integer, intent(in) :: year
     end function getCalendarMonth

     function convertToTime(this, year, month, date, hour) result(t)
       use BaseTime_mod, only: BaseTime
       import AbstractCalendar
       type (BaseTime) :: t
       class (AbstractCalendar), intent(in) :: this
       integer, intent(in) :: year
       integer, intent(in) :: month
       integer, intent(in) :: date
       integer, intent(in) :: hour
     end function convertToTime

     integer function getDaysInYear_year(this, year) result(daysInYear)
       import AbstractCalendar
       class (AbstractCalendar), intent(in) :: this
       integer, intent(in) :: year
     end function getDaysInYear_year

     integer function getDaysInMonth_monthAndYear(this, month, year) result(daysInYear)
       import AbstractCalendar
       class (AbstractCalendar), intent(in) :: this
       integer, intent(in) :: month
       integer, intent(in) :: year
     end function getDaysInMonth_monthAndYear


     function getSecondsInYear_year(this, year) result(secondsInYear)
       use TimeInterval_mod, only: TimeInterval
       import AbstractCalendar
       type (TimeInterval) :: secondsInYear
       class (AbstractCalendar), intent(in) :: this
       integer, intent(in) :: year

     end function getSecondsInYear_year

     
     integer function getMaxDaysInYear(this) result(daysInYear)
       import AbstractCalendar
       class (AbstractCalendar), intent(in) :: this
     end function getMaxDaysInYear


     function getConstInterval(this) result(interval)
       use TimeInterval_mod
       import AbstractCalendar
       type (TimeInterval) :: interval
       class (AbstractCalendar), intent(in) :: this
     end function getConstInterval


  end interface

contains

  function getSeconds(this, t) result (seconds)
    use BaseTime_mod, only: BaseTime
    use TimeInterval_mod, only: TimeInterval
    type (TimeInterval) :: seconds
    class (AbstractCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    type (BaseTime) :: t0
    integer :: year
    integer :: month
    integer :: date
    integer :: hour

    year = this%getYear(t)
    month = this%getMonth(t)
    date = this%getDate(t)
    hour = this%getHour(t)

    t0 = this%convertToTime(year, month, date, hour)
    seconds = TimeInterval(t - t0)
    
  end function getSeconds

  function getCalendarDate(this, t) result(cDate)
    use BaseTime_mod, only: BaseTime
    use CalendarDate_mod
    type (CalendarDate) :: cdate
    class (AbstractCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    cDate = CalendarDate(this%getMonth(t), this%getDate(t), this%getYear(t))

  end function getCalendarDate

  function getAnniversaryDate(this, t) result(cDate)
    use BaseTime_mod, only: BaseTime
    use AnniversaryDate_mod
    type (AnniversaryDate) :: cdate
    class (AbstractCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    cDate = AnniversaryDate(this%getMonth(t), this%getDate(t))

  end function getAnniversaryDate


  ! The procedures below should generally not be overridden by base classes.

  ! Determine year from time, then get days per year for that year.
  integer function getDaysInYear_time(this, t) result(daysInYear)
    use BaseTime_mod, only: BaseTime
    class (AbstractCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    integer :: year
    
    year = this%getYear(t)
    daysInYear = this%getDaysInYear(year)

  end function getDaysInYear_time


  ! Determine month and year from time. Pass those to per-calendar function.
  integer function getDaysInMonth_time(this, t) result(daysInMonth)
    use BaseTime_mod, only: BaseTime
    class (AbstractCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    integer :: month
    integer :: year

    month = this%getMonth(t)
    year = this%getYear(t)
    daysInMonth = this%getDaysInMonth(month, year)

  end function getDaysInMonth_time

  ! Determine year from time
  function getSecondsInYear_time(this, t) result(secondsInYear)
    use TimeInterval_mod, only: TimeInterval
    use BaseTime_mod, only: BaseTime
    type (TimeInterval) :: secondsInYear
    class (AbstractCalendar), intent(in) :: this
    class (BaseTime), intent(in) :: t

    integer :: year
    
    year = this%getYear(t)
    secondsInYear = this%getSecondsInYear_year(year)

  end function getSecondsInYear_time

  ! Static functions of calendar

  integer function getMonthsPerYear() result(monthsPerYear)
    monthsPerYear = MONTHS_PER_YEAR
  end function getMonthsPerYear

  integer function getHoursPerDay() result(hoursPerDay)
    hoursPerDay = HOURS_PER_DAY
  end function getHoursPerDay


  subroutine initTransitionDates(this)
    class (AbstractCalendar), intent(inout) :: this
    this%transitionDates = newAbstractTimeStampHashMap(10)
  end subroutine initTransitionDates


  subroutine addTransitionDate(this, name, timeStamp)
    use CalendarDate_mod
    class (AbstractCalendar), intent(inout) :: this
    character(len=*), intent(in) :: name
    class (AbstractTimeStamp), intent(in) :: timeStamp

    call this%transitionDates%insert(name, timeStamp)
    
  end subroutine addTransitionDate

  
  function findTransitionDate(this, name) result(date)
    class (AbstractTimeStamp), pointer :: date
    class (AbstractCalendar), intent(in) :: this
    character(len=*), intent(in) :: name

    date => this%transitionDates%getReference(name)

  end function findTransitionDate


  subroutine print(this, year, unit)
    use iso_fortran_env, only: OUTPUT_UNIT
    use CalendarMonth_mod, only: CalendarMonth
    class (AbstractCalendar), intent(in) :: this
    integer, intent(in) :: year
    integer, optional, intent(in) :: unit

    integer :: month
    class (CalendarMonth), allocatable :: m
    integer :: unit_
    character(len=32) :: fmt

    unit_ = OUTPUT_UNIT
    if (present(unit)) unit_ = unit

    write(unit_,'(58("-"),"|")')
    write(unit_,'(a20,i5,33(" "),"|")')'Calendar for Year:',year
    write(unit_,'(1x,a,a4,x,"|",4(x,a7,x,"|"))') &
         & 'Full Name   ','Abbr','# days ','1st day','mid day','lst day'
    write(unit_,'(58("-"),"|")')
    do month = 1, MONTHS_PER_YEAR
       allocate(m, source = this%getCalendarMonth(month, year))
       write(fmt,'("(1x,a,",i0,"x,a4,x,a1,4(x,i7,x,a1))")') 12 - len_trim(m%fullName)
       write(unit_,trim(fmt)) trim(m%fullName), m%abbreviation, "|",&
            & m%daysInMonth, "|",m%firstDayInMonth, "|",m%midDayInMonth, "|", m%lastDayInMonth, "|"
       deallocate(m)
    end do
    write(unit_,'(58("-"),"|")')
    write(unit_,*) ' '
    call this%printTransitionDates(unit_)

  end subroutine print

  subroutine printFancy(this, year, unit)
    use iso_fortran_env, only: OUTPUT_UNIT
    use CalendarMonth_mod, only: CalendarMonth
    class (AbstractCalendar), intent(in) :: this
    integer, intent(in) :: year
    integer, optional, intent(in) :: unit

    integer :: month
    class (CalendarMonth), allocatable :: m
    integer :: unit_

    unit_ = OUTPUT_UNIT
    if (present(unit)) unit_ = unit

    write(unit_,*)'Calendar for Year:',year

    do month = 1, MONTHS_PER_YEAR
       allocate(m, source = this%getCalendarMonth(month, year))
       call m%print(unit_, year)
       deallocate(m)
    end do

    call this%printTransitionDates(unit_)

  end subroutine printFancy

  subroutine printTransitionDates(this, unit)
    use AbstractTimeStampHashMap_mod
    class (AbstractCalendar), intent(in) :: this
    integer, intent(in) :: unit

    type (AbstractTimeStampHashMapIterator) :: iter
    class (AbstractTimeStamp), pointer :: p
    character(len=24) :: fmt

    iter = this%transitionDates%begin()
    do while(iter /= this%transitionDates%last())
       p => iter%value()
       write(fmt, '("(a,",i0,"("".""))")') 30-len_trim(iter%key())
       write(unit,fmt,advance='no') trim(iter%key())
       call p%print(unit)
       call iter%next()
    end do
       
    
  end subroutine printTransitionDates


  subroutine setVerbose(this, verbose)
     class (AbstractCalendar), intent(inout) :: this
     logical, intent(in) :: verbose

     this%verbose = verbose

  end subroutine setVerbose

  logical function getVerbose(this)
     class (AbstractCalendar), intent(in) :: this
     
     getVerbose = this%verbose

  end function getVerbose


end module AbstractCalendar_mod
