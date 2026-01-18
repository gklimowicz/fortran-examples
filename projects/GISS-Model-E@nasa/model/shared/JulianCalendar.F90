module JulianCalendar_mod
  use FixedCalendar_mod, only: FixedCalendar
  use AbstractCalendar_mod, only: MONTHS_PER_YEAR
  use CalendarMonth_mod, only: CalendarMonth
  use TimeInterval_mod, only: TimeInterval
  use TimeConstants_mod, only: INT_SECONDS_PER_DAY
  use TimeConstants_mod, only: SECONDS_PER_YEAR
  use TimeConstants_mod, only: INT_DAYS_PER_YEAR
  implicit none
  private

  public :: JulianCalendar
  public :: JULIAN_MONTHS ! for other calendars to see

  type, extends(FixedCalendar) :: JulianCalendar
   contains

  end type JulianCalendar

  type (CalendarMonth), parameter :: JULIAN_MONTHS(0:MONTHS_PER_YEAR+1) = &
       [ &
!                     Name        Abbr    #d  1st  Lst  Mid
       CalendarMonth('December',  'DEC ', 31, -30,   0, -15), &
       CalendarMonth('January',   'JAN ', 31, 001, 031, 016), &
       CalendarMonth('February',  'FEB ', 28, 032, 059, 045), &
       CalendarMonth('March',     'MAR ', 31, 060, 090, 075), &
       CalendarMonth('April',     'APR ', 30, 091, 120, 106), &
       CalendarMonth('May',       'MAY ', 31, 121, 151, 136), &
       CalendarMonth('June',      'JUN ', 30, 152, 181, 167), &
       CalendarMonth('July',      'JUL ', 31, 182, 212, 197), & 
       CalendarMonth('August',    'AUG ', 31, 213, 243, 228), &
       CalendarMonth('September', 'SEP ', 30, 244, 273, 259), &
       CalendarMonth('October',   'OCT ', 31, 274, 304, 289), &
       CalendarMonth('November',  'NOV ', 30, 305, 334, 320), &
       CalendarMonth('December',  'DEC ', 31, 335, 365, 350), &
       CalendarMonth('January',   'JAN ', 31, 366, 396, 381)  &
       ]

  ! Legacy support
  !         Legacy :  New
  public :: JDendOfM, LAST_JULIAN_DAY_IN_MONTH
  public :: JDmidOfM, MID_JULIAN_DAY_IN_MONTH

!@var LAST_JULIAN_DAY_IN_MONTH (JDendOfM, ) last Julian day in month
  integer, parameter :: LAST_JULIAN_DAY_IN_MONTH(0:MONTHS_PER_YEAR) = (/ &
       & 0,31,59,90,120,151,181,212,243,273,304,334,365 &
       & /)
  integer, parameter :: JDendOfM(0:MONTHS_PER_YEAR) = LAST_JULIAN_DAY_IN_MONTH
!@var MID_JULIAN_DAY_IN_MONTH(0:13) (JDmidOfM(0:13)) middle Julian day in month
  integer, parameter :: MID_JULIAN_DAY_IN_MONTH(0:MONTHS_PER_YEAR+1) = (/ &
       & -15,16,45,75,106,136,167,197,228,259,289,320,350,381 &
       & /)
  integer, parameter :: JDmidOfM(0:MONTHS_PER_YEAR+1) = MID_JULIAN_DAY_IN_MONTH

  ! Time is expressed as seconds since January 01 0h in BASE_YEAR
  integer, parameter :: BASE_YEAR = 1

  interface JulianCalendar              
     module procedure newJulianCalendar 
  end interface JulianCalendar

contains

  function newJulianCalendar() result(calendar)
    use CalendarDate_mod, only: CalendarDate
    type (JulianCalendar) :: calendar

    type (CalendarDate) :: birthday
    integer :: n

    call calendar%setDaysPerYear(INT_DAYS_PER_YEAR)
    call calendar%setSecondsPerDay(TimeInterval(INT_SECONDS_PER_DAY))
    
    do n = 0, MONTHS_PER_YEAR + 1
       call calendar%setNthCalendarMonth(n, JULIAN_MONTHS(n))
    end do

    
    call calendar%initTransitionDates()

    ! source <http://www.universetoday.com/12301/happy-birthday-johannes-kepler/>
    birthday = CalendarDate(month=12, date=27, year=1571)
    call calendar%addTransitionDate('Johannes Kepler Birthday', birthday)

  end function newJulianCalendar


end module JulianCalendar_mod
