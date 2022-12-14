module CalendarMonth_mod
  implicit none
  private

  public :: CalendarMonth
  public :: newCalendarMonth
  integer, parameter, public :: LEN_MONTH_ABBREVIATION = 4
  integer, parameter, public :: MAX_LEN_MONTH_NAME = 24

  ! Data struct - not a class
  type CalendarMonth
     character(len=MAX_LEN_MONTH_NAME) :: fullName
     character(len=LEN_MONTH_ABBREVIATION) :: abbreviation
     integer :: daysInMonth
     integer :: firstDayInMonth ! from start of year
     integer :: lastDayInMonth  ! from start of year
     integer :: midDayInMonth   ! from start of year
   contains
     procedure :: print
  end type CalendarMonth

  integer, parameter :: DAYS_PER_WEEK = 7

contains

  function newCalendarMonth(name, daysInMonth, firstDayInMonth, midDayInMonth) result(month)
    use StringUtilities_mod, only: toUpperCase
    type (CalendarMonth) :: month
    character(len=*), intent(in) :: name
    integer, intent(in) :: daysInMonth
    integer, intent(in) :: firstDayInMonth
    integer, intent(in) :: midDayInMonth

    month%fullName = trim(name)
    month%abbreviation = toUpperCase(name(1:3)) // ' '

    month%daysInMonth = daysInMonth
    month%firstDayInMonth = firstDayInMonth
    month%lastDayInMonth = firstDayInMonth + daysInMonth - 1
    month%midDayInMonth = midDayInMonth
    
  end function newCalendarMonth

  subroutine print(this, unit, year)
    class (CalendarMonth), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(in) :: year

    integer :: numDays
    integer :: numWeeks
    integer :: week, firstDay, lastDay
    integer :: offset
    integer :: day
    character(len=40) :: fmt

    write(unit,'()')
    write(unit,'(1x,a,i5.0)') trim(this%fullName), year
    write(unit,'(32("-"))')
    write(unit,'(3x,7(1x,a,1x))') 'Su ', 'M ', 'Tu', 'W ', 'Th', 'F ', 'Sa'

    numDays = this%lastDayInMonth - this%firstDayInMonth + 1

    offset = mod(this%firstDayInMonth-1, DAYS_PER_WEEK)
    week = 1
    firstDay = 1
    lastDay = DAYS_PER_WEEK - offset

    numWeeks = 1 + ((this%lastDayInMonth-1)/DAYS_PER_WEEK - (this%firstDayInMonth-1)/DAYS_PER_WEEK)


    write(fmt,'("((",i0,"x),",i0,"(i3,1x))")') (3+4*offset), DAYS_PER_WEEK - offset
    write(unit, trim(fmt)) [(day, day = firstDay, lastDay)]

    do week = 2, numWeeks
       firstDay = 1 + (week-1)*DAYS_PER_WEEK - offset
       lastDay = min(numDays, firstDay + DAYS_PER_WEEK-1)
       write(unit,'(3x,7(i3,1x))') [(day, day=firstDay, lastDay)]
    end do
    write(unit,'(32("-"))')
    
  end subroutine print

end module CalendarMonth_mod
