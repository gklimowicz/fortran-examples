! The Time class extends BaseTime by associating with a Calenar
! object.  Because the association is through a pointer, it is 
! crucial that clients do not destroy the calendar subsequent to
! creating the Time object.

module Time_mod
  use AbstractCalendar_mod, only: AbstractCalendar
  use Rational_mod
  use BaseTime_mod
  implicit none
  private

  public :: Time
  public :: newTime

  type, extends(BaseTime) :: Time
    class (AbstractCalendar), pointer :: calendar => null()
  contains
    procedure :: setByDate
    procedure :: setBaseTime

    procedure :: getYear
    procedure :: getMonth
    procedure :: getAbbreviation
    procedure :: getDate
    procedure :: getDayOfyear
    procedure :: getHour

    procedure :: add

    generic :: set => setByDate, setBaseTime
  end type Time

  interface tTime
     module procedure newTime
  end interface tTime

contains

  !----------------------------------------------------------------
  ! Creates a new time object - associating with a calendar.
  ! Initial value of time remains undefined until called with a 
  ! set() or setByDate() call.   
  !----------------------------------------------------------------
  function newTime(calendar) result(t)
    use AbstractCalendar_mod, only: AbstractCalendar
    type (Time) :: t
    class (AbstractCalendar), target, intent(in) :: calendar

    t%calendar => calendar
    
  end function newTime
  
  !----------------------------------------------------------------
  ! The following method sets the time using human (calendar) units
  ! of measurement.  Currently only times that are "on-the-hour" can be 
  ! set in this manner.   Finer-grained adjustments can be made with 
  ! arithmetic on times.
  !----------------------------------------------------------------
  subroutine setByDate(this, year, month, date, hour)
    class (Time), intent(inout) :: this
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: date
    integer, intent(in) :: hour

    this%BaseTime = this%calendar%convertToTime(year, month, date, hour)

  end subroutine setByDate

  subroutine setBaseTime(this, t)
    class (Time), intent(inout) :: this
    type (BaseTime), intent(in) :: t
    this%BaseTime = t
  end subroutine setBaseTime


  !----------------------------------------------------------------
  ! Access to year, month, date, etc is done by passing self to the
  ! contained calendar.
  !----------------------------------------------------------------

  integer function getYear(this) result(year)
    class (Time), intent(in) :: this
    year = this%calendar%getYear(this)
  end function getYear

  integer function getMonth(this) result(month)
    class (Time), intent(in) :: this
    month = this%calendar%getMonth(this)
  end function getMonth

  integer function getDate(this) result(date)
    class (Time), intent(in) :: this
    date = this%calendar%getDate(this)
  end function getDate

  integer function getDayOfYear(this) result(dayOfYear)
    class (Time), intent(in) :: this
    dayOfYear = this%calendar%getDayOfYear(this)
  end function getDayOfYear

  integer function getHour(this) result(hour)
    class (Time), intent(in) :: this
    hour = this%calendar%getHour(this)
  end function getHour

  function getAbbreviation(this) result(monthAbbrev)
    use CalendarMonth_mod, only: LEN_MONTH_ABBREVIATION
    character(len=LEN_MONTH_ABBREVIATION) :: monthAbbrev
    class (Time), intent(in) :: this
    monthAbbrev = this%calendar%getAbbrev(this)
  end function getAbbreviation


  subroutine add(this, dt)
     use BaseTime_mod, only: BaseTime
     class (Time), intent(inout) :: this
     class (Rational), intent(in) :: dt

     this%BaseTime = newBaseTime(this%BaseTime + dt)

  end subroutine add

end module Time_mod
