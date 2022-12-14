!--------------------------------------------------------------------
!
! CalendarDate is a simple data structure that extends AnniversaryDate
! by adding in an integer year component. I.e. it is used to represent
! non-recurring calendar events at the granularity of a day. (See
! CalendarDateAndTime for more precise time stamps.)  This class is an
! extension of the (empty) AbstractTimeStamp class, which provides
! print() methods.
!
! The structure is merely provided as a convenience for calendar
! diagnostics.
!
!--------------------------------------------------------------------

module CalendarDate_mod
  use AbstractTimeStamp_mod
  use StringUtilities_mod, only: toLowerCase
  use AnniversaryDate_mod
  implicit none
  private

  public :: CalendarDate

  type, extends(AnniversaryDate) :: CalendarDate
     integer :: year ! ignored if recurring
   contains
     procedure :: toString
  end type CalendarDate

  interface CalendarDate
     module procedure newCalendarDate
  end interface CalendarDate

contains


  ! Constructor
  function newCalendarDate(month, date, year) result(cDate)
    type (CalendarDate) :: cDate
    integer, intent(in) :: month
    integer, intent(in) :: date
    integer, optional, intent(in) :: year

    cDate%month = month
    cDate%date = date
    cDate%year = year

  end function newCalendarDate


  ! Convert data to human-readable string
  subroutine toString(this, string)
    class (CalendarDate), intent(in) :: this
    character(len=*), intent(out) :: string

    write(string,'(i2.2,"-",i2.2,"-",i4.4)') this%month, this%date, this%year

  end subroutine toString


end module CalendarDate_mod

#define TYPE CalendarDate
#define USE_MODULE CalendarDate_mod
#define TYPE_NAME CalendarDate
#define HAS_PRINT

#include "AssociativeArrayTemplate.h"

#define VALUE_TYPE CalendarDate
#define ASSOCIATIVE_ARRAY_TYPE CalendarDateAssociativeArray
#undef ITERATOR_TYPE
#define HASH_TYPE CalendarDateHashMap

#include "HashMapTemplate.h"

