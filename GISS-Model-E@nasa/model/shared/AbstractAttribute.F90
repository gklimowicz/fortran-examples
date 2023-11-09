module AbstractAttribute_mod
  implicit none
  private

  public :: AbstractAttribute
  public :: AttributePointer
  public :: DP
  public :: SP
  public :: MAX_LEN_ATTRIBUTE_STRING
  public :: MAX_LEN_LINE

  type, abstract :: AbstractAttribute
     integer :: placeholder
  contains
    procedure(equals), deferred :: equals
    ! TODO ifort stalls on higher modules if operator(==) is introduced
!!$    generic :: operator(==) => equals
    procedure(print), deferred :: print
    procedure(toString), deferred :: toString
    procedure(writeUnformatted), deferred :: writeUnformatted
    procedure(readUnformatted), deferred :: readUnformatted
    procedure(clean), deferred :: clean
  end type AbstractAttribute

  type AttributePointer
    class (AbstractAttribute), pointer :: p => null()
  end type AttributePointer

  integer, parameter :: DP = selected_real_kind(14)
  integer, parameter :: SP = selected_real_kind(6)
  integer, parameter :: MAX_LEN_ATTRIBUTE_STRING = 80
  integer, parameter :: MAX_LEN_LINE = 1000

  abstract interface

    logical function equals(this, b)
      import AbstractAttribute
      class (AbstractAttribute), intent(in) :: this
      class (AbstractAttribute), intent(in) :: b
    end function equals

    function toString(this) result(string)
      import AbstractAttribute, MAX_LEN_LINE
      class (AbstractAttribute), intent(in) :: this
      character(len=MAX_LEN_LINE) :: string
    end function toString

    subroutine print(this)
      import AbstractAttribute
      class (AbstractAttribute), intent(in) :: this
    end subroutine print

    subroutine writeUnformatted(this, unit)
      import AbstractAttribute
      class (AbstractAttribute), intent(in) :: this
      integer, intent(in) :: unit
    end subroutine writeUnformatted

    function readUnformatted(this, unit) result(new)
      import AbstractAttribute
      class (AbstractAttribute), intent(in) :: this
      integer, intent(in) :: unit
      class (AbstractAttribute), pointer :: new
    end function readUnformatted

    subroutine clean(this)
      import AbstractAttribute
      class (AbstractAttribute), intent(inout) :: this
    end subroutine clean

  end interface

end module AbstractAttribute_mod
