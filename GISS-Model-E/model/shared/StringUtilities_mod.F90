module StringUtilities_mod
!@sum Container module for various procedures which manipulate strings
!@auth T.Clune
  implicit none
  private

  public :: toLowerCase
  public :: toUpperCase
  public :: stringToInteger
  public :: stringToRealDP
  public :: stringToLogical

  public :: toString
  public :: forceTrim

  interface toLowerCase
    module procedure toLowerCase_scalar
    module procedure toLowerCase_array
  end interface toLowerCase

  interface toUpperCase
    module procedure toUpperCase_scalar
    module procedure toUpperCase_array
 end interface toUpperCase

  interface toString
    module procedure toString_integer
    module procedure toString_logical
    module procedure toString_realDP
    module procedure toString_string
  end interface toString

  integer, parameter :: MAX_LEN_STRING=32

contains

  pure function toLowerCase_scalar(string) result(newString)
!@auth I. Aleinov
!@sum Produces a copy of the string argument, but with any
!@+ capitalized letters replaced by their lower case equivalent.
!@+ Characters which are not alphabetic letters are copied without
!@+  modification.
    character(len=*), intent(in) :: string
    character(len=len(string)) :: newString

    integer n, i
    integer A, Z, shift, c

    A = iachar( 'A' )
    Z = iachar( 'Z' )
    shift = iachar( 'a' ) - iachar( 'A' )

    newString = trim(string)
    n = len(trim(newString))
    do i=1,n
      c = iachar( newString(i:i) )
      if ( c>=A .and. c<=Z ) newString(i:i) = achar( c + shift )
    enddo

  end function toLowerCase_scalar

  pure function toUpperCase_scalar(string) result(newString)
!@auth T. Clune
!@sum Produces a copy of the string argument, but with any 
!@+ lower case letters replaced by their upper case equivalent.
!@+ Characters which are not alphabetic letters are copied without
!@+  modification.
    character(len=*), intent(in) :: string
    character(len=len(string)) :: newString

    integer n, i
    integer A, Z, shift, c

    a = iachar( 'a' )
    z = iachar( 'z' )
    shift = iachar( 'a' ) - iachar( 'A' )

    newString = trim(string)
    n = len(trim(newString))
    do i=1,n
      c = iachar( newString(i:i) )
      if ( c>=a .and. c<=z ) newString(i:i) = achar( c - shift )
    enddo

 end function toUpperCase_scalar


  pure function toLowerCase_array(string) result(newString)
    character(len=*), intent(in) :: string(:)
    character(len=len(string(1))) :: newString(size(string))

    integer :: i
    do i = 1, size(string)
      newString(i) = trim(toLowerCase(string(i)))
    end do
  end function toLowerCase_array

  pure function toUpperCase_array(string) result(newString)
    character(len=*), intent(in) :: string(:)
    character(len=len(string(1))) :: newString(size(string))

    integer :: i
    do i = 1, size(string)
      newString(i) = trim(toUpperCase(string(i)))
    end do
  end function toUpperCase_array

  ! The following procedures are elemental and thus cannot throw exceptions
  ! Input values must be checked before calling.
  elemental integer function stringToInteger(string) result(value)
    character(len=*), intent(in) :: string

    read(string,'(i20)') value

  end function stringToInteger

  elemental real*8 function stringToRealDP(string) result(value)
    character(len=*), intent(in) :: string

      read(string,*) value

  end function stringToRealDP

  elemental logical function stringToLogical(string) result(flag)
    character(len=*), intent(in) :: string

    select case(trim(toLowerCase(string)))
    case ('true','t','.true.')
      flag = .true.
    case ('false','f','.false.')
      flag = .false.
    case default
      ! throw exception
    end select
      
  end function stringToLogical

  elemental function toString_logical(value) result(string)
    logical, intent(in) :: value
    character(len=MAX_LEN_STRING) :: string

    if (value) then
      write(string,'(a)') '.true.'
    else
      write(string,'(a)') '.false.'
    end if
  end function toString_logical

  elemental function toString_integer(value) result(string)
    integer, intent(in) :: value
    character(len=MAX_LEN_STRING) :: string

    write(string,'(i0)') value
    
  end function toString_integer

  elemental function toString_realDP(value) result(string)
    real(kind=8), intent(in) :: value
    character(len=MAX_LEN_STRING) :: string

    write(string,'(g25.17)') value
    
  end function toString_realDP
    
  elemental function toString_string(value) result(string)
    character(len=*), intent(in) :: value
    character(len=MAX_LEN_STRING) :: string

    string = trim(value)
  end function toString_string

  elemental subroutine forceTrim(string, length)
    character(len=*), intent(inout) :: string
    integer, intent(in) :: length
    string = string(1:length)
  end subroutine forceTrim

end module StringUtilities_mod
