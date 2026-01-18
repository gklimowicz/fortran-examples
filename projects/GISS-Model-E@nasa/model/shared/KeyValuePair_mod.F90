module KeyValuePair_mod
!@sum This module provides an interface which couples 
!@+ "keys" with an array of generics.   Such entites are
!@+ the building blocks of associative arrays (also called
!@+ "dictionaries", "maps", or "hash maps".
!@auth T. Clune
  use GenericType_mod
  implicit none
  private

  public :: KeyValuePair_type ! data type
  public :: KeyValuePair      ! constructor
  public :: clean             ! destructor
  ! accessors
  public :: getKey, getKeys
  public :: getValue, getValues
  public :: getNumValues  

  ! File access
  public :: writeUnformatted
  public :: readUnformatted

  ! Test equality of a pair.
  public :: operator(==)

  ! Parameters
  public :: MAX_LEN_KEY
  integer, parameter :: MAX_LEN_KEY = 32


  type KeyValuePair_type
    private
    character(len=MAX_LEN_KEY) :: key
    type (GenericType_type), pointer :: values(:) => null()
  contains
    procedure :: print
  end type KeyValuePair_type

  interface KeyValuePair
    module procedure KeyValuePair_scalar
    module procedure KeyValuePair_array
    module procedure KeyValuePair_copy
  end interface

  interface clean
    module procedure cleanKeyValuePair
  end interface

  interface getValue
    module procedure getValue_1
    module procedure getValue_i
  end interface

  interface operator(==)
    module procedure equals
  end interface

  interface readUnformatted
    module procedure readUnformatted_pair
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_pair
  end interface

  interface getKeys
    module procedure getKeys_array
  end interface

contains

  function KeyValuePair_scalar(key, value) result(pair)
!@sum Constructor for scalar value
    character(len=*), intent(in) :: key
    type (GenericType_type), intent(in) :: value
    type (KeyValuePair_type) :: pair
    
    pair%key = key
    allocate(pair%values(1))
    pair%values(1) = value
  end function KeyValuePair_scalar

  function KeyValuePair_array(key, values) result(pair)
!@sum Constructor for vector of values
    character(len=*), intent(in) :: key
    type (GenericType_type), intent(in) :: values(:)
    type (KeyValuePair_type) :: pair
    
    pair%key = key
    allocate(pair%values(size(values)))
    pair%values = values
  end function KeyValuePair_array

  function KeyValuePair_copy(original) result(pair)
!@sum Copy constructor
    type (KeyValuePair_type), intent(in) :: original
    type (KeyValuePair_type) :: pair
    
    pair%key = original%key
    allocate(pair%values(size(original%values)))
    pair%values = original%values
  end function KeyValuePair_copy

  ! Accessors
  function getKey(this) result(key)
    type (KeyValuePair_type) :: this
    character(len=MAX_LEN_KEY) :: key
    key = trim(this%key)
  end function getKey

  function getKeys_array(this) result(keys)
    type (KeyValuePair_type), target :: this(:)
    character(len=MAX_LEN_KEY), pointer :: keys(:)
#ifdef COMPILER_G95
!TODO GFortran breaks with a true pointer assignment.
    allocate(keys(size(this)))
    keys = this(:)%key
#else
    keys => this(:)%key
#endif
  end function getKeys_array

  function getValue_1(this) result(value)
    type (KeyValuePair_type), intent(in) :: this
    type (GenericType_type), pointer :: value
    value => this%values(1)
  end function getValue_1

  function getValues(this) result(values)
    type (KeyValuePair_type), intent(in) :: this
    type (GenericType_type), pointer :: values(:)
    values => this%values(:)
  end function getValues

  function getValue_i(this, ith) result(value)
!@ Get ith value
    type (KeyValuePair_type), intent(in) :: this
    integer, intent(in) :: ith
    type (GenericType_type), pointer :: value
    if (ith < 0 .or. ith > getNumValues(this)) then
      call stop_model('KeyValuePair_mod::getValue() - argument "ith" out of range.',14)
    end if
    value => this%values(ith)
  end function getValue_i

  function getNumValues(this) result(numValues)
    type (KeyValuePair_type), intent(in) :: this
    integer :: numValues
    numValues = size(this%values)
  end function getNumValues

  subroutine readUnformatted_pair(this, unit)
!@sum Read a KeyValuePair object from a unit attached to
!@+ an unformatted, sequential file.
    type (KeyValuePair_type), intent(out) :: this
    integer, intent(in) :: unit

    character (len=MAX_LEN_KEY) :: key
    type (GenericType_type), allocatable :: values(:)
    integer :: i, n

    read(unit) key, n
    allocate(values(n))
    do i = 1, n
      call readUnformatted(values(i), unit)
    end do
    this = KeyValuePair(key, values)
  end subroutine readUnformatted_pair

  subroutine writeUnformatted_pair(this, unit)
!@sum Write a KeyValuePair object to a unit attached to
!@+ an unformatted, sequential file.
    type (KeyValuePair_type), intent(in) :: this
    integer, intent(in) :: unit

    integer :: i, n

    n = getNumValues(this)
    write(unit) this%key, n
    do i = 1, n
      call writeUnformatted(this%values(i), unit)
    end do

  end subroutine writeUnformatted_pair

  logical function check(this, valueType, numValues)
!@sum Return true if object "this" has specified type
!@+ and number of values.   This utility function is used
!@+ to check conformance in later subroutines.
    type (KeyValuePair_type), intent(in) :: this
    integer, intent(in) :: valueType
    integer, intent(in) :: numValues

    check = all(getType(this%values) == valueType)
    if (.not. check) then
      call stop_model('Incorrect type for specified key: <' &
           & // trim(this%key) // '>', 14)
      return
    end if

    check = (numValues == size(this%values))
    if (.not. check) then
      call stop_model('Incorrect number of elements for specified key: <' &
           & // trim(this%key) // '>', 14)
      return
    end if

  end function check

  logical function equals(pairA, pairB) result(isEqual)
!@sum Return true if two pairs have identical contents.
    type (KeyValuePair_type), intent(in) :: pairA
    type (KeyValuePair_type), intent(in) :: pairB
    
    isEqual = .true.

    if (trim(getKey(pairA)) /= trim(getKey(pairB))) then
      isEqual = .false.
      return
    end if

    if (getNumValues(pairA) /= getNumValues(pairB)) then
      isEqual = .false.
      return
    end if

    if (.not. all(getValues(pairA)  == getValues(pairB))) then
      isEqual = .false.
      return
    end if

  end function equals

  subroutine cleanKeyValuePair(this)
    type (KeyValuePair_type), intent(inout) :: this
    deallocate(this%values)
  end subroutine cleanKeyValuePair

  subroutine print(this)
    class (KeyValuePair_type), intent(in) :: this
    integer :: i
    print*,__LINE__,__FILE__,' key: <',trim(this%key),'>'
    print*,'------------------------'
    do i = 1, size(this%values)
      call this%values(i)%print()
    end do
    print*,'------------------------'
  end subroutine print
    
end module KeyValuePair_mod
