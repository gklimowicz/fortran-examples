module GenericType_mod
!@sum This module provides a derived type which can contain one element
!@+ of any of the following intrinsic types: integer, real*8, logical, 
!@+ and string. This provides a flexible mechanism for organizing various
!@+ types of metadata.
!@auth T.Clune
  implicit none
  private

  public :: GenericType_type ! data type
  public :: GenericType      ! constructor
  public :: assignment(=)
  public :: operator(==)
  public :: toString
  public :: getType
  public :: readUnformatted
  public :: writeUnformatted

  ! parameters
  public :: MAX_LEN_VALUE
  public :: INTEGER_TYPE
  public :: REAL64_TYPE
  public :: LOGICAL_TYPE
  public :: STRING_TYPE

  integer, parameter :: MAX_LEN_VALUE = 128 ! longest permitted string value
  integer, parameter :: ILLEGAL_TYPE = -1
  integer, parameter :: INTEGER_TYPE = 1
  integer, parameter :: REAL64_TYPE  = 2
  integer, parameter :: LOGICAL_TYPE = 3
  integer, parameter :: STRING_TYPE  = 4


  type GenericType_type
    integer :: type = ILLEGAL_TYPE 
    ! Only one of the following has a legitimate value at any time.
    integer :: integerValue
    real*8  :: real64Value
    logical :: logicalvalue
    character(len=MAX_LEN_VALUE) :: stringValue
  contains
    procedure :: print
  end type GenericType_type

  ! Construct generic type from intrinsic data, read from
  ! a string, or copy from another generic object.
  interface GenericType ! constructors
    module procedure GenericType_copy
    module procedure GenericType_readString
    module procedure GenericType_integer
    module procedure GenericType_real64
    module procedure GenericType_logical
    module procedure GenericType_string
  end interface

  ! Overload assignment operator to facilitate extracting values of
  ! intrisic types.  Overload for 1D arrays as well.
  interface assignment(=)
    module procedure assignToInteger_sca_sca, assignToInteger_sca_arr
    module procedure assignToInteger_arr_sca, assignToInteger_arr_arr
    module procedure assignToReal64_sca_sca, assignToReal64_sca_arr
    module procedure assignToReal64_arr_sca, assignToReal64_arr_arr
    module procedure assignToLogical_sca_sca, assignToLogical_sca_arr
    module procedure assignToLogical_arr_sca, assignToLogical_arr_arr
    module procedure assignToString_sca_sca, assignToString_sca_arr
    module procedure assignToString_arr_sca, assignToString_arr_arr
  end interface

  ! This interface is provided mostly to facilitate testing.
  ! Overload operator "==" to facilitate testing.
  interface operator(==)
    module procedure equals_generic
    module procedure equals_integer
    module procedure equals_real64
    module procedure equals_logical
    module procedure equals_string
  end interface

  ! Return the integer parameter associated with generic.
  ! Overloaded to operate on intrinsic types.
  interface getType
    module procedure getType_generic
    module procedure getType_integer
    module procedure getType_real64
    module procedure getType_logical
    module procedure getType_string
  end interface

  ! Convert scalar (vector) of generics into scalar (vector) of
  ! strings.  Cannot use ELEMENTAL attribute because want to catch
  ! uninitialized variables which throws an exception.
  interface toString
    module procedure toString_single
    module procedure toString_multi
  end interface

  ! Read generic from a unit opened to a sequential, unformmatted file.
  interface readUnformatted
    module procedure readUnformatted_generic
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_generic
  end interface

  integer, parameter :: MAX_LEN_TYPE_STRING = 16

contains

  ! Begin constructors
  function GenericType_copy(original) result(generic)
    type (GenericType_type), intent(in) :: original
    type (GenericType_type) :: generic
    generic = original
  end function GenericType_copy

  elemental function GenericType_integer(value) result(generic)
    integer, intent(in) :: value
    type (GenericType_type) :: generic
    generic%type = INTEGER_TYPE
    generic%integerValue = value
  end function GenericType_integer

  elemental function GenericType_real64(value) result(generic)
    real*8, intent(in) :: value
    type (GenericType_type) :: generic
    generic%type = REAL64_TYPE
    generic%real64Value = value
  end function GenericType_real64

  elemental function GenericType_logical(value) result(generic)
    logical, intent(in) :: value
    type (GenericType_type) :: generic
    generic%type = LOGICAL_TYPE
    generic%logicalValue = value
  end function GenericType_logical

  elemental function GenericType_string(value) result(generic)
    character(len=*), intent(in) :: value
    type (GenericType_type) :: generic
    generic%type = STRING_TYPE
    generic%stringValue = value
  end function GenericType_string

  function GenericType_readString(string, type) result(generic)
!@sum Construct a generic from a string and a type specifier.
    use StringUtilities_mod, only: toLowerCase
    character(len=*), intent(in) :: string
    integer, intent(in) :: type
    type (GenericType_type) :: generic
    integer :: status

    generic%type = type

    select case (type)
    case (INTEGER_TYPE)
      read(string,'(i20)',iostat=status) generic%integerValue
    case (REAL64_TYPE)
      read(string,*,iostat=status) generic%real64Value
    case (LOGICAL_TYPE)
      generic%logicalValue = readLogical(toLowerCase(string), status)
    case (STRING_TYPE)
      generic%stringValue = string
      status = 0
    case default
      generic%type = ILLEGAL_TYPE
      call stop_model('GenericType::GenericType() - no such type.',14)
      return
    end select

    if (status /= 0) then
      call stop_model('GenericType::GenericType() - cannot convert string "' // &
           & trim(string) // '" to ' // trim(typeString(type)) // '.', 14)
    end if

  contains

    ! For convenience, allow multiple formats for logicals.
    logical function readLogical(string, status) result(flag)
      character(len=*), intent(in) :: string
      integer, intent(out) :: status
      status = 0
      select case(trim(toLowerCase(string)))
      case ('true','t','.true.')
        flag = .true.
      case ('false','f','.false.')
        flag = .false.
      case default
        status = -1 ! cannot convert
      end select
    end function readLogical

    function typeString(type)
      integer, intent(in) :: type

      character(len=MAX_LEN_TYPE_STRING) :: typeString

      select case (type)
      case (INTEGER_TYPE)
        typeString = 'integer'
      case (REAL64_TYPE)
        typeString = 'real*8'
      case (LOGICAL_TYPE)
        typeString = 'logical'
      case (STRING_TYPE)
        typeString = 'string'
      case default
        typeString = 'unknown type'
      end select
    end function typeString

  end function GenericType_readString
  ! End constructors

  ! Begin overload of assigment (=)
  subroutine assignToInteger_sca_sca(value, this)
    integer, intent(out) :: value
    type (GenericType_type), intent(in) :: this
    value = this%integerValue
  end subroutine assignToInteger_sca_sca

  subroutine assignToInteger_sca_arr(value, this)
    integer, intent(out) :: value
    type (GenericType_type), intent(in) :: this(:)

    if (nonconforming(1, size(this))) return
    value = this(1)%integerValue
  end subroutine assignToInteger_sca_arr

  subroutine assignToInteger_arr_arr(values, this)
    integer, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this(:)
    if (nonconforming(size(values), size(this))) return
    values(:) = this(:)%integerValue
  end subroutine assignToInteger_arr_arr

  subroutine assignToInteger_arr_sca(values, this)
    integer, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this
    if (nonconforming(size(values), 1)) return
    values(1) = this%integerValue
  end subroutine assignToInteger_arr_sca
  
  subroutine assignToReal64_sca_sca(value, this)
    real*8, intent(out) :: value
    type (GenericType_type), intent(in) :: this
    value = this%real64Value
  end subroutine assignToReal64_sca_sca

  subroutine assignToReal64_sca_arr(value, this)
    real*8, intent(out) :: value
    type (GenericType_type), intent(in) :: this(:)
    if (nonconforming(1, size(this))) return
    value = this(1)%real64Value
  end subroutine assignToReal64_sca_arr

  subroutine assignToReal64_arr_arr(values, this)
    real*8, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this(:)
    if (nonconforming(size(values), size(this))) return
    values = this(:)%real64Value
  end subroutine assignToReal64_arr_arr

  subroutine assignToReal64_arr_sca(values, this)
    real*8, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this
    if (nonconforming(size(values), 1)) return
    values(1) = this%real64Value
  end subroutine assignToReal64_arr_sca
  
  subroutine assignToLogical_sca_sca(value, this)
    logical, intent(out) :: value
    type (GenericType_type), intent(in) :: this
    value = this%logicalValue
  end subroutine assignToLogical_sca_sca

  subroutine assignToLogical_sca_arr(value, this)
    logical, intent(out) :: value
    type (GenericType_type), intent(in) :: this(:)
    if (nonconforming(1, size(this))) return
    value = this(1)%logicalValue
  end subroutine assignToLogical_sca_arr

  subroutine assignToLogical_arr_sca(values, this)
    logical, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this
    if (nonconforming(size(values), 1)) return
    values(1) = this%logicalValue
  end subroutine assignToLogical_arr_sca

  subroutine assignToLogical_arr_arr(values, this)
    logical, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this(:)
    if (nonconforming(size(values), size(this))) return
    values = this%logicalValue
  end subroutine assignToLogical_arr_arr

  subroutine assignToString_sca_arr(value, this)
    character(len=*), intent(out) :: value
    type (GenericType_type), intent(in) :: this(:)
    if (nonconforming(1, size(this))) return
    value = this(1)%stringValue
  end subroutine assignToString_sca_arr

  subroutine assignToString_sca_sca(value, this)
    character(len=*), intent(out) :: value
    type (GenericType_type), intent(in) :: this
    value = this%stringValue
  end subroutine assignToString_sca_sca

  subroutine assignToString_arr_sca(values, this)
    character(len=*), intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this
    if (nonconforming(size(values), 1)) return
    values(1) = this%stringValue
  end subroutine assignToString_arr_sca

  subroutine assignToString_arr_arr(values, this)
    character(len=*), intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this(:)
    if (nonconforming(size(values), size(this))) return
    values = this%stringValue
  end subroutine assignToString_arr_arr
  ! End of overload of assignment(=)

  logical function nonconforming(sizeA, sizeB)
!@sum Throws an exception if the sizes of two arrays/scalars are
!@+ not conforming.   
    integer, intent(in) :: sizeA
    integer, intent(in) :: sizeB

    if (sizeA /= sizeB) then
      call stop_model('GenericType_mod: nonconforming shapes.', 14)
      nonconforming = .true.
    else
      nonconforming = .false.
    end if

  end function nonconforming

  ! Begin overload of operator(==)

  elemental logical function equals_generic(genericA, genericB) result(same)
!@sum Test equality of two generics. Two generics are defined to be equal
!@+  iff they have the same type and the corresponding values are
!@+  the same by the criteria of that type.
    type (GenericType_type), intent(in) :: genericA
    type (GenericType_type), intent(in) :: genericB

    same = (genericA%type == genericB%type)

    if (same) then
      select case (genericB%type)
      case (INTEGER_TYPE)
        same = (genericA%integerValue == genericB%integerValue)
      case (REAL64_TYPE)
        same = (genericA%real64Value == genericB%real64Value)
      case (LOGICAL_TYPE)
        same = (genericA%logicalValue .eqv. genericB%logicalValue)
      case (STRING_TYPE)
        same = (genericA%stringValue == genericB%stringValue)
      end select
    end if

  end function equals_generic
  
  elemental logical function equals_integer(expected, generic) result(same)
    integer, intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (generic%type==INTEGER_TYPE) .and. (generic%integerValue==expected)
  end function equals_integer

  elemental logical function equals_real64(expected, generic) result(same)
    real*8, intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (generic%type == REAL64_TYPE) .and. (generic%real64Value==expected)
  end function equals_real64

  elemental logical function equals_logical(expected, generic) result(same)
    logical, intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (generic%type==LOGICAL_TYPE) .and. (generic%logicalValue .eqv. expected)
  end function equals_logical

  elemental logical function equals_string(expected, generic) result(same)
    character(len=*), intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (generic%type==STRING_TYPE) .and. (generic%stringValue == expected)
  end function equals_string
  ! End overload of operator(==)

  ! Overload of toString()

  ! Convert a generic into a human-readable string.  Useful for manipulating
  ! files containing metadata.
  function toString_single(this) result(string)
    type (GenericType_type), intent(in) :: this
    character(len=MAX_LEN_VALUE) :: string

    select case (this%type)
    case (INTEGER_TYPE)
      write(string, '(i0)') this%integerValue
    case (REAL64_TYPE)
      write(string, '(g25.17)') this%real64Value
    case (LOGICAL_TYPE)
      if (this%logicalValue .eqv. .true.) then
        write(string,'(a)') '.true.'
      else
        write(string,'(a)') '.false.'
      end if
    case (STRING_TYPE)
      string = trim(this%stringValue)
    case default
      call stop_model('GenericType::toString() - invalid type.',14)
    end select
  end function toString_single

  function toString_multi(this) result(string)
!@sum Convert an array of generics into an array of human-readable strings.
    type (GenericType_type), intent(in) :: this(:)
    character(len=MAX_LEN_VALUE) :: string(size(this))

    integer :: i

    do i = 1, size(this)
      string(i) = toString(this(i))
    end do

  end function toString_multi

  ! Various trivial accessor methods.
  elemental integer function getType_generic(this) result(type)
    type (GenericType_type), intent(in) :: this
    type = this%type
  end function getType_generic

  integer function getType_integer(value)
    integer, intent(in) :: value
    getType_integer = INTEGER_TYPE
  end function getType_integer

  integer function getType_real64(value)
    real*8, intent(in) :: value
    getType_real64 = REAL64_TYPE
  end function getType_real64

  integer function getType_logical(value)
    logical, intent(in) :: value
    getType_logical = LOGICAL_TYPE
  end function getType_logical

  integer function getType_string(value)
    character(len=*), intent(in) :: value
    getType_string = STRING_TYPE
  end function getType_string

  subroutine readUnformatted_generic(this, unit)
!@sum Read generic item from unit connected to an unformatted,
!@+ sequential file.
    type (GenericType_type), intent(out) :: this
    integer, intent(in) :: unit
    
    integer :: type
    integer :: integerValue
    real*8 :: real64Value
    logical :: logicalvalue
    character(len=MAX_LEN_VALUE) :: stringValue

    read(unit) type
    select case (type)
    case (INTEGER_TYPE)
      read(unit) integerValue
      this = GenericType(integerValue)
    case (REAL64_TYPE)
      read(unit) real64Value
      this = GenericType(real64Value)
    case (LOGICAL_TYPE)
      read(unit) logicalvalue
      this = GenericType(logicalvalue)
    case (STRING_TYPE)
      read(unit) stringValue
      this = GenericType(stringValue)
    case default
      call stop_model('GenericType_mod::readUnformatted() - unsupported type.', 14)
    end select

  end subroutine readUnformatted_generic

  subroutine writeUnformatted_generic(this, unit)
!@sum Write generic item to unit connected to an unformatted,
!@+ sequential file.
    type (GenericType_type), intent(in) :: this
    integer, intent(in) :: unit
    
    integer :: type

    write(unit) getType(this)
    select case (getType(this))
    case (INTEGER_TYPE)
      write(unit) this%integerValue
    case (REAL64_TYPE)
      write(unit) this%real64Value
    case (LOGICAL_TYPE)
      write(unit) this%logicalValue
    case (STRING_TYPE)
      write(unit) this%stringValue
    case default
      call stop_model('GenericType_mod::writeUnformatted() - unsupported type.', 14)
    end select

  end subroutine writeUnformatted_generic

  subroutine print(this)
    class (GenericType_type), intent(in):: this
    select case (this%type)
    case (INTEGER_TYPE)
      print*,__LINE__,__FILE__,' Integer: ', this%integerValue
    case (LOGICAL_TYPE)
      print*,__LINE__,__FILE__,' Logical: ', this%logicalValue
    case (REAL64_TYPE)
      print*,__LINE__,__FILE__,' Real64: ', this%real64Value
    case (STRING_TYPE)
      print*,__LINE__,__FILE__,' String: ', trim(this%stringValue)
    end select
  end subroutine print

end module GenericType_mod
