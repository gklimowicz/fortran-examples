!@sum This template produces a module that encapsulates a HashMap object.
!@+ An HashMap object stores elements formed by a combination of a key
!@+ value and a mapped value.
!@auth NCCS ASTG

#define IDENTITY(A) A
#define CONCAT(A,B) IDENTITY(A)IDENTITY(B)

#define DEFAULT_HASH_TABLE_SIZE 100

#define DEFAULT_ASSOCIATIVE_ARRAY_TYPE(type) CONCAT(type,AssociativeArray)

#define DEFAULT_HASH_TYPE(type) CONCAT(type,HashMap)
#define DEFAULT_ITERATOR_TYPE(type) CONCAT(type,Iterator)
#define DEFAULT_MODULE_NAME(type) CONCAT(type,_mod)
#define DEFAULT_CONSTRUCTOR(type) CONCAT(new,type)

#ifndef VALUE_MOD
#define VALUE_MOD DEFAULT_MODULE_NAME(VALUE_TYPE)
#endif

#ifndef ASSOCIATIVE_ARRAY_TYPE
#define ASSOCIATIVE_ARRAY_TYPE DEFAULT_ASSOCIATIVE_ARRAY_TYPE(VALUE_TYPE)
#endif

#ifndef ASSOCIATIVE_ARRAY_CONSTRUCTOR
#define ASSOCIATIVE_ARRAY_CONSTRUCTOR DEFAULT_CONSTRUCTOR(ASSOCIATIVE_ARRAY_TYPE)
#endif

#ifndef ASSOCIATIVE_ARRAY_MOD
#define ASSOCIATIVE_ARRAY_MOD DEFAULT_MODULE_NAME(ASSOCIATIVE_ARRAY_TYPE)
#endif

#ifndef ASSOCIATIVE_ARRAY_ITERATOR_TYPE
#define ASSOCIATIVE_ARRAY_ITERATOR_TYPE DEFAULT_ITERATOR_TYPE(ASSOCIATIVE_ARRAY_TYPE)
#endif

#ifndef HASH_TYPE
#define HASH_TYPE DEFAULT_HASH_TYPE(VALUE_TYPE)
#endif

#ifndef CONSTRUCTOR
#define CONSTRUCTOR DEFAULT_CONSTRUCTOR(HASH_TYPE)
#endif

#ifndef MODULE_NAME
#define MODULE_NAME DEFAULT_MODULE_NAME(HASH_TYPE)
#endif

#ifndef TYPE_NAME
#define TYPE_NAME VALUE_TYPE
#endif

#define STRINGIFY(str) '''str'''
  
#ifndef ITERATOR_TYPE
#define ITERATOR_TYPE DEFAULT_ITERATOR_TYPE(HASH_TYPE)
#endif

#ifndef REFERENCE_TYPE
#define REFERENCE_TYPE CONCAT(TYPE_NAME,Reference)
#endif

module MODULE_NAME
  use VALUE_MOD, only: TYPE_NAME
  use ASSOCIATIVE_ARRAY_MOD, only: Map => ASSOCIATIVE_ARRAY_TYPE
  use ASSOCIATIVE_ARRAY_MOD, only: MAX_LEN_KEY
  use ASSOCIATIVE_ARRAY_MOD, only: MapIterator => ASSOCIATIVE_ARRAY_ITERATOR_TYPE
  use ASSOCIATIVE_ARRAY_MOD, only: MapConstructor => ASSOCIATIVE_ARRAY_CONSTRUCTOR
  use ASSOCIATIVE_ARRAY_MOD, only: operator(==), operator(/=)

#ifdef WRAPPED_TYPE
  use ASSOCIATIVE_ARRAY_MOD, only: REFERENCE_TYPE
#endif
  implicit none
  private

  public :: HASH_TYPE
  public :: CONSTRUCTOR
  public :: ITERATOR_TYPE
!  public :: assignment(=)
  public :: operator(/=)
  public :: operator(==)
  public :: clean
#ifdef WRAPPED_TYPE
  public :: REFERENCE_TYPE
#endif
  public :: MAX_LEN_KEY

!!$  integer, parameter :: MAX_LEN_KEY = 32
  integer, parameter :: DONE = -1

! A hash map is an associative container that stores elements formed by a 
! combination of a key value and a mapped value. The (key,value) structure
! is managed by the AssociateArray object (Map type) and the actual mapped
! value is implemented in a hashFunction

  type HASH_TYPE
!!$    private
    integer :: tableSize = -1
    type (Map), allocatable :: table(:)
  contains
    procedure :: hashFunction
    procedure :: size => getSize
    procedure :: getReference
    procedure :: setValue
    procedure :: insertEntry
    generic :: insert => insertEntry
    procedure :: merge
    procedure :: has => hasIt
    procedure :: insertReference
    procedure :: print
    ! iterator operations
    procedure :: begin
    procedure :: last
  end type HASH_TYPE

! iterators are used to access the sequence of HashMap elements.
  type :: ITERATOR_TYPE
!!$    private
    class (HASH_TYPE), pointer :: reference => null()
    integer :: hashValue = 0
    type (MapIterator) :: subIterator
  contains
    procedure :: hasNext
    procedure :: next
    procedure :: key
    procedure :: value
    procedure :: copyIter
    generic, public :: assignment(=) => copyIter
  end type ITERATOR_TYPE

  interface clean
    module procedure clean_container
    module procedure clean_iterator
  end interface clean

  interface operator(/=)
    module procedure notEqual
  end interface operator(/=)

  interface operator(==)
    module procedure equal
  end interface operator(==)

contains

  function CONSTRUCTOR(hashTableSize) result(dictionary)
    integer, optional :: hashTableSize
    type (HASH_TYPE) :: dictionary

    integer :: hashTableSize_
    integer :: i

    hashTableSize_ = DEFAULT_HASH_TABLE_SIZE
    if (present(hashTableSize)) hashTableSize_ = hashTableSize
    
    allocate(dictionary%table(hashTableSize_))
    dictionary%tableSize = hashTableSize_

    do i = 1, hashTableSize_
      dictionary%table(i) = MapConstructor()
    end do

  end function CONSTRUCTOR

  function ITERATOR_CONSTRUCTOR (dictionary) result(iterator)
    type (HASH_TYPE), target :: dictionary
    type (ITERATOR_TYPE), pointer :: iterator

    allocate(iterator)
    iterator%reference => dictionary
    !The following makes a copy
    !allocate(iterator%reference, SOURCE=dictionary)
  end function ITERATOR_CONSTRUCTOR

  integer function getSize(this) 
    class (HASH_TYPE), intent(in) :: this

    integer :: i

    getSize = 0
    do i = 1, this%tableSize
      getSize = getSize + this%table(i)%size()
    end do

  end function getSize

  subroutine setValue(this, key, value)
    use StringUtilities_mod, only: toLowerCase
    class (HASH_TYPE), target, intent(inout) :: this
    character(len=*), intent(in) :: key
    class (TYPE_NAME) :: value

    integer :: hashValue

    hashValue = this%hashFunction(toLowerCase(key))
    call this%table(hashValue)%insert(key, value)

  end subroutine setValue

  subroutine insertEntry(this, key, value)
    use StringUtilities_mod, only: toLowerCase
    class (HASH_TYPE), target, intent(inout) :: this
    character(len=*), intent(in) :: key
    class (TYPE_NAME) :: value
    class (TYPE_NAME), pointer :: p

    integer :: hashValue

    hashValue = this%hashFunction(toLowerCase(key))
    call this%table(hashValue)%insert(key, value)
    
  end subroutine insertEntry

  subroutine insertReference(this, key, value)
    use StringUtilities_mod, only: toLowerCase
    class (HASH_TYPE), target, intent(inout) :: this
    character(len=*), intent(in) :: key
    class (TYPE_NAME), target :: value

    integer :: hashValue
    type (Map), pointer :: m

    hashValue = this%hashFunction(toLowerCase(key))
    m => this%table(hashValue)
    call m%insertReference(key, value)

  end subroutine insertReference

#ifdef WRAPPED_TYPE
  function getReference(this, key) result(ref)
    use StringUtilities_mod, only: toLowerCase
    class (HASH_TYPE), intent(in) :: this
    character(len=*), intent(in) :: key
    type (REFERENCE_TYPE) :: ref

    integer :: hashValue
    character(len=len(key)) lowerCaseKey

    lowerCaseKey = trim(toLowerCase(key))
    hashValue = this%hashFunction(lowerCaseKey)
    ref = this%table(hashValue)%getReference(lowerCaseKey)

  end function getReference
#else
  function getReference(this, key) result(ptr)
    use StringUtilities_mod, only: toLowerCase
    class (HASH_TYPE), intent(in) :: this
    character(len=*), intent(in) :: key
    class (TYPE_NAME), pointer :: ptr

    integer :: hashValue
    character(len=len(key)) lowerCaseKey

    lowerCaseKey = trim(toLowerCase(key))
    hashValue = this%hashFunction(lowerCaseKey)
    ptr => this%table(hashValue)%getReference(lowerCaseKey)

  end function getReference
#endif

  logical function hasIt(this, key)
    use StringUtilities_mod, only: toLowerCase
    class (HASH_TYPE), intent(in) :: this
    character(len=*), intent(in) :: key
    
    integer :: hashValue

    hashValue = this%hashFunction(key)
    hasIt = this%table(hashValue)%has(key)

  end function hasIt

 subroutine copyIter(a, b)
     class (ITERATOR_TYPE), intent(inout) :: a
     class (ITERATOR_TYPE), intent(in) :: b
!     type (ITERATOR_TYPE), intent(inout) :: a
!     type (ITERATOR_TYPE), intent(in) :: b

     a%reference => b%reference
     a%hashValue = b%hashValue
     a%subIterator = b%subIterator
  end subroutine copyIter

  subroutine merge(this, b)
    class (HASH_TYPE), intent(inout) :: this
    class (HASH_TYPE), intent(in) :: b

    type (ITERATOR_TYPE) :: iter
    class (TYPE_NAME), pointer :: t

    iter = b%begin()
    do while (iter /= b%last())
       if (.not. this%has(iter%key())) then
	 t => iter%value()
          call this%insert(iter%key(), t)
       end if
      call iter%next()
    end do

  end subroutine merge

  subroutine print(this)
    class (HASH_TYPE), intent(in) :: this
    type (ITERATOR_TYPE) :: iter
    class (TYPE_NAME), pointer :: t

    print*,'--------------------------'
    print*,' AssociativeArray: '
    print*,'--------------------------'

    iter = this%begin()
    do while (iter /= this%last())
      print*,'   key: <',trim(iter%key()),'>'
#ifdef HAS_PRINT
      t => iter%value()
      call t%print()
#endif
      call iter%next()
    end do
    print*,'--------------------------'
    print*,'--------------------------'
    print*,' '

  end subroutine print

  type (ITERATOR_TYPE) function begin(this) result(iterator)
    class (HASH_TYPE), target, intent(in) :: this

    integer :: i

    iterator%reference => this

    do i = 1, this%tableSize
      if (this%table(i)%size() > 0) then
        iterator%hashValue = i
        iterator%subIterator = this%table(i)%begin()
        return
      end if
    end do
    
    iterator%hashValue = DONE ! no entries

  end function begin

  type (ITERATOR_TYPE) function last(this) result(iterator)
    class (HASH_TYPE), target, intent(in) :: this

    iterator%reference => this
    iterator%hashValue = DONE ! no entries

  end function last

  logical function notEqual(a, b)
    class (ITERATOR_TYPE), intent(in) :: a
    class (ITERATOR_TYPE), intent(in) :: b
    ! TODO: throw exception if a and b do not have the same reference
    notEqual = .not. (a == b)
  end function notEqual

  logical function equal(a, b)
    class (ITERATOR_TYPE), intent(in) :: a
    class (ITERATOR_TYPE), intent(in) :: b

    ! TODO: throw exception if a and b do not have the same reference

    if (a%hashValue /= b%hashValue) then
      equal = .false.
      return
    end if
    if (a%hashValue == DONE) then
      equal = .true.
    else
      equal = (a%subIterator == b%subIterator)
    end if
  end function equal

  logical function hasNext(this)
    class (ITERATOR_TYPE), intent(in) :: this

    integer :: i

    if (this%hashValue == DONE) then
      hasNext = .false.
      return
    end if

    if (this%subIterator%hasNext()) then
      hasNext = .true.
      return
    else
      do i = this%hashValue + 1, this%reference%tableSize
        if (this%reference%table(i)%size() > 0) then
          hasNext = .true.
          return
        end if
      end do
    end if

    hasNext = .false.

  end function hasNext

  subroutine next(this)
    class (ITERATOR_TYPE), intent(inout) :: this

    integer :: i

    if (this%hashValue == DONE) then
      call stop_model('Cannot call next() when hash has no remaining iterations.',14)
      return
    end if

    if (this%subIterator%hasNext()) then
      call this%subIterator%next()
      return
    else
      do i = this%hashValue + 1, this%reference%tableSize
        if (this%reference%table(i)%size() > 0) then
          this%hashValue = i
          this%subIterator = this%reference%table(i)%begin()
          return
        end if
      end do
    end if

    this%hashValue = DONE

  end subroutine next

  function key(this)
    class (ITERATOR_TYPE), target, intent(in) :: this
    character(len=MAX_LEN_KEY), pointer :: key
    type (MapIterator), pointer :: iter
    iter => this%subIterator
    call iter%getKey(key)
  end function key

  function value(this)
    class (ITERATOR_TYPE), target, intent(in) :: this
    class (TYPE_NAME), pointer :: value
    value => this%subIterator%value()
  end function value

  integer function hashFunction(this, key) result(hashValue)
    use StringUtilities_mod, only: toLowerCase
    class (HASH_TYPE), intent(in) :: this
    character(len=*), intent(in) :: key

    integer :: i
    integer :: hashSum
    character(len=len(key)) :: lowerCaseKey

    lowerCaseKey = toLowerCase(key)

    hashSum = 0
    do i = 1, len_trim(lowerCaseKey)
      hashSum = hashSum + iachar(lowerCaseKey(i:i))
    end do

    hashValue = 1 + mod(hashSum-1, this%tableSize)

  end function hashFunction

  subroutine clean_container(this)
    type (HASH_TYPE), intent(inout) :: this
    deallocate(this%table)
  end subroutine clean_container

  subroutine clean_iterator(this)
    type (ITERATOR_TYPE), intent(inout) :: this
  end subroutine clean_iterator

end module MODULE_NAME

#undef VALUE_MOD
#undef ASSOCIATIVE_ARRAY_TYPE
#undef ASSOCIATIVE_ARRAY_CONSTRUCTOR
#undef ASSOCIATIVE_ARRAY_MOD
#undef ASSOCIATIVE_ARRAY_ITERATOR_TYPE
#undef HASH_TYPE
#undef CONSTRUCTOR
#undef MODULE_NAME
#undef TYPE_NAME
#undef TYPE
