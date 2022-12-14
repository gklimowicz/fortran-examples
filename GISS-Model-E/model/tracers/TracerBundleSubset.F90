module TracerPointer_mod
!@sum This module encapsulates the TracerPointer data type. This type
!@+ contains a reference (pointer) to a tracer data type as well as access
!@+ methods (get, set)
!@auth NCCS ASTG

  use Tracer_mod
  use Attributes_mod
  use AttributeHashMap_mod
  use TracerHashMap_mod
  implicit none
  private

  public :: TracerPointer
!!$  public :: newTracerPointer
!!$  public :: getTracerRefName

  type, extends(Tracer) :: TracerPointer
    private
    class(Tracer), pointer :: reference
  contains
    procedure :: set => setTracerPointer
    procedure :: get => getTracerPointer
    procedure :: print
    procedure :: size

!!$    procedure :: equals

!!$    procedure :: writeUnformatted
!!$    procedure :: readUnformatted
!!$    procedure :: clean

!!$    procedure :: merge
    procedure :: insertIntegerAttribute
    procedure :: insertInteger1dAttribute
    procedure :: insertLogical1dAttribute
    procedure :: insertLogicalAttribute
    procedure :: insertRealDP1dAttribute
    procedure :: insertRealDPAttribute
    procedure :: insertString1dAttribute
    procedure :: insertStringAttribute
  end type TracerPointer

  interface TracerPointer
    module procedure newTracerPointer
  end interface

  interface assignment(=)
    module procedure copyIt
  end interface assignment(=)

  ! TODO make type-bound when ifort is fixed - compilation slows to craws for higher-level modules
!!$  interface operator(==)
!!$    module procedure equals
!!$  end interface operator(==)

  interface clean
    module procedure clean_
  end interface clean

contains

  function newTracerPointer(tracerTarget)
    class (Tracer), target :: tracerTarget
    type (TracerPointer) :: newTracerPointer

    newTracerPointer%Tracer = newTracer()
    newTracerPointer%reference => tracerTarget

  end function newTracerPointer

  subroutine setTracerPointer(this, reference)
    class (TracerPointer), intent(inout) :: this
    class (Tracer), target :: reference

    this%reference => reference

  end subroutine setTracerPointer

  function getTracerPointer(this) result(reference)
    class (TracerPointer), intent(in) :: this
    class (Tracer), pointer :: reference
 
    reference  => this%reference

  end function getTracerPointer

  subroutine print(this)
    class (TracerPointer), intent(in) :: this
    type (AttributeHashMapIterator) :: iter
    class (Tracer), pointer :: t

    print*,'--------------------------'
    print*,' TracerPointer: '
    print*,'--------------------------'

    call this%reference%print()
!!$    iter = this%reference%begin()
!!$    do while (iter /= this%reference%last())
!!$      print*,'   key: <',trim(iter%key()),'>'
!!$      t => iter%value()
!!$      call t%print()
!!$      call iter%next()
!!$    end do
    print*,'--------------------------'
    print*,'--------------------------'
    print*,' '

  end subroutine print

  integer function size(this) 
    class (TracerPointer), intent(in) :: this
    size = this%reference%size()
  end function size

!!$  subroutine merge(this, b)
!!$    class (TracerPointer), intent(inout) :: this
!!$    class (TracerPointer), intent(in) :: b
!!$
!!$    type (AttributeHashMapIterator) :: iter
!!$
!!$    iter = b%begin()
!!$    do while (iter /= b%last())
!!$       if (.not. this%has(iter%key())) then
!!$          call this%AttributeDictionary%insert(iter%key(), iter%value())
!!$       else
!!$          call stop_model('AssociativeArray::merge() failed due to duplicate keys: <' &
!!$               & // trim(iter%key()) // '>.', 255)
!!$       end if
!!$      call iter%next()
!!$    end do
!!$
!!$  end subroutine merge

  subroutine insertIntegerAttribute(this, key, value)
    class (TracerPointer), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value 

    call this%reference%insert(key, newAttribute(value))
  end subroutine insertIntegerAttribute

  subroutine insertInteger1dAttribute(this, key, value)
    class (TracerPointer), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value (:)

    call this%reference%insert(key, newAttribute(value))
  end subroutine insertInteger1dAttribute

  subroutine insertLogical1dAttribute(this, key, value)
    class (TracerPointer), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value (:)

    call this%reference%insert(key, newAttribute(value))
  end subroutine insertLogical1dAttribute

  subroutine insertLogicalAttribute(this, key, value)
    class (TracerPointer), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value 

    call this%reference%insert(key, newAttribute(value))
  end subroutine insertLogicalAttribute

  subroutine insertRealDP1dAttribute(this, key, value)
    class (TracerPointer), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(kind=DP), intent(in) :: value (:)

    call this%reference%insert(key, newAttribute(value))
  end subroutine insertRealDP1dAttribute

  subroutine insertRealDPAttribute(this, key, value)
    class (TracerPointer), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(kind=DP), intent(in) :: value 

    call this%reference%insert(key, newAttribute(value))
  end subroutine insertRealDPAttribute

  subroutine insertString1dAttribute(this, key, value)
    class (TracerPointer), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value (:)

    call this%reference%insert(key, newAttribute(value))
  end subroutine insertString1dAttribute

  subroutine insertStringAttribute(this, key, value)
    class (TracerPointer), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value 

    call this%reference%insert(key, newAttribute(value))
  end subroutine insertStringAttribute


!!$  logical function equals(this, b)
!!$    class (TracerPointer), intent(in) :: this
!!$    class (TracerPointer), intent(in) :: b
!!$
!!$    type (AttributeHashMapIterator) :: iter
!!$    class (AbstractAttribute), pointer :: p1, p2
!!$
!!$    equals = .true.
!!$    if (this%reference%size() /= b%reference%size()) then
!!$      equals = .false. 
!!$      print*,'different size',this%reference%size(), b%reference%size()
!!$      return
!!$    end if
!!$
!!$    iter = this%reference%begin()
!!$    do while (iter /= this%reference%last())
!!$
!!$      if (.not. b%reference%has(iter%key())) then
!!$        equals = .false.
!!$        print*,'different key'
!!$        return
!!$      end if
!!$      
!!$      p1 => iter%value()
!!$      p2 => b%reference%getReference(iter%key())
!!$
!!$      if (.not. (p1%equals(p2))) then
!!$        equals = .false.
!!$        print*,'different value for key <',trim(iter%key()),'>'
!!$        call p1%print()
!!$        call p2%print()
!!$        return
!!$      end if
!!$
!!$      call iter%next()
!!$
!!$    end do
!!$
!!$  end function equals

  subroutine copyIt(a, b)
    type (TracerPointer), intent(inout) :: a
    type (TracerPointer), intent(in) :: b

    a%reference%AttributeHashMap = b%reference%AttributeHashMap
  end subroutine copyIt


!!$  subroutine writeUnformatted(this, unit)
!!$    class (TracerPointer), intent(in) :: this
!!$    integer, intent(in) :: unit
!!$
!!$    type (AttributeHashMapIterator) :: iter
!!$    class (AbstractAttribute), pointer :: p
!!$
!!$    write(unit) this%size()
!!$    iter = this%begin()
!!$    do while (iter /= this%last())
!!$
!!$      write(unit) iter%key()
!!$      p => iter%value()
!!$      write(unit) getAttributeType(p)
!!$      call p%writeUnformatted(unit)
!!$
!!$      call iter%next()
!!$    end do
!!$
!!$  contains
!!$
!!$    integer function getAttributeType(p) result(attributeType)
!!$      class (AbstractAttribute), pointer :: p
!!$      integer :: i
!!$
!!$      do i = 1, NUM_TYPES
!!$        if (same_type_as(p, prototypes(i)%p)) then
!!$          attributeType = i
!!$          return
!!$        end if
!!$      end do
!!$      
!!$      call stop_model('No prototype for attribute in AttributeDictionary writeUnformatted.',255)
!!$
!!$    end function getAttributeType
!!$
!!$  end subroutine writeUnformatted
!!$
!!$  subroutine readUnformatted(this, unit)
!!$    class (TracerPointer),intent(inout) :: this
!!$    integer, intent(in) :: unit
!!$
!!$    integer :: n
!!$    integer :: i
!!$    class (AbstractAttribute), pointer :: p, q
!!$    character(len=MAX_LEN_KEY) :: key
!!$    integer :: attributeType
!!$
!!$    read(unit) n
!!$    do i = 1, n
!!$      read(unit) key
!!$      read(unit) attributeType
!!$      p => prototypes(attributeType)%p
!!$      q => p%readUnformatted(unit)
!!$      call this%insert(trim(key), q)
!!$    end do
!!$
!!$  end subroutine readUnformatted

  subroutine clean_(this)
    type (TracerPointer), intent(inout) :: this
    call clean(this%reference%AttributeHashMap) ! parent
  end subroutine clean_

end module TracerPointer_mod

module TracerBundleSubset_mod
!@sum This module encapsulates the TracerBundleSubset data type. This type
!@+ contains a reference (pointer) to a TracerBundle data type and contains
!@+ methods that allow the subsetting a a parent bundle via filter functions.
!@auth NCCS ASTG

  use AttributeDictionary_mod, only: AttributeDictionary
  use Tracer_mod
  use TracerPointer_mod
  use TracerBundle_mod
  use TracerHashMap_mod
  implicit none
  private

  public :: TracerBundleSubset

  type, extends(TracerBundle) :: TracerBundleSubset
    private
    type (AttributeDictionary) :: privateValues
    class(TracerBundle), pointer :: reference
  contains
    procedure :: insertEntry 
    procedure :: getReference
    procedure :: getSubsetReference
    procedure :: print    
    !procedure :: merge        
    procedure :: getAttribute
    procedure :: findAttribute
    procedure :: setAttribute
    procedure :: hasAttribute
    ! iterator operations
    procedure :: begin        ! override base class method
    procedure :: last         ! override base class method
  end type TracerBundleSubset

  interface TracerBundleSubset
    module procedure newTracerBundleSubset
    module procedure newTracerBundleSubsetCopy
  end interface

  interface
    function filter (t)
      import
      class (Tracer), intent(in) :: t
      logical :: filter
    end function filter
  end interface

contains

  function newTracerBundleSubsetCopy(bundle) 
    class (TracerBundle), target :: bundle
    type (TracerBundleSubset) :: newTracerBundleSubsetCopy

    newTracerBundleSubsetCopy%TracerBundle = newTracerBundle() 
    newTracerBundleSubsetCopy%reference => bundle

  end function newTracerBundleSubsetCopy

  function newTracerBundleSubset(bundle, aFilter) 
    class (TracerBundle), target :: bundle
    procedure(filter) :: aFilter
    type (TracerBundleSubset) :: newTracerBundleSubset
!
    type (TracerIterator) :: iter
    class (Tracer), pointer :: tp
    type (TracerPointer) :: trp

    integer :: i

    i = 0

    newTracerBundleSubset%TracerBundle = newTracerBundle() 
    newTracerBundleSubset%reference => bundle
    iter = bundle%begin()
    do while (iter /= bundle%last())
       i = i + 1
      tp => iter%value()
      if (aFilter(tp)) then ! invoke parent insert method
        trp = TracerPointer(tp)
        call newTracerBundleSubset%TracerHashMap%insert(iter%key(), trp)
      end if
      call iter%next()
    end do

  end function newTracerBundleSubset

  subroutine print(this)
    class (TracerBundleSubset), intent(in) :: this
    type (TracerIterator) :: iter
    class (Tracer), pointer :: tp

    print*,'--------------------------'
    print*,' TracerBundleSubset: '
    print*,'--------------------------'

    iter = this%begin()
    do while (iter /= this%last())
      print*,'   key: <',trim(iter%key()),'>'
#ifdef HAS_PRINT
      tp => iter%value()
      call tp%print()
#endif
      call iter%next()
    end do
    print*,'--------------------------'
    print*,'--------------------------'
    print*,' '

  end subroutine print

  subroutine insertEntry(this, key, value)
    class (TracerBundleSubset), target, intent(inout) :: this
    character(len=*), intent(in) :: key ! name
    class (Tracer) :: value ! tracer
    class (Tracer), pointer :: p
    type (TracerPointer) :: tref
!    type (TracerReference) :: ref

    ! During a merge, reference bundle likely has
    ! the tracer already.  Attempting to insert creates
    ! an error, as the item getting inserted gets deleted
    ! to make room for itself ...   Working on a more robust
    ! solution to this issue, but for now, we just assume
    ! that identical key means identical tracer and skip
    ! the test.
    if (.not. this%reference%has(key)) then
       ! Only the bare tracer should go into reference bundle.
       select type (value)
       class is (TracerPointer)
          call this%reference%insert(key, value%get()) 
       class is (Tracer)
          call this%reference%insert(key, value) 
       end select
    end if

 !   ref = this%reference%getReference(key)
 !   p => ref%ptr
    p => this%reference%getReference(key)
    select type (p)
    class is (TracerPointer)
       call this%TracerHashMap%insert(key, p)
    class is (Tracer)
       tref = TracerPointer(p)
       call this%TracerHashMap%insert(key, tref)
    end select

  end subroutine insertEntry

  function getReference(this, key) result(ptr)
    use StringUtilities_mod, only: toLowerCase
    class (TracerBundleSubset), intent(in) :: this
    character(len=*), intent(in) :: key
    class (Tracer), pointer :: ptr
    class (Tracer), pointer :: ptmp
    
    ptmp => this%TracerBundle%getReference(key)
    select type (ptmp)
      class is (TracerPointer)   
      ptr => ptmp%get()
      do
        select type (q => ptr)

        class is (TracerPointer)
          ptr => q%get()
          class default
          exit
        end select
      end do
      class default
        call stop_model("Wrong type for contained object",1)
    end select
   
  end function getReference

  function getSubsetReference(this, key) result(ptr)
    class (TracerBundleSubset), intent(in) :: this
    character(len=*), intent(in) :: key
    class (Tracer), pointer :: ptr, p
!    type (TracerReference) :: ref

!    ref = this%TracerHashMap%getReference(key)
!    p => ref%ptr
    p => this%TracerHashMap%getReference(key)
    select type(p)
    type is(TracerPointer) 
      ptr => p%get()
    class default
    end select

  end function getSubsetReference

  !subroutine merge(this, b)
  !  class (TracerBundleSubset), intent(inout) :: this
  !  class (TracerBundleSubset), intent(in) :: b

  !  type (TracerIterator) :: iter
  !  class (Tracer), pointer :: tp

  !  print *, 'file->',__FILE__,'  merge!'
  !  iter = b%begin()
  !  do while (iter /= b%last())
!print *,iter%key()
  !     if (.not. this%has(iter%key())) then
  !        call this%TracerHashMap%insert(iter%key(), iter%value())
  !     end if
  !    call iter%next()
  !  end do

  !end subroutine merge

  type (TracerIterator) function begin(this) result(iterator)
    class (TracerBundleSubset), target, intent(in) :: this

    integer :: i

    iterator%reference => this

    do i = 1, this%tableSize
      if (this%table(i)%size() > 0) then
        iterator%hashValue = i
        iterator%subIterator = this%table(i)%begin()
        return
      end if
    end do
    
    iterator%hashValue = -1 ! no entries

  end function begin

  type (TracerIterator) function last(this) result(iterator)
    class (TracerBundleSubset), target, intent(in) :: this

    iterator%reference => this
    iterator%hashValue = -1 ! no entries

  end function last

  function findAttribute(this, species, attributeName) result(attribute)
    use AbstractAttribute_mod
    class (TracerBundleSubset), intent(in) :: this
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: attributeName

    class (AbstractAttribute), pointer :: attribute

    class (Tracer), pointer :: t
   
    t => this%getSubsetReference(species)
    attribute => t%getReference(attributeName) 
    
  end function findAttribute

  function getAttribute(this, species, attribute) result (attributeValue)
    use AbstractAttribute_mod
    class (TracerBundleSubset), intent(in) :: this
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: attribute
    class (AbstractAttribute), pointer :: attributeValue

    class (Tracer), pointer :: t

    t => this%getSubsetReference(trim(species))
    attributeValue => t%getReference(attribute)

  end function getAttribute

  subroutine setAttribute(this, species, attributeName, attributeValue)
    use AbstractAttribute_mod
    class (TracerBundleSubset), intent(inout) :: this
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: attributeName
    class (AbstractAttribute), intent(in) :: attributeValue

    class (Tracer), pointer :: t

    t => this%getSubsetReference(species)
    call t%insert(attributeName, attributeValue)

  end subroutine setAttribute

  function hasAttribute(this, attribute) result(has)
!@sum  This function returns a logical array of length equal to
!@+    the total number of tracers in the bundle.  Values are .true.
!@+    for those tracers that have the specified attribute.
!@+    Use should be limited, as order of tracers in bundle is meant to be
!@+    hidden from clients.
    class (TracerBundleSubset), intent(in) :: this
    character(len=*), intent(in) :: attribute
    logical, pointer :: has(:)

    type (TracerIterator) :: iter
    class (Tracer), pointer :: t, p
    integer :: i

    allocate(has(this%size()))
    has = .false.
    iter = this%begin()
    i = 0
    do while (iter /= this%last())
      i = i + 1
      t => iter%value()
      select type (t)
      type is (TracerPointer)
        p => t%get()
        has(i) = p%has(attribute)
      end select
      call iter%next()
    end do

  end function hasAttribute

end module TracerBundleSubset_mod
