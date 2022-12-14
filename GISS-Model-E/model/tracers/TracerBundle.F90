module TracerBundle_mod
  use AttributeDictionary_mod
  use Tracer_mod
  use TracerHashMap_mod
  implicit none
  private

  public :: TracerBundle ! derived type
  public :: newTracerBundle      ! constructor
  public :: readFromText ! constructor
  public :: readUnformattedBundle ! constructor
  public :: operator(==)
  public :: clean
!  public :: assignment(=) ! re-export from hash package

  public :: NOT_FOUND

  type, extends(TracerHashMap) :: TracerBundle
!!$    private
    type (AttributeDictionary) :: defaultValues
    character(len=MAX_LEN_KEY), allocatable :: mandatoryAttributes(:)
    logical :: locked = .false.
    type (AttributeDictionary) :: attributeVectorCache
  contains
    procedure :: delete
    procedure :: insertEntry ! override base class method
    procedure :: insertGetName ! extend generic
    generic   :: insert => insertGetName
    procedure :: findAttribute
    procedure :: getAttribute
    procedure :: setAttribute
    procedure :: hasAttribute
    procedure :: addDefault_integer
    procedure :: addDefault_logical
    procedure :: addDefault_real64
    procedure :: addDefault_string
    generic   :: addDefaultValue => addDefault_integer, &
                                    addDefault_logical, &
                                    addDefault_real64, &
                                    addDefault_string
    procedure :: countHaveAttribute
    procedure :: getAttributeVector
    procedure :: addMandatoryAttribute
    procedure :: writeFormatted
    procedure :: writeUnformatted => writeUnformatted_bundle
  end type TracerBundle

  interface clean
    module procedure cleanBundle
  end interface

  integer, parameter :: NOT_FOUND = -1

  interface operator(==)
    module procedure equals
  end interface

  integer, parameter :: LEN_HEADER = 80
  integer, parameter :: VERSION = 1
  character(len=*), parameter :: DESCRIPTION = 'TracerBundle'

contains

  function newTracerBundle()
!@sum Construct empty tracer    
    use AttributeDictionary_mod
    type (TracerBundle) :: newTracerBundle

    allocate(newTracerBundle%mandatoryAttributes(0))

    newTracerBundle%TracerHashMap = newTracerHashMap(100)
!!$    newTracerBundle%TracerHashMap = newTracerHashMap(1)
    newTracerBundle%defaultValues = newAttributeDictionary()
    
    newTracerBundle%attributeVectorCache = newAttributeDictionary()
    newTracerBundle%locked = .false.
    
  end function newTracerBundle

  subroutine setAttribute(this, species, attributeName, attributeValue)
    use AbstractAttribute_mod
    class (TracerBundle), intent(inout) :: this
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: attributeName
    class (AbstractAttribute), intent(in) :: attributeValue
!    type(TracerReference) :: ref
    class (Tracer), pointer :: t

    t => this%getReference(species)
!    ref = this%getReference(species)
!    t => ref%ptr
    call  t%insert(attributeName, attributeValue)

  end subroutine setAttribute
 
  function readFromText(unit, defaultValues) result(bundle)
!@sum Populate a TracerBundle from a unit attached to a formatted
!@+ file.  Optionally apply default values.
    use AttributeDictionary_mod
    use Parser_mod
    integer, intent(in) :: unit
    type (AttributeDictionary), optional, intent(in) :: defaultValues
    type (TracerBundle) :: bundle
    type (Tracer) :: aTracer
    type (Parser_type) :: parser

    integer :: status

    bundle = newTracerBundle()

    if (present(defaultValues)) then
      ! TODO might need a deep copy here?
      bundle%defaultValues = defaultValues
    else
      bundle%defaultValues = newAttributeDictionary()
    end if
    
    do
      aTracer = readOneTracer(unit, status)
      if (status /= 0) exit
      call bundle%insert(aTracer%getName(), aTracer)
    end do

  end function readFromText

  subroutine mergeDefaults(this, attributeName, attribute)
    use AbstractAttribute_mod
    type (TracerBundle), intent(inout) :: this
    character(len=*), intent(in) :: attributeName
    class (AbstractAttribute), intent(in) :: attribute

    type (TracerIterator) :: iter
    class (Tracer), pointer :: p

    iter = this%begin()
    do while (iter /= this%last())
      p => iter%value()
      if (.not. p%has(attributeName)) then
        call p%insert(attributeName, attribute)
      end if
      call iter%next()
    end do

  end subroutine mergeDefaults

  subroutine addDefault_integer(this, attribute, value)
!@sum Add a default attribute that applies to all tracers in bundle,
!@+ unless overridden with specific value for a given tracer.
    use Attributes_mod
    use AbstractAttribute_mod
    class (TracerBundle), target, intent(inout) :: this
    character(len=*), intent(in) :: attribute
    integer, intent(in) :: value

    call mergeDefaults(this, attribute, newAttribute(value))
    call this%defaultValues%insert(attribute, value)

  end subroutine addDefault_integer

  subroutine addDefault_logical(this, attribute, value)
!@sum Add a default attribute that applies to all tracers in bundle,
!@+ unless overridden with specific value for a given tracer.
    use Attributes_mod
    class (TracerBundle), target, intent(inout) :: this
    character(len=*), intent(in) :: attribute
    logical, intent(in) :: value

    call mergeDefaults(this, attribute, newAttribute(value))
    call this%defaultValues%insert(attribute, value)

  end subroutine addDefault_logical

  subroutine addDefault_real64(this, attribute, value)
!@sum Add a default attribute that applies to all tracers in bundle,
!@+ unless overridden with specific value for a given tracer.
    use Attributes_mod
    class (TracerBundle), target, intent(inout) :: this
    character(len=*), intent(in) :: attribute
    real*8, intent(in) :: value

    call mergeDefaults(this, attribute, newAttribute(value))
    call this%defaultValues%insert(attribute, value)

  end subroutine addDefault_real64

  subroutine addDefault_string(this, attribute, value)
!@sum Add a default attribute that applies to all tracers in bundle,
!@+ unless overridden with specific value for a given tracer.
    use Attributes_mod
    class (TracerBundle), target, intent(inout) :: this
    character(len=*), intent(in) :: attribute
    character(len=*), intent(in) :: value

    call mergeDefaults(this, attribute, newAttribute(value))
    call this%defaultValues%insert(attribute, value)

  end subroutine addDefault_string

  subroutine writeUnformatted_bundle(this, unit)
!@sum Write a bundle to a unit attached to an unformatted sequential file.
    class (TracerBundle), intent(in) :: this
    integer, intent(in) :: unit

    integer :: n
    type (TracerIterator) :: iter

    character(len=LEN_HEADER) :: header

    write(header, '(a," - version:",1x,i5.0)') DESCRIPTION, VERSION
    write(unit) header

    n = this%size()
    write(unit) n
    iter = this%begin()
    do while (iter /= this%last())
      call writeUnformatted(iter%value(), unit)
      call iter%next()
    end do

  end subroutine writeUnformatted_bundle

  function readUnformattedBundle(unit) result(this)
!@sum Read a bundle to a unit attached to an unformatted sequential file.
    type (TracerBundle) :: this
    integer, intent(in) :: unit

    integer :: i, n
    integer :: oldVersion
    character(len=len(DESCRIPTION)) :: tag
    character(len=LEN_HEADER) :: header
    type (Tracer) :: t

    read(unit) header
    read(header, '(a,11x,i10.0)') tag, oldVersion

    if (tag /= DESCRIPTION) then
      call stop_model(DESCRIPTION // '::readUnformatted() - incorrect header.', 14)
    end if

    if (oldVersion /= VERSION) then
      call stop_model(DESCRIPTION // '::readUnformatted() - unsupported format.', 14)
    end if
      
    read(unit) n
    this = newTracerBundle()
    do i = 1, n
      t = newTracer()
      call readUnformattedTracer(t, unit)
      call this%insert(t%getName(), t)
    end do

  end function readUnformattedBundle

  subroutine writeFormatted(this, unit)
    use Parser_mod, only: Parser_type, setBeginData, setEndData, setTokenSeparators
    use Parser_mod, only: parserWriteFormatted => writeFormatted
    class (TracerBundle), intent(in) :: this
    integer, intent(in) :: unit

    type (Parser_type) :: parser
    class (Tracer), pointer :: t
    type (TracerIterator) :: iter

    call setBeginData(parser, '{')
    call setEndData(parser, '}')
    call setTokenSeparators(parser, '=,')

    iter = this%begin()
    do while (iter /= this%last())
      t => iter%value()
      call parserWriteFormatted(parser, unit, t)
      call iter%next()
    end do

  end subroutine writeFormatted

  function hasAttribute(this, attribute) result(has)
!@sum  This function returns a logical array of length equal to
!@+    the total number of tracers in the bundle.  Values are .true.
!@+    for those tracers that have the specified attribute.
!@+    Use should be limited, as order of tracers in bundle is meant to be
!@+    hidden from clients.
    class (TracerBundle), intent(in) :: this
    character(len=*), intent(in) :: attribute
    logical, pointer :: has(:)

    type (TracerIterator) :: iter
    class (Tracer), pointer :: t
    integer :: i

    allocate(has(this%size()))
    iter = this%begin()
    i = 0
    do while (iter /= this%last())
      i = i + 1
      t => iter%value()
      has(i) = t%has(attribute)
      call iter%next()
    end do

  end function hasAttribute

  subroutine addMandatoryAttribute(this, attributeName)
    class (TracerBundle), intent(inout) :: this
    character(len=*), intent(in) :: attributeName

    character(len=MAX_LEN_KEY), allocatable :: oldAttributes(:)
    integer :: n
    type (TracerIterator) :: iter

    n = size(this%mandatoryAttributes)

    allocate(oldAttributes(n))
    oldAttributes = this%mandatoryAttributes
    deallocate(this%mandatoryAttributes)
    allocate(this%mandatoryAttributes(n+1))
    this%mandatoryAttributes(1:n) = oldAttributes
    deallocate(oldAttributes)

    this%mandatoryAttributes(n+1) = attributeName

    iter = this%begin()
    do while (iter /= this%last())
      if (.not. assertHasAttribute(iter%value(), attributeName)) return
      call iter%next()
    end do

  end subroutine addMandatoryAttribute

  subroutine assertHasAttributes(this, attributes)
    type (Tracer), intent(in) :: this
    character(len=MAX_LEN_KEY), intent(in) :: attributes(:)

    integer :: i

    do i = 1, size(attributes)
      if (.not. assertHasAttribute(this, attributes(i))) return
    end do

  end subroutine assertHasAttributes

  logical function assertHasAttribute(this, attribute)
    type (Tracer), intent(in) :: this
    character(len=*), intent(in) :: attribute

    character(len=MAX_LEN_KEY) :: name

    assertHasAttribute = .true.
    if (.not. this%has(attribute)) then
      name = this%getName()
      call stop_model("TracerBundle_mod - species '" // trim(name) // &
        & "' is missing mandatory attribute '" // trim(attribute) // "'.", 14)
      assertHasAttribute = .false.
    end if

  end function assertHasAttribute

  function getAttribute(this, species, attribute) result (attributeValue)
    use AbstractAttribute_mod
    class (TracerBundle), intent(in) :: this
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: attribute
    class (AbstractAttribute), pointer :: attributeValue
    class (AbstractAttribute), pointer :: attr
!    type (TRACERreference) :: ref

    class (Tracer), pointer :: t

    t => this%getReference(trim(species))
!    ref = this%getReference(trim(species))
!    t => ref%ptr
    attributeValue => t%getReference(attribute)

  end function getAttribute

  integer function countHaveAttribute(this, withAttribute) result(number)
    class (TracerBundle), intent(in) :: this
    character(len=*), intent(in) :: withAttribute

    number = count(hasAttribute(this, withAttribute))

  end function countHaveAttribute

  function makeSubset(this, withAttribute) result(subset)
    class (TracerBundle), target, intent(in) :: this
    character(len=*), intent(in) :: withAttribute
    type (TracerBundle) :: subset

    type (TracerIterator) :: iter
    class (Tracer), pointer :: t

    subset = newTracerBundle()
    subset%defaultValues = this%defaultValues

    iter = this%begin()
    do while (iter /= this%last())
      t => iter%value()
      if (t%has(withAttribute)) then
        call subset%insertReference(trim(iter%key()), iter%value())
      end if
      call iter%next()
    end do

    ! no further modifications are permitted
    subset%locked = .true. 

  end function makeSubset

  function getAttributeVector(this, attributeName) result(vector)
    use AbstractAttribute_mod
    use AttributeReference_mod
    class (TracerBundle), intent(inout) :: this
    character(len=*), intent(in) :: attributeName
    type (AttributeReference), pointer :: vector(:)

    class (Tracer), pointer :: t
    type (VectorAttribute) :: vecAttr
    type (TracerIterator) :: iter
    class (AbstractAttribute), pointer :: attribute

    integer :: i

    ! must be at least one tracer to determine type of result
    if (this%size() == 0) then
      vector => null()
      return
    end if

    if (this%attributeVectorCache%has(attributeName)) then
      ! Should be doable in 1 step, but compiler struggles ...
      attribute => this%attributeVectorCache%getReference(attributeName)
      vector => toPointer(attribute, vector)
      return
    end if

    this%locked = .true.  ! attributeVectorCache will be invalid if new tracers are added to bundle

    allocate(vector(this%size()))

    iter = this%begin()
    i = 1

    do while (iter /= this%last())
      t => iter%value()
      if (.not. t%has(attributeName)) then
        call stop_model('All tracers must have specified attribute to use getAttributeVector() method.',14)
        return
      end if
      attribute => t%getReference(attributeName)
      call vector(i)%set(attribute)
      i = i + 1
      call iter%next()
    end do

    vecAttr = newVectorAttribute(vector)
    call this%attributeVectorCache%insert(attributeName, vecAttr) ! save for efficient reference next time

  end function getAttributeVector

  logical function equals(bundleA, bundleB) result(isEqual)
    use Parser_mod, only: MAX_LEN_TOKEN
    type (TracerBundle), target, intent(in) :: bundleA
    type (TracerBundle), target, intent(in) :: bundleB

    character(len=MAX_LEN_TOKEN) :: name

    type (TracerIterator) :: iter
    class (Tracer), pointer :: t
!    type (TRACERreference) :: ref

    isEqual = .true.

    iter = bundleA%begin()
    do while (iter /= bundleA%last())
      name = trim(iter%key())
      t => iter%value()
!      ref = bundleB%getReference(name)
!      if (.not. (t%equals(ref%ptr))) then
      if (.not. (t%equals(bundleB%getReference(name)))) then
        isEqual = .false.
        exit
      end if
      call iter%next()
    end do

  end function equals

  subroutine cleanBundle(this)
    type (TracerBundle), intent(inout) :: this

    call clean(this%defaultValues)
    !    if (size(this%mandatoryAttributes)>0) then
    !       print *, 'SIZE = ',size(this%mandatoryAttributes)
    deallocate(this%mandatoryAttributes)
    !    end if
  end subroutine cleanBundle

  subroutine insertEntry(this, key, value)
    
    class (TracerBundle), target, intent(inout) :: this
    character(len=*), intent(in) :: key ! name
    class (Tracer) :: value ! tracer
    class (Tracer), pointer :: p

!    type(TRACERreference) :: ref

    if (this%locked) then
      call stop_model("TracerBundle_mod - cannot insert new tracer into subset. " // &
           & "Subsets are locked from modification.",14)
    end if

    call assertHasAttributes(value, this%mandatoryAttributes)

    call this%TracerHashMap%insertEntry(key, value) ! invoke parent method
!    ref = this%getReference(key)

!    call ref%ptr%merge(this%defaultValues)
    p => this%getReference(key)

    call p%merge(this%defaultValues)

  end subroutine insertEntry

  subroutine insertGetName(this, value)
    class (TracerBundle), intent(inout) :: this
    class (Tracer) :: value ! tracer

    call this%insertEntry(value%getName(), value)

  end subroutine insertGetName

  function findAttribute(this, species, attributeName) result(attribute)
    use AbstractAttribute_mod
    class (TracerBundle), intent(in) :: this
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: attributeName
    class (AbstractAttribute), pointer :: attribute
    class (Tracer), pointer :: t
!    type(TRACERreference) :: ref
!    type(AbstractAttributeReference) :: refAttr

    t => this%getReference(species)
! This is now the wrong interface:
!    attribute => t%getReference(attributeName)
    attribute => null()

!    ref = this%getReference(species)
!    refAttr = ref%ptr%getReference(attributeName)
!    attribute => refAttr%ptr

  end function findAttribute

  subroutine delete(this)
    class (TracerBundle), intent(inout) :: this
    integer :: i

    call clean(this%defaultValues)
    deallocate(this%mandatoryAttributes)

  end subroutine delete

end module TracerBundle_mod
