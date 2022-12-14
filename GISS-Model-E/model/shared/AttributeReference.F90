module AttributeReference_mod
  use AbstractAttribute_mod
  implicit none
  private

  public :: AttributeReference
  public :: VectorAttribute
  public :: newVectorAttribute
  public :: assignment(=)
  public :: toPointer
  
  type AttributeReference
!!$    private
    class (AbstractAttribute), pointer :: ptr => null()
  contains
    procedure :: set
    procedure :: get
  end type AttributeReference

  type, extends(AbstractAttribute) :: VectorAttribute
!!$    private
    type (AttributeReference), allocatable :: items(:)
  contains
    procedure :: equals
    procedure :: print
    procedure :: toString
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: clean
    procedure :: getReferenceScalar
    procedure :: getReferenceVector
    procedure :: size_
    generic :: size => size_
  end type VectorAttribute

  interface assignment(=)
    module procedure toType_
  end interface assignment(=)

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  subroutine set(this, reference)
    class (AttributeReference), intent(out) :: this
    class (AbstractAttribute), target :: reference

    this%ptr => reference

  end subroutine set

  function get(this) result(ptr)
    class (AttributeReference), intent(in) :: this
    class (AbstractAttribute), pointer :: ptr

    ptr => this%ptr

  end function get

  function newVectorAttribute(b) result(this)
    type (VectorAttribute) :: this
    type (AttributeReference), target :: b(:)

    integer :: i, n

    n = size(b)

    allocate(this%items(n))
    do i = 1, n
      this%items(i)%ptr => b(i)%ptr ! shallow copy - preserve references
    end do

  end function newVectorAttribute
  
  subroutine toType_(a, b)
    type (AttributeReference), allocatable, intent(out) :: a(:)
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    type is (VectorAttribute)
      a = p%items
    class default
      call stop_model('Illegal conversion of VectorAttribute.',255)
    end select

  end subroutine toType_


  logical function equals(this, b)
    class (VectorAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b
    equals = .true.
  end function Equals

  function toString(this) result(string)
    use AbstractAttribute_mod, only: MAX_LEN_LINE
    class (VectorAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string
  end function ToString

  subroutine print(this)
    class (VectorAttribute), intent(in) :: this
  end subroutine Print

  subroutine writeUnformatted(this, unit)
    class (VectorAttribute), intent(in) :: this
    integer, intent(in) :: unit
  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    class (VectorAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new
  end function readUnformatted

  subroutine clean(this)
    class (VectorAttribute), intent(inout) :: this
  end subroutine Clean

  subroutine getReferenceScalar(this, reference)
    class (VectorAttribute), target, intent(in) :: this
    class (*), pointer, intent(out) :: reference
    reference => null()
  end subroutine getReferenceScalar

  subroutine getReferenceVector(this, reference)
    class (VectorAttribute), target, intent(in) :: this
    class (*), pointer, intent(out) :: reference(:)
    reference => null()
  end subroutine getReferenceVector

  function toPointerType(this, vector) result(ptr)
     type (AttributeReference), pointer :: ptr(:)
     class (AbstractAttribute), target, intent(in) :: this
     type (AttributeReference), intent(in) :: vector(:)

     select type (this)
     class is (VectorAttribute)
        ptr => this%items
     end select

  end function toPointerType

  integer function size_(this)
     class (VectorAttribute), intent(in) :: this
     size_ = size(this%items)
  end function size_
  

end module AttributeReference_mod

