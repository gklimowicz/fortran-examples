module BaseTime_mod
  use Rational_mod
  implicit none
  private

  public :: BaseTime
  public :: newBaseTime
  public :: assignment(=)
  public :: real ! from Rational

  type, extends(Rational) :: BaseTime
  contains
    procedure, private :: set_rational
    procedure, private :: set_i4
    procedure, private :: set_i8
    procedure, private :: set_r8
    procedure, private :: add_dt

!!$    procedure, private :: multiply_by_int4
!!$    procedure, private :: divide_by_int4

!!$    generic :: operator(+) => add_dt
!!$    generic :: operator(*) => multiply_by_int4
!!$    generic :: operator(/) => divide_by_int4
!!$    generic :: assignment(=) => set_rational, set_i4, set_i8, set_r8

  end type BaseTime

  ! broken ifort 14.0.2 sigh (should be BaseTime, not newBaseTime)
  interface newBaseTime
    module procedure newBaseTime_int4
    module procedure newBaseTime_int8
    module procedure newBaseTime_r8
    module procedure newBaseTime_rational
    module procedure newBaseTime_string
  end interface NewBaseTime

  interface assignment(=)
     module procedure copyFromRational
  end interface

  integer, parameter :: DEFAULT_INT = kind(1)
  integer, parameter :: LONG_INT = selected_int_kind(10)
  integer, parameter :: LONG = max(DEFAULT_INT, LONG_INT) ! in case only 32 bit available

contains

  function newBaseTime_int4(i) result(t)
    type (BaseTime) :: t
    integer, intent(in) :: i
    t%Rational = Rational(i)
  end function newBaseTime_int4

  function newBaseTime_int8(i) result(t)
    type (BaseTime) :: t
    integer(kind=LONG), intent(in) :: i
    t%Rational = Rational(i)
  end function newBaseTime_int8

  function newBaseTime_r8(x) result(t)
    type (BaseTime) :: t
    real(kind=8), intent(in) :: x
    t%Rational = Rational(x, 1.d-6)
  end function newBaseTime_r8

  function newBaseTime_rational(r) result(t)
    type (BaseTime) :: t
    type (Rational), intent(in) :: r
    t%Rational = r
  end function newBaseTime_rational


  function newBaseTime_string(string) result(t)
    type (BaseTime) :: t
    character(len=*), intent(in) :: string
    t%Rational = Rational(string)
 end function newBaseTime_string


  subroutine set_rational(this, r)
    class (BaseTime), intent(inout) :: this
    type (Rational), intent(in) :: r
    this%Rational = r
  end subroutine set_rational

  subroutine set_i4(this, wholeSeconds)
    class (BaseTime), intent(inout) :: this
    integer, intent(in) :: wholeSeconds
    this%Rational = Rational(wholeSeconds)
  end subroutine set_i4

  subroutine set_i8(this, wholeSeconds)
    class (BaseTime), intent(inout) :: this
    integer*8, intent(in) :: wholeSeconds
    this%Rational = Rational(wholeSeconds)
  end subroutine set_i8

  subroutine set_r8(this, time)
    class (BaseTime), intent(inout) :: this
    real*8, intent(in) :: time
    this%Rational = Rational(time, tolerance=1.d-6) ! microseconds
  end subroutine set_r8

  function add_dt(this, dt) result(sum)
    class (BaseTime), intent(in) :: this
    type (BaseTime), intent(in) :: dt
    type (BaseTime) :: sum

    type (Rational) :: a, b

    a = this%Rational
    b = dt%Rational
    sum%Rational = a + b
    sum%Rational = this + dt

  end function add_dt
  
  function multiply_by_int4(this, n) result(product)
    class (BaseTime), intent(in) :: this
    integer, intent(in) :: n
    type (BaseTime) :: product

    type (Rational) :: a

    product = newBaseTime(this * n)
    
  end function multiply_by_int4
  
  function divide_by_int4(this, n) result(quotient)
    class (BaseTime), intent(in) :: this
    integer, intent(in) :: n
    type (BaseTime) :: quotient

    type (Rational) :: a

    a = this%Rational
    quotient%Rational = a / n
    
  end function divide_by_int4


  subroutine copyFromRational(t, r)
    type (BaseTime), intent(out) :: t
    type (Rational), intent(in) :: r
    t%Rational = r
  end subroutine copyFromRational

  
end module BaseTime_mod
