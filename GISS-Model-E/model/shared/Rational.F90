module Rational_mod
!@sum
!@+!---------------------------------------------------------------------------
!@+
!@+ This module provides a class which implements rational
!@+ numbers/arithmetic of the form q = w + n/d where {w,n,d} are
!@+ integers.  This representation, motivated by the treatment in ESMF,
!@+ provides a large range as well as sufficient accuracy to represent
!@+ most floating point values.
!@+
!@+ Times in the model are subclasses of Rational {BaseTime, Time,
!@+ TimeInterval}.  This ensures that certain ratios are exact integers
!@+ and prevents roundoff from inducing incorrect assertions about 
!@+ beginnings of cyclic periods (e.g. days).
!@+
!@+ NOTE: arithmetic among nontrivial Rational objects can easily lead to
!@+ integer overflow which is not detectable by most compilers.
!@+
!@auth T. Clune
!@+---------------------------------------------------------------------------

  use KindParameters_mod, only: SP, DP, DEFAULT_INT, LONG_INT
  implicit none
  private

!---------------------------------------------------------------------------
! Public entities

  ! Derived type and constructor
  public :: Rational

  ! Procedures for rounding to integer
  public :: nint
  public :: floor
  public :: ceiling
  public :: real

  ! modulo(q) returns   q - floor(q)
  public :: modulo
!---------------------------------------------------------------------------

  integer, parameter :: LONG = max(DEFAULT_INT, LONG_INT) ! in case only 32 bit available

  ! Numerical value is whole + numerator/denominator
  type Rational
    private
    integer(kind=LONG) :: whole       = 0
    integer(kind=LONG) :: numerator   = 0 
    integer(kind=LONG) :: denominator = 1 ! always positive
  contains
    procedure :: getWhole

    ! Arithmetic operations
    generic :: operator(+) => add_fraction
    generic :: operator(-) => subtract_fraction
    generic :: operator(*) => multiply_fraction, &
         & multiply_intLeft, multiply_intRight, &
         & multiply_longLeft, multiply_longRight
    generic :: operator(/) => divide_fraction, &
         &  divide_realDP, &
         &  divide_byInt4, divide_byInt8, &
         & divide_intoInt4, divide_intoInt8

    ! Comparison operators
    generic :: operator(==) => equals_fraction, equals_intLeft, equals_intRight
    generic :: operator(<) => lessThan_fraction
    generic :: operator(>) => greaterThan_fraction
    generic :: operator(>=) => greaterThanOrEqualTo_fraction

    procedure :: add_fraction       ! r = p + q
    procedure :: subtract_fraction  ! r = p - q

    procedure :: multiply_fraction  ! r = p * q
    procedure, pass(a) :: multiply_intLeft  ! r = i * p 
    procedure, pass(a) :: multiply_intRight ! r = p * i
    procedure, pass(a) :: multiply_longLeft ! r = i_8 * p
    procedure, pass(a) :: multiply_longRight ! r = p * i_8

    procedure :: divide_fraction   ! r = p / q
    procedure, pass(a) :: divide_realdp  ! r = p / x
    procedure, pass(a) :: divide_byInt4  ! r = p / i
    procedure, pass(a) :: divide_byInt8  ! r = p / i8
    procedure, pass(a) :: divide_intoInt4 ! r = i4 / p
    procedure, pass(a) :: divide_intoInt8 ! r = i8 / p

    ! Comparison operators
    procedure :: equals_fraction   ! r == p
    procedure, pass(a) :: equals_intRight  ! r == i
    procedure, pass(b) :: equals_intLeft   ! i == r
    procedure :: lessThan_fraction !  r < p
    procedure :: greaterThan_fraction ! r > p
    procedure :: greaterThanOrEqualTo_fraction ! r >= p

    ! conversion
    procedure :: toReal
    procedure :: toString

    ! internal use only
    procedure, private :: reduce ! put in canonical form

  end type Rational

  interface Rational
    module procedure newRational_defaultWhole
    module procedure newRational_default_n_over_d
    module procedure newRational_default
#ifndef LONG_UNSUPPORTED
    module procedure newRational_longWhole
    module procedure newRational_long_n_over_d
    module procedure newRational_long
#endif

    module procedure newRational_real_sp
    module procedure newRational_real_dp
    module procedure newRational_string
  end interface Rational

  interface nint
    module procedure nintRational
  end interface nint

  interface floor
    module procedure floorRational
  end interface floor

  interface ceiling
    module procedure ceilingRational
  end interface ceiling

  interface real
     module procedure real_rational
  end interface real

  interface modulo
    module procedure moduloRational
 end interface modulo

contains

  integer(kind=LONG) function getWhole(this) result(whole)
    class (Rational), intent(in) :: this
    whole = this%whole
  end function getWhole


  function toReal(this) result(x)
    real(kind=DP) :: x
    class(Rational), intent(in) :: this
    x = real(this%numerator,kind=DP) / real(this%denominator,kind=DP)
    x = x + this%whole
  end function toReal


  function real_rational(this) result(x)
    use iso_fortran_env, only: REAL64
    real(kind=REAL64) :: x
    class (Rational), intent(in) :: this

    x = this%toReal()
    
  end function real_rational

  ! Return the nearest integer
  ! Rounds to even for r = n + m/2
  integer function nintRational(q) result(n)
    class (Rational), intent(in) :: q

    n = abs(q%whole)

    if (2*abs(q%numerator) > q%denominator) then
      n = n + 1
    else if (2*abs(q%numerator) == q%denominator) then
      if (mod(n,2) == 1) n = n + 1
    end if
    
    if (q%whole < 0 .or. q%numerator < 0) n = - n

  end function nintRational

  integer function floorRational(q) result(floor)
    class (Rational), intent(in) :: q

    if (q%numerator >= 0) then ! exact integer
       floor = q%getWhole()
    else
       floor = q%getWhole() - 1
    end if

  end function floorRational


  integer function ceilingRational(q) result(ceiling)
    class (Rational), intent(in) :: q

    if (q%numerator <= 0) then ! exact integer
       ceiling = q%getWhole()
    else
       ceiling = q%getWhole() + 1
    end if

  end function ceilingRational


  function moduloRational(q, r) result(modulo)
    type (Rational) :: modulo
    class (Rational), intent(in) :: q
    class (Rational), intent(in) :: r

    type (Rational) :: s

    s = q / r
    modulo = q - (r * floor(s))

  end function moduloRational


  function newRational_defaultWhole(whole) result(r)
    type (Rational) :: r
    integer, intent(in) :: whole

    r%whole = whole
    r%numerator = 0
    r%denominator = 1
  end function newRational_defaultWhole


  function newRational_default_n_over_d(numerator, denominator) result(r)
    type (Rational) :: r
    integer, intent(in) :: numerator
    integer, intent(in) :: denominator

    r%whole = 0
    r%numerator = numerator
    r%denominator = denominator

    call r%reduce()

  end function newRational_default_n_over_d


  function newRational_default(whole, numerator, denominator) result(r)
    type (Rational) :: r
    integer, intent(in) :: whole
    integer, intent(in) :: numerator
    integer, intent(in) :: denominator

    ! TODO: Ifort - 14.0 does not use the correct interface for other part of constructor
    r = Rational(whole) + Rational(numerator, denominator)

 end function newRational_default

#ifndef LONG_UNSUPPORTED

  function newRational_longWhole(whole) result(r)
    type (Rational) :: r
    integer(kind=LONG), intent(in) :: whole

    r%whole = whole
    r%numerator = 0
    r%denominator = 1
  end function newRational_longWhole

  function newRational_long_n_over_d(numerator, denominator) result(r)
    type (Rational) :: r
    integer(kind=LONG), intent(in) :: numerator
    integer(kind=LONG), intent(in) :: denominator

    r%whole = 0
    r%numerator = numerator
    r%denominator = denominator
    call r%reduce()
  end function newRational_long_n_over_d

  function newRational_long(whole, numerator, denominator) result(r)
    type (Rational) :: r
    integer(kind=LONG), intent(in) :: whole
    integer(kind=LONG), intent(in) :: numerator
    integer(kind=LONG), intent(in) :: denominator

    r = Rational(whole) + Rational(numerator, denominator)
 end function newRational_long
#endif

! Use continued fractions to convert floating point to a fraction within
! a specified tolerance.
! http://en.wikipedia.org/wiki/Continued_fraction
! http://www-math.mit.edu/phase2/UJM/vol1/COLLIN~1.PDF

  function newRational_real_sp(x, tolerance) result(r)
     type (Rational) :: r
     real(kind=SP), intent(in) :: x
     real(kind=SP), intent(in) :: tolerance

     r = Rational(real(x,kind=DP), real(tolerance,kind=DP))
  end function newRational_real_sp

  function newRational_real_dp(x, tolerance) result(r)
     type (Rational) :: r
     real(kind=DP), intent(in) :: x
     real(kind=DP), intent(in) :: tolerance

     real(kind=DP) :: xx, absx
     integer(kind=LONG) :: a, w
     integer(kind=LONG) :: p_n, p_nm1, p_nm2
     integer(kind=LONG) :: q_n, q_nm1, q_nm2

     absx = abs(x)

     xx = absx
     w = floor(xx)
     xx = xx - w

     p_nm2 = 1
     q_nm2 = 0

     p_nm1 = 0
     q_nm1 = 1

     p_n = 0
     q_n = 1

     do while (xx /= 0.d0)
       xx = 1/xx
       a = floor(xx)
       xx = xx - a

        p_n = a * p_nm1 + p_nm2
        q_n = a * q_nm1 + q_nm2
        
        p_nm2 = p_nm1
        q_nm2 = q_nm1
        
        p_nm1 = p_n
        q_nm1 = q_n
        
        if (abs(q_n*(absx-w) - p_n) < q_n*tolerance) exit

     end do

     ! Restore sign of result
     if (x < 0) then
        w = -w
        p_n = -p_n
     end if

     r = newRational_long(w, p_n, q_n)

  end function newRational_real_dp

  ! Supports checkpoint/restart
  function newRational_string(string) result(r)
     type (Rational) :: r
     character(len=*), intent(in) :: string

     read(string,*) r%whole, r%numerator, r%denominator
  end function newRational_string

! Add two fractions and reduce to simplest form.
  function add_fraction(a, b) result(c)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b
    type (Rational) :: c

    integer(kind=LONG) :: gcf

    c%whole = a%whole + b%whole 

    gcf = greatestCommonFactor(a%denominator, b%denominator)
    
    c%numerator= a%numerator*(b%denominator/gcf) + b%numerator*(a%denominator/gcf)
    c%denominator = (a%denominator / gcf) * b%denominator

    call c%reduce()

  end function add_fraction

! Subtract two fractions and reduce to simplest form.
! a - b = a + (- b)
  function subtract_fraction(a, b) result(c)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b
    type (Rational) :: c

    c = a + b*(-1)

  end function subtract_fraction

! Multiply two fractions and reduce to simplest form.
  function multiply_fraction(a, b) result(c)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b
    type (Rational) :: c

    c%whole = a%whole * b%whole 
    c%numerator = b%whole*a%numerator*b%denominator + &
         & a%whole*b%numerator*a%denominator + &
         & a%numerator * b%numerator
    c%denominator = a%denominator*b%denominator
    call c%reduce()

  end function multiply_fraction

! Multiply by int on right
  function multiply_intRight(a, i) result(c)
    class (Rational), intent(in) :: a
    integer, intent(in) :: i
    type (Rational) :: c

    c%whole = a%whole * i
    c%numerator = a%numerator * i
    c%denominator = a%denominator
    call c%reduce()

  end function multiply_intRight

! Multiply by int on left
  function multiply_intLeft(i, a) result(c)
    integer, intent(in) :: i
    class (Rational), intent(in) :: a
    type (Rational) :: c

    c%whole = i * a%whole
    c%numerator = i * a%numerator
    c%denominator = a%denominator
    call c%reduce()

  end function multiply_intLeft
  
! Multiply by long int on right
  function multiply_longRight(a, i) result(c)
    class (Rational), intent(in) :: a
    integer(kind=LONG), intent(in) :: i
    type (Rational) :: c

    c%whole = a%whole * i
    c%numerator = a%numerator * i
    c%denominator = a%denominator
    call c%reduce()

  end function multiply_longRight

! Multiply by long int on left
  function multiply_longLeft(i, a) result(c)
    integer(kind=LONG), intent(in) :: i
    class (Rational), intent(in) :: a
    type (Rational) :: c

    c%whole = i * a%whole
    c%numerator = i * a%numerator
    c%denominator = a%denominator
    call c%reduce()

  end function multiply_longLeft
  

! Divide two fractions and reduce to simplest form.
  function divide_fraction(a, b) result(c)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b
    type (Rational) :: c
    
    c = Rational(a%whole*a%denominator + a%numerator, &
         & b%whole*b%denominator + b%numerator) * &
         & Rational(b%denominator, a%denominator)
    call c%reduce()

  end function divide_fraction

! Divide fractions by an double and reduce
  function divide_realDP(a, x) result(c)
    class (Rational), intent(in) :: a
    real(kind=DP), intent(in) :: x
    type (Rational) :: c

    c = a / Rational(x)
    call c%reduce()

  end function divide_realDP

! Divide fractions by an integer and reduce
  function divide_byInt4(a, i) result(c)
    class (Rational), intent(in) :: a
    integer(kind=DEFAULT_INT), intent(in) :: i
    type (Rational) :: c

    c = a / Rational(i)
    call c%reduce()

  end function divide_byInt4

! Divide integer by a fraction
  function divide_intoInt4(i, a) result(c)
    integer(kind=DEFAULT_INT), intent(in) :: i
    class (Rational), intent(in) :: a
    type (Rational) :: c

    c = Rational(i) / a
    call c%reduce()

  end function divide_intoInt4


! Divide integer by a fraction
  function divide_intoInt8(i, a) result(c)
    integer(kind=LONG), intent(in) :: i
    class (Rational), intent(in) :: a
    type (Rational) :: c

    c = Rational(i) / a
    call c%reduce()

  end function divide_intoInt8


! Divide fractions by an integer and reduce
  function divide_byInt8(a, i) result(c)
    class (Rational), intent(in) :: a
    integer(kind=LONG), intent(in) :: i
    type (Rational) :: c

    c = a / Rational(i)
    call c%reduce()

  end function divide_byInt8

  logical function equals_fraction(a, b) result(equals)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b

    equals = (a%whole == b%whole) .and. &
         & (a%numerator == b%numerator) .and. &
         & (a%denominator == b%denominator)

  end function equals_fraction

  logical function equals_intRight(a, b) result(equals)
    class (Rational), intent(in) :: a
    integer, intent(in) :: b

    equals = (a%whole == b) .and. (a%numerator == 0)

  end function equals_intRight

  logical function equals_intLeft(a, b) result(equals)
    integer, intent(in) :: a
    class (Rational), intent(in) :: b

    equals = (b%whole == a) .and. (b%numerator == 0)

  end function equals_intLeft

  logical function lessThan_fraction(r1, r2) result(lessThan)
    class (Rational), intent(in) :: r1
    class (Rational), intent(in) :: r2

    if (r1%whole < r2%whole) then
      lessThan = .true.
    else if (r1%whole > r2%whole) then
      lessThan = .false.
    else
      if (r1%numerator*r2%denominator < r2%numerator*r1%denominator) then
        lessThan = .true.
      else
        lessThan = .false.
      end if
    end if

  end function lessThan_fraction

  logical function greaterThan_fraction(r1, r2) result(greaterThan)
    class (Rational), intent(in) :: r1
    class (Rational), intent(in) :: r2

    if (r1%whole > r2%whole) then
      greaterThan = .true.
    else if (r1%whole < r2%whole) then
      greaterThan = .false.
    else
      if (r1%numerator*r2%denominator > r2%numerator*r1%denominator) then
        greaterThan = .true.
      else
        greaterThan = .false.
      end if
    end if

  end function greaterThan_fraction


  logical function greaterThanOrEqualTo_fraction(r1, r2) result(greaterThanOrEqualTo)
    class (Rational), intent(in) :: r1
    class (Rational), intent(in) :: r2

    if (r1%whole > r2%whole) then
      greaterThanOrEqualTo = .true.
    else if (r1%whole < r2%whole) then
      greaterThanOrEqualTo = .false.
    else
      if (r1%numerator*r2%denominator >= r2%numerator*r1%denominator) then
        greaterThanOrEqualTo = .true.
      else
        greaterThanOrEqualTo = .false.
      end if
    end if

 end function greaterThanOrEqualTo_fraction

  !
  ! Reduce rational number to standard form:
  !    -  |numerator| < denominator
  !    - sign(whole)*sign(numerator) >= 0
  !    - gcd(numerotor,denominator) = 1
  !
  subroutine reduce(this)
    class (Rational), intent(inout) :: this

    integer(kind=Long) :: whole, factor

    associate(w => this%whole, n => this%numerator, d => this%denominator)

      if (d < 0) then
         d = -d
         n = -n
      end if

      if (w > 0 .and. n < 0) then
         w = w - 1
         n = n + d
      else if (w < 0 .and. n > 0) then
         w =  w + 1
         n = n - d
      end if
         
      if (abs(n) >= d) then
         whole = n / d
         w = w + whole
         n = n - whole*d
      end if

      factor = greatestCommonFactor(abs(n), d)
      n = n / factor
      d = d / factor

    end associate

  end subroutine reduce

! Citation - Euclid
  function greatestCommonFactor(a, b) result(factor)
     integer(kind=LONG) :: factor
     integer(kind=LONG), intent(in) :: a
     integer(kind=LONG), intent(in) :: b

     integer(kind=LONG) :: ia, ib, f
     
     if (a==0 .and. b==0) then
        factor = 1
        return
     end if

     ia = a
     ib = b
     f = ia

     do while (f > 0)
        ia = mod(ib, f)
        ib = f
        f = ia
     end do

     factor = ib

  end function greatestCommonFactor

  ! Used for checkpointing.
  function toString(this) result(string)
     character(len=:), allocatable :: string
     class (Rational), intent(in) :: this

     allocate(character(len=60) :: string)

     write(string,'(3(I0,1x))') this%whole, this%numerator, this%denominator
     string = trim(string)

  end function toString

end module Rational_mod
