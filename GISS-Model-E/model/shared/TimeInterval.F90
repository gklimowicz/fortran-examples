module TimeInterval_mod
  use BaseTime_mod, only: BaseTime, real
  implicit none
  private

  public :: TimeInterval
  public :: real ! from BaseTime and Rational

  type, extends(BaseTime) :: TimeInterval
  end type TimeInterval

  interface TimeInterval
     module procedure newTimeInterval_rational
     module procedure newTimeInterval_integer
  end interface TimeInterval

contains

  function newTimeInterval_rational(r) result(interval)
    use Rational_mod
    type (TimeInterval) :: interval
    type (Rational), intent(in) :: r
    interval%Rational = r
  end function newTimeInterval_rational

  function newTimeInterval_integer(i) result(interval)
    use Rational_mod
    type (TimeInterval) :: interval
    integer, intent(in) :: i
    interval%Rational = Rational(i)
  end function newTimeInterval_integer

end module TimeInterval_mod
