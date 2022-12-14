!---------------------------------------------------------------------
! FixedOrbit extends AbstractOrbit for the common case of an ideal
! Keplerian orbit.  I.e. eccentricity, obliquity, and longitude at
! Periapsis can all be treated as constants.
!
! FixedOrbit is still abstract, as different subclasses will likely
! have alternative constructors, but more importantly will almost
! certainly have different methods for making a corresponding
! Calendar.  E.g. Earth has calendars that while based upon
! astronomical angles, are also imbedded with significant history,
! while one will want Planetary Calendars to be guided by the Earthly
! (Julian) counterparts
!
!---------------------------------------------------------------------


module FixedOrbit_mod
  use AbstractOrbit_mod, only: AbstractOrbit
  use KindParameters_mod, only: WP => DP, DP
  use BaseTime_mod, only: BaseTime
  use TimeInterval_mod, only: TimeInterval
  use Rational_mod
  implicit none
  private

  public :: FixedOrbit

  type, abstract, extends(AbstractOrbit) :: FixedOrbit
    private
  ! TODO: Does this go away?
    type (BaseTime) :: timeAtPeriapsis    ! seconds
    
    real(kind=WP) :: eccentricity
    real(kind=WP) :: obliquity
    real(kind=WP) :: longitudeAtPeriapsis

  contains

    procedure :: setYear

    procedure :: setEccentricity
    procedure :: getEccentricity

    procedure :: setObliquity
    procedure :: getObliquity

    procedure :: setLongitudeAtPeriapsis
    procedure :: getLongitudeAtPeriapsis

    procedure :: setTimeAtPeriapsis
    procedure :: getTimeAtPeriapsis

    procedure :: getMeanAnomaly
    procedure :: getTrueAnomaly

  end type FixedOrbit

  abstract interface

    function get(this, t) result (q)
      use KindParameters_mod, only: WP => DP
      import FixedOrbit
      class (FixedOrbit), intent(in) :: this

      real(kind=WP) :: q

    end function get

  end interface

  real (kind=DP), parameter :: PI = 2*asin(1.d0)

contains


  ! No-op
  subroutine setYear(this, year)
    class (FixedOrbit), intent(inout) :: this
    real(kind=WP), intent(in) :: year
  end subroutine setYear
  

  function getEccentricity(this) result(eccentricity)
    real(kind=WP) :: eccentricity
    class (FixedOrbit), intent(in) :: this
    eccentricity = this%eccentricity
  end function getEccentricity


  function getObliquity(this) result(obliquity)
    real(kind=WP) :: obliquity
    class (FixedOrbit), intent(in) :: this
    obliquity = this%obliquity
  end function getObliquity


  function getLongitudeAtPeriapsis(this) result(longitudeAtPeriapsis)
    real(kind=WP) :: longitudeAtPeriapsis
    class (FixedOrbit), intent(in) :: this
    longitudeAtPeriapsis = this%longitudeAtPeriapsis
  end function getLongitudeAtPeriapsis


  subroutine setEccentricity(this, eccentricity)
    class (FixedOrbit), intent(inout) :: this
    real(kind=WP), intent(in) :: eccentricity
    this%eccentricity = eccentricity
  end subroutine setEccentricity


  subroutine setObliquity(this, obliquity)
    class (FixedOrbit), intent(inout) :: this
    real(kind=WP), intent(in) :: obliquity
    this%obliquity = obliquity
  end subroutine setObliquity


  subroutine setLongitudeAtPeriapsis(this, longitudeAtPeriapsis)
    class (FixedOrbit), intent(inout) :: this
    real(kind=WP), intent(in) :: longitudeAtPeriapsis
    this%longitudeAtPeriapsis = longitudeAtPeriapsis
  end subroutine setLongitudeAtPeriapsis


  function getMeanAnomaly(this, t) result(meanAnomaly)
    use BaseTime_mod, only: BaseTime
    use Rational_mod
    real (kind=WP) :: meanAnomaly
    class (FixedOrbit), intent(in) :: this
    class (BaseTime), intent(in) :: t

    ! TODO: "fraction" is a bad name - Fortran intrinsic
    type (Rational) ::fraction

    type (TimeInterval) :: P

    P = this%getSiderealOrbitalPeriod()
    fraction = modulo(t,P) - modulo(this%timeAtPeriapsis,P)
    fraction = fraction / P

    meanAnomaly = real(fraction) * (2*PI)

  end function getMeanAnomaly


  function getTrueAnomaly(this, t) result(trueAnomaly)
    use BaseTime_mod
    use OrbitUtilities_mod, only: computeTrueAnomaly
    real (kind=WP) :: trueAnomaly
    class (FixedOrbit), intent(in) :: this
    class (BaseTime), intent(in) :: t

    real(kind=WP) :: meanAnomaly

    meanAnomaly = this%getMeanAnomaly(t)
    trueAnomaly = computeTrueAnomaly(meanAnomaly, this%eccentricity)

  end function getTrueAnomaly


  subroutine setTimeAtPeriapsis(this, timeAtPeriapsis)
    use BaseTime_mod, only: BaseTime
    class (FixedOrbit), intent(inout) :: this
    type (BaseTime), intent(in) :: timeAtPeriapsis
    this%timeAtPeriapsis = timeAtPeriapsis
  end subroutine setTimeAtPeriapsis


  function getTimeAtPeriapsis(this) result(timeAtPeriapsis)
    use BaseTime_mod, only: BaseTime
    class (FixedOrbit), intent(in) :: this
    type (BaseTime) :: timeAtPeriapsis
    timeAtPeriapsis = this%timeAtPeriapsis
  end function getTimeAtPeriapsis

end module FixedOrbit_mod
