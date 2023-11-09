!-----------------------------------------------------------------------
! PlanetaryOrbit extends FixedOrbit and is intended to be used for
! planets other than the Earth.  The only functionality added here is
! a constructor for specifying orbital parameters and a trivial
! makeCalendar() method which delegates to the PlanetaryCalendar class.
!-----------------------------------------------------------------------

module PlanetaryOrbit_mod
  use AbstractOrbit_mod
  use FixedOrbit_mod
  use KindParameters_Mod, only: WP => DP, DP
  implicit none
  private

  public :: PlanetaryOrbit

  type, extends(FixedOrbit) :: PlanetaryOrbit
    private
    integer :: foo
  contains
    procedure :: makeCalendar
    procedure :: print_unit
  end type PlanetaryOrbit


  interface PlanetaryOrbit
    module procedure newPlanetaryOrbit
    module procedure newPlanetaryOrbit_fromParams
  end interface PlanetaryOrbit


  real(kind=WP), parameter :: EARTH_LON_AT_PERIHELION = 282.9
  real (kind=DP), parameter :: PI = 2*asin(1.d0)

contains

  function newPlanetaryOrbit_fromParams(planetParams) result(orbit)
     use PlanetaryParams_mod
     type (PlanetaryOrbit) :: orbit
     type (PlanetaryParams), intent(in) :: planetParams

     associate(p => planetParams)
       orbit = PlanetaryOrbit( &
          & p%getObliquity(), &
          & p%getEccentricity(), &
          & p%getLongitudeAtPeriapsis(), &
          & p%getSiderealOrbitalPeriod(), &
          & p%getSiderealRotationPeriod(), &
          & p%getMeanDistance())
     end associate
     
  end function newPlanetaryOrbit_fromParams

  function newPlanetaryOrbit(obliquity, eccentricity, longitudeAtPeriapsis, &
       & siderealPeriod, siderealRotationPeriod, meanDistance) result(orbit)
    use Rational_mod
    use BaseTime_mod
    use TimeInterval_mod
    use OrbitUtilities_mod, only: computeMeanAnomaly
    type (PlanetaryOrbit) :: orbit
    real (kind=WP), intent(in) :: obliquity
    real (kind=WP), intent(in) :: eccentricity
    real (kind=WP), intent(in) :: longitudeAtPeriapsis
    real (kind=WP), intent(in) :: siderealPeriod
    real (kind=WP), intent(in) :: siderealRotationPeriod
    real (kind=WP), intent(in) :: meanDistance

    real (kind=WP) :: meanDay
    type (TimeInterval) :: meanDayInterval
    integer :: daysPerYear
    type (Rational) :: q
    real (kind=WP) :: MA0

    call orbit%setLongitudeAtPeriapsis(longitudeAtPeriapsis)
    call orbit%setObliquity(obliquity)
    call orbit%setEccentricity(eccentricity)

    call orbit%setMeanDistance(meanDistance)

    !--------------------------------------------------------------------------------------
    ! Note sidereal period and rotation period are adjusted to ensure integer days per year
    ! while preserving the length of the mean day.   Other conventions are possible.
    !--------------------------------------------------------------------------------------
    meanDay = 1/(1/siderealRotationPeriod - 1/siderealPeriod)
    daysPerYear = max(1, nint(siderealPeriod / meanDay))
    q=Rational(meanDay, tolerance=1.d-6)
    meanDayInterval = TimeInterval(q)
    call orbit%setMeanDay(meanDayInterval)
    call orbit%setSiderealOrbitalPeriod(TimeInterval(daysPerYear * meanDayInterval))
    call orbit%setSiderealRotationPeriod(TimeInterval(meanDayInterval * Rational(daysPerYear, daysPerYear+1)))
    
    MA0 = computeMeanAnomaly(PI/180*(longitudeAtPeriapsis - EARTH_LON_AT_PERIHELION), &
         & eccentricity)
    call orbit%setTimeAtPeriapsis(newBaseTime(MA0/(2*PI) * (daysPerYear*meanDay)))

  end function newPlanetaryOrbit


  ! Pass instance of self to constructor for PlanetaryCalendar.
  function makeCalendar(this) result(calendar)
    use PlanetaryCalendar_mod
    use AbstractCalendar_mod
    class (AbstractCalendar), allocatable :: calendar
    class (PlanetaryOrbit), intent(in) :: this

    allocate(calendar, source=PlanetaryCalendar(this))
  end function makeCalendar


  subroutine print_unit(this, unit)
    class (PlanetaryOrbit), intent(in) :: this
    integer, intent(in) :: unit

    write(unit,*) 'Fixed orbital parameters for planet.'
    write(unit,*) '  Eccentricity:', this%getEccentricity()
    write(unit,*) '  Obliquity (degs):',this%getObliquity()
    write(unit,*) '  Longitude at periapsis (degs from ve):', &
         & this%getLongitudeAtPeriapsis()

  end subroutine print_unit

end module PlanetaryOrbit_mod
