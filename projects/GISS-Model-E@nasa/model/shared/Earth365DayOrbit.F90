module Earth365DayOrbit_mod
  use KindParameters_mod, only: WP => DP
  use AbstractOrbit_mod, only: AbstractOrbit
  use FixedOrbit_mod, only: FixedOrbit
  use TimeConstants_mod, only: INT_SECONDS_PER_YEAR
  use Rational_mod
  use BaseTime_mod, only: BaseTime, newBaseTime
  use TimeInterval_mod, only: TimeInterval
  use PlanetaryParams_mod, only: DEFAULT_ECCENTRICITY
  use PlanetaryParams_mod, only: DEFAULT_OBLIQUITY
  use PlanetaryParams_mod, only: DEFAULT_LONGITUDE_AT_PERIAPSIS

  implicit none
  private

  public :: Earth365DayOrbit

  type, extends(FixedOrbit) :: Earth365DayOrbit
    private
    character(len=80) :: provenance
    type (BaseTime) :: timeAtVernalEquinox
   contains
     procedure :: makeCalendar
     procedure :: print_unit
  end type Earth365DayOrbit

  interface Earth365DayOrbit
     module procedure newEarth365DayOrbit_params  ! from orbpar()
     module procedure newEarth365DayOrbit_default ! hardwired values
     module procedure newEarth365DayOrbit_config  ! from rundeck
  end interface Earth365DayOrbit


  real(kind=WP), parameter :: PI = 2*asin(1.d0)
  real(kind=WP), parameter :: RADIANS_PER_DEGREE = PI/180

  integer :: EPOCH = 2000

contains

  
  ! Use long-term formulae for orbital parameters in orbpar()
  function newEarth365DayOrbit_config(referenceYear) result(orbit)
    type (Earth365DayOrbit) :: orbit
    real (kind=WP), intent(in) :: referenceYear

    real (kind=WP) :: eccen, obliq, omegvp

    call orbpar(referenceYear, eccen, obliq, omegvp)

    orbit = Earth365DayOrbit(eccen, obliq, omegvp)
    write(orbit%provenance,*)'Fixed orbital params from year = ', &
         & referenceYear, ' CE:'

  end function newEarth365DayOrbit_config


  ! Use default hardwired values for orbital parameters
  function newEarth365DayOrbit_default() result(orbit)
    type (Earth365DayOrbit) :: orbit

    orbit = Earth365DayOrbit(DEFAULT_ECCENTRICITY, DEFAULT_OBLIQUITY, &
         & DEFAULT_LONGITUDE_AT_PERIAPSIS)

    orbit%provenance = 'Fixed orbital params from hardwired defaults:'

  end function newEarth365DayOrbit_default


  function newEarth365DayOrbit_params(eccen, obliq, omegvp) result(orbit)
    use StringUtilities_mod, only: toLowerCase
    use OrbitUtilities_mod, only: computeMeanAnomaly, computeTrueAnomaly
    use JulianCalendar_mod, only: JulianCalendar
    use TimeConstants_mod, only: INT_SECONDS_PER_DAY
    use TimeInterval_mod, only: TimeInterval
    use Rational_mod

    type (Earth365DayOrbit) :: orbit
    real(kind=WP), intent(in) :: eccen
    real(kind=WP), intent(in) :: obliq
    real(kind=WP), intent(in) :: omegvp

    type (JulianCalendar) :: julian
    type (BaseTime) :: timeAtPeriapsis
    type (TimeInterval) :: siderealPeriod
    type (TimeInterval) :: siderealRotationPeriod
    real (kind=WP) :: meanAnomaly
    real (kind=WP) :: trueAnomaly

    call setFixed(eccen=eccen, obliq=obliq, omegvp=omegvp)
    call orbit%setMeanDistance(1.0_WP)
    julian = JulianCalendar()
    trueAnomaly = -orbit%getLongitudeAtPeriapsis() * RADIANS_PER_DEGREE
    meanAnomaly = computeMeanAnomaly(trueAnomaly, orbit%getEccentricity())

    ! Hardwired date for Vernal Equinox:   March 21 12:00
    orbit%timeAtVernalEquinox = julian%convertToTime(year=1, month=3, date=21, hour=12)
    timeAtPeriapsis = newBaseTime(orbit%timeAtVernalEquinox - &
         & Rational((meanAnomaly/(2*PI)) * INT_SECONDS_PER_YEAR, tolerance=1.d-15))

    call orbit%setTimeAtPeriapsis(timeAtPeriapsis)
    siderealPeriod = TimeInterval(Rational(INT_SECONDS_PER_YEAR))
    call orbit%setSiderealOrbitalPeriod(siderealPeriod)

    ! Around the world in 80 days ...
    siderealRotationPeriod = TimeInterval((Rational(INT_SECONDS_PER_DAY)*365)/366)
    call orbit%setSiderealRotationPeriod(siderealRotationPeriod)

    orbit%provenance = 'Orbital parameters set by rundeck:'

  contains
    
    subroutine setFixed(eccen, obliq, omegvp)
      real(kind=WP), intent(in) :: eccen
      real(kind=WP), intent(in) :: obliq
      real(kind=WP), intent(in) :: omegvp
      
      call orbit%setLongitudeAtPeriapsis(omegvp)
      call orbit%setObliquity(obliq)
      call orbit%setEccentricity(eccen)

      call orbit%setMeanDay(TimeInterval(INT_SECONDS_PER_DAY))

    end subroutine setFixed

  end function newEarth365DayOrbit_params

  function makeCalendar(this) result(calendar)
    use AbstractCalendar_mod, only: AbstractCalendar
    use JulianCalendar_mod, only: JulianCalendar
    class (AbstractCalendar), allocatable :: calendar
    class (Earth365DayOrbit), intent(in) :: this

    type (BaseTime) :: vernalEquinox
    type (BaseTime) :: autumnalEquinox
    type (BaseTime) :: winterSolstice
    type (BaseTime) :: summerSolstice

    allocate(calendar, source=JulianCalendar())

    ! Add orbital dates
    vernalEquinox  = this%timeAtVernalEquinox
    summerSolstice = this%rotate(vernalEquinox, PI/2)
    autumnalEquinox = this%rotate(vernalEquinox, PI)
    winterSolstice = this%rotate(vernalEquinox, 3*PI/2)

    call calendar%addTransitionDate('vernal equinox', &
         & calendar%getAnniversaryDate(vernalEquinox))
    call calendar%addTransitionDate('autumnal equinox', &
         & calendar%getAnniversaryDate(autumnalEquinox))
    call calendar%addTransitionDate('winter solstice', &
         & calendar%getAnniversaryDate(winterSolstice))
    call calendar%addTransitionDate('summer solstice', &
         & calendar%getAnniversaryDate(summerSolstice))

    
  end function makeCalendar
    
  subroutine print_unit(this, unit)
    class (Earth365DayOrbit), intent(in) :: this
    integer, intent(in) :: unit

    write(unit,'(a)') trim(this%provenance)

    write(unit,*) '  Eccentricity:', this%getEccentricity()
    write(unit,*) '  Obliquity (degs):',this%getObliquity()
    write(unit,*) '  Longitude at periapsis (degs from ve):', &
         & this%getLongitudeAtPeriapsis()

  end subroutine print_unit

end module Earth365DayOrbit_mod
