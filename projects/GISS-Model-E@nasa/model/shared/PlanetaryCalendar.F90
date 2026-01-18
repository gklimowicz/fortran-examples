
module PlanetaryCalendar_mod
  use AbstractCalendar_mod
  use FixedCalendar_mod
  use FixedOrbit_mod
  use AbstractOrbit_mod
  use CalendarMonth_mod, only: CalendarMonth
  use KindParameters_mod, only: WP => DP
  use TimeInterval_mod
  use Rational_mod
  implicit none
  private

  public :: PlanetaryCalendar

  type, extends(FixedCalendar) :: PlanetaryCalendar
     private
     integer :: foo
   contains
  end type PlanetaryCalendar

  interface PlanetaryCalendar
     module procedure newPlanetaryCalendar_likeJulian
     module procedure newPlanetaryCalendar_keplerian
     module procedure newPlanetaryCalendar_longitudes
  end interface PlanetaryCalendar

  real (kind=WP), parameter :: PI = 2*asin(1.d0)

  integer, parameter :: MIN_DAYS_PER_YEAR = 120

contains


  ! Create a calendar that preserves the longitudes of the month
  ! starts relative to the vernal equinox.  This choice results in a
  ! correspondence between months and seasons that will be familiar to
  ! climate scientists.  Of course if the eccentricity is large and
  ! the periapsis is in the Northern summer, the results might be a
  ! bit counterintuitive

  ! Note that if the rotational period of a planet is sufficiently
  ! slow this approach may result in months that have zero days.  The
  ! code will issue a message and terminate under such conditions.
  ! The most obvious case for this would be for a tidally locked
  ! planet.

  function newPlanetaryCalendar_likeJulian(orbit) result(calendar)
    use AbstractOrbit_mod
    use FixedOrbit_mod
    use Earth365DayOrbit_mod, only: Earth365DayOrbit
    type (PlanetaryCalendar) :: calendar
    class (FixedOrbit), intent(in) :: orbit

    calendar = PlanetaryCalendar(orbit, Earth365DayOrbit())
    
  end function newPlanetaryCalendar_likeJulian


  function newPlanetaryCalendar_keplerian(orbit, referenceOrbit) result(calendar)
    use AbstractOrbit_mod
    use FixedOrbit_mod
    use AbstractCalendar_mod, only: AbstractCalendar
    use Rational_mod
    type (PlanetaryCalendar) :: calendar
    class (FixedOrbit), intent(in) :: orbit
    class (FixedOrbit), intent(in) :: referenceOrbit

    real (kind=WP) :: TA, TA0
    class (AbstractCalendar), allocatable :: referenceCalendar
    integer :: i
    real (kind=WP) :: longitudes(MONTHS_PER_YEAR)

    allocate(referenceCalendar, source=referenceOrbit%makeCalendar())

    TA0 = -referenceOrbit%getLongitudeAtPeriapsis()*(PI/180)

    do i = 1, MONTHS_PER_YEAR
      TA = referenceOrbit%getTrueAnomaly(referenceCalendar%convertToTime(1,i,1,0))
      longitudes(i) = modulo(TA-TA0, 2*PI)
    end do

    calendar = PlanetaryCalendar(orbit, longitudes)

    deallocate(referenceCalendar)
      
  end function newPlanetaryCalendar_keplerian


  ! This intermediate constructor uses the longitudes of month starts as 
  ! determiners of month lengths for the specified orbit.
  function newPlanetaryCalendar_longitudes(orbit, monthLongitudes) result(calendar)
    use OrbitUtilities_mod, only: computeMeanAnomaly
    use AbstractOrbit_mod
    use FixedOrbit_mod
    use JulianCalendar_mod, only: JULIAN_MONTHS
    use OrbitUtilities_mod, only: computeMeanAnomaly
    use Rational_mod
    use TimeInterval_mod
    use BaseTime_mod
    type (PlanetaryCalendar) :: calendar
    class (FixedOrbit), intent(in) :: orbit
    real(kind=WP), intent(in) :: monthLongitudes(:) ! in radians

    real(kind=WP) :: M     ! mean anomaly
    real(kind=WP) :: M0    ! mean anomaly for Jan 01
    real(kind=WP) :: trueAnomaly
    integer :: i
    real(kind=WP) :: Lp
    type (CalendarMonth) :: months(0:MONTHS_PER_YEAR+1)
    type (BaseTime) :: tVE
    integer :: daysPerYear
    
    call calendar%setSecondsPerDay(orbit%getMeanDay())

    daysPerYear = nint(orbit%getSiderealOrbitalPeriod()/orbit%getMeanDay())
    ! Calendar does not care whether an orbit is retrograde.
    daysPerYear = abs(daysPerYear)

    ! Require a minimum number of days per year. (Suggested by G. Schmidt.)
    if (daysPerYear < MIN_DAYS_PER_YEAR) then
       if (calendar%getVerbose()) then
          write(*,*) '***********************************************************'
          write(*,*) '* Warning calendar days do not correspond to solar days.  *'
          write(*,*) '* Hourly diagnostics should not be used.                  *'
          write(*,*) '***********************************************************'
       end if
       daysPerYear = MIN_DAYS_PER_YEAR
       call calendar%setSecondsPerDay( &
           & TimeInterval( orbit%getSiderealOrbitalPeriod() / daysPerYear ))
    end if
    
    call calendar%setDaysPerYear(daysPerYear)

    if (size(monthLongitudes) /= MONTHS_PER_YEAR) then
       call stop_model('PlanetaryCalendar assumes 12 month years.')
       return
    end if

    Lp = orbit%getLongitudeAtPeriapsis()*(PI/180)
    trueAnomaly = monthLongitudes(1) - Lp
    M0 = computeMeanAnomaly(trueAnomaly, orbit%getEccentricity())

    months%fullName = JULIAN_MONTHS%fullName
    months%abbreviation = JULIAN_MONTHS%abbreviation

    ! First pass - just determine first day in month
    months(1)%firstDayInMonth = 1
    months(MONTHS_PER_YEAR+1)%firstDayInMonth = 1 + calendar%getDaysInYear()
    do i = 2, MONTHS_PER_YEAR

       trueAnomaly = monthLongitudes(i) - Lp
       M     = modulo(computeMeanAnomaly(trueAnomaly, orbit%getEccentricity()) - M0, 2*PI)

       months(i)%firstDayInMonth = 1 + nint(M * calendar%getDaysInYear() / (2*PI))

    end do
    months(0)%firstDayInMonth = months(MONTHS_PER_YEAR)%firstDayInMonth - calendar%getDaysInYear()


    ! second pass - remaining data
    months(0)%lastDayInMonth = 0
    months(MONTHS_PER_YEAR)%lastDayInMonth = calendar%getDaysInYear()

    do i = 1, MONTHS_PER_YEAR - 1
       months(i)%lastDayInMonth = months(i+1)%firstDayInMonth - 1
    end do
    months(MONTHS_PER_YEAR+1)%lastDayInMonth = months(1)%lastDayInMonth + calendar%getDaysInYear()

    ! And finally mid day, days per,  ...
    months(:)%midDayInMonth = (months(:)%firstDayInMonth + months(:)%lastDayInMonth)/2
    months(:)%daysInMonth = 1 + months(:)%lastDayInMonth - months(:)%firstDayInMonth

    do i = 1, MONTHS_PER_YEAR
      call calendar%setNthCalendarMonth(i, months(i))
    end do


    ! Special dates

    call calendar%initTransitionDates()

    tVE = orbit%rotate(calendar%convertToTime(1,1,1,0), -monthLongitudes(1))

    call calendar%addTransitionDate('periapsis', calendar%getAnniversaryDate(orbit%rotate(tVE, Lp)))
    call calendar%addTransitionDate('apsis', calendar%getAnniversaryDate(orbit%rotate(tVE, Lp+PI)))


    call calendar%addTransitionDate('vernal equinox', calendar%getAnniversaryDate(orbit%rotate(tVE, 0.d0)))
    call calendar%addTransitionDate('summer solstice', calendar%getAnniversaryDate(orbit%rotate(tVE, PI/2)))
    call calendar%addTransitionDate('autumnal equinox', calendar%getAnniversaryDate(orbit%rotate(tVE, PI)))
    call calendar%addTransitionDate('winter solstice', calendar%getAnniversaryDate(orbit%rotate(tVE, 3*PI/2)))

 end function newPlanetaryCalendar_longitudes

end module PlanetaryCalendar_mod
