!----------------------------------------------------------------------
!
! This representation of orbital mechanics is a bit biased by the
! requirements of supporting climate modeling.  In particular these
! interfaces support computation of zenith angle of the primary which
! in turn depends upon obliquity and rotation period in addition to
! traditional Keplerian orbital parameters.
!
! Note that the interfaces provided here do not assume an ideal
! Keplerian orbit.  (See FixedOrbit for that.)  Keplerian parameters
! can evolve slowly over time, which unfortunately means that a Time
! object must be passed to many of the methods that would seemingly be
! returning constants.
!
! The implementation does provide an interface for returning the
! "osculating" orbit - a fancy term for the orbit in the absence of
! perturbations).  Due to the requirements of Fortran, the FixedOrbit
! class must be implemented in the same module.
!
! One interesting feature of this design is that concrete
! implementations must implement the makeCalendar() method.  This
! allows the architecture to enforce that the Calendar and Orbit used
! by the model are consistent with one another.
!
! One awkward aspect of implementing orbits is that the constructors
! will generally have many parameters.  In many intances it may be
! better to use an AttributeDictionary than a long list of parameters.
!
! The implementation of getDeclinationAngle() and getHourAngle() are 
! based upon the original implementation by Gary Russel.  His implementation
! was based upon   V.M.Blanco and  S.W.McCuskey, 1961, 
! "Basic Physics of the Solar System", pages 135 - 151. 
! Note - Existence of Moon and heavenly bodies other than
! Earth and Sun are ignored.  Earth is assumed to be spherical.
!
!----------------------------------------------------------------------

module AbstractOrbit_mod
  use AbstractCalendar_mod
  use KindParameters_mod, only: WP => DP
  use Rational_mod
  use TimeInterval_mod
  implicit none
  private

  public :: AbstractOrbit

  type, abstract :: AbstractOrbit
    private

    type (TimeInterval) :: siderealOrbitalPeriod ! seconds
    type (TimeInterval) :: siderealRotationPeriod ! seconds
    type (TimeInterval) :: meanDay        ! seconds
    real(kind=WP) :: meanDistance  ! in astronomical units (AU's)
    logical :: verbose = .false.

  contains

    ! Primary methods provide information necessary for 
    ! determining incoming radiation from primary.
    procedure :: getSinDeclinationAngle
    procedure :: getDeclinationAngle
    procedure :: getHourAngle
    procedure :: getDistance
    
    ! Secondary methods provide other orbital parameters
    procedure(setYear), deferred :: setYear
    procedure(getSlow), deferred :: getEccentricity
    procedure(getSlow), deferred :: getObliquity
    procedure(getSlow), deferred :: getLongitudeAtPeriapsis
    
    procedure :: getSiderealOrbitalPeriod
    procedure :: getSiderealRotationPeriod
    procedure :: getMeanDay
    procedure :: getMeanDistance
    
    procedure :: setSiderealRotationPeriod
    procedure :: setSiderealOrbitalPeriod
    procedure :: setMeanDay
    procedure :: setMeanDistance
    
    procedure(get), deferred :: getMeanAnomaly
    procedure(get), deferred :: getTrueAnomaly
    
    procedure(makeCalendar), deferred :: makeCalendar

    procedure(print_unit), deferred :: print_unit
    procedure :: print_stdout
    generic :: print => print_unit, print_stdout

    procedure :: rotate

    procedure :: setVerbose
    procedure :: getVerbose

  end type AbstractOrbit


  abstract interface
     

    ! Constant - or nearly constant parameters 
    function getSlow(this) result(q)
      use KindParameters_mod, only: WP => DP
      import AbstractOrbit
      real (kind=WP) :: q
      class (AbstractOrbit), intent(in) :: this
    end function getSlow

    
    ! Modify slowly evolving values
    subroutine setYear(this, year)
      use KindParameters_mod, only: WP => DP
      import AbstractOrbit
      class (AbstractOrbit), intent(inout) :: this
      real(kind=WP), intent(in) :: year
    end subroutine setYear


    function get(this, t) result(q)
      use KindParameters_mod, only: WP => DP
      use BaseTime_mod
      import AbstractOrbit
      real (kind=WP) :: q
      class (AbstractOrbit), intent(in) :: this
      class (BaseTime), intent(in) :: t
    end function get


    function getInterval(this) result(q)
      use KindParameters_mod, only: WP => DP
      use BaseTime_mod
      use TimeInterval_mod, only: TimeInterval
      import AbstractOrbit
      type (TimeInterval) :: q
      class (AbstractOrbit), intent(in) :: this
    end function getInterval


    function makeCalendar(this) result(calendar)
      use AbstractCalendar_mod, only: AbstractCalendar
      import AbstractOrbit
      class (AbstractCalendar), allocatable :: calendar
      class (AbstractOrbit), intent(in) :: this
    end function makeCalendar


    subroutine print_unit(this, unit) 
      import AbstractOrbit
      class (AbstractOrbit), intent(in) :: this
      integer, intent(in) :: unit
    end subroutine print_unit


  end interface

  real(kind=WP), parameter :: PI = 2*asin(1.d0)

contains


  function getSinDeclinationAngle(this, t) result(sinDeclinationAngle)
    use BaseTime_mod, only: BaseTime
    real (kind=WP) :: sinDeclinationAngle
    class (AbstractOrbit), intent(in) :: this
    class (BaseTime), intent(in) :: t

    real (kind=WP) :: trueAnomaly
    real (kind=WP) :: deltaAnomaly
    real (kind=WP) :: sinObliquity
    real (kind=WP), parameter :: pi = 2*asin(1.0d0)

    real (KIND=WP), parameter :: RADIANS_PER_DEGREE = pi/180

    trueAnomaly = this%getTrueAnomaly(t)
    deltaAnomaly = trueAnomaly + this%getLongitudeAtPeriapsis()*RADIANS_PER_DEGREE

    sinObliquity = sin(this%getObliquity() * RADIANS_PER_DEGREE)
    sinDeclinationAngle = sinObliquity * sin(deltaAnomaly)

  end function getSinDeclinationAngle

  function getDeclinationAngle(this, t) result(declinationAngle)
    use BaseTime_mod, only: BaseTime
    real (kind=WP) :: declinationAngle
    class (AbstractOrbit), intent(in) :: this
    class (BaseTime), intent(in) :: t

    declinationAngle = asin(this%getSinDeclinationAngle(t))

  end function getDeclinationAngle


  function getHourAngle(this, t) result(hourAngle)
    use Rational_mod
    use BaseTime_mod
    use TimeInterval_mod
    real(kind=WP) :: hourAngle
    class (AbstractOrbit), intent(in) :: this
    class (BaseTime), intent(in) :: t


    real (kind=WP), parameter :: pi = 2*asin(1.0d0)
    real (kind=WP) :: TA
    real (kind=WP) :: MA

    type (BaseTime) :: t0 ! T=0
    type (Rational) :: days
    real (kind=WP) :: meanHourAngle

    t0 = newBaseTime(0)

    days = (t-t0) / this%getMeanDay()

    meanHourAngle = 2*PI*(real(days) - floor(days))

    TA = this%getTrueAnomaly(t)
    MA = this%getMeanAnomaly(t)

    hourAngle = modulo(meanHourAngle + TA - MA, 2*PI)

  end function getHourAngle

  function getDistance(this, t) result(distance)
    use BaseTime_mod, only: BaseTime
    use OrbitUtilities_mod, only: computeDistance
    real(kind=WP) :: distance
    class (AbstractOrbit), intent(in) :: this
    class (BaseTime), intent(in) :: t

    associate(a => this%meanDistance, e => this%getEccentricity(), &
         & ma => this%getMeanAnomaly(t))
      distance = computeDistance(a, e, ma)
    end associate

  end function getDistance

  subroutine setMeanDistance(this, meanDistance)
    class (AbstractOrbit), intent(inout) :: this
    real(kind=WP), intent(in) :: meanDistance
    this%meanDistance = meanDistance
  end subroutine setMeanDistance


  function getMeanDistance(this) result(meanDistance)
    real(kind=WP) :: MeanDistance
    class (AbstractOrbit), intent(in) :: this
    meanDistance = this%meanDistance
  end function getMeanDistance


  subroutine setSiderealRotationPeriod(this, siderealRotationPeriod)
    use TimeInterval_mod, only: TimeInterval
    class (AbstractOrbit), intent(inout) :: this
    type (TimeInterval) :: siderealRotationPeriod
    this%siderealRotationPeriod = siderealRotationPeriod
  end subroutine setSiderealRotationPeriod


  function getSiderealRotationPeriod(this) result(siderealRotationPeriod)
    use TimeInterval_mod, only: TimeInterval
    type (TimeInterval) :: siderealRotationPeriod
    class (AbstractOrbit), intent(in) :: this
    siderealRotationPeriod = this%siderealRotationPeriod
  end function getSiderealRotationPeriod


  subroutine setMeanDay(this, meanDay)
    use TimeInterval_mod, only: TimeInterval
    class (AbstractOrbit), intent(inout) :: this
    type (TimeInterval) :: meanDay
    this%meanDay = meanDay
  end subroutine setMeanDay


  function getMeanDay(this) result(meanDay)
    use TimeInterval_mod, only: TimeInterval
    type (TimeInterval) :: meanDay
    class (AbstractOrbit), intent(in) :: this
    meanDay = this%meanDay
  end function getMeanDay


  subroutine setSiderealOrbitalPeriod(this, siderealOrbitalPeriod)
    use TimeInterval_mod, only: TimeInterval
    class (AbstractOrbit), intent(inout) :: this
    type (TimeInterval), intent(in) :: siderealOrbitalPeriod
    this%siderealOrbitalPeriod = siderealOrbitalPeriod
 end subroutine setSiderealOrbitalPeriod


  function getSiderealOrbitalPeriod(this) result(siderealOrbitalPeriod)
    use TimeInterval_mod, only: TimeInterval
    type (TimeInterval) :: siderealOrbitalPeriod
    class (AbstractOrbit), intent(in) :: this
    siderealOrbitalPeriod = this%siderealOrbitalPeriod
  end function getSiderealOrbitalPeriod


  ! For diagnsotic purposes - default to stdout
  subroutine print_stdout(this)
    use iso_fortran_env, only: OUTPUT_UNIT
    class (AbstractOrbit), intent(in) :: this

    call this%print(OUTPUT_UNIT)

  end subroutine print_stdout


  ! Useful method for determining times of orbit events such as
  ! solstices and equinoctes.
  ! Returns time after orbit has subtended angle (in radians) from
  ! time t.  
  function rotate(this, t, angle) result(newT)
    use OrbitUtilities_mod, only: computeMeanAnomaly
    use TimeInterval_mod
    use BaseTime_mod
    type (BaseTime) :: newT

    class (AbstractOrbit), intent(in) :: this
    type (BaseTime), intent(in) :: t
    real (kind=WP), intent(in) :: angle
    
    real (kind=WP) :: trueAnomaly
    real (kind=WP) :: M0, M1
    type (TimeInterval) :: tOrbit
    
    trueAnomaly =  this%getTrueAnomaly(t)
    M0 = this%getMeanAnomaly(t)
    M1 = computeMeanAnomaly(trueAnomaly + angle, this%getEccentricity())
    
    tOrbit = this%getSiderealOrbitalPeriod()
    newT = newBaseTime(t + Rational((M1-M0)/(2*PI) * real(tOrbit), 1.d-6))
  end function rotate


  subroutine setVerbose(this, verbose)
     class (AbstractOrbit), intent(inout) :: this
     logical, intent(in) :: verbose

     this%verbose = verbose

  end subroutine setVerbose

  logical function getVerbose(this)
     class (AbstractOrbit), intent(in) :: this
     
     getVerbose = this%verbose

  end function getVerbose


end module AbstractOrbit_mod

