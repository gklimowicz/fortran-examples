module PlanetaryParams_mod
!@type PlanetaryParams a derived type that manages physical parameters
!@+    that are directly dependent on the planet being simulated.
!@+    Default values are those for modern Earth as specified in AR5.
!@+    
!@+    Note that data components here are encapsulated to ensure that 
!@+    they are not accidentally modified after being set. 
!@+    It is expected that a higher container (e.g. Constants_mod) will
!@+    provide these parameters and various related quantities using
!@+    PROTECTED to prevent accidental change during a simulation.

   use KindParameters_mod, only: WP => DP
   use TimeConstants_mod, only: SECONDS_PER_DAY
   use TimeConstants_mod, only: SECONDS_PER_YEAR
   use TimeConstants_mod, only: DAYS_PER_YEAR
   implicit none
   private
   
   public :: PlanetaryParams

   ! Earth defaults

   public :: DEFAULT_ECCENTRICITY
   public :: DEFAULT_OBLIQUITY
   public :: DEFAULT_LONGITUDE_AT_PERIAPSIS
   public :: DEFAULT_SIDEREAL_ORBITAL_PERIOD
   public :: DEFAULT_SIDEREAL_ROTATION_PERIOD

!@par DEFAULT_ECCENTRICITY
   real (kind=WP), parameter :: DEFAULT_ECCENTRICITY = 0.0167D0
!@par DEFAULT_OBLIQUITY - degrees
   real (kind=WP), parameter :: DEFAULT_OBLIQUITY = 23.44D0
!@par DEFAULT_LONGITUDE_AT_PERIAPSIS - degrees
   real (kind=WP), parameter :: DEFAULT_LONGITUDE_AT_PERIAPSIS = 282.9D0
!@par DEFAULT_SIDEREAL_ORBITAL_PERIOD - seconds
   real (kind=WP), parameter :: DEFAULT_SIDEREAL_ORBITAL_PERIOD = &
        & SECONDS_PER_YEAR ! seconds
!@par DEFAULT_SIDEREAL_ROTATION_PERIOD - sidereal rotation (seconds)
   real (kind=WP), parameter :: DEFAULT_SIDEREAL_ROTATION_PERIOD = &
        & SECONDS_PER_DAY * DAYS_PER_YEAR/(DAYS_PER_YEAR+1) ! seconds

   type PlanetaryParams
      private

      ! Orbital parameters
      real (kind=WP) :: obliquity = DEFAULT_OBLIQUITY
      real (kind=WP) :: eccentricity = DEFAULT_ECCENTRICITY
      real (kind=WP) :: longitudeAtPeriapsis = DEFAULT_LONGITUDE_AT_PERIAPSIS
      real (kind=WP) :: siderealOrbitalPeriod = DEFAULT_SIDEREAL_ORBITAL_PERIOD
      real (kind=WP) :: meanDistance  = 1.0 ! A.U.

      ! Intrinsic parameters
      real (kind=WP) :: siderealRotationPeriod = SECONDS_PER_DAY * &
           &  (DAYS_PER_YEAR/(1+DAYS_PER_YEAR))  ! seconds
!!$      real (kind=WP) :: radius          ! meters
!!$      real (kind=WP) :: gravity         ! m s^-2

   contains

      ! Orbital parameters
      procedure :: getObliquity
      procedure :: getEccentricity
      procedure :: getLongitudeAtPeriapsis
      procedure :: getSiderealOrbitalPeriod
      procedure :: getMeanDistance

      ! Intrinsic parameters
      procedure :: getSiderealRotationPeriod
!!$      procedure :: getRadius
!!$      procedure :: getGravity

   end type PlanetaryParams


   interface PlanetaryParams
      module procedure newPlanetaryParams_fromRundeck
   end interface PlanetaryParams

contains

   function newPlanetaryParams_fromRundeck()  result(params)
      use Dictionary_mod
      type (PlanetaryParams) :: params

      call sync_param('obliquity', params%obliquity)
      call sync_param('eccentricity', params%eccentricity)
      call sync_param('longitudeAtPeriapsis', params%longitudeAtPeriapsis)
      call sync_param('siderealOrbitalPeriod', params%siderealOrbitalPeriod)
      call sync_param('siderealRotationPeriod', params%siderealRotationPeriod)
      call sync_param('meanDistance', params%meanDistance)
      
   end function newPlanetaryParams_fromRundeck


   real (kind=WP) function getObliquity(this) result(obliquity)
      class (PlanetaryParams), intent(in) :: this
      obliquity = this%obliquity
   end function getObliquity


   real (kind=WP) function getEccentricity(this) result(eccentricity)
      class (PlanetaryParams), intent(in) :: this
      eccentricity = this%eccentricity
   end function getEccentricity


   real (kind=WP) function getLongitudeAtPeriapsis(this) result(longitudeAtPeriapsis)
      class (PlanetaryParams), intent(in) :: this
      longitudeAtPeriapsis = this%longitudeAtPeriapsis
   end function getLongitudeAtPeriapsis


   real (kind=WP) function getSiderealOrbitalPeriod(this) result(siderealOrbitalPeriod)
      class (PlanetaryParams), intent(in) :: this
      siderealOrbitalPeriod = this%siderealOrbitalPeriod
   end function getSiderealOrbitalPeriod
   

   real (kind=WP) function getMeanDistance(this) result(meanDistance)
      class (PlanetaryParams), intent(in) :: this
      meanDistance = this%meanDistance
   end function getMeanDistance


   real (kind=WP) function getSiderealRotationPeriod(this) result(siderealRotationPeriod)
      class (PlanetaryParams), intent(in) :: this
      siderealRotationPeriod = this%siderealRotationPeriod
   end function getSiderealRotationPeriod


!!$   real (kind=WP) function getRadius(this) result(radius)
!!$      class (PlanetaryParams), intent(in) :: this
!!$      radius = this%radius
!!$   end function radius
!!$
!!$
!!$   real (kind=WP) function getGravity(this) result(gravity)
!!$      class (PlanetaryParams), intent(in) :: this
!!$      gravity = this%gravity
!!$   end function gravity
!!$

end module PlanetaryParams_mod
