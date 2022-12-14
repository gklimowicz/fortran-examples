#include "rundeck_opts.h"

module PlanetParams_mod
#ifdef PLANET_PARAMS
!@type PlanetParams_t a derived type for the collection of
!@+    parameters users may want to change when simulating
!@+    environments other than Earth.  When this module is
!@+    in use, the corresponding parameters in the constants
!@+    module take their values from one of the instances of
!@+    PlanetParams_t declared below.
  type PlanetParams_t
!@param name name of planet being simulated.
!@+   Not intended to be widely used to control physics settings etc.
!@+   May be removed in future.
     character(len=16) :: name
!@param grav gravitational acceleration (m/s2)
     real*8 :: grav
!@param radius planetary radius (m)
     real*8 :: radius
!@param mair molar mass of "dry" (fixed-composition) air (g/mol)
     real*8 :: mair
!@param srat ratio of specific heats at const. press. and vol. (for fixed-composition air)
     real*8 :: srat
!@param psf (hPa) global mean suface pressure
     real*8 :: psf
!@param ptop (hPa) the pressure level at which the GCM vertical coordinate should
!@+     transition from terrain-following to constant-pressure.  While the details
!@+     of a GCM vertical coordinate are not a fundamental planetary property, the
!@+     ratio ptop/psf is one nondimensional measure of the height of the terrain
!@+     relative to the depth of the atmosphere.   Once the GCM grid information
!@+     is set at run time rather than compile time, this parameter may be determined
!@+     from the topography file, but it must be set here for the moment.
     real*8 :: ptop
  end type PlanetParams_t

!
! Declare the values of parameters for various planets.
!
  type(PlanetParams_t), parameter, private :: &
       notEarth = &
       PlanetParams_t( &
       name = 'notEarth', &
       grav = 9.80665d0, &
       radius = 6371000d0, &
       mair = 28.9655d0, &
       srat = 1.401d0, &
       psf = 984d0, &
       ptop = 150d0 &
       )
!
  type(PlanetParams_t), parameter, private :: &
       likeMars = &
       PlanetParams_t( &
       name = 'likeMars', &
       grav = 3.711d0, &
       radius = 6371000d0*.531d0, &
       mair = 44.01d0, &  ! pure CO2
       srat = 1.31d0, &   ! 0 C value at en.wikipedia.org/wiki/Carbon_dioxide_(data_page)
       psf = 6.36d0, &
       ptop = .3d0 & ! Olympus Mons forces this small ptop/psf ratio.
       )

!@param PlanetParams the instance of PlanetParams_t to be used. CPP selects
!@+     which instance of PlanetParams_t is copied into PlanetParams.
  type(PlanetParams_t), parameter :: PlanetParams = PLANET_PARAMS
#endif /* PLANET_PARAMS defined */
end module PlanetParams_mod
