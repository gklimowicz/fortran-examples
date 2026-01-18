module Geometry_mod
!@sum Provides resolution dependend geometric quantities.
!@+   replaces functionality in various GEOM*.f files.
!@auth Tom Clune <thomas.l.clune@nasa.gov> 
  use constant, only: pi, twopi, radian
  implicit none
  private

  public :: Geometry_type
  public :: newGeometry
  public :: newGeometry4x5
  public :: newGeometry8x10
  public :: getIndicesFromLatLon

  type GridSpacing_type
    private
    real*8 :: radians
    real*8 :: degrees
    real*8 :: minutes
    integer :: numCells
  end type GridSpacing_type

  type Geometry_type
    private
    integer :: numLongitudes
    integer :: numLatitudes
    integer :: numLevels

    type (GridSpacing_type) :: latSpacing
    type (GridSpacing_type) :: lonSpacing

  end type Geometry_type

  interface newGeometry
    module procedure newGeometryDefault
    module procedure newGeometryPoleFraction
!!$    module procedure newGeometry8x10Default
  end interface

  integer, parameter :: degreesToMinutes = 60
  real, parameter :: degreesToRadians = radian
  real, parameter :: half = 1/2.d+0

contains

  function newGeometryDefault(numLongitudes, numLatitudes, numLevels) result(geometry)
!@sum Construct new geometry object.
    type (Geometry_type) :: geometry
    
    integer, intent(in) :: numLongitudes
    integer, intent(in) :: numLatitudes
    integer, intent(in) :: numLevels

    ! full cell at poles
    geometry = newGeometry(numLongitudes, numLatitudes, numLevels, 1.0d0) 

  end function newGeometryDefault

  function newGeometryPoleFraction(numLongitudes, numLatitudes, numLevels, &
       & poleFraction) result(geometry)
!@sum Construct new geometry object with special allowance for small
!@+   cell at poles. This special treatment at the poles is required
!@+   to support 1/2 cell for 4x5 and 1/4 cell for 8x10.
    type (Geometry_type) :: geometry
    
    integer, intent(in) :: numLongitudes
    integer, intent(in) :: numLatitudes
    integer, intent(in) :: numLevels
    real*8, intent(in) :: poleFraction ! in fraction of grid cell
    
    geometry%numLongitudes = numLongitudes
    geometry%numLatitudes  = numLatitudes
    geometry%numLevels     = numLevels
    
    call setLongitudeSpacing(geometry%lonSpacing, numLongitudes)
    call setLatitudeSpacing(geometry%latSpacing, numLatitudes, poleFraction)
    
  contains

    subroutine setLongitudeSpacing(this, numLongitudes)
      type (GridSpacing_type), intent(inout) :: this
      integer, intent(in) :: numLongitudes
      
      this%numCells = numLongitudes
      this%degrees = 360.d+0/numLongitudes
      this%radians = this%degrees * degreesToRadians
      this%minutes = this%degrees * degreesToMinutes
      
    end subroutine setLongitudeSpacing
    
    subroutine setLatitudeSpacing(this, numLatitudes, poleFraction)
      type (GridSpacing_type), intent(inout) :: this
      integer, intent(in) :: numLatitudes
      real*8, intent(in) :: poleFraction
      
      this%numCells = numLatitudes
      this%degrees = 180.d+0/(numLatitudes - 2 + 2*poleFraction)
      this%radians = this%degrees * degreesToRadians
      this%minutes = this%minutes * degreesToMinutes

    end subroutine setLatitudeSpacing
    
  end function newGeometryPoleFraction

  function newGeometry4x5(numLevels) result(geometry)
!@sum Construct new geometry object.
    type (Geometry_type) :: geometry
    integer, intent(in) :: numLevels

    integer, parameter :: numLongitudes = 72
    integer, parameter :: numLatitudes = 46

    ! half cell at poles
    geometry = newGeometry(numLongitudes, numLatitudes, numLevels, &
         & poleFraction=0.5d+0)

  end function newGeometry4x5

  function newGeometry8x10(numLevels) result(geometry)
!@sum Construct new geometry object.
    type (Geometry_type) :: geometry
    integer, intent(in) :: numLevels

    integer, parameter :: numLongitudes = 36
    integer, parameter :: numLatitudes = 24

    ! half cell at poles
    geometry = newGeometry(numLongitudes, numLatitudes, numLevels, &
         & poleFraction=0.25d+0)

  end function newGeometry8x10

  function getIndicesFromLatLon(this, degreesLongitude, degreesLatitude) &
       & result(ij)
!@sum Retrieve the grid cell indices for a given longitude and latitude 
!@+expressed in degrees.
    type (Geometry_type) :: this
    real*8, intent(in) :: degreesLongitude ! degrees (not radians)
    real*8, intent(in) :: degreesLatitude  ! degrees (not radians)
    integer :: ij(2)

    ij(1) = getIndex(this%lonSpacing, degreesLongitude)
    ij(2) = getIndex(this%latSpacing, degreesLatitude)

  contains

    integer function getIndex(this, coordinateInDegrees) result(idx)
      type (GridSpacing_type) :: this
      real*8, intent(in) :: coordinateInDegrees ! degrees (not radians)

      real*8 :: offset

      offset = HALF*(1 + this%numCells)
      idx = nint(offset + coordinateInDegrees/this%degrees)

      idx = max(idx, 1)
      idx = min(idx, this%numCells)

    end function getIndex

  end function getIndicesFromLatLon


end module Geometry_mod
