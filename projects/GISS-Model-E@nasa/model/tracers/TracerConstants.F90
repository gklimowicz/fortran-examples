module TracerConstants_mod
  use CONSTANT, only: molMassWater => mwat
  implicit none
  private

  ! Usage:
  !    mm_CO2 = TRACER_CONSTANTS%CO2%molMass

  public :: TRACER_CONSTANTS
  public :: H2O18

  type Properties
    real*8 :: molMass    ! g
    real*8 :: solubility = 0 ! units?
  end type Properties

  type (Properties), parameter ::          &
       CO2   = Properties(44.d0, 0.d0),    &
       O3    = Properties(48.d0, 0.d0),    &
       water = Properties(18.015d0, 0.d0), &
       H2O18 = Properties(20.d0, 0.d0),    &
       HDO   = Properties(19d0, 0.d0),     &
       H2O17 = Properties(19d0, 0.d0)

  type TracerMetadata
    type (Properties) :: CO2
    type (Properties) :: O3
    type (Properties) :: water
    type (Properties) :: H2O18
    type (Properties) :: HDO
    type (Properties) :: H2O17
  end type TracerMetadata


  ! This syntax must wait until F2003 features are available in the default ifort
  ! and gfortran compilers
!!$$  type (TracerMetadata), parameter :: TRACER_CONSTANTS = TracerMetadata( &
!!$$       & CO2   = Properties(molMass = 44.d0),    &
!!$$       & O3    = Properties(molMass = 48.d0),    &
!!$$       & water = Properties(molMass = 18.015d0), &
!!$$       & H2O18 = H2O18, &
!!$$       & HDO   = Properties(molMass = 19d0),     &
!!$$       & H2O17 = Properties(molMass = 19d0)      &
!!$$       & )
!!$$

  type (TracerMetadata), parameter :: TRACER_CONSTANTS = &
       TracerMetadata( CO2, O3, water, H2O18, HDO, H2O17)


end module TracerConstants_mod

