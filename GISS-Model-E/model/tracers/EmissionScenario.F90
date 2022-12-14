module EmissionScenario_mod
  use EmissionRegion_mod
  use domain_decomp_atm, only: dist_grid
  implicit none
  private

  public :: EmissionScenario
  public :: newEmissionScenario
  public :: clean

  integer, parameter :: MAX_LEN_SECTOR_NAME = 10

  type EmissionSector ! internal data structure
    character(len=MAX_LEN_SECTOR_NAME) :: name ! set, but not currently used
    real*8, pointer :: scalingFactors(:)
  end type EmissionSector

  type EmissionScenario
    type (EmissionRegion), pointer :: regions(:) => null()
    type (EmissionSector), pointer :: sectors(:) => null()
    ! these store references for later computation
    type (dist_grid), pointer :: grid => null()
    real*8, pointer :: latitudes(:,:) => null()
    real*8, pointer :: longitudes(:,:) => null()
  contains
    procedure :: appendSector
    procedure :: extendSectorList
    procedure :: scaleSource
  end type EmissionScenario

contains

  function newEmissionSector(name, scalingFactors) result(sector)
    character(len=*), intent(in) :: name
    real*8, intent(in) :: scalingFactors(:)
    type (EmissionSector), pointer :: sector

    allocate(sector)

    sector%name = name
    allocate(sector%scalingFactors(size(scalingFactors)))
    sector%scalingFactors = scalingFactors

  end function newEmissionSector

  function newEmissionScenario(regions, grid, latitudes, longitudes) result(scenario)
    use domain_decomp_atm, only: dist_grid
    type (EmissionScenario) :: scenario
    type (EmissionRegion), intent(in) :: regions(:)
    type (dist_grid), target :: grid
!TODO - these should be part of some structure that is passed in; possibly dist_grid
    real*8, target, intent(in) :: latitudes(:,:)
    real*8, target, intent(in) :: longitudes(:,:)

    allocate(scenario%regions(size(regions)))
    scenario%regions = regions
    scenario%grid => grid
    scenario%latitudes => latitudes
    scenario%longitudes => longitudes

    allocate(scenario%sectors(0))

  end function newEmissionScenario

  subroutine appendSector(this, sectorName, scalingFactors)
    class (EmissionScenario), intent(inout) :: this
    character(len=*), intent(in) :: sectorName
    real*8, intent(in) :: scalingFactors(:)

    type (EmissionSector), pointer :: sector

    if (size(scalingFactors) /= size(this%regions)) then
      call stop_model('Incorrect number of scaling factors for this emission scenario.',999)
    end if

    sector => newEmissionSector(sectorName, scalingFactors)
    call this%extendSectorList(sector)
    deallocate(sector)

  end subroutine appendSector

  subroutine extendSectorList(this, sector)
    class (EmissionScenario), intent(inout) :: this
    type (EmissionSector), intent(in) :: sector

    integer :: n
    type (EmissionSector), allocatable :: tmp(:)
    
! TODO - use F2003 reallocate facility
    n = size(this%sectors)
    allocate(tmp(n))

    tmp(1:n) = this%sectors(1:n)
    deallocate(this%sectors)

    allocate(this%sectors(n+1))
    this%sectors(1:n) = tmp(1:n)
    deallocate(tmp)

    this%sectors(n+1) = sector
    
  end subroutine extendSectorList

  subroutine scaleSource(this, trSource, sectorIndices)
    use domain_decomp_atm, only: getDomainBounds
    class (EmissionScenario), intent(in) :: this
    real*8, intent(inout) :: trSource(this%grid%i_strt_halo:,this%grid%j_strt_halo:)
    integer, intent(in) :: sectorIndices(:)

    integer :: i, j
    integer :: iSector, iRegion
    real*8 :: scale

    integer :: j_0, j_1
    integer :: i_0, i_1
    real*8 :: lat, lon


    call getDomainBounds(this%grid, j_strt=j_0, j_stop=j_1, i_strt = i_0, i_stop = i_1)

    do iSector = 1, size(sectorIndices)
      do iRegion = 1, size(this%regions)

        scale = this%sectors(iSector)%scalingFactors(iRegion)
        if (scale == 1.d0) cycle ! 

        do j = j_0, j_1
          do i = i_0, i_1
            lat = this%latitudes(i,j)
            lon = this%longitudes(i,j)
            if (this%regions(iRegion)%hasLatLon(lat,lon)) then
              trSource(i,j) = trSource(i,j) * scale
            end if
          end do
        end do

      end do
    end do

  end subroutine scaleSource

  subroutine clean(this)
    type (EmissionScenario), intent(inout) :: this
    integer :: i

    do i = 1, size(this%sectors)
      deallocate(this%sectors(i)%scalingFactors)
    end do
    deallocate(this%sectors)

    deallocate(this%regions) ! possibly need to clean inside of regions?
  end subroutine clean

end module EmissionScenario_mod
