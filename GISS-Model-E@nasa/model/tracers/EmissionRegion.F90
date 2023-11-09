module EmissionRegion_mod
  implicit none
  private

  public :: initializeEmissionsRegions
  public :: EmissionRegion
  public :: newEmissionRegion

!TODO - encapsulate these
  public :: regions
  public :: numRegions

  integer, parameter :: NUM_MAX_REGIONS = 10

  type EmissionRegion
    real*8 :: northernEdge
    real*8 :: southernEdge
    real*8 :: easternEdge
    real*8 :: westernEdge
  contains
    procedure :: hasLatLon
  end type EmissionRegion

  type (EmissionRegion), save :: regions(NUM_MAX_REGIONS)
  integer, save :: numRegions

contains

  ! see how many regions there are, save names in array:
  subroutine initializeEmissionsRegions()
    use Dictionary_mod, only : sync_param
    integer, parameter :: MAX_LEN_REGIONS=124
    character(len=MAX_LEN_REGIONS) :: regions_are
    real*8 :: tmp(NUM_MAX_REGIONS)

    regions_are=' '
    call sync_param("regions_are",regions_are)

    print*,__FILE__,__LINE__,trim(regions_are)
    numRegions = countRegions(regions_are)

    if (numRegions > NUM_MAX_REGIONS) call stop_model("n_max_reg must be increased",255)

    if (numRegions > 0) then
      regions(1:numRegions)%northernEdge = getEdges("REG_N", numRegions)
      regions(1:numRegions)%southernEdge = getEdges("REG_S", numRegions)
      regions(1:numRegions)%easternEdge = getEdges("REG_E", numRegions)
      regions(1:numRegions)%westernEdge = getEdges("REG_W", numRegions)
    end if

    contains

      ! Count the number of space-delimited region names in a string
      integer function countRegions(regions) result(numRegions)
        character(len=*), intent(in) :: regions
        integer :: i, j
        
        numRegions = 0
        i = 1
        do while(i < len(regions))
          j = index(regions_are(i:len(regions))," ")
          if (j > 1) then
            numRegions = numRegions + 1
            i=i+j
          else
            i=i+1
          end if
        enddo
    end function countRegions

      function getEdges(name, n) result(edges)
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real*8 :: edges(n)

        real*8 :: tmp(n)
        real*8, parameter :: UNDEFINED = -1.e30

        tmp = UNDEFINED
        call sync_param(trim(name), edges, n)
      end function getEdges

  end subroutine initializeEmissionsRegions

  function newEmissionRegion(lowerLatLon, upperLatLon) result(region)
    real*8, intent(in) :: lowerLatLon(2)
    real*8, intent(in) :: upperLatLon(2)
    type (EmissionRegion) :: region

    region%northernEdge = upperLatLon(1)
    region%easternEdge = upperLatLon(2)

    region%southernEdge = lowerLatLon(1)
    region%westernEdge = lowerLatLon(2)
    
  end function newEmissionRegion


  logical function hasLatLon(this, lat, lon)
    class (EmissionRegion), intent(in) :: this
    real*8, intent(in) :: lat
    real*8, intent(in) :: lon
    
    hasLatLon = &
         & (lat >= this%southernEdge .and. lat <= this%northernEdge) .and. &
         & (lon >= this%westernEdge .and. lon < this%easternEdge) ! East should be <=, but preserving original behavior
  end function hasLatLon

end module EmissionRegion_mod
