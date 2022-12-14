module Domain_mod
  implicit none
  private

  public :: Domain_type
  public :: newDomain

  type Domain1D_type
    private
    integer :: interior(2) ! exclude halo cells
    integer :: exterior(2) ! include halo cells
  end type Domain1D_type

  type Domain_type
    private
    type (Domain1D_type) :: axes(2)
    integer :: tile  = 1 ! support for cubed-sphere
    integer :: globalExtents(2) 
  end type Domain_type

contains

  function newDomain1d(interior, exterior) result(domain1d)
    type (Domain1d_type) :: domain1d
    integer :: interior(2) ! exclude halo cells
    integer :: exterior(2) ! include halo cells
    
    domain1d%interior = interior
    domain1d%exterior = exterior

  end function newDomain1d

  function newDomain(im, jm) result(domain)
    type (Domain_type) :: domain
    integer, intent(in) :: im
    integer, intent(in) :: jm

    domain%globalExtents = (/im, jm/)

  end function newDomain

  subroutine setTile(this, tile)
    type (Domain_type), intent(inout) :: this
    integer, intent(in) :: tile
    this%tile = tile
  end subroutine setTile

  subroutine setAxis(this, ith, axis)
    type (Domain_type), intent(inout) :: this
    integer, intent(in) :: ith
    type (Domain1D_type), intent(in) :: axis

    this%axes(ith) = axis

  end subroutine setAxis

end module Domain_mod
