module class_LatLonGrid

	implicit none

	type, public :: LatLonGrid
		integer :: ims		! # of longitude grid points
		integer :: jms		! # of latitude grid points

		real*8 :: dlat_dg	! Latitude grid spacing (degrees)
		real*8 :: dlon_dg	! Longitutde grid spacing (degrees)
	contains
		procedure :: init => LatLonGrid_init
		procedure :: ij2latlon => LatLonGrid_ij2latlon
		procedure :: i2lon => LatLonGrid_i2lon
		procedure :: j2lat => LatLonGrid_j2lat
	end type LatLonGrid

contains

	subroutine LatLonGrid_init(this, ims, jms)
		class(LatLonGrid), intent(inout) :: this
		integer, intent(in) :: ims
		integer, intent(in) :: jms
	! begin

		this%ims = ims
		this%jms = jms

        this%dlat_dg=180./real(jms) ! even spacing (default)
		! WARNING: HACK@
        IF (jms.eq.46) this%dlat_dg=180./real(jms-1) ! 1/2 box at pole for 4x5
        this%dlon_dg = 360.d0/real(ims)
	end subroutine LatLonGrid_init
	! ---------------------------------------------------------
	subroutine LatLonGrid_ij2latlon(this, i, j, lat, lon)
		class(LatLonGrid), intent(in) :: this
		integer, intent(in) :: i, j
		real*8, intent(out) :: lat, lon

		real*8 fjeq
	! begin
		lon = -180.+(i-0.5)*this%dlon_dg

        if (j .ne. 1 .and. j .ne. this%jms) then
           fjeq=0.5*(1+this%jms)
           lat = this%dlat_dg*(j - fjeq)
        elseif (j .eq. 1) then
           lat = -90.
        elseif (j .eq. this%jms) then
           lat = 90.
        endif
	end subroutine LatLonGrid_ij2latlon
	! ---------------------------------------------------------
	function LatLonGrid_i2lon(this, i)
		real*8 :: LatLonGrid_i2lon
		class(LatLonGrid), intent(in) :: this
		integer, intent(in) :: i
	! begin
		LatLonGrid_i2lon = -180.+(i-0.5)*this%dlon_dg
	end function LatLonGrid_i2lon
	! ---------------------------------------------------------
	function LatLonGrid_j2lat(this, j)
		real*8 :: LatLonGrid_j2lat
		class(LatLonGrid), intent(in) :: this
		integer, intent(in) :: j

		real*8 fjeq
	! begin
        if (j .ne. 1 .and. j .ne. this%jms) then
           fjeq=0.5*(1+this%jms)
           LatLonGrid_j2lat = this%dlat_dg*(j - fjeq)
        elseif (j .eq. 1) then
           LatLonGrid_j2lat = -90.
        elseif (j .eq. this%jms) then
           LatLonGrid_j2lat = 90.
        endif
	end function LatLonGrid_j2lat
	! ---------------------------------------------------------

end module class_LatLonGrid

! ==================================================================
! 
! PROGRAM test
! 
! use class_LatLonGrid
! 
! type(LatLonGrid) :: llg
! real*8 :: lat, lon
! 
! integer, parameter :: numi = 72
! integer, parameter :: numj = 46
! 
! call llg%init(numi, numj)
! 
! !do i=1,numi
! !	do j=1,numj
! !		call llg%ij2latlon(i,j, lat, lon)
! !		write(*,*) i,j,lat, lon
! !	end do
! !end do
! 
! do i=1,numi
! 	lon = llg%i2lon(i)
! 	write(*,*) '(i, lon) ==> ', i, lon
! end do
! 
! do j=1,numj
! 	lat = llg%j2lat(j)
! 	write(*,*) '(j, lat) ==> ', j, lat
! end do
! 
! 
! END program test
! 
