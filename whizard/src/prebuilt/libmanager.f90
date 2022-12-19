! WHIZARD library manager
!
! This file handles the prebuilt process libraries
! (currently: empty)

subroutine dispatch_prclib_static (driver, basename, modellibs_ldflags)
  use iso_varying_string, string_t => varying_string
  use prclib_interfaces
  implicit none
  class(prclib_driver_t), intent(inout), allocatable :: driver
  type(string_t), intent(in) :: basename
  logical, intent(in), optional :: modellibs_ldflags
end subroutine dispatch_prclib_static

subroutine get_prclib_static (libname)
  use iso_varying_string, string_t => varying_string
  implicit none
  type(string_t), dimension(:), intent(inout), allocatable :: libname
  allocate (libname (0))
end subroutine get_prclib_static
