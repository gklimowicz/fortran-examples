! --------------------------------------------------------------------------
! This file enables overloading for type/kind/rank for generic I/O
! actions.   CPP is used to include this file multiple times
! with different values for the _DECLARE_ macro and the _SUFFIX_ token.
! --------------------------------------------------------------------------
   subroutine _PROC_NAME_(handle, action, arr, varinfo, &
        & unusable, &
        & r4_on_disk, & ! defvar params
        & jdim, & ! extra dist params
        & no_xdim, & ! extra write params
        & record, record1 & ! extra read params
        )
      use iso_fortran_env, only: real64
      class (ParallelIo), intent(in) :: handle
      character(len=*), intent(in) :: action
      _DECLARE_(arr)
      character(len=*), intent(in) :: varinfo
      ! 
      type (UnusableArgument), optional :: unusable
      ! remaining arguments must use keyword
      logical, optional, intent(in) :: r4_on_disk
      integer, optional, intent(in) :: jdim
      logical, optional, intent(in) :: no_xdim
      integer, optional, intent(in) :: record, record1

      ! local variables
      character(len=:), allocatable :: variableName

      variableName = getVariableNameFrom(varInfo)

      select case (action)
      case ('define')
         if (allocated(handle%compareString)) then
            call defvar(handle%grid, handle%fileId, arr, varinfo, &
                 & r4_on_disk, defby=handle%compareString)
         else
            call defvar(handle%grid, handle%fileId, arr, varinfo, &
                 & r4_on_disk)
         end if
      case ('read_dist')
         call read_dist_data(handle%grid, handle%fileId, variableName, arr, &
              & jdim=jdim, no_xdim=no_xdim, record=record, record1=record1)
      case ('write_dist')
         call write_dist_data(handle%grid, handle%fileId, variableName, &
              & arr, jdim=jdim, no_xdim=no_xdim)
      case ('default')
         call stop_model('NetCdFHandle::do() unsupported option for iAction.')
      end select

   end subroutine _PROC_NAME_
