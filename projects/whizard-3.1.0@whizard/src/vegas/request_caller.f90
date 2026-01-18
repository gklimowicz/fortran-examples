! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'
module request_caller

  use kinds, only: default
  use diagnostics

  use request_base
  use balancer_base
  use balancer_channel
  use request_state
  use request_callback

  use mpi_f08 !NODEP!

  implicit none
  private

  public :: request_caller_t

  type, extends (request_base_t):: request_caller_t
     private
     integer :: n_channels = 0
     integer :: n_workers = 0
     type(request_state_t) :: state
   contains
     procedure :: init => request_caller_init
     procedure :: write => request_caller_write
     procedure :: has_workers => request_caller_has_workers
     procedure :: update_balancer => request_caller_update_balancer
     procedure :: handle_workload => request_caller_handle_workload
     procedure :: request_workload => request_caller_request_workload
     procedure :: release_workload => request_caller_release_workload
     procedure :: handle_and_release_workload => &
          request_caller_handle_and_release_workload
     procedure :: request_terminate => request_caller_request_terminate
  end type request_caller_t


  interface
    module subroutine request_caller_write (req, unit)
      class(request_caller_t), intent(in) :: req
      integer, intent(in), optional :: unit
    end subroutine request_caller_write
    module function request_caller_has_workers (req) result (flag)
      class(request_caller_t), intent(in) :: req
      logical :: flag
    end function request_caller_has_workers
    module subroutine request_caller_update_balancer (req, weight, parallel_grid)
      class(request_caller_t), intent(inout) :: req
      real(default), dimension(:), intent(in) :: weight
      logical, dimension(:), intent(in) :: parallel_grid
    end subroutine request_caller_update_balancer
    module subroutine request_caller_handle_workload (req)
      class(request_caller_t), intent(inout) :: req
    end subroutine request_caller_handle_workload
    module subroutine request_caller_request_workload (req, request)
      class(request_caller_t), intent(inout) :: req
      type(request_t), intent(out) :: request
    end subroutine request_caller_request_workload
    module subroutine request_caller_release_workload (req, request)
      class(request_caller_t), intent(inout) :: req
      type(request_t), intent(in) :: request
    end subroutine request_caller_release_workload
    module subroutine request_caller_handle_and_release_workload (req, request)
      class(request_caller_t), intent(inout) :: req
      type(request_t), intent(in) :: request
    end subroutine request_caller_handle_and_release_workload
    module subroutine request_caller_request_terminate (req)
      class(request_caller_t), intent(inout) :: req
    end subroutine request_caller_request_terminate
  end interface

contains

  subroutine request_caller_init (req, comm, n_channels)
    class(request_caller_t), intent(out) :: req
    type(MPI_COMM), intent(in) :: comm
    integer, intent(in) :: n_channels
    call req%base_init (comm)
    call MPI_COMM_SIZE (req%comm, req%n_workers)
    !! Exclude master rank (0) from set of workers.
    req%n_workers = req%n_workers - 1
    if (.not. req%has_workers ()) then
       call msg_warning ("Must not handle less than 3 ranks in a master/slave global queue.")
       call MPI_ABORT (req%comm, 1)
    end if
    req%n_channels = n_channels
    call req%state%init (comm, req%n_workers)
    if (req%is_master ()) then
       call allocate_balancer ()
    end if
  contains
    subroutine allocate_balancer ()
      class(balancer_base_t), allocatable :: balancer
      allocate (balancer_channel_t :: balancer)
      select type (balancer)
      type is (balancer_channel_t)
         call balancer%init (n_workers = req%n_workers, n_resources = req%n_channels)
      end select
      call req%add_balancer (balancer)
    end subroutine allocate_balancer
  end subroutine request_caller_init


end module request_caller
