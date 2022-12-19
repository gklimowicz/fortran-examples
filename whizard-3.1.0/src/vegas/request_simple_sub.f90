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

submodule (request_simple) request_simple_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine request_simple_update (req, parallel_grid)
    class(request_simple_t), intent(inout) :: req
    logical, dimension(:), intent(in) :: parallel_grid
    integer :: me, worker
    call req%reset ()
    call MPI_COMM_RANK (req%comm, me)
    worker = SHIFT_RANK_TO_WORKER(me)
    select type (balancer => req%balancer)
    type is (balancer_simple_t)
       call balancer%update_state (worker, parallel_grid)
    end select
    req%parallel_grid = parallel_grid
  end subroutine request_simple_update

  module subroutine request_simple_write (req, unit)
    class(request_simple_t), intent(in) :: req
    integer, intent(in), optional :: unit
    integer :: u, n_size
    u = given_output_unit (unit)
    write (u, "(A)") "[REQUEST_SIMPLE]"
    write (u, "(A,1X,I0)") "N_CHANNELS", req%n_channels
    write (u, "(A,1X,I0)") "N_WORKERS", req%n_workers
    n_size = min (25, req%n_channels)
    write (u, "(A,25(1X,L1))") "PARALLEL_GRID", req%parallel_grid(:n_size)
    call req%base_write (u)
  end subroutine request_simple_write

  module function request_simple_has_workers (req) result (flag)
    class(request_simple_t), intent(in) :: req
    logical :: flag
    flag = (req%n_workers > 1)
  end function request_simple_has_workers

  module function request_simple_get_request_master &
       (req, channel) result (rank)
    class(request_simple_t), intent(in) :: req
    integer, intent(in) :: channel
    integer :: rank
    if (.not. allocated (req%balancer)) then
       call msg_bug ("Error: Balancer is not allocated.")
    end if
    rank = shift_worker_to_rank (req%balancer%get_resource_master (channel))
    !! "Caveat emptor" hits here:
    !! The balancer returns either a valid worker id or (-1) depending on the associated resource (it must be active...)
    !! We have to check whether returned worker index is plausible.
  end function request_simple_get_request_master

  !> Request workload.
  !!
  !! Depending on parallel_grid, we fill the request object differently.
  !! First, we do not set commnuicator for .not. parallel_grid (group and group master are set to .false., also).
  !! And the callback needs to be executed.
  !! Second, for parallel_grid, we set req%comm to the associated communicator and set group to .true..
  !! However, the overall master has the grid's result, therefore, only the master needs to the callback.
  !! Remark: We can actually intercept the callback for the master to himself; the results are already in the current position.
  module subroutine request_simple_request_workload (req, request)
    class(request_simple_t), intent(inout) :: req
    type(request_t), intent(out) :: request
    integer :: rank, worker_id
    call MPI_COMM_RANK (req%comm, rank)
    worker_id = shift_rank_to_worker (rank)
    if (.not. req%balancer%is_pending () &
         .or. .not. req%balancer%is_assignable (worker_id)) then
       request%terminate = .true.
       return
    end if
    call req%balancer%assign_worker (worker_id, request%handler_id)
    associate (channel => request%handler_id)
      if (req%parallel_grid (channel)) then
         request%comm = req%external_comm
         request%group = .true.
         !! The object communicator is master.
         request%group_master = &
              (req%get_request_master (channel) == rank)
         request%callback = req%is_master ()
      else
         request%comm = req%external_comm
         request%group = .false.
         request%group_master = .true.
         request%callback = .true.
      end if
    end associate
  end subroutine request_simple_request_workload

  module subroutine request_simple_release_workload (req, request)
    class(request_simple_t), intent(inout) :: req
    type(request_t), intent(in) :: request
    integer :: rank, worker_id
    call MPI_COMM_RANK (req%comm, rank)
    worker_id = shift_rank_to_worker (rank)
    call req%balancer%free_worker (worker_id, request%handler_id)
  end subroutine request_simple_release_workload

  module subroutine request_simple_handle_and_release_workload (req, request)
    class(request_simple_t), intent(inout) :: req
    type(request_t), intent(in) :: request
    call req%call_client_handler (request%handler_id)
    call req%release_workload (request)
  end subroutine request_simple_handle_and_release_workload


end submodule request_simple_s

