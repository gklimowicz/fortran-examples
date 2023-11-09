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

submodule (request_caller) request_caller_s

  use io_units

  implicit none

contains

  module subroutine request_caller_write (req, unit)
    class(request_caller_t), intent(in) :: req
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(A)") "[REQUEST_CALLER]"
    call req%base_write (u)
    if (req%is_master ()) &
       call req%state%write (u)
  end subroutine request_caller_write

  module function request_caller_has_workers (req) result (flag)
    class(request_caller_t), intent(in) :: req
    logical :: flag
    !! n_workers excludes the master rank.
    flag = (req%n_workers > 1)
  end function request_caller_has_workers

  module subroutine request_caller_update_balancer (req, weight, parallel_grid)
    class(request_caller_t), intent(inout) :: req
    real(default), dimension(:), intent(in) :: weight
    logical, dimension(:), intent(in) :: parallel_grid
    !! \note bug if not allocated?
    if (.not. allocated (req%balancer)) return
    call req%state%reset ()
    call req%reset ()
    select type (balancer => req%balancer)
    type is (balancer_channel_t)
       call balancer%update_state(weight, parallel_grid)
    end select
  end subroutine request_caller_update_balancer

  module subroutine request_caller_handle_workload (req)
    class(request_caller_t), intent(inout) :: req
    integer :: handler, tag, source, worker_id
    if (.not. allocated (req%balancer)) then
       call msg_warning ("Request: Error occured, load balancer not allocated.&
          & Terminate all workers.")
       !! We postpone to stop the program here so we can terminate all workers gracefully.
       !! First, we receive their requests, then we overwrite their "original" tag to MPI_TAG_TERMINATE.
       !! Second, we iterate this, until all workers are terminated and return without doing any besides.
    end if
    call req%state%init_request ()
    call req%state%receive_request ()
    do while (.not. req%state%is_terminated ())
       call req%state%await_request ()
       do while (req%state%has_request ())
          call req%state%get_request (source, tag, handler)
          !! Formally differentiate between worker_id and source.
          worker_id = source
          if (.not. allocated (req%balancer)) tag = MPI_TAG_TERMINATE
          select case (tag)
          case (MPI_TAG_REQUEST)
             if (req%balancer%is_assignable (worker_id)) then
                call req%balancer%assign_worker (worker_id, handler)
                if (.not. req%balancer%has_resource_group (handler)) then
                   call req%state%update_request (source, MPI_TAG_ASSIGN_SINGLE, handler)
                else
                   call req%state%update_request (source, MPI_TAG_ASSIGN_GROUP, handler)
                   call provide_request_group (handler, source)
                end if
             else
                call req%state%terminate (source)
             end if
          case (MPI_TAG_HANDLER_AND_RELEASE)
             call req%call_handler (handler, source_rank = source)
             call req%balancer%free_worker (worker_id, handler)
          case (MPI_TAG_RELEASE)
             call req%balancer%free_worker (worker_id, handler)
          case (MPI_TAG_TERMINATE)
             call req%state%terminate (source)
          case (MPI_TAG_CLIENT_TERMINATE)
             !! Allow workers to request their own termination.
             call req%state%terminate (source)
          case default
             call msg_warning ()
          end select
       end do
       call req%state%receive_request ()
    end do
    !! If we are here, there should be no leftover communnication.
    !! Hence, we must check whether there is no left-over communication call (from server-side).
    call req%state%free_request ()
  contains
    subroutine provide_request_group (handler_id, dest_rank)
      integer, intent(in) :: handler_id
      integer, intent(in) :: dest_rank
      integer, dimension(:), allocatable :: rank
      !! Rank indices and worker indices are identical, as we skip the master worker deliberately,
      !! we can reuse the worker indices as rank indices.
      call req%balancer%get_resource_group (handler_id, rank)
      call req%state%provide_request_group (dest_rank, rank)
    end subroutine provide_request_group
  end subroutine request_caller_handle_workload

  module subroutine request_caller_request_workload (req, request)
    class(request_caller_t), intent(inout) :: req
    type(request_t), intent(out) :: request
    type(MPI_STATUS) :: status
    call req%state%client_serve (request%handler_id, status)
    request%terminate = .false.
    request%group = .false.
    request%callback = .false.
    request%comm = MPI_COMM_NULL
    select case (status%MPI_TAG)
    case (MPI_TAG_ASSIGN_SINGLE)
       !! Default to req's communicator.
       request%comm = req%external_comm
       request%group_master = .true.
       request%callback = .true.
    case (MPI_TAG_ASSIGN_GROUP)
       request%group = .true.
       call retrieve_request_group (request%handler_id)
       call req%cache%get_comm (request%comm)
       request%group_master = req%cache%is_master ()
       request%callback = req%cache%is_master ()
    case (MPI_TAG_TERMINATE)
       request%terminate = status%MPI_TAG == MPI_TAG_TERMINATE
    end select
  contains
    subroutine retrieve_request_group (handler_id)
      integer, intent(in) :: handler_id
      integer, dimension(:), allocatable :: rank
      !! Here, worker and rank indices are interchangeable.
      call req%state%retrieve_request_group (rank)
      call req%cache%update (handler_id, rank)
    end subroutine retrieve_request_group
  end subroutine request_caller_request_workload

  module subroutine request_caller_release_workload (req, request)
    class(request_caller_t), intent(inout) :: req
    type(request_t), intent(in) :: request
    call req%state%client_free (request%handler_id, &
         has_callback = request%group_master)
  end subroutine request_caller_release_workload

  module subroutine request_caller_handle_and_release_workload (req, request)
    class(request_caller_t), intent(inout) :: req
    type(request_t), intent(in) :: request
    if (.not. req%handler%has_handler (request%handler_id)) then
       call msg_bug ("Request: Handler is not registered for this worker.")
    end if
    call req%release_workload (request)
    call req%call_client_handler (request%handler_id)
  end subroutine request_caller_handle_and_release_workload

  module subroutine request_caller_request_terminate (req)
    class(request_caller_t), intent(inout) :: req
    if (req%is_master ()) return
    call req%state%client_terminate ()
  end subroutine request_caller_request_terminate


end submodule request_caller_s

