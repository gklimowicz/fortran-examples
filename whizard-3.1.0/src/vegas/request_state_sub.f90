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

submodule (request_state) request_state_s

  use diagnostics

  implicit none

contains

  module subroutine request_state_init (state, comm, n_workers)
    class(request_state_t), intent(out) :: state
    type(MPI_COMM), intent(in) :: comm
    integer, intent(in) :: n_workers
    integer :: rank
    call MPI_COMM_DUP (comm, state%comm)
    state%n_workers = n_workers
    state%n_workers_done = n_workers
    call state%request_iterator%init (1, n_workers)
    allocate (state%request(state%n_workers), source = MPI_REQUEST_NULL)
    allocate (state%status(state%n_workers), source = MPI_STATUS_IGNORE)
    allocate (state%handler(state%n_workers), source = MPI_EMPTY_HANDLER)
    allocate (state%indices(state%n_workers), source = 0)
    allocate (state%terminated(state%n_workers), source = .false.)
    state%indices = [(rank, rank = 1, n_workers)]
  end subroutine request_state_init

  module subroutine request_state_write (state, unit)
    class(request_state_t), intent(in) :: state
    integer, intent(in), optional :: unit
    integer :: u, i
    u = ERROR_UNIT; if (present (unit)) u = unit
    write (u, "(A)") "[REQUEST_STATE]"
    write (u, "(A,1X,I0)") "N_WORKERS", state%n_workers
    write (u, "(A,1X,I0)") "N_WORKERS_DONE", state%n_workers_done
    write (u, "(A)") "RANK | SOURCE | TAG | ERROR | REQUEST_NULL"
    do i = 1, state%n_workers_done
       write (u, "(A,4(1X,I0),1X,L1)") "REQUEST", state%indices(i), &
            state%status(i)%MPI_SOURCE, &
            state%status(i)%MPI_TAG, &
            state%status(i)%MPI_ERROR, &
            (state%request(i) == MPI_REQUEST_NULL)
    end do
    write (u, "(A,999(1X,I0))") "HANDLER", state%handler
    write (u, "(A,999(1X,L1))") "TERMINATED", state%terminated
  end subroutine request_state_write

  module subroutine request_state_reset (state)
    class(request_state_t), intent(inout) :: state
    integer :: rank
    state%n_workers_done = state%n_workers
    call state%request_iterator%init (1, state%n_workers)
    state%handler = MPI_EMPTY_HANDLER
    state%indices = [(rank, rank = 1, state%n_workers)]
    state%terminated = .false.
  end subroutine request_state_reset

  ! pure module function request_state_is_terminated (state) result (flag)
  module function request_state_is_terminated (state) result (flag)
    class(request_state_t), intent(in) :: state
    logical :: flag
    flag = all (state%terminated)
  end function request_state_is_terminated

  module subroutine request_state_set_terminated (state, rank)
    class(request_state_t), intent(inout) :: state
    integer, intent(in) :: rank
    state%terminated(rank) = .true.
  end subroutine request_state_set_terminated

  module subroutine request_state_terminate (state, rank)
    class(request_state_t), intent(inout) :: state
    integer, intent(in) :: rank
    integer :: error
    call MPI_SEND (MPI_EMPTY_HANDLER, 1, MPI_INTEGER, &
         rank, MPI_TAG_TERMINATE, state%comm, error)
    if (error /= 0) then
       write (msg_buffer, "(A,1X,I3)") "Request: Error occured " // &
            "during terminate, RANK", rank
       call msg_bug ()
    end if
    call state%set_terminated (rank)
  end subroutine request_state_terminate

  module subroutine request_state_client_terminate (state)
    class(request_state_t), intent(in) :: state
    integer :: error
    call MPI_SEND (MPI_EMPTY_HANDLER, 1, MPI_INTEGER, &
         0, MPI_TAG_CLIENT_TERMINATE, state%comm, error)
    if (error /= 0) then
       write (msg_buffer, "(A,1X,I3)") "Request: Error occured " // &
            "during client-sided terminate"
       call msg_bug ()
    end if
  end subroutine request_state_client_terminate

  module subroutine request_state_init_request (state)
    class(request_state_t), intent(inout) :: state
    integer :: i, rank, error
    do i = 1, state%n_workers_done
       rank = state%indices(i)
       call MPI_RECV_INIT (state%handler(rank), 1, MPI_INTEGER, &
            rank, MPI_ANY_TAG, state%comm, state%request(rank), error)
       if (error /= 0) then
          write (msg_buffer, "(A,2(A,1X,I0))") "Request: Error occured during receive init, &
               & RANK", rank, "HANDLER", state%handler(rank)
          call msg_message ()
          call MPI_ABORT (state%comm, MPI_STATE_ERR)
       end if
    end do
  end subroutine request_state_init_request

  module subroutine request_state_receive_request (state)
    class(request_state_t), intent(inout) :: state
    integer :: i, rank
    integer :: error
    if (state%is_terminated ()) return
    call sanitize_from_terminated_ranks ()
    !! Receive new requests from (still active) workers.
    do i = 1, state%n_workers_done
       rank = state%indices(i)
       call MPI_START (state%request(rank), error)
       if (error /= 0) then
          write (msg_buffer, "(A,2(A,1X,I6))") "Request: Error occured during receive request, &
             & RANK", rank, "HANDLER", state%handler(rank)
          call msg_message ()
          call MPI_ABORT (state%comm, MPI_STATE_ERR)
       end if
    end do
  contains
    subroutine sanitize_from_terminated_ranks ()
      integer :: n_workers_done
      integer, dimension(:), allocatable :: indices
      !! Remove terminated ranks from done workers.
      indices = pack(state%indices(:state%n_workers_done), &
           .not. state%terminated(state%indices(:state%n_workers_done)))
      state%n_workers_done = size (indices)
      state%indices(:state%n_workers_done) = indices
    end subroutine sanitize_from_terminated_ranks
  end subroutine request_state_receive_request

  module subroutine request_state_await_request (state)
    class(request_state_t), intent(inout) :: state
    integer :: error
    if (state%is_terminated ()) return
    !! We verify that we have active handles associated with request state.
    call MPI_TESTSOME (state%n_workers, state%request, state%n_workers_done, &
         state%indices, state%status, error)
    if (error /= 0) then
       write (ERROR_UNIT, "(A)") "Error occured during await (testing) request..."
       call state%write (ERROR_UNIT)
       call MPI_ABORT (state%comm, MPI_STATE_ERR)
    else if (state%n_workers_done == MPI_UNDEFINED) then
       write (ERROR_UNIT, "(A)") "TEST_WAITSOME returned with MPI_UNDEFINED."
       call state%write (ERROR_UNIT)
       call MPI_ABORT (state%comm, MPI_STATE_ERR)
    end if
    !! Wait a little bit...
    if (state%n_workers_done == 0) then
       !! Proof: REQUEST(i), i \in {1, N_workers}, i is equivalent to rank.
       !! Proof: INDICES(j), STATUS(j), j \in {1, N_workers_done}
       !! Proof: INDICES(j) -> i, injectiv.
       call MPI_WAITSOME (state%n_workers, state%request, state%n_workers_done, &
            state%indices, state%status, error)
       if (error /= 0) then
          write (ERROR_UNIT, "(A)") "Error occured during await request..."
          call state%write (ERROR_UNIT)
          call MPI_ABORT (state%comm, MPI_STATE_ERR)
       end if
    endif
    call state%request_iterator%init (1, state%n_workers_done)
  end subroutine request_state_await_request

  pure module function request_state_has_request (state) result (flag)
    class(request_state_t), intent(in) :: state
    logical :: flag
    flag = state%request_iterator%is_iterable ()
  end function request_state_has_request

  module subroutine request_state_get_request (state, rank, tag, handler)
    class(request_state_t), intent(inout) :: state
    integer, intent(out) :: rank
    integer, intent(out) :: tag
    integer, intent(out) :: handler
    integer :: ndx
    if (.not. state%has_request ()) then
       call msg_bug ("Request: Cannot access missing request.")
    end if
    ndx = state%request_iterator%next ()
    rank = state%indices(ndx)
    if (rank /= state%status(ndx)%MPI_SOURCE) then
       write (msg_buffer, "(A,2(1X,I3))") &
          "Request: RANK and SOURCE mismatch", rank, &
          state%status(ndx)%MPI_SOURCE
       call msg_bug ()
    end if
    tag = state%status(ndx)%MPI_TAG
    handler = state%handler(rank)
  end subroutine request_state_get_request

  module subroutine request_state_update_request (state, rank, tag, handler)
    class(request_state_t), intent(inout) :: state
    integer, intent(in) :: rank
    integer, intent(in) :: tag
    integer, intent(in) :: handler
    integer :: error
    state%handler(rank) = handler
    call MPI_SEND (handler, 1, MPI_INTEGER, &
         rank, tag, state%comm, error)
    if (error /= 0) then
       write (msg_buffer, "(A,3(A,1X,I3))") "Request: Error occured during update, &
          &RANK", rank, "TAG", tag, "HANDLER", handler
       call msg_bug ()
    end if
  end subroutine request_state_update_request

  module subroutine request_state_free_request (state)
    class(request_state_t), intent(inout) :: state
    integer :: rank
    do rank = 1, state%n_workers
       if (state%request(rank) == MPI_REQUEST_NULL) cycle
       call MPI_REQUEST_FREE (state%request(rank))
    end do
  end subroutine request_state_free_request

  module subroutine request_state_provide_request_group &
       (state, dest_rank, worker)
    class(request_state_t), intent(in) :: state
    integer, intent(in) :: dest_rank
    integer, dimension(:), intent(in) :: worker
    call MPI_SEND (worker, size (worker), MPI_INTEGER, &
         dest_rank, MPI_TAG_COMMUNICATOR_GROUP, state%comm)
  end subroutine request_state_provide_request_group

  module subroutine request_state_retrieve_request_group (state, worker)
    class(request_state_t), intent(inout) :: state
    integer, dimension(:), allocatable, intent(out) :: worker
    type(MPI_STATUS) :: status
    integer :: n_workers
    call MPI_PROBE (0, MPI_TAG_COMMUNICATOR_GROUP, state%comm, status)
    call MPI_GET_COUNT(status, MPI_INTEGER, n_workers)
    allocate (worker (n_workers), source = 0)
    call MPI_RECV (worker, n_workers, MPI_INTEGER, &
         0, MPI_TAG_COMMUNICATOR_GROUP, state%comm, status)
  end subroutine request_state_retrieve_request_group

  module subroutine request_state_client_serve (state, handler_id, status)
    class(request_state_t), intent(in) :: state
    integer, intent(out) :: handler_id
    type(MPI_STATUS), intent(out) :: status
    call MPI_SEND (MPI_EMPTY_HANDLER, 1, MPI_INTEGER, &
         0, MPI_TAG_REQUEST, state%comm)
    call MPI_RECV (handler_id, 1, MPI_INTEGER, &
         0, MPI_ANY_TAG, state%comm, status)
  end subroutine request_state_client_serve

  module subroutine request_state_client_free (state, handler_id, has_callback)
    class(request_state_t), intent(in) :: state
    integer, intent(in) :: handler_id
    logical, intent(in) :: has_callback
    integer :: tag
    tag = merge (MPI_TAG_HANDLER_AND_RELEASE, MPI_TAG_RELEASE, has_callback)
    call MPI_SEND (handler_id, 1, MPI_INTEGER, &
         0, tag, state%comm)
  end subroutine request_state_client_free


end submodule request_state_s

