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

submodule (request_callback) request_callback_s

  use diagnostics

  implicit none

contains

  module subroutine request_handler_base_write (handler, unit)
    class(request_handler_t), intent(in) :: handler
    integer, intent(in), optional :: unit
    integer :: u, i
    u = ERROR_UNIT; if (present (unit)) u = unit
    write (u, "(A,1X,I0)") "N_REQUESTS", handler%n_requests
    write (u, "(A,1X,I0)") "TAG_OFFSET", handler%tag_offset
    write (u, "(A,1X,L1)") "FINISHED", handler%finished
    write (u, "(A,1X,L1)") "ACTIVATED", handler%activated
    write (u, "(A)") "I | SOURCE | TAG | ERROR | REQUEST_NULL"
    do i = 1, handler%n_requests
       write (u, "(A,4(1X,I0),1X,L1)") "REQUEST", i, &
            handler%status(i)%MPI_SOURCE, &
            handler%status(i)%MPI_TAG, &
            handler%status(i)%MPI_ERROR, &
            (handler%request(i) == MPI_REQUEST_NULL)
    end do
  end subroutine request_handler_base_write

  !! \param[inout] handler Handler must be intent inout, as the calling function may already manipulated the extended object.
  !! \param[in] n_requests Number of MPI requests the objects needs to be able to handle.
  !! \param[in] tag_offset First tag to be used, all other must follow in an increasing manner until tag_offset + (N_r + 1).
  !! Proof: tag \in {tag_offset, tag_offset + n_requests}.
  module subroutine request_handler_allocate (handler, n_requests, tag_offset)
    class(request_handler_t), intent(inout) :: handler
    integer, intent(in) :: n_requests
    integer, intent(in) :: tag_offset
    allocate (handler%request(n_requests), source = MPI_REQUEST_NULL)
    allocate (handler%status(n_requests))
    handler%n_requests = n_requests
    if (mod (tag_offset, n_requests) /= 0) &
         call msg_bug ("Error during handler allocate, tag_offset " // &
         "is not a multiple of n_requests.")
    !! What is the max.-allowed MPI_TAG?
    handler%tag_offset = tag_offset
    handler%activated = .false.
    handler%finished = .false.
  end subroutine request_handler_allocate

  module subroutine request_handler_get_status (handler)
    class(request_handler_t), intent(inout) :: handler
    integer :: i
    logical :: flag
    if (.not. handler%activated) return
    handler%finished = .true.
    do i = 1, handler%n_requests
       call MPI_REQUEST_GET_STATUS (handler%request(i), flag, &
            handler%status(i))
       handler%finished = handler%finished .and. flag
    end do
  end subroutine request_handler_get_status

  module subroutine request_handler_waitall (handler)
    class(request_handler_t), intent(inout) :: handler
    integer :: error
    if (.not. handler%activated .or. handler%finished) return
    call MPI_WAITALL (handler%n_requests, handler%request, &
         handler%status, error)
    if (error /= 0) then
       call msg_bug ("Request: Error occured during waitall on handler.")
    end if
    handler%finished = .true.
  end subroutine request_handler_waitall

  module function request_handler_testall (handler) result (flag)
    class(request_handler_t), intent(inout) :: handler
    logical :: flag
    integer :: error
    if (.not. handler%activated .or. .not. handler%finished) then
       call MPI_TESTALL (handler%n_requests, handler%request, &
            handler%finished, handler%status, error)
       ! call print_status ()
       if (error /= 0) then
          call msg_bug ("Request: Error occured during testall on handler.")
       end if
    end if
    flag = handler%finished
  contains
    subroutine print_status ()
      integer :: i
      do i = 1, handler%n_requests
         associate (status => handler%status(i))
           write (ERROR_UNIT, *) status%MPI_SOURCE, status%MPI_TAG, &
                status%MPI_ERROR
         end associate
      end do
    end subroutine print_status
  end function request_handler_testall

  module subroutine request_handler_free (handler)
    class(request_handler_t), intent(inout) :: handler
    integer :: i, error
    do i = 1, handler%n_requests
       if (handler%request(i) == MPI_REQUEST_NULL) cycle
       call MPI_REQUEST_FREE (handler%request(i), error)
       if (error /= 0) then
          call msg_bug ("Request: Error occured during free " // &
               "request on handler.")
       end if
    end do
  end subroutine request_handler_free

  module subroutine request_handler_manager_init (rhm, comm)
    class(request_handler_manager_t), intent(out) :: rhm
    type(MPI_COMM), intent(in) :: comm
    call MPI_COMM_DUP (comm, rhm%comm)
  end subroutine request_handler_manager_init

  module subroutine request_handler_manager_write (rhm, unit)
    class(request_handler_manager_t), intent(in) :: rhm
    integer, intent(in), optional :: unit
    integer :: u
    u = ERROR_UNIT; if (present (unit)) u = unit
    write (u, "(A)") "[REQUEST_CALLBACK_MANAGER]"
    call rhm%tree%write (u)
  end subroutine request_handler_manager_write

  module subroutine request_handler_manager_add (rhm, handler_id, handler)
    class(request_handler_manager_t), intent(inout) :: rhm
    integer, intent(in) :: handler_id
    class(request_handler_t), pointer, intent(in) :: handler
    class(*), pointer :: obj
    obj => handler
    call rhm%tree%insert (handler_id, obj)
  end subroutine request_handler_manager_add

  module subroutine request_handler_manager_clear (rhm)
    class(request_handler_manager_t), intent(inout) :: rhm
    call rhm%tree%clear ()
  end subroutine request_handler_manager_clear

  module subroutine request_handler_manager_fill_status (rhm)
    class(request_handler_manager_t), intent(inout) :: rhm
    type(binary_tree_iterator_t) :: iterator
    integer :: handler_id
    class(request_handler_t), pointer :: handler
    call iterator%init (rhm%tree)
    do while (iterator%is_iterable ())
       call iterator%next (handler_id)
       call rhm%handler_at (handler_id, handler)
       call handler%get_status ()
    end do
  end subroutine request_handler_manager_fill_status

  module function request_handler_manager_test &
       (rhm, handler_id) result (flag)
    class(request_handler_manager_t), intent(inout) :: rhm
    integer, intent(in) :: handler_id
    logical :: flag
    class(request_handler_t), pointer :: handler
    call rhm%handler_at (handler_id, handler)
    flag = handler%testall ()
  end function request_handler_manager_test

  module subroutine request_handler_manager_wait (rhm, handler_id)
    class(request_handler_manager_t), intent(inout) :: rhm
    integer, intent(in) :: handler_id
    class(request_handler_t), pointer :: handler
    call rhm%handler_at (handler_id, handler)
    call handler%waitall ()
  end subroutine request_handler_manager_wait

  module subroutine request_handler_manager_waitall (rhm)
    class(request_handler_manager_t), intent(inout) :: rhm
    type(binary_tree_iterator_t) :: iterator
    integer :: handler_id
    call iterator%init (rhm%tree)
    do while (iterator%is_iterable ())
       call iterator%next (handler_id)
       !! Test handler (destructive test on request handler).
       if (.not. rhm%test (handler_id)) &
            call rhm%wait (handler_id)
    end do
  end subroutine request_handler_manager_waitall

  module subroutine request_handler_manager_handler_at &
       (rhm, handler_id, handler)
    class(request_handler_manager_t), intent(in) :: rhm
    integer, intent(in) :: handler_id
    class(request_handler_t), pointer, intent(out) :: handler
    class(*), pointer :: obj
    call rhm%tree%search (handler_id, obj)
    select type (obj)
    class is (request_handler_t)
       handler => obj
    class default
       call msg_bug ("Object is not derived from request_handler_t.")
    end select
  end subroutine request_handler_manager_handler_at

  module function request_handler_manager_has_handler &
       (rhm, handler_id) result (flag)
    class(request_handler_manager_t), intent(inout) :: rhm
    integer, intent(in) :: handler_id
    logical :: flag
    flag = rhm%tree%has_key (handler_id)
  end function request_handler_manager_has_handler

  module subroutine request_handler_manager_callback &
       (rhm, handler_id, source_rank)
    class(request_handler_manager_t), intent(inout) :: rhm
    integer, intent(in) :: handler_id
    integer, intent(in) :: source_rank
    class(request_handler_t), pointer :: handler
    if (.not. rhm%tree%has_key (handler_id)) return
    call rhm%handler_at (handler_id, handler)
    call handler%handle (source_rank = source_rank, comm = rhm%comm)
  end subroutine request_handler_manager_callback

  module subroutine request_handler_manager_client_callback &
       (rhm, handler_id, dest_rank)
    class(request_handler_manager_t), intent(inout) :: rhm
    integer, intent(in) :: handler_id
    integer, intent(in) :: dest_rank
    class(request_handler_t), pointer :: handler
    if (.not. rhm%tree%has_key (handler_id)) return
    call rhm%handler_at (handler_id, handler)
    call handler%client_handle (dest_rank = dest_rank, comm = rhm%comm)
  end subroutine request_handler_manager_client_callback


end submodule request_callback_s

