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

submodule (request_base) request_base_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine request_group_cache_init (cache, comm)
    class(request_group_cache_t), intent(inout) :: cache
    type(MPI_COMM), intent(in) :: comm
    call MPI_COMM_DUP (comm, cache%parent_comm)
    !! Local operation.
    call MPI_COMM_GROUP (cache%parent_comm, cache%parent_group)
    cache%group = MPI_GROUP_EMPTY
    cache%comm = MPI_COMM_NULL
  end subroutine request_group_cache_init

  module subroutine request_group_cache_reset (cache)
    class(request_group_cache_t), intent(inout) :: cache
    cache%group = MPI_GROUP_EMPTY
    cache%comm = MPI_COMM_NULL
  end subroutine request_group_cache_reset

  module subroutine request_group_cache_update (cache, tag, rank)
    class(request_group_cache_t), intent(inout) :: cache
    integer, intent(in) :: tag
    integer, dimension(:), allocatable, intent(inout) :: rank
    type(MPI_GROUP) :: group
    integer :: result, error
    call move_alloc (rank, cache%rank)
    call MPI_GROUP_INCL (cache%parent_group, size (cache%rank), cache%rank, group)
    call MPI_GROUP_COMPARE (cache%group, group, result)
    if (result /= MPI_IDENT) then
       cache%group = group
       if (cache%comm /= MPI_COMM_NULL) call MPI_COMM_FREE (cache%comm)
       !! Group-local operation. However, time consuming.
       call MPI_COMM_CREATE_GROUP (cache%parent_comm, cache%group, tag, &
            cache%comm, error)
       if (error /= 0) then
          call msg_bug ("Error occured during communicator creation...")
       end if
    ! else
       ! call msg_message ("CACHE UPDATE: GROUPS ARE (NEARLY) IDENTICAL")
    end if
  end subroutine request_group_cache_update

  module subroutine request_group_cache_get_comm (cache, comm)
    class(request_group_cache_t), intent(in) :: cache
    type(MPI_COMM), intent(out) :: comm
    comm = cache%comm
  end subroutine request_group_cache_get_comm

  module function request_group_cache_is_master (cache) result (flag)
    class(request_group_cache_t), intent(in) :: cache
    integer :: rank, error
    logical :: flag
    call MPI_COMM_RANK (cache%comm, rank, error)
    if (error /= 0) then
       call msg_bug ("Error: Could not retrieve group rank.")
    end if
    flag = (rank == 0)
  end function request_group_cache_is_master

  !! =================================================
  !! request_base_t
  !! =================================================

  !> Initialize request base with parent communicator.
  !!
  !! In order to separate the communication between different parts of the request library,
  !! duplicate the parent communicator using MPI_COMM_DUP, also done by cache and handler objects.
  !!
  !! \param[in] comm Parent MPI communicator for overall library.
  module subroutine request_base_init (req, comm)
    class(request_base_t), intent(out) :: req
    type(MPI_COMM), intent(in) :: comm
    call MPI_COMM_DUP (comm, req%comm)
    call MPI_COMM_DUP (comm, req%external_comm)
    call req%cache%init (comm)
    call req%handler%init (comm)
  end subroutine request_base_init

  module subroutine request_base_write (req, unit)
    class(request_base_t), intent(in) :: req
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (allocated (req%balancer)) then
       call req%balancer%write (u)
    else
       write (u, "(A)") "[BALANCER]"
       write (u, "(A)") "=> Not allocated"
    end if
    call req%handler%write (u)
  end subroutine request_base_write

  !> Check whether current worker is master rank in object communicator.
  !!
  !! Do not confuse with a group's master !!!
  !! Proof: rank == 0
  module function request_base_is_master (req) result (flag)
    class(request_base_t), intent(in) :: req
    integer :: rank, ierr
    logical :: flag
    call MPI_COMM_RANK (req%comm, rank, ierr)
    if (ierr /= 0) then
       write (*, "(A,1X,I0)") "MPI Error: request_base_is_master", ierr
       stop 1
    end if
    flag = (rank == 0)
  end function request_base_is_master

  !> Provide external communicator.
  !!
  !! The external communicator is just a duplicate of request%comm,
  !! in order to provide the same group of workers to external communication,
  !! however, in a different context, such that communication from outside does not interfere with request.
  module subroutine request_base_get_external_comm (req, comm)
    class(request_base_t), intent(in) :: req
    type(MPI_COMM), intent(out) :: comm
    comm = req%external_comm
  end subroutine request_base_get_external_comm

  !> Add balancer to request.
  !!
  !! \param[inout] balancer
  module subroutine request_base_add_balancer (req, balancer)
    class(request_base_t), intent(inout) :: req
    class(balancer_base_t), allocatable, intent(inout) :: balancer
    if (allocated (req%balancer)) deallocate (req%balancer)
    call move_alloc (balancer, req%balancer)
  end subroutine request_base_add_balancer

  !> Add request handler with handler_id.
  !!
  !! \param[in] handler_id
  !! \param[in] handler Pointer to handler object.
  module subroutine request_base_add_handler (req, handler_id, handler)
    class(request_base_t), intent(inout) :: req
    integer, intent(in) :: handler_id
    class(request_handler_t), pointer, intent(in) :: handler
    call req%handler%add (handler_id, handler)
  end subroutine request_base_add_handler

  !> Reset request.
  !! Clear handler manager from associated callbacks,
  !! deallocate balancer, iff allocated, and reset communicator cache.
  module subroutine request_base_reset (req, deallocate_balancer)
    class(request_base_t), intent(inout) :: req
    logical, intent(in), optional :: deallocate_balancer
    logical :: flag
    flag = .false.; if (present (deallocate_balancer)) &
         flag = deallocate_balancer
    if (flag .and. allocated (req%balancer)) then
       deallocate (req%balancer)
    end if
    call req%handler%clear ()
    call req%cache%reset ()
  end subroutine request_base_reset

  !> Call handler for master communication for handler_id.
  !!
  !! \param[in] handler_id The associated key of the callback object.
  !! \param[in] source_rank The rank of the result's source.
  module subroutine request_base_call_handler &
       (req, handler_id, source_rank)
    class(request_base_t), intent(inout) :: req
    integer, intent(in) :: handler_id
    integer, intent(in) :: source_rank
    call req%handler%callback (handler_id, source_rank)
  end subroutine request_base_call_handler

  !> Call handler for slave communication for handler_id.
  !!
  !! \param[in] handler_id The associated key of the callback object.
  module subroutine request_base_call_client_handler (req, handler_id)
    class(request_base_t), intent(inout) :: req
    integer, intent(in) :: handler_id
    call req%handler%client_callback (handler_id, 0)
  end subroutine request_base_call_client_handler

  !> Wait on all handler in request handler manager to finish communication.
  module subroutine request_base_await_handler (req)
    class(request_base_t), intent(inout) :: req
    call req%handler%waitall ()
  end subroutine request_base_await_handler

  module subroutine request_base_barrier (req)
    class(request_base_t), intent(in) :: req
    integer :: error
    call MPI_BARRIER (req%comm, error)
    if (error /= MPI_SUCCESS) then
       call msg_fatal ("Request: Error occured during MPI_BARRIER synchronisation.")
    end if
  end subroutine request_base_barrier


end submodule request_base_s

