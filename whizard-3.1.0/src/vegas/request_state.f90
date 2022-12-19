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
module request_state

  use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
  use iterator

  use mpi_f08 !NODEP!

  implicit none
  private

  public :: request_state_t

  integer, parameter, public :: MPI_EMPTY_HANDLER = 0

  integer, parameter, public :: MPI_TAG_NULL = 0, &
       MPI_TAG_REQUEST = 1, &
       MPI_TAG_RELEASE = 2, &
       MPI_TAG_HANDLER_AND_RELEASE = 4, &
       MPI_TAG_TERMINATE = 8, &
       MPI_TAG_CLIENT_TERMINATE = 16, &
       MPI_TAG_ASSIGN_SINGLE = 32, &
       MPI_TAG_ASSIGN_GROUP = 64, &
       MPI_TAG_COMMUNICATOR_GROUP = 128

  integer, parameter :: MPI_STATE_ERR = 1


  type :: request_state_t
     private
     type(MPI_COMM) :: comm
     integer :: n_workers = 0
     integer :: n_workers_done = 0
     !! From MPI-3.1 book
     !! i \in {1, N_workes_done}, max size = N_workers
     type(MPI_Request), dimension(:), allocatable :: request
     type(MPI_Status), dimension(:), allocatable :: status
     integer, dimension(:), allocatable :: indices
     !! i \in {1, N_workers}
     integer, dimension(:), allocatable :: handler
     logical, dimension(:), allocatable :: terminated
     type(iterator_t) :: request_iterator
   contains
     procedure :: init => request_state_init
     procedure :: write => request_state_write
     procedure :: reset => request_state_reset
     procedure :: is_terminated => request_state_is_terminated
     procedure, private :: set_terminated => request_state_set_terminated
     procedure :: terminate => request_state_terminate
     procedure :: client_terminate => request_state_client_terminate
     procedure :: init_request => request_state_init_request
     procedure :: receive_request => request_state_receive_request
     procedure :: await_request => request_state_await_request
     procedure :: has_request => request_state_has_request
     procedure :: get_request => request_state_get_request
     procedure :: update_request => request_state_update_request
     procedure :: free_request => request_state_free_request
     procedure :: provide_request_group => &
          request_state_provide_request_group
     procedure :: retrieve_request_group => &
          request_state_retrieve_request_group
     procedure :: client_serve => request_state_client_serve
     procedure :: client_free => request_state_client_free
  end type request_state_t


  interface
    module subroutine request_state_init (state, comm, n_workers)
      class(request_state_t), intent(out) :: state
      type(MPI_COMM), intent(in) :: comm
      integer, intent(in) :: n_workers
    end subroutine request_state_init
    module subroutine request_state_write (state, unit)
      class(request_state_t), intent(in) :: state
      integer, intent(in), optional :: unit
    end subroutine request_state_write
    module subroutine request_state_reset (state)
      class(request_state_t), intent(inout) :: state
    end subroutine request_state_reset
    module function request_state_is_terminated (state) result (flag)
      class(request_state_t), intent(in) :: state
      logical :: flag
    end function request_state_is_terminated
    module subroutine request_state_set_terminated (state, rank)
      class(request_state_t), intent(inout) :: state
      integer, intent(in) :: rank
    end subroutine request_state_set_terminated
    module subroutine request_state_terminate (state, rank)
      class(request_state_t), intent(inout) :: state
      integer, intent(in) :: rank
    end subroutine request_state_terminate
    module subroutine request_state_client_terminate (state)
      class(request_state_t), intent(in) :: state
    end subroutine request_state_client_terminate
    module subroutine request_state_init_request (state)
      class(request_state_t), intent(inout) :: state
    end subroutine request_state_init_request
    module subroutine request_state_receive_request (state)
      class(request_state_t), intent(inout) :: state
    end subroutine request_state_receive_request
    module subroutine request_state_await_request (state)
      class(request_state_t), intent(inout) :: state
    end subroutine request_state_await_request
    pure module function request_state_has_request (state) result (flag)
      class(request_state_t), intent(in) :: state
      logical :: flag
    end function request_state_has_request
    module subroutine request_state_get_request (state, rank, tag, handler)
      class(request_state_t), intent(inout) :: state
      integer, intent(out) :: rank
      integer, intent(out) :: tag
      integer, intent(out) :: handler
    end subroutine request_state_get_request
    module subroutine request_state_update_request (state, rank, tag, handler)
      class(request_state_t), intent(inout) :: state
      integer, intent(in) :: rank
      integer, intent(in) :: tag
      integer, intent(in) :: handler
    end subroutine request_state_update_request
    module subroutine request_state_free_request (state)
      class(request_state_t), intent(inout) :: state
    end subroutine request_state_free_request
    module subroutine request_state_provide_request_group &
         (state, dest_rank, worker)
      class(request_state_t), intent(in) :: state
      integer, intent(in) :: dest_rank
      integer, dimension(:), intent(in) :: worker
    end subroutine request_state_provide_request_group
    module subroutine request_state_retrieve_request_group (state, worker)
      class(request_state_t), intent(inout) :: state
      integer, dimension(:), allocatable, intent(out) :: worker
    end subroutine request_state_retrieve_request_group
    module subroutine request_state_client_serve (state, handler_id, status)
      class(request_state_t), intent(in) :: state
      integer, intent(out) :: handler_id
      type(MPI_STATUS), intent(out) :: status
    end subroutine request_state_client_serve
    module subroutine request_state_client_free (state, handler_id, has_callback)
      class(request_state_t), intent(in) :: state
      integer, intent(in) :: handler_id
      logical, intent(in) :: has_callback
    end subroutine request_state_client_free
  end interface

end module request_state
