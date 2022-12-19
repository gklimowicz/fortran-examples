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
module request_base

  use balancer_base
  use request_callback
  use mpi_f08 !NODEP!

  implicit none
  private

  public :: request_t
  public :: request_base_t

  type :: request_t
     integer :: handler_id = 0
     logical :: terminate = .false.
     logical :: group = .false.
     logical :: group_master = .false.
     logical :: callback = .false.
     type(MPI_COMM) :: comm
  end type request_t

  type :: request_group_cache_t
     private
     type(MPI_COMM) :: parent_comm
     type(MPI_GROUP) :: parent_group
     type(MPI_COMM) :: comm
     type(MPI_GROUP) :: group
     integer, dimension(:), allocatable :: rank
   contains
     procedure :: init => request_group_cache_init
     procedure :: reset => request_group_cache_reset
     procedure :: update => request_group_cache_update
     procedure :: get_comm => request_group_cache_get_comm
     procedure :: is_master => request_group_cache_is_master
  end type request_group_cache_t

  type, abstract :: request_base_t
     type(MPI_COMM) :: comm
     type(MPI_COMM) :: external_comm !! communicator for use outside of request, however, just a duplicate of comm.
     class(balancer_base_t), allocatable :: balancer
     type(request_group_cache_t) :: cache
     type(request_handler_manager_t) :: handler
   contains
     procedure :: base_init => request_base_init
     procedure :: base_write => request_base_write
     procedure(request_base_deferred_write), deferred :: write
     procedure :: is_master => request_base_is_master
     procedure(request_base_has_workers), deferred :: has_workers
     procedure :: get_external_comm => request_base_get_external_comm
     procedure :: add_balancer => request_base_add_balancer
     procedure :: add_handler => request_base_add_handler
     procedure :: reset => request_base_reset
     procedure :: call_handler => request_base_call_handler
     procedure :: call_client_handler => request_base_call_client_handler
     procedure :: await_handler => request_base_await_handler
     procedure :: barrier => request_base_barrier
     procedure(request_base_request_workload), deferred :: request_workload
     procedure(request_base_release_workload), deferred :: release_workload
     procedure(request_base_handle_and_release_workload), deferred :: handle_and_release_workload
  end type request_base_t


  !> The basic idea behind the request mechanism is that each associated worker can request a workload either from a predefined local stack, from local stealing or from a global queue.
  !!
  !! We note that it is not necessary to differentiate between master and worker on this level of abstraction.
  !! Hence, the request interface ignores any notion regarding a possible parallelization concept.
  abstract interface
     subroutine request_base_deferred_write (req, unit)
       import :: request_base_t
       class(request_base_t), intent(in) :: req
       integer, intent(in), optional :: unit
     end subroutine request_base_deferred_write

     !> Verify if request object has workers.
     !!
     !! An implementation shall return if there at least two workers, or otherwisely stated,
     !! one master and one slave at least, when both are used as computing ranks.
     logical function request_base_has_workers (req) result (flag)
       import :: request_base_t
       class(request_base_t), intent(in) :: req
     end function request_base_has_workers
  end interface

     !> Request workload and returns an request_t object.
     !!
     !! The request_t object has an associated handler_id and provide several ways
     !! to indicate whether the execution is to be terminated, or the request has an associated communictor.
     !! Finally, whether we expect that the handler id will be connected to an callback.
     !!
     !! \param[out] request Request container.
  abstract interface
     subroutine request_base_request_workload (req, request)
       import :: request_base_t, request_t
       class(request_base_t), intent(inout) :: req
       type(request_t), intent(out) :: request
     end subroutine request_base_request_workload
  end interface

     !> Release workload with the information from the request container.
     !!
     !! The release procedure may notify the master about the finishing of the workload associated with the handler_id.
     !! Or, it may just bookkeep whether the workload has finished.
     !! Additionally, if request%callback was true, it could handle the callback (from client side.)
  abstract interface
     subroutine request_base_release_workload (req, request)
       import :: request_base_t, request_t
       class(request_base_t), intent(inout) :: req
       type(request_t), intent(in) :: request
     end subroutine request_base_release_workload
  end interface

     !> Handle associated callback and release workload with the information from the request container.
     !!
     !! The procedure must call the associated callback handler using the handler_id.
     !! Remark: The callback manager is quite squishy regarding a missing handler (silent failure).
     !!         The procedure has to take care whether the callback was actually successful.
     !! The further release of the workload can then be deferred to the release_workload procedure.
     !! \param[in] request.
  abstract interface
     subroutine request_base_handle_and_release_workload (req, request)
       import :: request_base_t, request_t
       class(request_base_t), intent(inout) :: req
       type(request_t), intent(in) :: request
     end subroutine request_base_handle_and_release_workload
  end interface


  interface
    module subroutine request_group_cache_init (cache, comm)
      class(request_group_cache_t), intent(inout) :: cache
      type(MPI_COMM), intent(in) :: comm
    end subroutine request_group_cache_init
    module subroutine request_group_cache_reset (cache)
      class(request_group_cache_t), intent(inout) :: cache
    end subroutine request_group_cache_reset
    module subroutine request_group_cache_update (cache, tag, rank)
      class(request_group_cache_t), intent(inout) :: cache
      integer, intent(in) :: tag
      integer, dimension(:), allocatable, intent(inout) :: rank
    end subroutine request_group_cache_update
    module subroutine request_group_cache_get_comm (cache, comm)
      class(request_group_cache_t), intent(in) :: cache
      type(MPI_COMM), intent(out) :: comm
    end subroutine request_group_cache_get_comm
    module function request_group_cache_is_master (cache) result (flag)
      class(request_group_cache_t), intent(in) :: cache
      integer :: rank, error
      logical :: flag
    end function request_group_cache_is_master
    module subroutine request_base_init (req, comm)
      class(request_base_t), intent(out) :: req
      type(MPI_COMM), intent(in) :: comm
    end subroutine request_base_init
    module subroutine request_base_write (req, unit)
      class(request_base_t), intent(in) :: req
      integer, intent(in), optional :: unit
    end subroutine request_base_write
    module function request_base_is_master (req) result (flag)
      class(request_base_t), intent(in) :: req
      integer :: rank, ierr
      logical :: flag
    end function request_base_is_master
    module subroutine request_base_get_external_comm (req, comm)
      class(request_base_t), intent(in) :: req
      type(MPI_COMM), intent(out) :: comm
    end subroutine request_base_get_external_comm
    module subroutine request_base_add_balancer (req, balancer)
      class(request_base_t), intent(inout) :: req
      class(balancer_base_t), allocatable, intent(inout) :: balancer
    end subroutine request_base_add_balancer
    module subroutine request_base_add_handler (req, handler_id, handler)
      class(request_base_t), intent(inout) :: req
      integer, intent(in) :: handler_id
      class(request_handler_t), pointer, intent(in) :: handler
    end subroutine request_base_add_handler
    module subroutine request_base_reset (req, deallocate_balancer)
      class(request_base_t), intent(inout) :: req
      logical, intent(in), optional :: deallocate_balancer
      logical :: flag
    end subroutine request_base_reset
    module subroutine request_base_call_handler &
         (req, handler_id, source_rank)
      class(request_base_t), intent(inout) :: req
      integer, intent(in) :: handler_id
      integer, intent(in) :: source_rank
    end subroutine request_base_call_handler
    module subroutine request_base_call_client_handler (req, handler_id)
      class(request_base_t), intent(inout) :: req
      integer, intent(in) :: handler_id
    end subroutine request_base_call_client_handler
    module subroutine request_base_await_handler (req)
      class(request_base_t), intent(inout) :: req
    end subroutine request_base_await_handler
    module subroutine request_base_barrier (req)
      class(request_base_t), intent(in) :: req
      integer :: error
    end subroutine request_base_barrier
  end interface

end module request_base
