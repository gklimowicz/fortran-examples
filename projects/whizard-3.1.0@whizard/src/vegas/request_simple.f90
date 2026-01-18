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
module request_simple

  use array_list

  use balancer_base
  use balancer_simple
  use request_base

  use mpi_f08 !NODEP!

  implicit none
  private

  public :: request_simple_t

  type, extends (request_base_t) :: request_simple_t
     integer :: n_workers = 0
     integer :: n_channels = 0
     logical, dimension(:), allocatable :: parallel_grid
   contains
     procedure :: init => request_simple_init
     procedure :: update => request_simple_update
     procedure :: write => request_simple_write
     procedure :: has_workers => request_simple_has_workers
     procedure :: get_request_master => request_simple_get_request_master
     !! deferred.
     procedure :: request_workload => request_simple_request_workload
     procedure :: release_workload => request_simple_release_workload
     procedure :: handle_and_release_workload => &
          request_simple_handle_and_release_workload
  end type request_simple_t


  interface
    module subroutine request_simple_update (req, parallel_grid)
      class(request_simple_t), intent(inout) :: req
      logical, dimension(:), intent(in) :: parallel_grid
    end subroutine request_simple_update
    module subroutine request_simple_write (req, unit)
      class(request_simple_t), intent(in) :: req
      integer, intent(in), optional :: unit
    end subroutine request_simple_write
    module function request_simple_has_workers (req) result (flag)
      class(request_simple_t), intent(in) :: req
      logical :: flag
    end function request_simple_has_workers
    module function request_simple_get_request_master &
         (req, channel) result (rank)
      class(request_simple_t), intent(in) :: req
      integer, intent(in) :: channel
      integer :: rank
    end function request_simple_get_request_master
    module subroutine request_simple_request_workload (req, request)
      class(request_simple_t), intent(inout) :: req
      type(request_t), intent(out) :: request
    end subroutine request_simple_request_workload
    module subroutine request_simple_release_workload (req, request)
      class(request_simple_t), intent(inout) :: req
      type(request_t), intent(in) :: request
    end subroutine request_simple_release_workload
    module subroutine request_simple_handle_and_release_workload (req, request)
      class(request_simple_t), intent(inout) :: req
      type(request_t), intent(in) :: request
    end subroutine request_simple_handle_and_release_workload
  end interface

contains

  module subroutine request_simple_init (req, comm, n_channels)
    class(request_simple_t), intent(out) :: req
    type(MPI_COMM), intent(in) :: comm
    integer, intent(in) :: n_channels
    integer :: n_workers
    call req%base_init (comm)
    call MPI_COMM_SIZE (req%comm, req%n_workers)
    req%n_channels = n_channels
    allocate (req%parallel_grid (n_channels), source = .false.)
    call allocate_balancer ()
  contains
    subroutine allocate_balancer ()
      class(balancer_base_t), allocatable :: balancer
      allocate (balancer_simple_t :: balancer)
      select type (balancer)
      type is (balancer_simple_t)
         call balancer%init (n_workers = req%n_workers, n_resources = req%n_channels)
      end select
      call req%add_balancer (balancer)
    end subroutine allocate_balancer
  end subroutine request_simple_init


end module request_simple
