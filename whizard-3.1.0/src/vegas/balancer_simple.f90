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
module balancer_simple

  use balancer_base

  implicit none
  private

  integer, parameter :: N_BALANCER_SIMPLE_STATES = 1, &
       BALANCER_SIMPLE_CHANNEL = 1

  public :: balancer_simple_t

  type, extends (balancer_base_t) :: balancer_simple_t
     logical, dimension(:), allocatable :: parallel_grid
   contains
     procedure :: init => balancer_simple_init
     procedure :: write => balancer_simple_write
     procedure :: update_state => balancer_simple_update_state
     procedure :: has_resource_group => balancer_simple_has_resource_group
     procedure :: get_resource_group => balancer_simple_get_resource_group
     procedure :: get_resource_master => balancer_simple_get_resource_master
     procedure, private :: map_channel_to_worker => &
          balancer_simple_map_channel_to_worker
     procedure :: assign_worker => balancer_simple_assign_worker
     procedure :: free_worker => balancer_simple_free_worker
  end type balancer_simple_t


  interface
    module subroutine balancer_simple_init (balancer, n_workers, n_resources)
      class(balancer_simple_t), intent(out) :: balancer
      integer, intent(in) :: n_workers
      integer, intent(in) :: n_resources
    end subroutine balancer_simple_init
    module subroutine balancer_simple_write (balancer, unit)
      class(balancer_simple_t), intent(in) :: balancer
      integer, intent(in), optional :: unit
    end subroutine balancer_simple_write
    module subroutine balancer_simple_update_state &
         (balancer, worker_id, parallel_grid)
      class(balancer_simple_t), intent(inout) :: balancer
      integer, intent(in) :: worker_id
      logical, dimension(:), intent(in) :: parallel_grid
    end subroutine balancer_simple_update_state
    pure module function balancer_simple_has_resource_group &
         (balancer, resource_id) result (flag)
      class(balancer_simple_t), intent(in) :: balancer
      integer, intent(in) :: resource_id
      logical :: flag
    end function balancer_simple_has_resource_group
    pure module subroutine balancer_simple_get_resource_group &
         (balancer, resource_id, group)
      class(balancer_simple_t), intent(in) :: balancer
      integer, intent(in) :: resource_id
      integer, dimension(:), allocatable, intent(out) :: group
    end subroutine balancer_simple_get_resource_group
    pure module function balancer_simple_get_resource_master &
         (balancer, resource_id) result (worker_id)
      class(balancer_simple_t), intent(in) :: balancer
      integer, intent(in) :: resource_id
      integer :: worker_id
    end function balancer_simple_get_resource_master
    pure module function balancer_simple_map_channel_to_worker &
         (balancer, channel) result (worker)
      class(balancer_simple_t), intent(in) :: balancer
      integer, intent(in) :: channel
      integer :: worker
    end function balancer_simple_map_channel_to_worker
    module subroutine balancer_simple_assign_worker (balancer, worker_id, resource_id)
      class(balancer_simple_t), intent(inout) :: balancer
      integer, intent(in) :: worker_id
      integer, intent(out) :: resource_id
    end subroutine balancer_simple_assign_worker
    module subroutine balancer_simple_free_worker &
         (balancer, worker_id, resource_id)
      class(balancer_simple_t), intent(inout) :: balancer
      integer, intent(in) :: worker_id
      integer, intent(in) :: resource_id
    end subroutine balancer_simple_free_worker
  end interface

end module balancer_simple
