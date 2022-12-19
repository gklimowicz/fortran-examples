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
module balancer_channel

  use kinds, only: default
  use balancer_base

  implicit none
  private

  public :: balancer_channel_t

  real(default), parameter :: BETA = 1.5_default

  integer, parameter :: N_BALANCER_CHANNEL_STATE = 2, &
       CHANNEL_STATE = 1, &
       GRID_STATE = 2


  type, extends(balancer_base_t) :: balancer_channel_t
     private
     integer :: n_parallel_grids = 0
     integer :: n_parallel_channels = 0
     integer :: n_grid_workers = 0
     integer :: n_channel_workers = 0
     logical, dimension(:), allocatable :: parallel_grid
   contains
     procedure :: init => balancer_channel_init
     procedure :: write => balancer_channel_write
     procedure :: update_state => balancer_channel_update_state
     procedure :: has_resource_group => balancer_channel_has_resource_group
     procedure :: get_resource_group => balancer_channel_get_resource_group
     procedure :: get_resource_master => balancer_channel_get_resource_master
     procedure :: assign_worker => balancer_channel_assign_worker
     procedure :: free_worker => balancer_channel_free_worker
  end type balancer_channel_t


  interface
    module subroutine balancer_channel_init (balancer, n_workers, n_resources)
      class(balancer_channel_t), intent(out), target :: balancer
      integer, intent(in) :: n_workers
      integer, intent(in) :: n_resources
    end subroutine balancer_channel_init
    module subroutine balancer_channel_write (balancer, unit)
      class(balancer_channel_t), intent(in) :: balancer
      integer, intent(in), optional :: unit
    end subroutine balancer_channel_write
    module subroutine balancer_channel_update_state &
         (balancer, weight, parallel_grid)
      class(balancer_channel_t), intent(inout) :: balancer
      real(default), dimension(:), intent(in) :: weight
      logical, dimension(:), intent(in) :: parallel_grid
    end subroutine balancer_channel_update_state
    pure module function balancer_channel_has_resource_group &
         (balancer, resource_id) result (flag)
      class(balancer_channel_t), intent(in) :: balancer
      integer, intent(in) :: resource_id
      logical :: flag
    end function balancer_channel_has_resource_group
    pure module subroutine balancer_channel_get_resource_group &
         (balancer, resource_id, group)
      class(balancer_channel_t), intent(in) :: balancer
      integer, intent(in) :: resource_id
      integer, dimension(:), allocatable, intent(out) :: group
    end subroutine balancer_channel_get_resource_group
    pure module function balancer_channel_get_resource_master &
         (balancer, resource_id) result (worker_id)
      class(balancer_channel_t), intent(in) :: balancer
      integer, intent(in) :: resource_id
      integer :: worker_id
    end function balancer_channel_get_resource_master
    module subroutine balancer_channel_assign_worker &
         (balancer, worker_id, resource_id)
      class(balancer_channel_t), intent(inout) :: balancer
      integer, intent(in) :: worker_id
      integer, intent(out) :: resource_id
    end subroutine balancer_channel_assign_worker
    module subroutine balancer_channel_free_worker &
         (balancer, worker_id, resource_id)
      class(balancer_channel_t), intent(inout) :: balancer
      integer, intent(in) :: worker_id
      integer, intent(in) :: resource_id
    end subroutine balancer_channel_free_worker
  end interface

end module balancer_channel
