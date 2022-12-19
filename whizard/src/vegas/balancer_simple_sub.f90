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

submodule (balancer_simple) balancer_simple_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine balancer_simple_init (balancer, n_workers, n_resources)
    class(balancer_simple_t), intent(out) :: balancer
    integer, intent(in) :: n_workers
    integer, intent(in) :: n_resources
    type(resource_state_t), dimension(:), allocatable :: state
    call balancer%base_init (n_workers, n_resources)
    allocate (balancer%parallel_grid(n_resources), source = .false.)
    allocate (state (N_BALANCER_SIMPLE_STATES))
    call state(BALANCER_SIMPLE_CHANNEL)%init ( &
         mode = STATE_SINGLE, &
         n_workers = balancer%n_workers)
    call balancer%add_state (state)
  end subroutine balancer_simple_init

  module subroutine balancer_simple_write (balancer, unit)
    class(balancer_simple_t), intent(in) :: balancer
    integer, intent(in), optional :: unit
    integer :: u, n_size
    u = given_output_unit (unit)
    call balancer%base_write (u)
    n_size = min (25, size (balancer%parallel_grid))
    write (u, "(A,25(1X,L1))") "Parallel Grids:", balancer%parallel_grid(:n_size)
  end subroutine balancer_simple_write

  module subroutine balancer_simple_update_state &
       (balancer, worker_id, parallel_grid)
    class(balancer_simple_t), intent(inout) :: balancer
    integer, intent(in) :: worker_id
    logical, dimension(:), intent(in) :: parallel_grid
    integer :: ch, worker
    balancer%parallel_grid = parallel_grid
    if (.not. allocated (balancer%state)) then
       call msg_bug ("Error: balancer state not allocated.")
    end if
    associate (state => balancer%state(BALANCER_SIMPLE_CHANNEL))
      call state%clear ()
      do ch = 1, balancer%n_resources
         if (parallel_grid(ch)) then
            call state%add_resource (ch)
         else
            worker = balancer%map_channel_to_worker (ch)
            if (worker == worker_id) then
               call state%add_resource (ch)
            end if
         end if
      end do
      call state%freeze ()
    end associate
  end subroutine balancer_simple_update_state

  pure module function balancer_simple_has_resource_group &
       (balancer, resource_id) result (flag)
    class(balancer_simple_t), intent(in) :: balancer
    integer, intent(in) :: resource_id
    logical :: flag
    if (.not. balancer%resource(resource_id)%is_active ()) then
       flag = .false.
       return
    end if
    flag = balancer%parallel_grid (resource_id)
  end function balancer_simple_has_resource_group

  pure module subroutine balancer_simple_get_resource_group &
       (balancer, resource_id, group)
    class(balancer_simple_t), intent(in) :: balancer
    integer, intent(in) :: resource_id
    integer, dimension(:), allocatable, intent(out) :: group
    integer :: i
    if (.not. balancer%has_resource_group (resource_id)) return
    group = pack ([(i, i=1,balancer%n_workers)], &
         mask = balancer%worker%get_resource () == resource_id)
  end subroutine balancer_simple_get_resource_group

  pure module function balancer_simple_get_resource_master &
       (balancer, resource_id) result (worker_id)
    class(balancer_simple_t), intent(in) :: balancer
    integer, intent(in) :: resource_id
    integer :: worker_id
    !! \note Do NOT check on resource activation (see interface prescription).
    if (balancer%parallel_grid(resource_id)) then
       worker_id = 1
    else
       worker_id = balancer%map_channel_to_worker (resource_id)
    end if
  end function balancer_simple_get_resource_master

  pure module function balancer_simple_map_channel_to_worker &
       (balancer, channel) result (worker)
    class(balancer_simple_t), intent(in) :: balancer
    integer, intent(in) :: channel
    integer :: worker
    !! Proof: channel \in {1, N_c}, number of workers N, rank \in {0, …, N - 1}
    !! Proof: worker \in {1, …, N}
    !! a = b mod c, then 0 <= a < c
    worker = mod (channel - 1, balancer%n_workers) + 1
  end function balancer_simple_map_channel_to_worker

  module subroutine balancer_simple_assign_worker (balancer, worker_id, resource_id)
    class(balancer_simple_t), intent(inout) :: balancer
    integer, intent(in) :: worker_id
    integer, intent(out) :: resource_id
    integer :: i
    if (.not. balancer%is_assignable (worker_id)) then
       resource_id = -1
       RETURN
    end if
    if (balancer%worker(worker_id)%is_assigned ()) then
       resource_id = balancer%worker(worker_id)%get_resource ()
       RETURN
    end if
    associate (state => balancer%state(BALANCER_SIMPLE_CHANNEL))
      if (.not. state%has_resource ()) then
         resource_id = 0
         return
      end if
      resource_id = state%assign_resource ()
      if (balancer%parallel_grid(resource_id)) then
         do i = 1, balancer%n_workers
            if (balancer%is_worker_pending (i)) then
               write (msg_buffer, "(A,1X,I0,1X,A,1X,I0,1X,A)") &
                    "WORKER", i, "ASSIGNED"
               call msg_bug ()
            end if
            call balancer%worker(i)%add_resource (resource_id)
         end do
         call balancer%resource(resource_id)%set_active &
              (n_workers = balancer%n_workers)
      else
         call balancer%worker(worker_id)%add_resource (resource_id)
         call balancer%resource(resource_id)%set_active (n_workers = 1)
      end if
    end associate
  end subroutine balancer_simple_assign_worker

  module subroutine balancer_simple_free_worker &
       (balancer, worker_id, resource_id)
    class(balancer_simple_t), intent(inout) :: balancer
    integer, intent(in) :: worker_id
    integer, intent(in) :: resource_id
    integer :: i
    if (.not. balancer%worker(worker_id)%is_assigned ()) return
    if (.not. resource_id == &
         balancer%worker(worker_id)%get_resource ()) then
       call msg_bug ("Balancer simple: resource and " // &
            "associated resource do not match.")
    end if
    associate (state => balancer%state(BALANCER_SIMPLE_CHANNEL))
      call balancer%resource(resource_id)%set_inactive ()
      call state%free_resource (resource_id)
      if (balancer%parallel_grid(resource_id)) then
         do i = 1, balancer%n_workers
            call balancer%worker(i)%free ()
         end do
      else
         call balancer%worker(worker_id)%free ()
      end if
    end associate
  end subroutine balancer_simple_free_worker


end submodule balancer_simple_s

