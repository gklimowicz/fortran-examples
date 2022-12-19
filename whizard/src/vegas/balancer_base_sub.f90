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

submodule (balancer_base) balancer_base_s

  use io_units
  use diagnostics

  implicit none

contains

  !> Shift rank index to worker index.
  !! Proof: rank \in {0, …, N - 1}, worker \in {1, …, N}
  elemental module function shift_rank_to_worker (rank) result (worker)
    integer, intent(in) :: rank
    integer :: worker
    worker = rank + 1
  end function shift_rank_to_worker

  !> Shift worker index to rank index.
  !! Proof: rank \in {0, …, N - 1}, worker \in {1, …, N}
  elemental module function shift_worker_to_rank (worker) result (rank)
    integer, intent(in) :: worker
    integer :: rank
    rank = worker - 1
  end function shift_worker_to_rank

  module subroutine worker_write (worker, unit)
    class(worker_t), intent(in) :: worker
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3(A,1X,I3,1X),A,1X,L1)") "RESOURCE", worker%resource, &
         "STATE", worker%state, &
         "N_RESOURCES", worker%n_resources, &
         "ASSIGNED", worker%assigned
  end subroutine worker_write

  elemental module function worker_is_assigned (worker) result (flag)
    class(worker_t), intent(in) :: worker
    logical :: flag
    flag = worker%assigned
  end function worker_is_assigned

  elemental module function worker_get_resource (worker) result (resource_id)
    class(worker_t), intent(in) :: worker
    integer :: resource_id
    resource_id = worker%resource
  end function worker_get_resource

  elemental module function worker_get_state (worker) result (i_state)
    class(worker_t), intent(in) :: worker
    integer :: i_state
    i_state = worker%state
  end function worker_get_state

  elemental module subroutine worker_add_resource (worker, resource_id)
    class(worker_t), intent(inout) :: worker
    integer, intent(in) :: resource_id
    worker%n_resources = worker%n_resources + 1
    worker%assigned = .true.
    worker%resource = resource_id
  end subroutine worker_add_resource

  elemental module subroutine worker_free (worker)
    class(worker_t), intent(inout) :: worker
    worker%assigned = .false.
    worker%resource = 0
  end subroutine worker_free

  module subroutine resource_write (resource, unit)
    class(resource_t), intent(in) :: resource
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(A,1X,I3,1X,A,1X,L1,1X,A,1X,I3)") &
         "RESOURCE_ID", resource%resource_id, &
         "ACTIVE", resource%active, &
         "N_ASSIGNED_WORKERS", resource%n_assigned_workers
  end subroutine resource_write

  elemental module function resource_is_active (resource) result (flag)
    class(resource_t), intent(in) :: resource
    logical :: flag
    flag = resource%active
  end function resource_is_active

  module subroutine resource_set_active (resource, n_workers)
    class(resource_t), intent(inout) :: resource
    integer, intent(in) :: n_workers
    resource%active = .true.
    resource%n_assigned_workers = n_workers
  end subroutine resource_set_active

  module subroutine resource_set_inactive (resource)
    class(resource_t), intent(inout) :: resource
    resource%active = .false.
  end subroutine resource_set_inactive

  module subroutine resource_state_write (state, unit)
    class(resource_state_t), intent(in) :: state
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(A,1X,I0)") "N_STATE_WORKERS", state%n_workers
    select case (state%mode)
    case (STATE_SINGLE)
       write (u, "(A)") "MODE ONE-TO-ONE"
    case (STATE_ALL)
       write (u, "(A)") "MODE ALL-TO-ONE"
    case default
       write (u, "(A)") "UNSUPPORTED MODE"
    end select
    write (u, "(A)") "RESOURCE"
    call state%resource_stack%write (u)
    write (u, "(A)") "FINISHED"
    call state%finished_stack%write (u)
  end subroutine resource_state_write

  module subroutine resource_state_init (state, mode, n_workers)
    class(resource_state_t), intent(out) :: state
    integer, intent(in) :: mode
    integer, intent(in) :: n_workers
    state%mode = mode
    state%n_workers = n_workers
    call state%resource_stack%init ()
    call state%finished_stack%init ()
  end subroutine resource_state_init

  module subroutine resource_state_add_resource (state, i_resource)
    class(resource_state_t), intent(inout) :: state
    integer, intent(in) :: i_resource
    call state%resource_stack%add (i_resource)
  end subroutine resource_state_add_resource

  module subroutine resource_state_freeze (state)
    class(resource_state_t), intent(inout) :: state
    call state%resource_stack%sort ()
    call state%resource_stack %reverse_order ()
  end subroutine resource_state_freeze

  module subroutine resource_state_clear (state)
    class(resource_state_t), intent(inout) :: state
    call state%resource_stack%clear ()
    call state%finished_stack%clear ()
  end subroutine resource_state_clear

  elemental module function resource_state_has_resource (state) result (flag)
    class(resource_state_t), intent(in) :: state
    logical :: flag
    flag = .not. state%resource_stack%is_empty ()
  end function resource_state_has_resource

  module function resource_state_assign_resource (state) result (i_resource)
    class(resource_state_t), intent(inout) :: state
    integer :: i_resource
    if (state%resource_stack%is_empty ()) then
       i_resource = 0
       call msg_bug ("Error: No leftover resource on stack.")
       return
    end if
    i_resource = state%resource_stack%remove () !! Pop last element from stack.
  end function resource_state_assign_resource

  module subroutine resource_state_free_resource (state, i_resource)
    class(resource_state_t), intent(inout) :: state
    integer, intent(in) :: i_resource
    if (state%resource_stack%is_element (i_resource)) then
       call msg_bug &
            ("Error: Cannot free resource, still on resource stack.")
    end if
    call state%finished_stack%add (i_resource)
  end subroutine resource_state_free_resource

  module subroutine balancer_base_write (balancer, unit)
    class(balancer_base_t), intent(in) :: balancer
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(A)") "[REQUEST BALANCER]"
    write (u, "(3(A,1X,I3,1X))") "N_WORKERS", balancer%n_workers, &
         "N_RESOURCES", balancer%n_resources, &
         "N_STATES", balancer%n_states
    write (u, "(A)") "[WORKER]"
    do i = 1, balancer%n_workers
       call balancer%worker(i)%write (u)
    end do
    write (u, "(A)") "[RESOURCE]"
    do i = 1, balancer%n_resources
       call balancer%resource(i)%write (u)
    end do
    write (u, "(A)") "[STATES]"
    do i = 1, balancer%n_states
       call balancer%state(i)%write (u)
    end do
  end subroutine balancer_base_write

  module subroutine balancer_base_base_init (balancer, n_workers, n_resources)
    class(balancer_base_t), intent(out) :: balancer
    integer, intent(in) :: n_workers
    integer, intent(in) :: n_resources
    balancer%n_workers = n_workers
    balancer%n_resources = n_resources
    allocate (balancer%worker (n_workers))
    allocate (balancer%resource (n_resources))
    call init_resource ()
  contains
    subroutine init_resource ()
      integer :: i
      do i = 1, balancer%n_resources
         balancer%resource(i)%resource_id = i
      end do
    end subroutine init_resource
  end subroutine balancer_base_base_init

  module subroutine balancer_base_add_state (balancer, state)
    class(balancer_base_t), intent(inout) :: balancer
    type(resource_state_t), dimension(:), allocatable, intent(inout) :: state
    balancer%n_states = size (state)
    call move_alloc (state, balancer%state)
    call balancer%link_worker_and_state ()
  end subroutine balancer_base_add_state

  module subroutine balancer_base_link_worker_and_state (balancer)
    class(balancer_base_t), intent(inout) :: balancer
    integer :: i, j, i_worker
    if (.not. allocated (balancer%state)) &
         call msg_bug ("Error: resource state not allocated.")
    !! Link worker to a state.
    i_worker = 1
    do i = 1, balancer%n_states
       do j = 1, balancer%state(i)%n_workers
          if (i_worker > balancer%n_workers) then
             call msg_bug ("Balancer: Number of state workers&
                  & exceeding global number of workers")
          end if
          associate (worker => balancer%worker(i_worker))
            worker%state = i
            !! Reset worker attributes.
            worker%resource = 0
            worker%n_resources = 0
            worker%assigned = .false.
          end associate
          i_worker = i_worker + 1
       end do
    end do
  end subroutine balancer_base_link_worker_and_state

  pure module function balancer_base_is_assignable &
       (balancer, worker_id) result (flag)
    class(balancer_base_t), intent(in) :: balancer
    integer, intent(in) :: worker_id
    integer :: i_state, resource_id
    logical :: flag
    flag = .false.
    if (balancer%worker(worker_id)%assigned) then
       resource_id = balancer%worker(worker_id)%resource
       flag = balancer%resource(resource_id)%is_active ()
    else
       i_state = balancer%worker(worker_id)%get_state ()
       flag = balancer%state(i_state)%has_resource ()
    end if
  end function balancer_base_is_assignable

  pure module function balancer_base_is_worker_pending &
       (balancer, worker_id) result (flag)
    class(balancer_base_t), intent(in) :: balancer
    integer, intent(in) :: worker_id
    integer :: resource_id
    logical :: flag
    flag = balancer%worker(worker_id)%assigned
    if (flag) then
       resource_id = balancer%worker(worker_id)%get_resource ()
       flag = balancer%resource(resource_id)%is_active ()
    end if
  end function balancer_base_is_worker_pending

  module function balancer_base_is_pending (balancer) result (flag)
    class(balancer_base_t), intent(in) :: balancer
    logical :: flag
    flag = all (balancer%state%has_resource ())
  end function balancer_base_is_pending


end submodule balancer_base_s

