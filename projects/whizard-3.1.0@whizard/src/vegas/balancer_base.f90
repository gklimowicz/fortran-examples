module balancer_base

  use array_list

  implicit none
  private

  public :: resource_state_t
  public :: balancer_base_t
  public :: shift_rank_to_worker
  public :: shift_worker_to_rank

  integer, parameter, public :: STATE_SINGLE = 1, &
       STATE_ALL = 2


  type :: worker_t
     private
     integer :: resource = 0
     integer :: state = 0
     integer :: n_resources = 0
     logical :: assigned = .false.
   contains
     procedure :: write => worker_write
     procedure :: is_assigned => worker_is_assigned
     procedure :: get_resource => worker_get_resource
     procedure :: get_state => worker_get_state
     procedure :: add_resource => worker_add_resource
     procedure :: free => worker_free
  end type worker_t

  type :: resource_t
     private
     integer :: resource_id = 0
     logical :: active = .false.
     integer :: n_assigned_workers = 0
   contains
     procedure :: write => resource_write
     procedure :: is_active => resource_is_active
     procedure :: set_active => resource_set_active
     procedure :: set_inactive => resource_set_inactive
  end type resource_t

  type :: resource_state_t
     integer :: n_workers = 0
     integer :: mode = 0
     type(array_list_t) :: resource_stack
     type(array_list_t) :: finished_stack
   contains
     procedure :: write => resource_state_write
     procedure :: init => resource_state_init
     procedure :: add_resource => resource_state_add_resource
     procedure :: freeze => resource_state_freeze
     procedure :: clear => resource_state_clear
     procedure :: has_resource => resource_state_has_resource
     procedure :: assign_resource => resource_state_assign_resource
     procedure :: free_resource => resource_state_free_resource
  end type resource_state_t

  type, abstract :: balancer_base_t
     integer :: n_workers = 0
     integer :: n_resources = 0
     integer :: n_states = 0
     type(worker_t), dimension(:), allocatable :: worker
     type(resource_t), dimension(:), allocatable :: resource
     type(resource_state_t), dimension(:), allocatable :: state
   contains
     procedure :: base_write => balancer_base_write
     procedure(balancer_base_deferred_write), deferred :: write
     procedure :: base_init => balancer_base_base_init
     procedure :: add_state => balancer_base_add_state
     procedure, private :: link_worker_and_state => &
          balancer_base_link_worker_and_state
     procedure :: is_assignable => balancer_base_is_assignable
     procedure :: is_worker_pending => balancer_base_is_worker_pending
     procedure :: is_pending => balancer_base_is_pending
     procedure(balancer_base_has_resource_group), deferred :: has_resource_group
     procedure(balancer_base_get_resource_group), deferred :: get_resource_group
     procedure(balancer_base_get_resource_master), deferred :: get_resource_master
     procedure(balancer_base_assign_worker), deferred :: assign_worker
     procedure(balancer_base_free_worker), deferred :: free_worker
  end type balancer_base_t


  abstract interface
     subroutine balancer_base_deferred_write (balancer, unit)
       import :: balancer_base_t
       class(balancer_base_t), intent(in) :: balancer
       integer, intent(in), optional :: unit
     end subroutine balancer_base_deferred_write
   end interface

     !> Has resource an associated resource group.
     !!
     !! \note .true. only on an active resource, else .false.
  abstract interface
     pure logical function balancer_base_has_resource_group (balancer, resource_id) &
          result (flag)
       import :: balancer_base_t
       class(balancer_base_t), intent(in) :: balancer
       integer, intent(in) :: resource_id
     end function balancer_base_has_resource_group
   end interface

     !> Get resource group.
     !!
     !! \note Implementation must check against group existence.
     !! \return group (allocated|NOT allocated for (inactive|non-group) resource)
  abstract interface
     pure subroutine balancer_base_get_resource_group (balancer, resource_id, group)
       import :: balancer_base_t
       class(balancer_base_t), intent(in) :: balancer
       integer, intent(in) :: resource_id
       integer, dimension(:), allocatable, intent(out) :: group
     end subroutine balancer_base_get_resource_group
   end interface

     !> Get resource master (worker).
     !!
     !! Return worker as given, however, if extended type is used in a non-local
     !! or in combination with a commnuicative request type, then check on activation status of associated resource.
     !!
     !! \return worker Valid worker index (\in {1, …, N}) only on active resource*, else -1.
  abstract interface
     pure integer function balancer_base_get_resource_master (balancer, resource_id) &
          result (worker)
       import :: balancer_base_t
       class(balancer_base_t), intent(in) :: balancer
       integer, intent(in) :: resource_id
     end function balancer_base_get_resource_master
  end interface

     !> Assign resource to a given worker or retrieve current assigned resource.
     !!
     !! If worker has already a resource assigned, return resource.
     !! If worker has not been assigned a resource, retrieve new resource from state.
     !!
     !! \note Each call must check if a worker is assignable, if not, the procedure must return resource_id = -1.
  abstract interface
     subroutine balancer_base_assign_worker (balancer, worker_id, resource_id)
       import :: balancer_base_t
       class(balancer_base_t), intent(inout) :: balancer
       integer, intent(in) :: worker_id
       integer, intent(out) :: resource_id
     end subroutine balancer_base_assign_worker
   end interface

  abstract interface
     subroutine balancer_base_free_worker (balancer, worker_id, resource_id)
       import :: balancer_base_t
       class(balancer_base_t), intent(inout) :: balancer
       integer, intent(in) :: worker_id
       integer, intent(in) :: resource_id
     end subroutine balancer_base_free_worker
  end interface


  interface
    elemental module function shift_rank_to_worker (rank) result (worker)
      integer, intent(in) :: rank
      integer :: worker
    end function shift_rank_to_worker
    elemental module function shift_worker_to_rank (worker) result (rank)
      integer, intent(in) :: worker
      integer :: rank
    end function shift_worker_to_rank
    module subroutine worker_write (worker, unit)
      class(worker_t), intent(in) :: worker
      integer, intent(in), optional :: unit
    end subroutine worker_write
    elemental module function worker_is_assigned (worker) result (flag)
      class(worker_t), intent(in) :: worker
      logical :: flag
    end function worker_is_assigned
    elemental module function worker_get_resource (worker) result (resource_id)
      class(worker_t), intent(in) :: worker
      integer :: resource_id
    end function worker_get_resource
    elemental module function worker_get_state (worker) result (i_state)
      class(worker_t), intent(in) :: worker
      integer :: i_state
    end function worker_get_state
    elemental module subroutine worker_add_resource (worker, resource_id)
      class(worker_t), intent(inout) :: worker
      integer, intent(in) :: resource_id
    end subroutine worker_add_resource
    elemental module subroutine worker_free (worker)
      class(worker_t), intent(inout) :: worker
    end subroutine worker_free
    module subroutine resource_write (resource, unit)
      class(resource_t), intent(in) :: resource
      integer, intent(in), optional :: unit
    end subroutine resource_write
    elemental module function resource_is_active (resource) result (flag)
     class(resource_t), intent(in) :: resource
      logical :: flag
    end function resource_is_active
    module subroutine resource_set_active (resource, n_workers)
      class(resource_t), intent(inout) :: resource
      integer, intent(in) :: n_workers
    end subroutine resource_set_active
    module subroutine resource_set_inactive (resource)
      class(resource_t), intent(inout) :: resource
    end subroutine resource_set_inactive
    module subroutine resource_state_write (state, unit)
      class(resource_state_t), intent(in) :: state
      integer, intent(in), optional :: unit
    end subroutine resource_state_write
    module subroutine resource_state_init (state, mode, n_workers)
      class(resource_state_t), intent(out) :: state
      integer, intent(in) :: mode
      integer, intent(in) :: n_workers
    end subroutine resource_state_init
    module subroutine resource_state_add_resource (state, i_resource)
      class(resource_state_t), intent(inout) :: state
      integer, intent(in) :: i_resource
    end subroutine resource_state_add_resource
    module subroutine resource_state_freeze (state)
      class(resource_state_t), intent(inout) :: state
    end subroutine resource_state_freeze
    module subroutine resource_state_clear (state)
      class(resource_state_t), intent(inout) :: state
    end subroutine resource_state_clear
    elemental module function resource_state_has_resource (state) result (flag)
      class(resource_state_t), intent(in) :: state
      logical :: flag
    end function resource_state_has_resource
    module function resource_state_assign_resource (state) result (i_resource)
      class(resource_state_t), intent(inout) :: state
      integer :: i_resource
    end function resource_state_assign_resource
    module subroutine resource_state_free_resource (state, i_resource)
      class(resource_state_t), intent(inout) :: state
      integer, intent(in) :: i_resource
    end subroutine resource_state_free_resource
    module subroutine balancer_base_write (balancer, unit)
      class(balancer_base_t), intent(in) :: balancer
      integer, intent(in), optional :: unit
    end subroutine balancer_base_write
    module subroutine balancer_base_base_init (balancer, n_workers, n_resources)
      class(balancer_base_t), intent(out) :: balancer
      integer, intent(in) :: n_workers
      integer, intent(in) :: n_resources
    end subroutine balancer_base_base_init
    module subroutine balancer_base_add_state (balancer, state)
      class(balancer_base_t), intent(inout) :: balancer
      type(resource_state_t), dimension(:), allocatable, intent(inout) :: state
    end subroutine balancer_base_add_state
    module subroutine balancer_base_link_worker_and_state (balancer)
      class(balancer_base_t), intent(inout) :: balancer
    end subroutine balancer_base_link_worker_and_state
    pure module function balancer_base_is_assignable &
         (balancer, worker_id) result (flag)
      class(balancer_base_t), intent(in) :: balancer
      integer, intent(in) :: worker_id
      integer :: i_state, resource_id
      logical :: flag
    end function balancer_base_is_assignable
    pure module function balancer_base_is_worker_pending &
         (balancer, worker_id) result (flag)
      class(balancer_base_t), intent(in) :: balancer
      integer, intent(in) :: worker_id
      integer :: resource_id
      logical :: flag
    end function balancer_base_is_worker_pending
    module function balancer_base_is_pending (balancer) result (flag)
      class(balancer_base_t), intent(in) :: balancer
      logical :: flag
    end function balancer_base_is_pending
  end interface

end module balancer_base
