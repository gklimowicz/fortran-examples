module request_callback

  use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
  use binary_tree

  use mpi_f08 !NODEP!

  implicit none
  private

  public :: request_handler_t
  public :: request_handler_manager_t

  type, abstract :: request_handler_t
     integer :: n_requests = 0
     integer :: tag_offset = 0
     type(MPI_REQUEST), dimension(:), allocatable :: request
     type(MPI_STATUS), dimension(:), allocatable :: status
     logical :: activated = .false.
     logical :: finished = .false.
   contains
     procedure :: base_write => request_handler_base_write
     procedure(request_handler_write), deferred :: write
     !! \todo{sbrass} implement initialization procedure.
     procedure(request_handler_handle), deferred :: handle
     procedure(request_handler_client_handle), deferred :: client_handle
     procedure :: allocate => request_handler_allocate
     procedure :: get_status => request_handler_get_status
     procedure :: testall => request_handler_testall
     procedure :: waitall => request_handler_waitall
     procedure :: free => request_handler_free
  end type request_handler_t

  type :: request_handler_manager_t
     private
     type(MPI_COMM) :: comm
     type(binary_tree_t) :: tree
   contains
     procedure :: init => request_handler_manager_init
     procedure :: write => request_handler_manager_write
     procedure :: add => request_handler_manager_add
     procedure :: clear => request_handler_manager_clear
     procedure, private :: fill_status => request_handler_manager_fill_status
     procedure :: test => request_handler_manager_test
     procedure :: wait => request_handler_manager_wait
     procedure :: waitall => request_handler_manager_waitall
     procedure, private :: handler_at => request_handler_manager_handler_at
     procedure :: has_handler => request_handler_manager_has_handler
     procedure :: callback => request_handler_manager_callback
     procedure :: client_callback => request_handler_manager_client_callback
  end type request_handler_manager_t


  abstract interface
     subroutine request_handler_write (handler, unit)
       import :: request_handler_t
       class(request_handler_t), intent(in) :: handler
       integer, intent(in), optional :: unit
     end subroutine request_handler_write
  end interface

  abstract interface
     !> Handle a request from server side.
     !!
     !! The message tag can be used in order to uniquify the respective messages between master and slave.
     !! E.g. by explicitly setting it, or by using it in a computation i * N_R + j, i \in {1, …, N} and j \in {1, …, N_R}.
     !!
     !! Must set *activated* to .true. when called.
     !! \param[in] source Integer rank of the source in comm.
     !! \param[in] tag Specify the message tag.
     !! \param[in] comm MPI communicator.
     subroutine request_handler_handle (handler, source_rank, comm)
       import :: request_handler_t, MPI_COMM
       class(request_handler_t), intent(inout) :: handler
       integer, intent(in) :: source_rank
       type(MPI_COMM), intent(in) :: comm
     end subroutine request_handler_handle
  end interface

  abstract interface
     !> Handle a request from client side.
     !!
     !! Must set *activated* to .true. when called.
     !! \param[in] rank Integer of the receiver in comm.
     !! \param[in] tag Specify the message tag.
     !! \param[in] comm MPI communicator.
     subroutine request_handler_client_handle (handler, dest_rank, comm)
       import :: request_handler_t, MPI_COMM
       class(request_handler_t), intent(inout) :: handler
       integer, intent(in) :: dest_rank
       type(MPI_COMM), intent(in) :: comm
     end subroutine request_handler_client_handle
  end interface


  interface
    module subroutine request_handler_base_write (handler, unit)
      class(request_handler_t), intent(in) :: handler
      integer, intent(in), optional :: unit
    end subroutine request_handler_base_write
    module subroutine request_handler_allocate (handler, n_requests, tag_offset)
      class(request_handler_t), intent(inout) :: handler
      integer, intent(in) :: n_requests
      integer, intent(in) :: tag_offset
    end subroutine request_handler_allocate
    module subroutine request_handler_get_status (handler)
      class(request_handler_t), intent(inout) :: handler
    end subroutine request_handler_get_status
    module subroutine request_handler_waitall (handler)
      class(request_handler_t), intent(inout) :: handler
    end subroutine request_handler_waitall
    module function request_handler_testall (handler) result (flag)
      class(request_handler_t), intent(inout) :: handler
      logical :: flag
    end function request_handler_testall
    module subroutine request_handler_free (handler)
      class(request_handler_t), intent(inout) :: handler
    end subroutine request_handler_free
    module subroutine request_handler_manager_init (rhm, comm)
      class(request_handler_manager_t), intent(out) :: rhm
      type(MPI_COMM), intent(in) :: comm
    end subroutine request_handler_manager_init
    module subroutine request_handler_manager_write (rhm, unit)
      class(request_handler_manager_t), intent(in) :: rhm
      integer, intent(in), optional :: unit
    end subroutine request_handler_manager_write
    module subroutine request_handler_manager_add (rhm, handler_id, handler)
      class(request_handler_manager_t), intent(inout) :: rhm
      integer, intent(in) :: handler_id
      class(request_handler_t), pointer, intent(in) :: handler
    end subroutine request_handler_manager_add
    module subroutine request_handler_manager_clear (rhm)
      class(request_handler_manager_t), intent(inout) :: rhm
    end subroutine request_handler_manager_clear
    module subroutine request_handler_manager_fill_status (rhm)
      class(request_handler_manager_t), intent(inout) :: rhm
    end subroutine request_handler_manager_fill_status
    module function request_handler_manager_test &
         (rhm, handler_id) result (flag)
      class(request_handler_manager_t), intent(inout) :: rhm
      integer, intent(in) :: handler_id
      logical :: flag
    end function request_handler_manager_test
    module subroutine request_handler_manager_wait (rhm, handler_id)
      class(request_handler_manager_t), intent(inout) :: rhm
      integer, intent(in) :: handler_id
    end subroutine request_handler_manager_wait
    module subroutine request_handler_manager_waitall (rhm)
      class(request_handler_manager_t), intent(inout) :: rhm
    end subroutine request_handler_manager_waitall
    module subroutine request_handler_manager_handler_at &
         (rhm, handler_id, handler)
      class(request_handler_manager_t), intent(in) :: rhm
      integer, intent(in) :: handler_id
      class(request_handler_t), pointer, intent(out) :: handler
    end subroutine request_handler_manager_handler_at
    module function request_handler_manager_has_handler &
         (rhm, handler_id) result (flag)
      class(request_handler_manager_t), intent(inout) :: rhm
      integer, intent(in) :: handler_id
      logical :: flag
    end function request_handler_manager_has_handler
    module subroutine request_handler_manager_callback &
         (rhm, handler_id, source_rank)
      class(request_handler_manager_t), intent(inout) :: rhm
      integer, intent(in) :: handler_id
      integer, intent(in) :: source_rank
    end subroutine request_handler_manager_callback
    module subroutine request_handler_manager_client_callback &
         (rhm, handler_id, dest_rank)
      class(request_handler_manager_t), intent(inout) :: rhm
      integer, intent(in) :: handler_id
      integer, intent(in) :: dest_rank
    end subroutine request_handler_manager_client_callback
  end interface

end module request_callback
