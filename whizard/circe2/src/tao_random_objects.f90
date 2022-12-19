module tao_random_objects
  use kinds
  use tao_random_numbers
  use circe2
  implicit none
  private
  public :: rng_tao
  type, extends (rng_type) :: rng_tao
     integer :: seed = 0
     integer :: n_calls = 0
     type(tao_random_state) :: state
   contains
     procedure :: generate => rng_tao_generate
     procedure :: init => rng_tao_init
  end type rng_tao
contains
  subroutine rng_tao_generate (rng_obj, u)
    class(rng_tao), intent(inout) :: rng_obj
    real(default), intent(out) :: u
    call tao_random_number (rng_obj%state, u)
    rng_obj%n_calls = rng_obj%n_calls + 1
  end subroutine rng_tao_generate
  subroutine rng_tao_init (rng_obj, seed)
    class(rng_tao), intent(inout) :: rng_obj
    integer, intent(in) :: seed
    rng_obj%seed = seed
    call tao_random_create (rng_obj%state, seed)
  end subroutine rng_tao_init
end module tao_random_objects
