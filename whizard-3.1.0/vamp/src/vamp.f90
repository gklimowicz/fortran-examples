! vamp.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vamp_grid_type
  use kinds
  use divisions
  private
  type, public :: vamp_grid
     ! private !: forced by \texttt{use} association in interface
     type(division_t), dimension(:), pointer :: div => null ()
     real(kind=default), dimension(:,:), pointer :: map => null ()
     real(kind=default), dimension(:), pointer :: mu_x => null ()
     real(kind=default), dimension(:), pointer :: sum_mu_x => null ()
     real(kind=default), dimension(:,:), pointer :: mu_xx => null ()
     real(kind=default), dimension(:,:), pointer :: sum_mu_xx => null ()
     real(kind=default), dimension(2) :: mu
     real(kind=default), dimension(2) :: mu_plus, mu_minus
     real(kind=default) :: sum_integral, sum_weights, sum_chi2
     real(kind=default) :: calls, dv2g, jacobi
     real(kind=default) :: f_min, f_max
     real(kind=default) :: mu_gi, sum_mu_gi
     integer, dimension(:), pointer :: num_div => null ()
     integer :: num_calls, calls_per_cell
     logical :: stratified = .true.
     logical :: all_stratified = .true.
     logical :: quadrupole = .false.
     logical :: independent
     integer :: equivalent_to_ch, multiplicity
  end type vamp_grid
end module vamp_grid_type
module vamp_equivalences
  use kinds
  use divisions
  use vamp_grid_type !NODEP!
  implicit none
  private
  public :: vamp_equivalences_init
  public :: vamp_equivalences_final
  public :: vamp_equivalences_write
  public :: vamp_equivalence_set
  public :: vamp_equivalences_complete
  integer, parameter, public :: &
       VEQ_IDENTITY = 0, VEQ_INVERT = 1, VEQ_SYMMETRIC = 2, VEQ_INVARIANT = 3
  type, public :: vamp_equivalence_t
     integer :: left, right
     integer, dimension(:), allocatable :: permutation
     integer, dimension(:), allocatable :: mode
  end type vamp_equivalence_t
  type, public :: vamp_equivalences_t
     type(vamp_equivalence_t), dimension(:), allocatable :: eq
     integer :: n_eq, n_ch
     integer, dimension(:), allocatable :: pointer
     logical, dimension(:), allocatable :: independent
     integer, dimension(:), allocatable :: equivalent_to_ch
     integer, dimension(:), allocatable :: multiplicity
     integer, dimension(:), allocatable :: symmetry
     logical, dimension(:,:), allocatable :: div_is_invariant
  end type vamp_equivalences_t
contains
  subroutine vamp_equivalence_init (eq, n_dim)
    type(vamp_equivalence_t), intent(inout) :: eq
    integer, intent(in) :: n_dim
    allocate (eq%permutation(n_dim), eq%mode(n_dim))
  end subroutine vamp_equivalence_init
  subroutine vamp_equivalences_init (eq, n_eq, n_ch, n_dim)
    type(vamp_equivalences_t), intent(inout) :: eq
    integer, intent(in) :: n_eq, n_ch, n_dim
    integer :: i
    eq%n_eq = n_eq
    eq%n_ch = n_ch
    allocate (eq%eq(n_eq))
    allocate (eq%pointer(n_ch+1))
    do i=1, n_eq
       call vamp_equivalence_init (eq%eq(i), n_dim)
    end do
    allocate (eq%independent(n_ch), eq%equivalent_to_ch(n_ch))
    allocate (eq%multiplicity(n_ch), eq%symmetry(n_ch))
    allocate (eq%div_is_invariant(n_ch, n_dim))
    eq%independent = .true.
    eq%equivalent_to_ch = 0
    eq%multiplicity = 0
    eq%symmetry = 0
    eq%div_is_invariant = .false.
  end subroutine vamp_equivalences_init
  subroutine vamp_equivalence_final (eq)
    type(vamp_equivalence_t), intent(inout) :: eq
    deallocate (eq%permutation, eq%mode)
  end subroutine vamp_equivalence_final
  subroutine vamp_equivalences_final (eq)
    type(vamp_equivalences_t), intent(inout) :: eq
  ! integer :: i
  ! do i=1, eq%n_eq
  !    call vamp_equivalence_final (eq%eq(i))
  ! end do
    if (allocated (eq%eq))  deallocate (eq%eq)
    if (allocated (eq%pointer))  deallocate (eq%pointer)
    if (allocated (eq%multiplicity))  deallocate (eq%multiplicity)
    if (allocated (eq%symmetry))  deallocate (eq%symmetry)
    if (allocated (eq%independent))  deallocate (eq%independent)
    if (allocated (eq%equivalent_to_ch))  deallocate (eq%equivalent_to_ch)
    if (allocated (eq%div_is_invariant))  deallocate (eq%div_is_invariant)
    eq%n_eq = 0
    eq%n_ch = 0
  end subroutine vamp_equivalences_final
  subroutine vamp_equivalence_write (eq, unit)
    integer, intent(in), optional :: unit
    integer :: u
    type(vamp_equivalence_t), intent(in) :: eq
    u = 6;  if (present (unit))  u = unit
    write (u, "(3x,A,2(1x,I0))") "Equivalent channels:", eq%left, eq%right
    write (u, "(5x,A,99(1x,I0))") "Permutation:", eq%permutation
    write (u, "(5x,A,99(1x,I0))") "Mode:       ", eq%mode
  end subroutine vamp_equivalence_write
  subroutine vamp_equivalences_write (eq, unit)
    type(vamp_equivalences_t), intent(in) :: eq
    integer, intent(in), optional :: unit
    integer :: u
    integer :: ch, i
    u = 6;  if (present (unit))  u = unit
    write (u, "(1x,A)") "Inequivalent channels:"
    if (allocated (eq%independent)) then
       do ch=1, eq%n_ch
          if (eq%independent(ch)) then
             write (u, "(3x,A,1x,I0,A,4x,A,I0,4x,A,I0,4x,A,999(L1))") &
                  "Channel", ch, ":", &
                  "Mult. = ", eq%multiplicity(ch), &
                  "Symm. = ", eq%symmetry(ch), &
                  "Invar.: ", eq%div_is_invariant(ch,:)
          end if
       end do
    else
       write (u, "(3x,A)") "[not allocated]"
    end if
    write (u, "(1x,A)") "Equivalence list:"
    if (allocated (eq%eq)) then
       do i=1, size (eq%eq)
          call vamp_equivalence_write (eq%eq(i), u)
       end do
    else
       write (u, "(3x,A)") "[not allocated]"
    end if
  end subroutine vamp_equivalences_write
  subroutine vamp_equivalence_set (eq, i, left, right, perm, mode)
    type(vamp_equivalences_t), intent(inout) :: eq
    integer, intent(in) :: i
    integer, intent(in) :: left, right
    integer, dimension(:), intent(in) :: perm, mode
    eq%eq(i)%left = left
    eq%eq(i)%right = right
    eq%eq(i)%permutation = perm
    eq%eq(i)%mode = mode
  end subroutine vamp_equivalence_set
  subroutine vamp_equivalences_complete (eq)
    type(vamp_equivalences_t), intent(inout) :: eq
    integer :: i, ch
    ch = 0
    do i=1, eq%n_eq
       if (ch /= eq%eq(i)%left) then
          ch = eq%eq(i)%left
          eq%pointer(ch) = i
       end if
    end do
    eq%pointer(ch+1) = eq%n_eq + 1
    do ch=1, eq%n_ch
       call set_multiplicities (eq%eq(eq%pointer(ch):eq%pointer(ch+1)-1))
    end do
  ! call write (6, eq)
  contains
    subroutine set_multiplicities (eq_ch)
      type(vamp_equivalence_t), dimension(:), intent(in) :: eq_ch
      integer :: i
      if (.not. all(eq_ch%left == ch) .or. eq_ch(1)%right > ch) then
         do i = 1, size (eq_ch)
            call vamp_equivalence_write (eq_ch(i))
         end do
         stop "VAMP: Equivalences: Something's wrong with equivalence ordering"
      end if
      eq%symmetry(ch) = count (eq_ch%right == ch)
      if (mod (size(eq_ch), eq%symmetry(ch)) /= 0) then
         do i = 1, size (eq_ch)
            call vamp_equivalence_write (eq_ch(i))
         end do
         stop "VAMP: Equivalences: Something's wrong with permutation count"
      end if
      eq%multiplicity(ch) = size (eq_ch) / eq%symmetry(ch)
      eq%independent(ch) = all (eq_ch%right >= ch)
      eq%equivalent_to_ch(ch) = eq_ch(1)%right
      eq%div_is_invariant(ch,:) = eq_ch(1)%mode == VEQ_INVARIANT
    end subroutine set_multiplicities
  end subroutine vamp_equivalences_complete
end module vamp_equivalences
module vamp_rest
  use kinds
  use utils
  use exceptions
  use divisions
  use tao_random_numbers
  use vamp_stat
  use linalg
  use iso_fortran_env
  use vamp_grid_type !NODEP!
  use vamp_equivalences !NODEP!
  implicit none
  private
  public :: vamp_copy_grid, vamp_delete_grid
  public :: vamp_create_grid, vamp_create_empty_grid
  private :: map_domain
  public :: vamp_discard_integral
  private :: set_grid_options
  private :: vamp_reshape_grid_internal
  public :: vamp_reshape_grid
  public :: vamp_nullify_f_limits
  public :: vamp_rigid_divisions
  public :: vamp_get_covariance, vamp_nullify_covariance
  public :: vamp_get_variance, vamp_nullify_variance
  public :: vamp_sample_grid
  public :: vamp_sample_grid0
  public :: vamp_refine_grid
  public :: vamp_refine_grids
  public :: vamp_probability
  public :: vamp_average_iterations
  private :: vamp_average_iterations_grid
  public :: vamp_fork_grid
  private :: vamp_fork_grid_single, vamp_fork_grid_multi
  public :: vamp_join_grid
  private :: vamp_join_grid_single, vamp_join_grid_multi
  public :: vamp_fork_grid_joints
  public :: vamp_sample_grid_parallel
  public :: vamp_distribute_work
  public :: vamp_create_history, vamp_copy_history, vamp_delete_history
  public :: vamp_terminate_history
  public :: vamp_get_history, vamp_get_history_single
  public :: vamp_print_history, vamp_write_history
  private :: vamp_print_one_history, vamp_print_histories
  ! private :: vamp_write_one_history, vamp_write_histories
  public :: vamp_multi_channel, vamp_multi_channel0
  public :: vamp_jacobian, vamp_check_jacobian
  public :: vamp_create_grids, vamp_create_empty_grids
  public :: vamp_copy_grids, vamp_delete_grids
  public :: vamp_discard_integrals
  public :: vamp_update_weights
  public :: vamp_reshape_grids
  public :: vamp_sample_grids
  public :: vamp_reduce_channels
  public :: vamp_refine_weights
  private :: vamp_average_iterations_grids
  private :: vamp_get_history_multi
  public :: vamp_sum_channels
  public :: select_rotation_axis
  public :: select_rotation_subspace
  private :: more_pancake_than_cigar
  private :: select_subspace_explicit
  private :: select_subspace_guess
  public :: vamp_print_covariance
  ! private :: condense
  public :: condense
  ! private :: condense_action
  public :: condense_action
  public :: vamp_next_event
  private :: vamp_next_event_single, vamp_next_event_multi
  public :: vamp_warmup_grid, vamp_warmup_grids
  public :: vamp_integrate
  private :: vamp_integrate_grid, vamp_integrate_region
  public :: vamp_integratex
  private :: vamp_integratex_region
  public :: vamp_write_grid
  private :: write_grid_unit, write_grid_name
  public :: vamp_read_grid
  private :: read_grid_unit, read_grid_name
  public :: vamp_write_grids
  private :: write_grids_unit, write_grids_name
  public :: vamp_read_grids
  private :: read_grids_unit, read_grids_name
  public :: vamp_read_grids_raw
  private :: read_grids_raw_unit, read_grids_raw_name
  public :: vamp_read_grid_raw
  private :: read_grid_raw_unit, read_grid_raw_name
  public :: vamp_write_grids_raw
  private :: write_grids_raw_unit, write_grids_raw_name
  public :: vamp_write_grid_raw
  private :: write_grid_raw_unit, write_grid_raw_name
  public :: vamp_marshal_grid_size, vamp_marshal_grid, vamp_unmarshal_grid
  public :: vamp_marshal_history_size, vamp_marshal_history
  public :: vamp_unmarshal_history
  interface vamp_average_iterations
     module procedure vamp_average_iterations_grid
  end interface
  interface vamp_fork_grid
     module procedure vamp_fork_grid_single, vamp_fork_grid_multi
  end interface
  interface vamp_join_grid
     module procedure vamp_join_grid_single, vamp_join_grid_multi
  end interface
  interface vamp_get_history
     module procedure vamp_get_history_single
  end interface
  interface vamp_print_history
     module procedure vamp_print_one_history, vamp_print_histories
  end interface
  interface vamp_write_history
     module procedure vamp_write_one_history_unit, vamp_write_histories_unit
  end interface
  interface vamp_average_iterations
     module procedure  vamp_average_iterations_grids
  end interface
  interface vamp_get_history
     module procedure vamp_get_history_multi
  end interface
  interface select_rotation_subspace
     module procedure select_subspace_explicit, select_subspace_guess
  end interface
  interface vamp_next_event
     module procedure vamp_next_event_single, vamp_next_event_multi
  end interface
  interface vamp_integrate
     module procedure vamp_integrate_grid, vamp_integrate_region
  end interface
  interface vamp_integratex
     module procedure vamp_integratex_region
  end interface
  interface vamp_write_grid
     module procedure write_grid_unit, write_grid_name
  end interface
  interface vamp_read_grid
     module procedure read_grid_unit, read_grid_name
  end interface
  interface vamp_write_grids
     module procedure write_grids_unit, write_grids_name
  end interface
  interface vamp_read_grids
     module procedure read_grids_unit, read_grids_name
  end interface
  interface vamp_write_grid_raw
     module procedure write_grid_raw_unit, write_grid_raw_name
  end interface
  interface vamp_read_grid_raw
     module procedure read_grid_raw_unit, read_grid_raw_name
  end interface
  interface vamp_write_grids_raw
     module procedure write_grids_raw_unit, write_grids_raw_name
  end interface
  interface vamp_read_grids_raw
     module procedure read_grids_raw_unit, read_grids_raw_name
  end interface
  integer, parameter, private :: MAGIC_GRID = 22222222
  integer, parameter, private :: MAGIC_GRID_BEGIN = MAGIC_GRID + 1
  integer, parameter, private :: MAGIC_GRID_END = MAGIC_GRID + 2
  integer, parameter, private :: MAGIC_GRID_EMPTY = MAGIC_GRID + 3
  integer, parameter, private :: MAGIC_GRID_MAP = MAGIC_GRID + 4
  integer, parameter, private :: MAGIC_GRID_MU_X = MAGIC_GRID + 5
  integer, parameter, private :: MAGIC_GRIDS = 33333333
  integer, parameter, private :: MAGIC_GRIDS_BEGIN = MAGIC_GRIDS + 1
  integer, parameter, private :: MAGIC_GRIDS_END = MAGIC_GRIDS + 2
  type, public :: vamp_data_t
  end type vamp_data_t

  type(vamp_data_t), parameter, public :: NO_DATA = vamp_data_t ()

  type, public :: vamp_history
     private
     real(kind=default) :: &
          integral, std_dev, avg_integral, avg_std_dev, avg_chi2, f_min, f_max
     integer :: calls
     logical :: stratified
     logical :: verbose
     type(div_history), dimension(:), pointer :: div => null ()
  end type vamp_history
  type, public :: vamp_grids
     !!! private !: \emph{used by \texttt{vampi}}
     real(kind=default), dimension(:), pointer :: weights => null ()
     type(vamp_grid), dimension(:), pointer :: grids => null ()
     integer, dimension(:), pointer :: num_calls => null ()
     real(kind=default) :: sum_chi2, sum_integral, sum_weights
  end type vamp_grids
  integer, private, parameter :: NUM_DIV_DEFAULT = 20
  real(kind=default), private, parameter :: QUAD_POWER = 0.5_default
  integer, private, parameter :: BUFFER_SIZE = 50
  character(len=*), parameter, private :: &
       descr_fmt =         "(1x,a)", &
       integer_fmt =       "(1x,a17,1x,i15)", &
       integer_array_fmt = "(1x,i17,1x,i15)", &
       logical_fmt =       "(1x,a17,1x,l1)", &
       double_fmt =        "(1x,a17,1x,e30.22e4)", &
       double_array_fmt =  "(1x,i17,1x,e30.22e4)", &
       double_array2_fmt =  "(2(1x,i8),1x,e30.22e4)"
contains
  pure subroutine vamp_create_grid &
       (g, domain, num_calls, num_div, &
        stratified, quadrupole, covariance, map, exc)
    type(vamp_grid), intent(inout) :: g
    real(kind=default), dimension(:,:), intent(in) :: domain
    integer, intent(in) :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole, covariance
    real(kind=default), dimension(:,:), intent(in), optional :: map
    type(exception), intent(inout), optional :: exc
    character(len=*), parameter :: FN = "vamp_create_grid"
    real(kind=default), dimension(size(domain,dim=2)) :: &
         x_min, x_max, x_min_true, x_max_true
    integer :: ndim
    ndim = size (domain, dim=2)
    allocate (g%div(ndim), g%num_div(ndim))
    x_min = domain(1,:)
    x_max = domain(2,:)
    if (present (map)) then
       allocate (g%map(ndim,ndim))
       g%map = map
       x_min_true = x_min
       x_max_true = x_max
       call map_domain (g%map, x_min_true, x_max_true, x_min, x_max)
       call create_division (g%div, x_min, x_max, x_min_true, x_max_true)
    else
       nullify (g%map)
       call create_division (g%div, x_min, x_max)
    end if
    g%num_calls = num_calls
    if (present (num_div)) then
       g%num_div = num_div
    else
       g%num_div = NUM_DIV_DEFAULT
    end if
    g%stratified = .true.
    g%quadrupole = .false.
    g%independent = .true.
    g%equivalent_to_ch = 0
    g%multiplicity = 1
    nullify (g%mu_x, g%mu_xx, g%sum_mu_x, g%sum_mu_xx)
    call vamp_discard_integral &
         (g, num_calls, num_div, stratified, quadrupole, covariance, exc)
  end subroutine vamp_create_grid
  pure subroutine map_domain (map, true_xmin, true_xmax, xmin, xmax)
    real(kind=default), dimension(:,:), intent(in) :: map
    real(kind=default), dimension(:), intent(in) :: true_xmin, true_xmax
    real(kind=default), dimension(:), intent(out) :: xmin, xmax
    real(kind=default), dimension(2**size(xmin),size(xmin)) :: corners
    integer, dimension(size(xmin)) :: zero_to_n
    integer :: j, ndim, perm
    ndim = size (xmin)
    zero_to_n = (/ (j, j=0,ndim-1) /)
    do perm = 1, 2**ndim
       corners (perm,:) = &
            merge (true_xmin, true_xmax, btest (perm-1, zero_to_n))
    end do
    corners = matmul (corners, map)
    xmin = minval (corners, dim=1)
    xmax = maxval (corners, dim=1)
  end subroutine map_domain
  elemental subroutine vamp_create_empty_grid (g)
    type(vamp_grid), intent(inout) :: g
    nullify (g%div, g%num_div, g%map, g%mu_x, g%mu_xx, g%sum_mu_x, g%sum_mu_xx)
  end subroutine vamp_create_empty_grid
  pure subroutine vamp_discard_integral &
       (g, num_calls, num_div, stratified, quadrupole, covariance, exc, &
      & independent, equivalent_to_ch, multiplicity)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole, covariance
    type(exception), intent(inout), optional :: exc
    logical, intent(in), optional :: independent
    integer, intent(in), optional :: equivalent_to_ch, multiplicity
    character(len=*), parameter :: FN = "vamp_discard_integral"
    g%mu = 0.0
    g%mu_plus = 0.0
    g%mu_minus = 0.0
    g%mu_gi = 0.0
    g%sum_integral = 0.0
    g%sum_weights = 0.0
    g%sum_chi2 = 0.0
    g%sum_mu_gi = 0.0
    if (associated (g%sum_mu_x)) then
       g%sum_mu_x = 0.0
       g%sum_mu_xx = 0.0
    end if
    call set_grid_options (g, num_calls, num_div, stratified, quadrupole, &
                           independent, equivalent_to_ch, multiplicity)
    if ((present (num_calls)) &
        .or. (present (num_div)) &
        .or. (present (stratified)) &
        .or. (present (quadrupole)) &
        .or. (present (covariance))) then
       call vamp_reshape_grid &
            (g, g%num_calls, g%num_div, &
             g%stratified, g%quadrupole, covariance, exc)
    end if
  end subroutine vamp_discard_integral
  pure subroutine set_grid_options &
       (g, num_calls, num_div, stratified, quadrupole, &
        independent, equivalent_to_ch, multiplicity)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    logical, intent(in), optional :: independent
    integer, intent(in), optional :: equivalent_to_ch, multiplicity
    if (present (num_calls)) then
       g%num_calls = num_calls
    end if
    if (present (num_div)) then
       g%num_div = num_div
    end if
    if (present (stratified)) then
       g%stratified = stratified
    end if
    if (present (quadrupole)) then
       g%quadrupole = quadrupole
    end if
    if (present (independent)) then
       g%independent = independent
    end if
    if (present (equivalent_to_ch)) then
       g%equivalent_to_ch = equivalent_to_ch
    end if
    if (present (multiplicity)) then
       g%multiplicity = multiplicity
    end if
  end subroutine set_grid_options
  pure subroutine vamp_reshape_grid_internal &
       (g, num_calls, num_div, &
        stratified, quadrupole, covariance, exc, use_variance, &
        independent, equivalent_to_ch, multiplicity)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole, covariance
    type(exception), intent(inout), optional :: exc
    logical, intent(in), optional :: use_variance
    logical, intent(in), optional :: independent
    integer, intent(in), optional :: equivalent_to_ch, multiplicity
    integer :: ndim, num_cells
    integer, dimension(size(g%div)) :: ng
    character(len=*), parameter :: FN = "vamp_reshape_grid_internal"
    ndim = size (g%div)
    call set_grid_options &
         (g, num_calls, num_div, stratified, quadrupole, &
        & independent, equivalent_to_ch, multiplicity)
    if (g%stratified) then
       ng = (g%num_calls / 2.0 + 0.25)**(1.0/ndim)
    !  ng = ng * real (g%num_div, kind=default) &
    !          / (product (real (g%num_div, kind=default)))**(1.0/ndim)
    else
       ng = 1
    end if
    call reshape_division (g%div, g%num_div, ng, use_variance)
    call clear_integral_and_variance (g%div)
    num_cells = product (rigid_division (g%div))
    g%calls_per_cell = max (g%num_calls / num_cells, 2)
    g%calls = real (g%calls_per_cell) * real (num_cells)
    g%jacobi = product (volume_division (g%div)) / g%calls
    g%dv2g = (g%calls / num_cells)**2 &
         / g%calls_per_cell / g%calls_per_cell / (g%calls_per_cell - 1.0)
    call vamp_nullify_f_limits (g)
    g%all_stratified = all (stratified_division (g%div))
    if (present (covariance)) then
       ndim = size (g%div)
       if (covariance .and. (.not. associated (g%mu_x))) then
          allocate (g%mu_x(ndim), g%mu_xx(ndim,ndim))
          allocate (g%sum_mu_x(ndim), g%sum_mu_xx(ndim,ndim))
          g%sum_mu_x = 0.0
          g%sum_mu_xx = 0.0
       else if ((.not. covariance) .and. (associated (g%mu_x))) then
          deallocate (g%mu_x, g%mu_xx, g%sum_mu_x, g%sum_mu_xx)
       end if
    end if
  end subroutine vamp_reshape_grid_internal
  pure subroutine vamp_reshape_grid &
       (g, num_calls, num_div, stratified, quadrupole, covariance, exc, &
        independent, equivalent_to_ch, multiplicity)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole, covariance
    type(exception), intent(inout), optional :: exc
    logical, intent(in), optional :: independent
    integer, intent(in), optional :: equivalent_to_ch, multiplicity
    call vamp_reshape_grid_internal &
         (g, num_calls, num_div, stratified, quadrupole, covariance, &
          exc, use_variance = .false., &
          independent=independent, equivalent_to_ch=equivalent_to_ch, &
          multiplicity=multiplicity)
  end subroutine vamp_reshape_grid
  elemental subroutine vamp_nullify_f_limits (g)
    type(vamp_grid), intent(inout) :: g
    g%f_min = 1.0
    g%f_max = 0.0
  end subroutine vamp_nullify_f_limits
  pure function vamp_rigid_divisions (g) result (ng)
    type(vamp_grid), intent(in) :: g
    integer, dimension(size(g%div)) :: ng
    ng = rigid_division (g%div)
  end function vamp_rigid_divisions
  pure function vamp_get_covariance (g) result (cov)
    type(vamp_grid), intent(in) :: g
    real(kind=default), dimension(size(g%div),size(g%div)) :: cov
    if (associated (g%mu_x)) then
       if (abs (g%sum_weights) <= tiny (cov(1,1))) then
          where (g%sum_mu_xx == 0.0_default)
             cov = 0.0
          elsewhere
             cov = huge (cov(1,1))
          endwhere
       else
          cov = g%sum_mu_xx / g%sum_weights &
                  - outer_product (g%sum_mu_x, g%sum_mu_x) / g%sum_weights**2
       end if
    else
       cov = 0.0
    end if
  end function vamp_get_covariance
  elemental subroutine vamp_nullify_covariance (g)
    type(vamp_grid), intent(inout) :: g
    if (associated (g%mu_x)) then
       g%sum_mu_x = 0
       g%sum_mu_xx = 0
    end if
  end subroutine vamp_nullify_covariance
  elemental function vamp_get_variance (g) result (v)
    type(vamp_grid), intent(in) :: g
    real(kind=default) :: v
    if (abs (g%sum_weights) <= tiny (v)) then
       if (g%sum_mu_gi == 0.0_default) then
          v = 0.0
       else
          v = huge (v)
       end if
    else
       v = g%sum_mu_gi / g%sum_weights
    end if
  end function vamp_get_variance
  elemental subroutine vamp_nullify_variance (g)
    type(vamp_grid), intent(inout) :: g
    g%sum_mu_gi = 0
  end subroutine vamp_nullify_variance
    subroutine vamp_sample_grid0 &
         (rng, g, func, data, channel, weights, grids, exc, &
          negative_weights)
      type(tao_random_state), intent(inout) :: rng
      type(vamp_grid), intent(inout) :: g
      class(vamp_data_t), intent(in) :: data
      integer, intent(in), optional :: channel
      real(kind=default), dimension(:), intent(in), optional :: weights
      type(vamp_grid), dimension(:), intent(in), optional :: grids
      type(exception), intent(inout), optional :: exc
      interface
         function func (xi, data, weights, channel, grids) result (f)
           use kinds
           use vamp_grid_type !NODEP!
           import vamp_data_t
           real(kind=default), dimension(:), intent(in) :: xi
           class(vamp_data_t), intent(in) :: data
           real(kind=default), dimension(:), intent(in), optional :: weights
           integer, intent(in), optional :: channel
           type(vamp_grid), dimension(:), intent(in), optional :: grids
           real(kind=default) :: f
         end function func
      end interface
      character(len=*), parameter :: FN = "vamp_sample_grid0"
      logical, intent(in), optional :: negative_weights
      real(kind=default), parameter :: &
           eps =  tiny (1._default) / epsilon (1._default)
      character(len=6) :: buffer
      integer :: j, k
      integer, dimension(size(g%div)) :: cell
      real(kind=default) :: wgt, f, f2
      real(kind=default) :: sum_f, sum_f2, var_f
      real(kind=default) :: sum_f_plus, sum_f2_plus, var_f_plus
      real(kind=default) :: sum_f_minus, sum_f2_minus, var_f_minus
      real(kind=default), dimension(size(g%div)):: x, x_mid, wgts
      real(kind=default), dimension(size(g%div)):: r
      integer, dimension(size(g%div)) :: ia
      integer :: ndim
      logical :: neg_w
      ndim = size (g%div)
      neg_w = .false.
      if (present (negative_weights)) neg_w = negative_weights
      if (present (channel) .neqv. present (weights)) then
         call raise_exception (exc, EXC_FATAL, FN, &
              "channel and weights required together")
         return
      end if
      g%mu = 0.0
      g%mu_plus = 0.0
      g%mu_minus = 0.0
      cell = 1
      call clear_integral_and_variance (g%div)
      if (associated (g%mu_x)) then
         g%mu_x = 0.0
         g%mu_xx = 0.0
      end if
      if (present (channel)) then
         g%mu_gi = 0.0
      end if
      loop_over_cells: do
         sum_f = 0.0
         sum_f_plus = 0.0
         sum_f_minus = 0.0
         sum_f2 = 0.0
         sum_f2_plus = 0.0
         sum_f2_minus = 0.0
         do k = 1, g%calls_per_cell
            call tao_random_number (rng, r)
            call inject_division (g%div, real (r, kind=default), &
                                  cell, x, x_mid, ia, wgts)
            wgt = g%jacobi * product (wgts)
            if (associated (g%map)) then
               x = matmul (g%map, x)
            end if
            if (associated (g%map)) then
               if (all (inside_division (g%div, x))) then
                  f = wgt * func (x, data, weights, channel, grids)
               else
                  f = 0.0
               end if
            else
               f = wgt * func (x, data, weights, channel, grids)
            end if
            if (g%f_min > g%f_max) then
               g%f_min = abs (f) * g%calls
               g%f_max = abs (f) * g%calls
            else if (abs (f) * g%calls < g%f_min) then
               g%f_min = abs (f) * g%calls
            else if (abs (f) * g%calls > g%f_max) then
               g%f_max = abs (f) * g%calls
            end if
            f2 = f * f
            sum_f = sum_f + f
            sum_f2 = sum_f2 + f2
            if (f > 0) then
               sum_f_plus = sum_f_plus + f
               sum_f2_plus = sum_f2_plus + f * f
            else if (f < 0) then
               sum_f_minus = sum_f_minus + f
               sum_f2_minus = sum_f2_minus + f * f
            end if
            call record_integral (g%div, ia, f)
            ! call record_efficiency (g%div, ia, f/g%f_max)
            if ((associated (g%mu_x)) .and. (.not. g%all_stratified)) then
               g%mu_x = g%mu_x + x * f
               g%mu_xx = g%mu_xx + outer_product (x, x) * f
            end if
            if (present (channel)) then
               g%mu_gi = g%mu_gi + f2
            end if
         end do
         var_f = sum_f2 * g%calls_per_cell - sum_f**2
         var_f_plus = sum_f2_plus * g%calls_per_cell - sum_f_plus**2
         var_f_minus = sum_f2_minus * g%calls_per_cell - sum_f_minus**2
         if (var_f <= 0.0) then
            var_f = tiny (1.0_default)
         end if
         if (sum_f_plus /= 0 .and. var_f_plus <= 0) then
            var_f_plus = tiny (1.0_default)
         end if
         if (sum_f_minus /= 0 .and. var_f_minus <= 0) then
            var_f_minus = tiny (1.0_default)
         end if
         g%mu = g%mu + (/ sum_f, var_f /)
         g%mu_plus = g%mu_plus + (/ sum_f_plus, var_f_plus /)
         g%mu_minus = g%mu_minus + (/ sum_f_minus, var_f_minus /)
         call record_variance (g%div, ia, var_f)
         if ((associated (g%mu_x)) .and. g%all_stratified) then
            if (associated (g%map)) then
               x_mid = matmul (g%map, x_mid)
            end if
            g%mu_x = g%mu_x + x_mid * var_f
            g%mu_xx = g%mu_xx + outer_product (x_mid, x_mid) * var_f
         end if
         do j = ndim, 1, -1
            cell(j) = modulo (cell(j), rigid_division (g%div(j))) + 1
            if (cell(j) /= 1) then
               cycle loop_over_cells
            end if
         end do
         exit loop_over_cells
      end do loop_over_cells
      g%mu(2) = g%mu(2) * g%dv2g
      if (g%mu(2) < eps * max (g%mu(1)**2, 1._default)) then
         g%mu(2) = eps * max (g%mu(1)**2, 1._default)
      end if
      if (neg_w) then
         g%mu_plus(2) = g%mu_plus(2) * g%dv2g
         if (g%mu_plus(2) < eps * max (g%mu_plus(1)**2, 1._default)) then
            g%mu_plus(2) = eps * max (g%mu_plus(1)**2, 1._default)
         end if
         g%mu_minus(2) = g%mu_minus(2) * g%dv2g
         if (g%mu_minus(2) < eps * max (g%mu_minus(1)**2, 1._default)) then
            g%mu_minus(2) = eps * max (g%mu_minus(1)**2, 1._default)
         end if
      end if
      if (g%mu(1)>0) then
         g%sum_integral = g%sum_integral + g%mu(1) / g%mu(2)
         g%sum_weights = g%sum_weights + 1.0 / g%mu(2)
         g%sum_chi2 = g%sum_chi2 + g%mu(1)**2 / g%mu(2)
         if (associated (g%mu_x)) then
            if (g%all_stratified) then
               g%mu_x = g%mu_x / g%mu(2)
               g%mu_xx = g%mu_xx / g%mu(2)
            else
               g%mu_x = g%mu_x / g%mu(1)
               g%mu_xx = g%mu_xx / g%mu(1)
            end if
            g%sum_mu_x = g%sum_mu_x + g%mu_x / g%mu(2)
            g%sum_mu_xx = g%sum_mu_xx + g%mu_xx / g%mu(2)
         end if
         if (present (channel)) then
            g%sum_mu_gi = g%sum_mu_gi + g%mu_gi / g%mu(2)
         end if
      else if (neg_w) then
         g%sum_integral = g%sum_integral + g%mu(1) / g%mu(2)
         g%sum_weights = g%sum_weights + 1.0 / g%mu(2)
         g%sum_chi2 = g%sum_chi2 + g%mu(1)**2 / g%mu(2)
         if (associated (g%mu_x)) then
            if (g%all_stratified) then
               g%mu_x = g%mu_x / g%mu(2)
               g%mu_xx = g%mu_xx / g%mu(2)
            else
               g%mu_x = g%mu_x / g%mu(1)
               g%mu_xx = g%mu_xx / g%mu(1)
            end if
            g%sum_mu_x = g%sum_mu_x + g%mu_x / g%mu(2)
            g%sum_mu_xx = g%sum_mu_xx + g%mu_xx / g%mu(2)
         end if
         if (present (channel)) then
            g%sum_mu_gi = g%sum_mu_gi + g%mu_gi / g%mu(2)
         end if
         else
         if (present(channel) .and. g%mu(1)==0) then
            write (buffer, "(I6)")  channel
            call raise_exception (exc, EXC_WARN, "! vamp", &
                 "Function identically zero in channel " // buffer)
         else if (present(channel) .and. g%mu(1)<0) then
            write (buffer, "(I6)")  channel
            call raise_exception (exc, EXC_ERROR, "! vamp", &
                 "Negative integral in channel " // buffer)
         end if
         g%sum_integral = 0
         g%sum_chi2 = 0
         g%sum_weights = 0
      end if
    end subroutine vamp_sample_grid0

  pure function vamp_probability (g, x) result (p)
    type(vamp_grid), intent(in) :: g
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default) :: p
    p = product (probability (g%div, x))
  end function vamp_probability
  subroutine vamp_apply_equivalences (g, eq)
    type(vamp_grids), intent(inout) :: g
    type(vamp_equivalences_t), intent(in) :: eq
    integer :: n_ch, n_dim, nb, i, ch, ch_src, dim, dim_src
    integer, dimension(:,:), allocatable :: n_bin
    real(kind=default), dimension(:,:,:), allocatable :: var_tmp
    n_ch = size (g%grids)
    if (n_ch == 0)  return
    n_dim = size (g%grids(1)%div)
    allocate (n_bin(n_ch, n_dim))
    do ch = 1, n_ch
       do dim = 1, n_dim
          n_bin(ch, dim) = size (g%grids(ch)%div(dim)%variance)
       end do
    end do
    allocate (var_tmp (maxval(n_bin), n_dim, n_ch))
    var_tmp = 0
    do i=1, eq%n_eq
       ch = eq%eq(i)%left
       ch_src = eq%eq(i)%right
       do dim=1, n_dim
          nb = n_bin(ch_src, dim)
          dim_src = eq%eq(i)%permutation(dim)
          select case (eq%eq(i)%mode(dim))
          case (VEQ_IDENTITY)
             var_tmp(:nb,dim,ch) = var_tmp(:nb,dim,ch) &
                  & + g%grids(ch_src)%div(dim_src)%variance
          case (VEQ_INVERT)
             var_tmp(:nb,dim,ch) = var_tmp(:nb,dim,ch) &
                  & + g%grids(ch_src)%div(dim_src)%variance(nb:1:-1)
          case (VEQ_SYMMETRIC)
             var_tmp(:nb,dim,ch) = var_tmp(:nb,dim,ch) &
                  & + g%grids(ch_src)%div(dim_src)%variance / 2 &
                  & + g%grids(ch_src)%div(dim_src)%variance(nb:1:-1)/2
          case (VEQ_INVARIANT)
             var_tmp(:nb,dim,ch) = 1
          end select
       end do
    end do
    do ch=1, n_ch
       do dim=1, n_dim
          g%grids(ch)%div(dim)%variance = var_tmp(:n_bin(ch, dim),dim,ch)
       end do
    end do
    deallocate (var_tmp)
    deallocate (n_bin)
  end subroutine vamp_apply_equivalences
  pure subroutine vamp_refine_grid (g, exc)
    type(vamp_grid), intent(inout) :: g
    type(exception), intent(inout), optional :: exc
    real(kind=default), dimension(size(g%div)) :: quad
    integer :: ndim
    if (g%quadrupole) then
       ndim = size (g%div)
       quad = (quadrupole_division (g%div))**QUAD_POWER
       call vamp_reshape_grid_internal &
            (g, use_variance = .true., exc = exc, &
             num_div = int (quad / product (quad)**(1.0/ndim) * g%num_div))
    else
       call refine_division (g%div)
       call vamp_nullify_f_limits (g)
    end if
  end subroutine vamp_refine_grid
  subroutine vamp_refine_grids (g)
    type(vamp_grids), intent(inout) :: g
    integer :: ch
    do ch=1, size(g%grids)
       call refine_division (g%grids(ch)%div)
       call vamp_nullify_f_limits (g%grids(ch))
    end do
  end subroutine vamp_refine_grids
  subroutine vamp_sample_grid &
       (rng, g, func, data, iterations, &
        integral, std_dev, avg_chi2, accuracy, &
        channel, weights, grids, exc, history)
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grid), intent(inout) :: g
    class(vamp_data_t), intent(in) :: data
    integer, intent(in) :: iterations
    real(kind=default), intent(out), optional :: integral, std_dev, avg_chi2
    real(kind=default), intent(in), optional :: accuracy
    integer, intent(in), optional :: channel
    real(kind=default), dimension(:), intent(in), optional :: weights
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    type(exception), intent(inout), optional :: exc
    type(vamp_history), dimension(:), intent(inout), optional :: history
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    character(len=*), parameter :: FN = "vamp_sample_grid"
    real(kind=default) :: local_integral, local_std_dev, local_avg_chi2
    integer :: iteration, ndim
    ndim = size (g%div)
    iterate: do iteration = 1, iterations
       call vamp_sample_grid0 &
            (rng, g, func, data, channel, weights, grids, exc)
       call vamp_average_iterations &
            (g, iteration, local_integral, local_std_dev, local_avg_chi2)
       if (present (history)) then
          if (iteration <= size (history)) then
             call vamp_get_history &
                  (history(iteration), g, local_integral, local_std_dev, &
                   local_avg_chi2)
          else
             call raise_exception (exc, EXC_WARN, FN, "history too short")
          end if
          call vamp_terminate_history (history(iteration+1:))
       end if
       if (present (accuracy)) then
          if (local_std_dev <= accuracy * local_integral) then
             call raise_exception (exc, EXC_INFO, FN, &
                  "requested accuracy reached")
             exit iterate
          end if
       end if
       if (iteration < iterations) call vamp_refine_grid (g)
    end do iterate
    if (present (integral)) then
       integral = local_integral
    end if
    if (present (std_dev)) then
       std_dev = local_std_dev
    end if
    if (present (avg_chi2)) then
       avg_chi2 = local_avg_chi2
    end if
  end subroutine vamp_sample_grid
  elemental subroutine vamp_average_iterations_grid &
       (g, iteration, integral, std_dev, avg_chi2)
    type(vamp_grid), intent(in) :: g
    integer, intent(in) :: iteration
    real(kind=default), intent(out) :: integral, std_dev, avg_chi2
    real(kind=default), parameter :: eps = 1000 * epsilon (1._default)
    if (g%sum_weights>0) then
       integral = g%sum_integral / g%sum_weights
       std_dev = sqrt (1.0 / g%sum_weights)
       avg_chi2 = &
            max ((g%sum_chi2 - g%sum_integral * integral) / (iteration-0.99), &
                 0.0_default)
       if (avg_chi2 < eps * g%sum_chi2)  avg_chi2 = 0
    else
       integral = 0
       std_dev = 0
       avg_chi2 = 0
    end if
  end subroutine vamp_average_iterations_grid
  pure subroutine vamp_fork_grid_single (g, gs, d, exc)
    type(vamp_grid), intent(in) :: g
    type(vamp_grid), dimension(:), intent(inout) :: gs
    integer, intent(in) :: d
    type(exception), intent(inout), optional :: exc
    character(len=*), parameter :: FN = "vamp_fork_grid_single"
    type(division_t), dimension(:), allocatable :: d_tmp
    integer :: i, j, num_grids, num_div, ndim, num_cells
    num_grids = size (gs)
    ndim = size (g%div)
    num_div = size (g%div)
    do i = 1, size (gs)
       if (associated (gs(i)%div)) then
          if (size (gs(i)%div) /= num_div) then
             allocate (gs(i)%div(num_div))
             call create_empty_division (gs(i)%div)
          end if
       else
          allocate (gs(i)%div(num_div))
          call create_empty_division (gs(i)%div)
       end if
    end do
    do j = 1, ndim
       if (j == d) then
          allocate (d_tmp(num_grids))
          do i = 1, num_grids
             d_tmp(i) = gs(i)%div(j)
          end do
          call fork_division (g%div(j), d_tmp, g%calls_per_cell, gs%calls_per_cell, exc)
          do i = 1, num_grids
             gs(i)%div(j) = d_tmp(i)
          end do
          deallocate (d_tmp)
          if (present (exc)) then
             if (exc%level > EXC_WARN) then
                return
             end if
          end if
       else
          do i = 1, num_grids
             call copy_division (gs(i)%div(j), g%div(j))
          end do
       end if
    end do
    if (d == 0) then
       if (any (stratified_division (g%div))) then
          call raise_exception (exc, EXC_FATAL, FN, &
                                "d == 0 incompatiple w/ stratification")
       else
          gs(2:)%calls_per_cell = ceiling (real (g%calls_per_cell) / num_grids)
          gs(1)%calls_per_cell = g%calls_per_cell - sum (gs(2:)%calls_per_cell)
       end if
    end if
    do i = 1, num_grids
       call copy_array_pointer (gs(i)%num_div, g%num_div)
       if (associated (g%map)) then
          call copy_array_pointer (gs(i)%map, g%map)
       end if
       if (associated (g%mu_x)) then
          call create_array_pointer (gs(i)%mu_x, ndim)
          call create_array_pointer (gs(i)%sum_mu_x, ndim)
          call create_array_pointer (gs(i)%mu_xx, (/ ndim, ndim /))
          call create_array_pointer (gs(i)%sum_mu_xx, (/ ndim, ndim /))
       end if
    end do
    gs%mu(1) = 0.0
    gs%mu(2) = 0.0
    gs%mu_plus(1) = 0.0
    gs%mu_plus(2) = 0.0
    gs%mu_minus(1) = 0.0
    gs%mu_minus(2) = 0.0
    gs%sum_integral = 0.0
    gs%sum_weights = 0.0
    gs%sum_chi2 = 0.0
    gs%mu_gi = 0.0
    gs%sum_mu_gi = 0.0
    gs%stratified = g%stratified
    gs%all_stratified = g%all_stratified
    gs%quadrupole = g%quadrupole
    do i = 1, num_grids
       num_cells = product (rigid_division (gs(i)%div))
       gs(i)%calls = gs(i)%calls_per_cell * num_cells
       gs(i)%num_calls = gs(i)%calls
       gs(i)%jacobi = product (volume_division (gs(i)%div)) / gs(i)%calls
       gs(i)%dv2g = (gs(i)%calls / num_cells)**2 &
            / gs(i)%calls_per_cell / gs(i)%calls_per_cell / (gs(i)%calls_per_cell - 1.0)
    end do
    gs%f_min = g%f_min * (gs%jacobi * gs%calls) / (g%jacobi * g%calls)
    gs%f_max = g%f_max * (gs%jacobi * gs%calls) / (g%jacobi * g%calls)
  end subroutine vamp_fork_grid_single
  pure subroutine vamp_join_grid_single (g, gs, d, exc)
    type(vamp_grid), intent(inout) :: g
    type(vamp_grid), dimension(:), intent(inout) :: gs
    integer, intent(in) :: d
    type(exception), intent(inout), optional :: exc
    type(division_t), dimension(:), allocatable :: d_tmp
    integer :: i, j, num_grids
    num_grids = size (gs)
    do j = 1, size (g%div)
       if (j == d) then
          allocate (d_tmp(num_grids))
          do i = 1, num_grids
             d_tmp(i) = gs(i)%div(j)
          end do
          call join_division (g%div(j), d_tmp, exc)
          deallocate (d_tmp)
          if (present (exc)) then
             if (exc%level > EXC_WARN) then
                return
             end if
          end if
       else
          allocate (d_tmp(num_grids))
          do i = 1, num_grids
             d_tmp(i) = gs(i)%div(j)
          end do
          call sum_division (g%div(j), d_tmp)
          deallocate (d_tmp)
       end if
    end do
    g%f_min = minval (gs%f_min * (g%jacobi * g%calls) / (gs%jacobi * gs%calls))
    g%f_max = maxval (gs%f_max * (g%jacobi * g%calls) / (gs%jacobi * gs%calls))
    g%mu(1) = sum (gs%mu(1))
    g%mu(2) = sum (gs%mu(2))
    g%mu_plus(1) = sum (gs%mu_plus(1))
    g%mu_plus(2) = sum (gs%mu_plus(2))
    g%mu_minus(1) = sum (gs%mu_minus(1))
    g%mu_minus(2) = sum (gs%mu_minus(2))
    g%mu_gi = sum (gs%mu_gi)
    g%sum_mu_gi = g%sum_mu_gi + g%mu_gi / g%mu(2)
    g%sum_integral = g%sum_integral + g%mu(1) / g%mu(2)
    g%sum_chi2 = g%sum_chi2 + g%mu(1)**2 / g%mu(2)
    g%sum_weights = g%sum_weights + 1.0 / g%mu(2)
    if (associated (g%mu_x)) then
       do i = 1, num_grids
          g%mu_x = g%mu_x + gs(i)%mu_x
          g%mu_xx = g%mu_xx + gs(i)%mu_xx
       end do
       g%sum_mu_x = g%sum_mu_x + g%mu_x / g%mu(2)
       g%sum_mu_xx = g%sum_mu_xx + g%mu_xx / g%mu(2)
    end if
  end subroutine vamp_join_grid_single
  pure recursive subroutine vamp_fork_grid_multi (g, gs, gx, d, exc)
    type(vamp_grid), intent(in) :: g
    type(vamp_grid), dimension(:), intent(inout) :: gs, gx
    integer, dimension(:,:), intent(in) :: d
    type(exception), intent(inout), optional :: exc
    character(len=*), parameter :: FN = "vamp_fork_grid_multi"
    integer :: i, offset, stride, joints_offset, joints_stride
    select case (size (d, dim=2))
       case (0)
          return
       case (1)
          call vamp_fork_grid_single (g, gs, d(1,1), exc)
       case default
          offset = 1
          stride = product (d(2,2:))
          joints_offset = 1 + d(2,1)
          joints_stride = vamp_fork_grid_joints (d(:,2:))
          call vamp_create_empty_grid (gx(1:d(2,1)))
          call vamp_fork_grid_single (g, gx(1:d(2,1)), d(1,1), exc)
          do i = 1, d(2,1)
             call vamp_fork_grid_multi &
                  (gx(i), gs(offset:offset+stride-1), &
                   gx(joints_offset:joints_offset+joints_stride-1), &
                   d(:,2:), exc)
             offset = offset + stride
             joints_offset = joints_offset + joints_stride
          end do
    end select
  end subroutine vamp_fork_grid_multi
  pure function vamp_fork_grid_joints (d) result (s)
    integer, dimension(:,:), intent(in) :: d
    integer :: s
    integer :: i
    s = 0
    do i = size (d, dim=2) - 1, 1, -1
       s = (s + 1) * d(2,i)
    end do
  end function vamp_fork_grid_joints
  pure recursive subroutine vamp_join_grid_multi (g, gs, gx, d, exc)
    type(vamp_grid), intent(inout) :: g
    type(vamp_grid), dimension(:), intent(inout) :: gs, gx
    integer, dimension(:,:), intent(in) :: d
    type(exception), intent(inout), optional :: exc
    character(len=*), parameter :: FN = "vamp_join_grid_multi"
    integer :: i, offset, stride, joints_offset, joints_stride
    select case (size (d, dim=2))
       case (0)
          return
       case (1)
          call vamp_join_grid_single (g, gs, d(1,1), exc)
       case default
          offset = 1
          stride = product (d(2,2:))
          joints_offset = 1 + d(2,1)
          joints_stride = vamp_fork_grid_joints (d(:,2:))
          do i = 1, d(2,1)
             call vamp_join_grid_multi &
                  (gx(i), gs(offset:offset+stride-1), &
                   gx(joints_offset:joints_offset+joints_stride-1), &
                   d(:,2:), exc)
             offset = offset + stride
             joints_offset = joints_offset + joints_stride
          end do
          call vamp_join_grid_single (g, gx(1:d(2,1)), d(1,1), exc)
          call vamp_delete_grid (gx(1:d(2,1)))
    end select
  end subroutine vamp_join_grid_multi
  subroutine vamp_sample_grid_parallel &
       (rng, g, func, data, iterations, &
        integral, std_dev, avg_chi2, accuracy, &
        channel, weights, grids, exc, history)
    type(tao_random_state), dimension(:), intent(inout) :: rng
    type(vamp_grid), intent(inout) :: g
    class(vamp_data_t), intent(in) :: data
    integer, intent(in) :: iterations
    real(kind=default), intent(out), optional :: integral, std_dev, avg_chi2
    real(kind=default), intent(in), optional :: accuracy
    integer, intent(in), optional :: channel
    real(kind=default), dimension(:), intent(in), optional :: weights
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    type(exception), intent(inout), optional :: exc
    type(vamp_history), dimension(:), intent(inout), optional :: history
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    character(len=*), parameter :: FN = "vamp_sample_grid_parallel"
    real(kind=default) :: local_integral, local_std_dev, local_avg_chi2
    type(exception), dimension(size(rng)) :: excs
    type(vamp_grid), dimension(:), allocatable :: gs, gx
    !hpf$ processors p(number_of_processors())
    !hpf$ distribute gs(cyclic(1)) onto p
    integer, dimension(:,:), pointer :: d
    integer :: iteration, i
    integer :: num_workers
    nullify (d)
    call clear_exception (excs)
    iterate: do iteration = 1, iterations
       call vamp_distribute_work (size (rng), vamp_rigid_divisions (g), d)
       num_workers = max (1, product (d(2,:)))
       if (num_workers > 1) then
          allocate (gs(num_workers), gx(vamp_fork_grid_joints (d)))
          call vamp_create_empty_grid (gs)
          !: \texttt{vamp\_fork\_grid} is certainly not local.  Speed freaks might
          !: want to tune it to the processor topology, but the gain will be small.
          call vamp_fork_grid (g, gs, gx, d, exc)
          !hpf$ independent
          do i = 1, num_workers
             call vamp_sample_grid0 &
                  (rng(i), gs(i), func, data, &
                   channel, weights, grids, exc)
          end do
          if ((present (exc)) .and. (any (excs(1:num_workers)%level > 0))) then
             call gather_exceptions (exc, excs(1:num_workers))
          end if
          call vamp_join_grid (g, gs, gx, d, exc)
          call vamp_delete_grid (gs)
          deallocate (gs, gx)
       else
          call vamp_sample_grid0 &
               (rng(1), g, func, data, channel, weights, grids, exc)
       end if
       if (present (exc)) then
          if (exc%level > EXC_WARN) then
             return
          end if
       end if
       call vamp_average_iterations &
            (g, iteration, local_integral, local_std_dev, local_avg_chi2)
       if (present (history)) then
          if (iteration <= size (history)) then
             call vamp_get_history &
                  (history(iteration), g, local_integral, local_std_dev, &
                   local_avg_chi2)
          else
             call raise_exception (exc, EXC_WARN, FN, "history too short")
          end if
          call vamp_terminate_history (history(iteration+1:))
       end if
       if (present (accuracy)) then
          if (local_std_dev <= accuracy * local_integral) then
             call raise_exception (exc, EXC_INFO, FN, &
                  "requested accuracy reached")
             exit iterate
          end if
       end if
       if (iteration < iterations) call vamp_refine_grid (g)
    end do iterate
    deallocate (d)
    if (present (integral)) then
       integral = local_integral
    end if
    if (present (std_dev)) then
       std_dev = local_std_dev
    end if
    if (present (avg_chi2)) then
       avg_chi2 = local_avg_chi2
    end if
  end subroutine vamp_sample_grid_parallel
  pure subroutine vamp_distribute_work (num_workers, ng, d)
    integer, intent(in) :: num_workers
    integer, dimension(:), intent(in) :: ng
    integer, dimension(:,:), pointer :: d
    integer, dimension(32) :: factors
    integer :: n, num_factors, i, j
    integer, dimension(size(ng)) :: num_forks
    integer :: nfork
    try: do n = num_workers, 1, -1
       call factorize (n, factors, num_factors)
       num_forks = 1
       do i = num_factors, 1, -1
          j = sum (maxloc (ng / num_forks))
          nfork = num_forks(j) * factors(i)
          if (nfork <= ng(j)) then
             num_forks(j) = nfork
          else
             cycle try
          end if
       end do
       j = count (num_forks > 1)
       if (associated (d)) then
          if (size (d, dim = 2) /= j) then
             deallocate (d)
             allocate (d(2,j))
          end if
       else
          allocate (d(2,j))
       end if
       j = 1
       do i = 1, size (ng)
          if (num_forks(i) > 1) then
             d(:,j) = (/ i, num_forks(i) /)
             j = j + 1
          end if
       end do
       return
    end do try
  end subroutine vamp_distribute_work
  elemental subroutine vamp_create_history (h, ndim, verbose)
    type(vamp_history), intent(out) :: h
    integer, intent(in), optional :: ndim
    logical, intent(in), optional :: verbose
    if (present (verbose)) then
       h%verbose = verbose
    else
       h%verbose = .false.
    end if
    h%calls = 0.0
    if (h%verbose .and. (present (ndim))) then
       if (associated (h%div)) then
          deallocate (h%div)
       end if
       allocate (h%div(ndim))
    end if
  end subroutine vamp_create_history
  elemental subroutine vamp_terminate_history (h)
    type(vamp_history), intent(inout) :: h
    h%calls = 0.0
  end subroutine vamp_terminate_history
  pure subroutine vamp_get_history_single (h, g, integral, std_dev, avg_chi2)
    type(vamp_history), intent(inout) :: h
    type(vamp_grid), intent(in) :: g
    real(kind=default), intent(in) :: integral, std_dev, avg_chi2
    h%calls = g%calls
    h%stratified = g%all_stratified
    h%integral = g%mu(1)
    h%std_dev = sqrt (g%mu(2))
    h%avg_integral = integral
    h%avg_std_dev = std_dev
    h%avg_chi2 = avg_chi2
    h%f_min = g%f_min
    h%f_max = g%f_max
    if (h%verbose) then
       if (associated (h%div)) then
          if (size (h%div) /= size (g%div)) then
             deallocate (h%div)
             allocate (h%div(size(g%div)))
          end if
       else
          allocate (h%div(size(g%div)))
       end if
       call copy_history (h%div, summarize_division (g%div))
    end if
  end subroutine vamp_get_history_single
  subroutine vamp_print_one_history (h, tag)
    type(vamp_history), dimension(:), intent(in) :: h
    character(len=*), intent(in), optional :: tag
    type(div_history), dimension(:), allocatable :: h_tmp
    character(len=BUFFER_SIZE) :: pfx
    character(len=1) :: s
    integer :: i, imax, j
    if (present (tag)) then
       pfx = tag
    else
       pfx = "[vamp]"
    end if
    print "(1X,A78)", repeat ("-", 78)
    print "(1X,A8,1X,A2,A9,A1,1X,A11,1X,8X,1X," &
                          // "1X,A13,1X,8X,1X,A5,1X,A5)", &
         pfx, "it", "#calls", "", "integral", "average", "chi2", "eff."
    imax = size (h)
    iterations: do i = 1, imax
       if (h(i)%calls <= 0) then
          imax = i - 1
          exit iterations
       end if
       ! *JR: Skip zero channel
       if (h(i)%f_max==0) cycle
       if (h(i)%stratified) then
          s = "*"
       else
          s = ""
       end if
       print "(1X,A8,1X,I2,I9,A1,1X,E11.4,A1,E8.2,A1," &
                             // "1X,E13.6,A1,E8.2,A1,F5.1,1X,F5.3)", pfx, &
            i, h(i)%calls, s, h(i)%integral, "(", h(i)%std_dev, ")", &
            h(i)%avg_integral, "(", h(i)%avg_std_dev, ")", h(i)%avg_chi2, &
            h(i)%integral / h(i)%f_max
    end do iterations
    print "(1X,A78)", repeat ("-", 78)
    if (all (h%verbose) .and. (imax >= 1)) then
       if (associated (h(1)%div)) then
          allocate (h_tmp(imax))
          dimensions: do j = 1, size (h(1)%div)
             do i = 1, imax
                call copy_history (h_tmp(i), h(i)%div(j))
             end do
             if (present (tag)) then
                write (unit = pfx, fmt = "(A,A1,I2.2)") &
                     trim (tag(1:min(len_trim(tag),8))), "#", j
             else
                write (unit = pfx, fmt = "(A,A1,I2.2)") "[vamp]", "#", j
             end if
             call print_history (h_tmp, tag = pfx)
             print "(1X,A78)", repeat ("-", 78)
          end do dimensions
          deallocate (h_tmp)
       end if
    end if
    flush (output_unit)
  end subroutine vamp_print_one_history
  subroutine vamp_print_histories (h, tag)
    type(vamp_history), dimension(:,:), intent(in) :: h
    character(len=*), intent(in), optional :: tag
    character(len=BUFFER_SIZE) :: pfx
    integer :: i
    print "(1X,A78)", repeat ("=", 78)
    channels: do i = 1, size (h, dim=2)
       if (present (tag)) then
          write (unit = pfx, fmt = "(A4,A1,I3.3)") tag, "#", i
       else
          write (unit = pfx, fmt = "(A4,A1,I3.3)") "chan", "#", i
       end if
       call vamp_print_one_history (h(:,i), pfx)
    end do channels
    print "(1X,A78)", repeat ("=", 78)
    flush (output_unit)
  end subroutine vamp_print_histories
  subroutine vamp_write_one_history_unit (u, h, tag)
    integer, intent(in) :: u
    type(vamp_history), dimension(:), intent(in) :: h
    character(len=*), intent(in), optional :: tag
    type(div_history), dimension(:), allocatable :: h_tmp
    character(len=BUFFER_SIZE) :: pfx
    character(len=1) :: s
    integer :: i, imax, j
    if (present (tag)) then
       pfx = tag
    else
       pfx = "[vamp]"
    end if
    write (u, "(1X,A78)") repeat ("-", 78)
    write (u, "(1X,A8,1X,A2,A9,A1,1X,A11,1X,8X,1X," &
         // "1X,A13,1X,8X,1X,A5,1X,A5)") &
         pfx, "it", "#calls", "", "integral", "average", "chi2", "eff."
    imax = size (h)
    iterations: do i = 1, imax
       if (h(i)%calls <= 0) then
          imax = i - 1
          exit iterations
       end if
       ! *WK: Skip zero channel
       if (h(i)%f_max==0) cycle
       if (h(i)%stratified) then
          s = "*"
       else
          s = ""
       end if
       write (u, "(1X,A8,1X,I2,I9,A1,1X,ES11.4,A1,ES8.2,A1," &
            // "1X,ES13.6,A1,ES8.2,A1,F5.1,1X,F5.3)") pfx, &
            i, h(i)%calls, s, h(i)%integral, "(", h(i)%std_dev, ")", &
            h(i)%avg_integral, "(", h(i)%avg_std_dev, ")", h(i)%avg_chi2, &
            h(i)%integral / h(i)%f_max
    end do iterations
    write (u, "(1X,A78)") repeat ("-", 78)
    if (all (h%verbose) .and. (imax >= 1)) then
       if (associated (h(1)%div)) then
          allocate (h_tmp(imax))
          dimensions: do j = 1, size (h(1)%div)
             do i = 1, imax
                call copy_history (h_tmp(i), h(i)%div(j))
             end do
             if (present (tag)) then
                write (unit = pfx, fmt = "(A,A1,I2.2)") &
                     trim (tag(1:min(len_trim(tag),8))), "#", j
             else
                write (unit = pfx, fmt = "(A,A1,I2.2)") "[vamp]", "#", j
             end if
             call write_history (u, h_tmp, tag = pfx)
             print "(1X,A78)", repeat ("-", 78)
          end do dimensions
          deallocate (h_tmp)
       end if
    end if
    flush (u)
  end subroutine vamp_write_one_history_unit
  subroutine vamp_write_histories_unit (u, h, tag)
    integer, intent(in) :: u
    type(vamp_history), dimension(:,:), intent(in) :: h
    character(len=*), intent(in), optional :: tag
    character(len=BUFFER_SIZE) :: pfx
    integer :: i
    write (u, "(1X,A78)") repeat ("=", 78)
    channels: do i = 1, size (h, dim=2)
       if (present (tag)) then
          write (unit = pfx, fmt = "(A4,A1,I3.3)") tag, "#", i
       else
          write (unit = pfx, fmt = "(A4,A1,I3.3)") "chan", "#", i
       end if
       call vamp_write_one_history_unit (u, h(:,i), pfx)
    end do channels
    write (u, "(1X,A78)") repeat ("=", 78)
    flush (u)
  end subroutine vamp_write_histories_unit
  function vamp_multi_channel &
       (func, data, phi, ihp, jacobian, x, weights, channel, grids) result (w_x)
    class(vamp_data_t), intent(in) :: data
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), dimension(:), intent(in) :: weights
    integer, intent(in) :: channel
    type(vamp_grid), dimension(:), intent(in) :: grids
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    interface
       pure function phi (xi, channel) result (x)
         use kinds
         real(kind=default), dimension(:), intent(in) :: xi
         integer, intent(in) :: channel
         real(kind=default), dimension(size(xi)) :: x
       end function phi
    end interface
    interface
       pure function ihp (x, channel) result (xi)
         use kinds
         real(kind=default), dimension(:), intent(in) :: x
         integer, intent(in) :: channel
         real(kind=default), dimension(size(x)) :: xi
       end function ihp
    end interface
    interface
       pure function jacobian (x, data, channel) result (j)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: x
         class(vamp_data_t), intent(in) :: data
         integer, intent(in) :: channel
         real(kind=default) :: j
       end function jacobian
    end interface
    real(kind=default) :: w_x
    integer :: i
    real(kind=default), dimension(size(x)) :: phi_x
    real(kind=default), dimension(size(weights)) :: g_phi_x, g_pi_x
    phi_x = phi (x, channel)
    do i = 1, size (weights)
       if (i == channel) then
          g_pi_x(i) = vamp_probability (grids(i), x)
       else
          g_pi_x(i) = vamp_probability (grids(i), ihp (phi_x, i))
       end if
    end do
    do i = 1, size (weights)
       g_phi_x(i) = g_pi_x(i) / g_pi_x(channel) * jacobian (phi_x, data, i)
    end do
    w_x = func (phi_x, data, weights, channel, grids) &
         / dot_product (weights, g_phi_x)
  end function vamp_multi_channel
  function vamp_multi_channel0 &
       (func, data, phi, jacobian, x, weights, channel) result (w_x)
    class(vamp_data_t), intent(in) :: data
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), dimension(:), intent(in) :: weights
    integer, intent(in) :: channel
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    interface
       pure function phi (xi, channel) result (x)
         use kinds
         real(kind=default), dimension(:), intent(in) :: xi
         integer, intent(in) :: channel
         real(kind=default), dimension(size(xi)) :: x
       end function phi
    end interface
    interface
       pure function jacobian (x, data, channel) result (j)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: x
         class(vamp_data_t), intent(in) :: data
         integer, intent(in) :: channel
         real(kind=default) :: j
       end function jacobian
    end interface
    real(kind=default) :: w_x
    real(kind=default), dimension(size(x)) :: x_prime
    real(kind=default), dimension(size(weights)) :: g_phi_x
    integer :: i
    x_prime = phi (x, channel)
    do i = 1, size (weights)
       g_phi_x(i) = jacobian (x_prime, data, i)
    end do
    w_x = func (x_prime, data) / dot_product (weights, g_phi_x)
  end function vamp_multi_channel0
  pure subroutine vamp_jacobian (phi, channel, x, region, jacobian, delta_x)
    integer, intent(in) :: channel
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), dimension(:,:), intent(in) :: region
    real(kind=default), intent(out) :: jacobian
    real(kind=default), intent(in), optional :: delta_x
    interface
       pure function phi (xi, channel) result (x)
         use kinds
         real(kind=default), dimension(:), intent(in) :: xi
         integer, intent(in) :: channel
         real(kind=default), dimension(size(xi)) :: x
       end function phi
    end interface
    real(kind=default), dimension(size(x)) :: x_min, x_max
    real(kind=default), dimension(size(x)) :: x_plus, x_minus
    real(kind=default), dimension(size(x),size(x)) :: d_phi
    real(kind=default), parameter :: &
         dx_default = 10.0_default**(-precision(jacobian)/3)
    real(kind=default) :: dx
    integer :: j
    if (present (delta_x)) then
       dx = delta_x
    else
       dx = dx_default
    end if
    x_min = region(1,:)
    x_max = region(2,:)
    x_minus = max (x_min, x)
    x_plus = min (x_max, x)
    do j = 1, size (x)
       x_minus(j) = max (x_min(j), x(j) - dx)
       x_plus(j) = min (x_max(j), x(j) + dx)
       d_phi(:,j) = (phi (x_plus, channel) - phi (x_minus, channel)) &
            / (x_plus(j) - x_minus(j))
       x_minus(j) = max (x_min(j), x(j))
       x_plus(j) = min (x_max(j), x(j))
    end do
    call determinant (d_phi, jacobian)
    jacobian = abs (jacobian)
  end subroutine vamp_jacobian
  subroutine vamp_check_jacobian &
          (rng, n, func, data, phi, channel, region, delta, x_delta)
    type(tao_random_state), intent(inout) :: rng
    integer, intent(in) :: n
    class(vamp_data_t), intent(in) :: data
    integer, intent(in) :: channel
    real(kind=default), dimension(:,:), intent(in) :: region
    real(kind=default), intent(out) :: delta
    real(kind=default), dimension(:), intent(out), optional :: x_delta
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    interface
       pure function phi (xi, channel) result (x)
         use kinds
         real(kind=default), dimension(:), intent(in) :: xi
         integer, intent(in) :: channel
         real(kind=default), dimension(size(xi)) :: x
       end function phi
    end interface
    real(kind=default), dimension(size(region,dim=2)) :: x, r
    real(kind=default) :: jac, d
    real(kind=default), dimension(0) :: wgts
    integer :: i
    delta = 0.0
    do i = 1, max (1, n)
       call tao_random_number (rng, r)
       x = region(1,:) + (region(2,:) - region(1,:)) * r
       call vamp_jacobian (phi, channel, x, region, jac)
       d = func (phi (x, channel), data, wgts, channel) * jac &
            - 1.0_default
       if (abs (d) >= abs (delta)) then
          delta = d
          if (present (x_delta)) then
             x_delta = x
          end if
       end if
     end do
  end subroutine vamp_check_jacobian
  pure subroutine vamp_create_grids &
       (g, domain, num_calls, weights, maps, num_div, &
        stratified, quadrupole, exc)
    type(vamp_grids), intent(inout) :: g
    real(kind=default), dimension(:,:), intent(in) :: domain
    integer, intent(in) :: num_calls
    real(kind=default), dimension(:), intent(in) :: weights
    real(kind=default), dimension(:,:,:), intent(in), optional :: maps
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    type(exception), intent(inout), optional :: exc
    character(len=*), parameter :: FN = "vamp_create_grids"
    integer :: ch, nch
    nch = size (weights)
    allocate (g%grids(nch), g%weights(nch), g%num_calls(nch))
    g%weights = weights / sum (weights)
    g%num_calls = g%weights * num_calls
    do ch = 1, size (g%grids)
       if (present (maps)) then
          call vamp_create_grid &
               (g%grids(ch), domain, g%num_calls(ch), num_div, &
                stratified, quadrupole, map = maps(:,:,ch), exc = exc)
       else
          call vamp_create_grid &
               (g%grids(ch), domain, g%num_calls(ch), num_div, &
                stratified, quadrupole, exc = exc)
       end if
    end do
    g%sum_integral = 0.0
    g%sum_chi2 = 0.0
    g%sum_weights = 0.0
  end subroutine vamp_create_grids
  pure subroutine vamp_create_empty_grids (g)
    type(vamp_grids), intent(inout) :: g
    nullify (g%grids, g%weights, g%num_calls)
  end subroutine vamp_create_empty_grids
  pure subroutine vamp_discard_integrals &
       (g, num_calls, num_div, stratified, quadrupole, exc, eq)
    type(vamp_grids), intent(inout) :: g
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    type(exception), intent(inout), optional :: exc
    type(vamp_equivalences_t), intent(in), optional :: eq
    integer :: ch
    character(len=*), parameter :: FN = "vamp_discard_integrals"
    g%sum_integral = 0.0
    g%sum_weights = 0.0
    g%sum_chi2 = 0.0
    do ch = 1, size (g%grids)
       call vamp_discard_integral (g%grids(ch))
    end do
    if (present (num_calls)) then
       call vamp_reshape_grids &
            (g, num_calls, num_div, stratified, quadrupole, exc, eq)
    end if
  end subroutine vamp_discard_integrals
  pure subroutine vamp_update_weights &
       (g, weights, num_calls, num_div, stratified, quadrupole, exc)
    type(vamp_grids), intent(inout) :: g
    real(kind=default), dimension(:), intent(in) :: weights
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    type(exception), intent(inout), optional :: exc
    character(len=*), parameter :: FN = "vamp_update_weights"
    if (sum (weights) > 0) then
       g%weights = weights / sum (weights)
    else
       g%weights = 1._default / size(g%weights)
    end if
    if (present (num_calls)) then
       call vamp_discard_integrals (g, num_calls, num_div, &
                                   stratified, quadrupole, exc)
    else
       call vamp_discard_integrals (g, sum (g%num_calls), num_div, &
                                   stratified, quadrupole, exc)
    end if
  end subroutine vamp_update_weights
  pure subroutine vamp_reshape_grids &
       (g, num_calls, num_div, stratified, quadrupole, exc, eq)
    type(vamp_grids), intent(inout) :: g
    integer, intent(in) :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    type(exception), intent(inout), optional :: exc
    type(vamp_equivalences_t), intent(in), optional :: eq
    integer, dimension(size(g%grids(1)%num_div)) :: num_div_new
    integer :: ch
    character(len=*), parameter :: FN = "vamp_reshape_grids"
    g%num_calls = g%weights * num_calls
    do ch = 1, size (g%grids)
       if (g%num_calls(ch) >= 2) then
          if (present (eq)) then
             if (present (num_div)) then
                num_div_new = num_div
             else
                num_div_new = g%grids(ch)%num_div
             end if
             where (eq%div_is_invariant(ch,:))
                num_div_new = 1
             end where
             call vamp_reshape_grid (g%grids(ch), g%num_calls(ch), &
                     num_div_new, stratified, quadrupole, exc = exc, &
                     independent = eq%independent(ch), &
                     equivalent_to_ch = eq%equivalent_to_ch(ch), &
                     multiplicity = eq%multiplicity(ch))
          else
             call vamp_reshape_grid (g%grids(ch), g%num_calls(ch), &
                     num_div, stratified, quadrupole, exc = exc)
          end if
       else
          g%num_calls(ch) = 0
       end if
    end do
  end subroutine vamp_reshape_grids
    subroutine vamp_sample_grids &
         (rng, g, func, data, iterations, integral, std_dev, avg_chi2, &
          accuracy, history, histories, exc, eq, warn_error, negative_weights)
      type(tao_random_state), intent(inout) :: rng
      type(vamp_grids), intent(inout) :: g
      class(vamp_data_t), intent(in) :: data
      integer, intent(in) :: iterations
      real(kind=default), intent(out), optional :: integral, std_dev, avg_chi2
      real(kind=default), intent(in), optional :: accuracy
      type(vamp_history), dimension(:), intent(inout), optional :: history
      type(vamp_history), dimension(:,:), intent(inout), optional :: histories
      type(exception), intent(inout), optional :: exc
      type(vamp_equivalences_t), intent(in), optional :: eq
      logical, intent(in), optional :: warn_error, negative_weights
      interface
         function func (xi, data, weights, channel, grids) result (f)
           use kinds
           use vamp_grid_type !NODEP!
           import vamp_data_t
           real(kind=default), dimension(:), intent(in) :: xi
           class(vamp_data_t), intent(in) :: data
           real(kind=default), dimension(:), intent(in), optional :: weights
           integer, intent(in), optional :: channel
           type(vamp_grid), dimension(:), intent(in), optional :: grids
           real(kind=default) :: f
         end function func
      end interface
      integer :: ch, iteration
      logical :: neg_w
      type(exception), dimension(size(g%grids)) :: excs
      logical, dimension(size(g%grids)) :: active
      real(kind=default), dimension(size(g%grids)) :: weights, integrals, std_devs
      real(kind=default) :: local_integral, local_std_dev, local_avg_chi2
      character(len=*), parameter :: FN = "vamp_sample_grids"
      integrals = 0
      std_devs = 0
      neg_w = .false.
      if (present (negative_weights)) neg_w = negative_weights
      active = (g%num_calls >= 2)
      where (active)
         weights = g%num_calls
      elsewhere
         weights = 0.0
      endwhere
      if (sum (weights) /= 0)  weights = weights / sum (weights)
      call clear_exception (excs)
      iterate: do iteration = 1, iterations
         do ch = 1, size (g%grids)
            if (active(ch)) then
               call vamp_discard_integral (g%grids(ch))
               call vamp_sample_grid0 &
                    (rng, g%grids(ch), func, data, &
                     ch, weights, g%grids, excs(ch), neg_w)
               if (present (exc) .and. present (warn_error)) then
                  if (warn_error) call handle_exception (excs(ch))
               end if
               call vamp_average_iterations &
                    (g%grids(ch), iteration, integrals(ch), std_devs(ch), local_avg_chi2)
               if (present (histories)) then
                  if (iteration <= ubound (histories, dim=1)) then
                     call vamp_get_history &
                          (histories(iteration,ch), g%grids(ch), &
                           integrals(ch), std_devs(ch), local_avg_chi2)
                  else
                     call raise_exception (exc, EXC_WARN, FN, "history too short")
                  end if
                  call vamp_terminate_history (histories(iteration+1:,ch))
               end if
            else
               call vamp_nullify_variance (g%grids(ch))
               call vamp_nullify_covariance (g%grids(ch))
            end if
         end do
         if (present(eq))  call vamp_apply_equivalences (g, eq)
         if (iteration < iterations) then
            do ch = 1, size (g%grids)
               active(ch) = (integrals(ch) /= 0)
               if (active(ch)) then
                  call vamp_refine_grid (g%grids(ch))
                end if
            end do
         end if
         if (present (exc) .and. (any (excs%level > 0))) then
            call gather_exceptions (exc, excs)
    !       return
         end if
         call vamp_reduce_channels (g, integrals, std_devs, active)
         call vamp_average_iterations &
              (g, iteration, local_integral, local_std_dev, local_avg_chi2)
         if (present (history)) then
            if (iteration <= size (history)) then
               call vamp_get_history &
                    (history(iteration), g, local_integral, local_std_dev, &
                     local_avg_chi2)
            else
               call raise_exception (exc, EXC_WARN, FN, "history too short")
            end if
            call vamp_terminate_history (history(iteration+1:))
         end if
         if (present (accuracy)) then
            if (local_std_dev <= accuracy * local_integral) then
               call raise_exception (exc, EXC_INFO, FN, &
                    "requested accuracy reached")
               exit iterate
            end if
         end if
      end do iterate
      if (present (integral)) then
         integral = local_integral
      end if
      if (present (std_dev)) then
         std_dev = local_std_dev
      end if
      if (present (avg_chi2)) then
         avg_chi2 = local_avg_chi2
      end if
    end subroutine vamp_sample_grids

  pure subroutine vamp_reduce_channels (g, integrals, std_devs, active)
    type(vamp_grids), intent(inout) :: g
    real(kind=default), dimension(:), intent(in) :: integrals, std_devs
    logical, dimension(:), intent(in) :: active
    real(kind=default) :: this_integral, this_weight, total_calls
    real(kind=default) :: total_variance
    if (.not.any(active)) return
    total_calls = sum (g%num_calls, mask=active)
    if (total_calls > 0) then
       this_integral = sum (g%num_calls * integrals, mask=active) / total_calls
    else
       this_integral = 0
    end if
    total_variance = sum ((g%num_calls*std_devs)**2, mask=active)
    if (total_variance > 0) then
       this_weight = total_calls**2 / total_variance
    else
       this_weight = 0
    end if
    g%sum_weights = g%sum_weights + this_weight
    g%sum_integral = g%sum_integral + this_weight * this_integral
    g%sum_chi2 = g%sum_chi2 + this_weight * this_integral**2
  end subroutine vamp_reduce_channels
  elemental subroutine vamp_average_iterations_grids &
       (g, iteration, integral, std_dev, avg_chi2)
    type(vamp_grids), intent(in) :: g
    integer, intent(in) :: iteration
    real(kind=default), intent(out) :: integral, std_dev, avg_chi2
    real(kind=default), parameter :: eps = 1000 * epsilon (1._default)
    if (g%sum_weights>0) then
       integral = g%sum_integral / g%sum_weights
       std_dev = sqrt (1.0 / g%sum_weights)
       avg_chi2 = &
            max ((g%sum_chi2 - g%sum_integral * integral) / (iteration-0.99), &
                 0.0_default)
       if (avg_chi2 < eps * g%sum_chi2)  avg_chi2 = 0
    else
       integral = 0
       std_dev = 0
       avg_chi2 = 0
    end if
  end subroutine vamp_average_iterations_grids
  pure subroutine vamp_refine_weights (g, power)
    type(vamp_grids), intent(inout) :: g
    real(kind=default), intent(in), optional :: power
    real(kind=default) :: local_power
    real(kind=default), parameter :: DEFAULT_POWER = 0.5_default
    if (present (power)) then
       local_power = power
    else
       local_power = DEFAULT_POWER
    end if
    call vamp_update_weights &
         (g, g%weights * vamp_get_variance (g%grids) ** local_power)
  end subroutine vamp_refine_weights
  pure subroutine vamp_get_history_multi (h, g, integral, std_dev, avg_chi2)
    type(vamp_history), intent(inout) :: h
    type(vamp_grids), intent(in) :: g
    real(kind=default), intent(in) :: integral, std_dev, avg_chi2
    h%calls = sum (g%grids%calls)
    h%stratified = all (g%grids%all_stratified)
    h%integral = 0.0
    h%std_dev = 0.0
    h%avg_integral = integral
    h%avg_std_dev = std_dev
    h%avg_chi2 = avg_chi2
    h%f_min = 0.0
    h%f_max = huge (h%f_max)
    if (h%verbose) then
       h%verbose = .false.
       if (associated (h%div)) then
          deallocate (h%div)
       end if
    end if
  end subroutine vamp_get_history_multi
  function vamp_sum_channels (x, weights, func, data, grids) result (g)
    real(kind=default), dimension(:), intent(in) :: x, weights
    class(vamp_data_t), intent(in) :: data
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    real(kind=default) :: g
    integer :: ch
    g = 0.0
    do ch = 1, size (weights)
       g = g + weights(ch) * func (x, data, weights, ch, grids)
    end do
  end function vamp_sum_channels
  subroutine more_pancake_than_cigar (eval, yes_or_no)
    real(kind=default), dimension(:), intent(in) :: eval
    logical, intent(out) :: yes_or_no
    integer, parameter :: N_CL = 2
    real(kind=default), dimension(size(eval)) :: evals
    real(kind=default), dimension(N_CL) :: cluster_pos
    integer, dimension(N_CL,2) :: clusters
    evals = eval
    call sort (evals)
    call condense (evals, cluster_pos, clusters)
    print *, clusters(1,2) - clusters(1,1) + 1, "small EVs: ", &
         evals(clusters(1,1):clusters(1,2))
    print *, clusters(2,2) - clusters(2,1) + 1, "large EVs: ", &
         evals(clusters(2,1):clusters(2,2))
    if ((clusters(1,2) - clusters(1,1)) &
         < (clusters(2,2) - clusters(2,1))) then
       print *, " => PANCAKE!"
       yes_or_no = .true.
    else
       print *, " => CIGAR!"
       yes_or_no = .false.
    end if
  end subroutine more_pancake_than_cigar
  subroutine select_rotation_axis (cov, r, pancake, cigar)
    real(kind=default), dimension(:,:), intent(in) :: cov
    real(kind=default), dimension(:,:), intent(out) :: r
    integer, intent(in), optional :: pancake, cigar
    integer :: num_pancake, num_cigar
    logical :: like_pancake
    real(kind=default), dimension(size(cov,dim=1),size(cov,dim=2)) :: evecs
    real(kind=default), dimension(size(cov,dim=1)) :: evals, abs_evec
    integer :: iv
    integer, dimension(2) :: i
    real(kind=default) :: cos_theta, sin_theta, norm
    if (present (pancake)) then
       num_pancake = pancake
    else
       num_pancake = -1
    endif
    if (present (cigar)) then
       num_cigar = cigar
    else
       num_cigar = -1
    endif
    call diagonalize_real_symmetric (cov, evals, evecs)
    if (num_pancake > 0) then
       print *, "FORCED PANCAKE: ", num_pancake
       iv = sum (minloc (evals))
    else if (num_cigar > 0) then
       print *, "FORCED CIGAR: ", num_cigar
       iv = sum (maxloc (evals))
    else
       call more_pancake_than_cigar (evals, like_pancake)
       if (like_pancake) then
          iv = sum (minloc (evals))
       else
          iv = sum (maxloc (evals))
       end if
    end if
    abs_evec = abs (evecs(:,iv))
    i(1) = sum (maxloc (abs_evec))
    abs_evec(i(1)) = -1.0
    i(2) = sum (maxloc (abs_evec))
    print *, iv, evals(iv), " => ", evecs(:,iv)
    print *, i(1), abs_evec(i(1)), ", ", i(2), abs_evec(i(2))
    print *, i(1), evecs(i(1),iv), ", ", i(2), evecs(i(2),iv)
    cos_theta = evecs(i(1),iv)
    sin_theta = evecs(i(2),iv)
    norm = 1.0 / sqrt (cos_theta**2 + sin_theta**2)
    cos_theta = cos_theta * norm
    sin_theta = sin_theta * norm
    call unit (r)
    r(i(1),i) =  (/   cos_theta, - sin_theta /)
    r(i(2),i) =  (/   sin_theta,   cos_theta /)
  end subroutine select_rotation_axis
  subroutine select_subspace_explicit (cov, r, subspace)
    real(kind=default), dimension(:,:), intent(in) :: cov
    real(kind=default), dimension(:,:), intent(out) :: r
    integer, dimension(:), intent(in) :: subspace
    real(kind=default), dimension(size(subspace)) :: eval_sub
    real(kind=default), dimension(size(subspace),size(subspace)) :: &
         cov_sub, evec_sub
    cov_sub = cov(subspace,subspace)
    call diagonalize_real_symmetric (cov_sub, eval_sub, evec_sub)
    call unit (r)
    r(subspace,subspace) = evec_sub
  end subroutine select_subspace_explicit
  subroutine select_subspace_guess (cov, r, ndim, pancake, cigar)
    real(kind=default), dimension(:,:), intent(in) :: cov
    real(kind=default), dimension(:,:), intent(out) :: r
    integer, intent(in) :: ndim
    integer, intent(in), optional :: pancake, cigar
    integer :: num_pancake, num_cigar
    logical :: like_pancake
    real(kind=default), dimension(size(cov,dim=1),size(cov,dim=2)) :: evecs
    real(kind=default), dimension(size(cov,dim=1)) :: evals, abs_evec
    integer :: iv, i
    integer, dimension(ndim) :: subspace
    if (present (pancake)) then
       num_pancake = pancake
    else
       num_pancake = -1
    endif
    if (present (cigar)) then
       num_cigar = cigar
    else
       num_cigar = -1
    endif
    call diagonalize_real_symmetric (cov, evals, evecs)
    if (num_pancake > 0) then
       print *, "FORCED PANCAKE: ", num_pancake
       iv = sum (minloc (evals))
    else if (num_cigar > 0) then
       print *, "FORCED CIGAR: ", num_cigar
       iv = sum (maxloc (evals))
    else
       call more_pancake_than_cigar (evals, like_pancake)
       if (like_pancake) then
          iv = sum (minloc (evals))
       else
          iv = sum (maxloc (evals))
       end if
    end if
    abs_evec = abs (evecs(:,iv))
    subspace(1) = sum (maxloc (abs_evec))
    do i = 2, ndim
       abs_evec(subspace(i-1)) = -1.0
       subspace(i) = sum (maxloc (abs_evec))
    end do
    call select_subspace_explicit (cov, r, subspace)
  end subroutine select_subspace_guess
  subroutine vamp_print_covariance (cov)
    real(kind=default), dimension(:,:), intent(in) :: cov
    real(kind=default), dimension(size(cov,dim=1)) :: &
         evals, abs_evals, tmp
    real(kind=default), dimension(size(cov,dim=1),size(cov,dim=2)) :: &
         evecs, abs_evecs
    integer, dimension(size(cov,dim=1)) :: idx
    integer :: i, i_max, j
    i_max = size (evals)
    call diagonalize_real_symmetric (cov, evals, evecs)
    call sort (evals, evecs)
    abs_evals = abs (evals)
    abs_evecs = abs (evecs)
    print "(1X,A78)", repeat ("-", 78)
    print "(1X,A)", "Eigenvalues and eigenvectors:"
    print "(1X,A78)", repeat ("-", 78)
    do i = 1, i_max
       print "(1X,I2,A1,1X,E11.4,1X,A1,10(10(1X,F5.2)/,18X))", &
            i, ":", evals(i), "|", evecs(:,i)
    end do
    print "(1X,A78)", repeat ("-", 78)
    print "(1X,A)", "Approximate subspaces:"
    print "(1X,A78)", repeat ("-", 78)
    do i = 1, i_max
       idx = (/ (j, j=1,i_max) /)
       tmp = abs_evecs(:,i)
       call sort (tmp, idx, reverse = .true.)
       print "(1X,I2,A1,1X,E11.4,1X,A1,10(1X,I5))", &
            i, ":", evals(i), "|", idx(1:min(10,size(idx)))
       print "(17X,A1,10(1X,F5.2))", &
                              "|", evecs(idx(1:min(10,size(idx))),i)
    end do
    print "(1X,A78)", repeat ("-", 78)
  end subroutine vamp_print_covariance
  subroutine condense (x, cluster_pos, clusters, linear)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), dimension(:), intent(out) :: cluster_pos
    integer, dimension(:,:), intent(out) :: clusters
    logical, intent(in), optional :: linear
    logical :: linear_metric
    real(kind=default), dimension(size(x)) :: cl_pos
    real(kind=default) :: wgt0, wgt1
    integer :: cl_num
    integer, dimension(size(x),2) :: cl
    integer :: i, gap
    linear_metric = .false.
    if (present (linear)) then
       linear_metric = linear
    end if
    cl_pos = x
    cl_num = size (cl_pos)
    cl = spread ((/ (i, i=1,cl_num) /), dim = 2, ncopies = 2)
    do cl_num = size (cl_pos), size (cluster_pos) + 1, -1
       if (linear_metric) then
          gap = sum (minloc (cl_pos(2:cl_num) - cl_pos(1:cl_num-1)))
       else
          gap = sum (minloc (cl_pos(2:cl_num) / cl_pos(1:cl_num-1)))
       end if
       wgt0 = cl(gap,2) - cl(gap,1) + 1
       wgt1 = cl(gap+1,2) - cl(gap+1,1) + 1
       cl_pos(gap) = (wgt0 * cl_pos(gap) + wgt1 * cl_pos(gap+1)) / (wgt0 + wgt1)
       cl(gap,2) = cl(gap+1,2)
       cl_pos(gap+1:cl_num-1) = cl_pos(gap+2:cl_num)
       cl(gap+1:cl_num-1,:) = cl(gap+2:cl_num,:)
       print *, cl_num, ": action = ", condense_action (x, cl)
    end do
    cluster_pos = cl_pos(1:cl_num)
    clusters = cl(1:cl_num,:)
  end subroutine condense
  function condense_action (positions, clusters) result (s)
    real(kind=default), dimension(:), intent(in) :: positions
    integer, dimension(:,:), intent(in) :: clusters
    real(kind=default) :: s
    integer :: i
    integer, parameter :: POWER = 2
    s = 0
    do i = 1, size (clusters, dim = 1)
       s = s + standard_deviation (positions(clusters(i,1) &
                                             :clusters(i,2))) ** POWER
    end do
  end function condense_action
  subroutine vamp_next_event_single &
       (x, rng, g, func, data, &
        weight, channel, weights, grids, exc)
    real(kind=default), dimension(:), intent(out) :: x
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grid), intent(inout) :: g
    real(kind=default), intent(out), optional :: weight
    class(vamp_data_t), intent(in) :: data
    integer, intent(in), optional :: channel
    real(kind=default), dimension(:), intent(in), optional :: weights
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    type(exception), intent(inout), optional :: exc
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    character(len=*), parameter :: FN = "vamp_next_event_single"
    real(kind=default), dimension(size(g%div)):: wgts
    real(kind=default), dimension(size(g%div)):: r
    integer, dimension(size(g%div)):: ia
    real(kind=default) :: f, wgt
    real(kind=default) :: r0
    rejection: do
       call tao_random_number (rng, r)
       call inject_division_short (g%div, real(r, kind=default), x, ia, wgts)
       wgt = g%jacobi * product (wgts)
       wgt = g%calls * wgt !: the calling procedure will divide by \#calls
       if (associated (g%map)) then
          x = matmul (g%map, x)
       end if
       if (associated (g%map)) then
          if (all (inside_division (g%div, x))) then
             f = wgt * func (x, data, weights, channel, grids)
          else
             f = 0.0
          end if
       else
          f = wgt * func (x, data, weights, channel, grids)
       end if
       ! call record_efficiency (g%div, ia, f/g%f_max)
       if (present (weight)) then
          weight = f
          exit rejection
       else
          if (abs(f) > g%f_max) then
             g%f_max = f
             call raise_exception (exc, EXC_WARN, FN, "weight > 1")
             exit rejection
          end if
          call tao_random_number (rng, r0)
          if (r0 * g%f_max <= abs(f)) then
             exit rejection
          end if
       end if
    end do rejection
  end subroutine vamp_next_event_single
  subroutine vamp_next_event_multi &
       (x, rng, g, func, data, phi, weight, excess, positive, exc)
    real(kind=default), dimension(:), intent(out) :: x
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grids), intent(inout) :: g
    class(vamp_data_t), intent(in) :: data
    real(kind=default), intent(out), optional :: weight
    real(kind=default), intent(out), optional :: excess
    logical, intent(out), optional :: positive
    type(exception), intent(inout), optional :: exc
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    interface
       pure function phi (xi, channel) result (x)
         use kinds
         real(kind=default), dimension(:), intent(in) :: xi
         integer, intent(in) :: channel
         real(kind=default), dimension(size(xi)) :: x
       end function phi
    end interface
    character(len=*), parameter :: FN = "vamp_next_event_multi"
    real(kind=default), dimension(size(x)) :: xi
    real(kind=default) :: r, wgt
    real(kind=default), dimension(size(g%weights)) :: weights
    integer :: channel
    if (any (g%grids%f_max > 0)) then
       weights = g%weights * g%grids%f_max
    else
       weights = g%weights
    end if
    weights = weights / sum (weights)
    rejection: do
       call tao_random_number (rng, r)
       select_channel: do channel = 1, size (g%weights)
          r = r - weights(channel)
          if (r <= 0.0) then
             exit select_channel
          end if
       end do select_channel
       channel = min (channel, size (g%weights)) !: for $r=1$ and rounding errors
       call vamp_next_event_single &
            (xi, rng, g%grids(channel), func, data, wgt, &
             channel, g%weights, g%grids, exc)
       if (present (weight)) then
          weight = wgt * g%weights(channel) / weights(channel)
          exit rejection
       else
          if (abs (wgt) > g%grids(channel)%f_max) then
             if (present(excess)) then
                excess = abs (wgt) / g%grids(channel)%f_max - 1
             else
               call raise_exception (exc, EXC_WARN, FN, "weight > 1")
          !      print *, "weight > 1 (", wgt/g%grids(channel)%f_max, &
          !           & ") in channel ", channel

             end if
          !  exit rejection
          else
             if (present(excess)) excess = 0
          end if
          call tao_random_number (rng, r)
          if (r * g%grids(channel)%f_max <= abs (wgt)) then
             if (present (positive))  positive = wgt >= 0
             exit rejection
          end if
       end if
    end do rejection
    x = phi (xi, channel)
  end subroutine vamp_next_event_multi
  subroutine vamp_warmup_grid &
       (rng, g, func, data, iterations, exc, history)
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grid), intent(inout) :: g
    class(vamp_data_t), intent(in) :: data
    integer, intent(in) :: iterations
    type(exception), intent(inout), optional :: exc
    type(vamp_history), dimension(:), intent(inout), optional :: history
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    call vamp_sample_grid &
       (rng, g, func, data, &
        iterations - 1, exc = exc, history = history)
    call vamp_sample_grid0 (rng, g, func, data, exc = exc)
  end subroutine vamp_warmup_grid
  subroutine vamp_warmup_grids &
       (rng, g, func, data, iterations, history, histories, exc)
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grids), intent(inout) :: g
    class(vamp_data_t), intent(in) :: data
    integer, intent(in) :: iterations
    type(vamp_history), dimension(:), intent(inout), optional :: history
    type(vamp_history), dimension(:,:), intent(inout), optional :: histories
    type(exception), intent(inout), optional :: exc
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    integer :: ch
    logical, dimension(size(g%grids)) :: active
    real(kind=default), dimension(size(g%grids)) :: weights
    active = (g%num_calls >= 2)
    where (active)
       weights = g%num_calls
    elsewhere
       weights = 0.0
    end where
    weights = weights / sum (weights)
    call vamp_sample_grids (rng, g, func, data, iterations - 1, &
                            exc = exc, history = history, histories = histories)
    do ch = 1, size (g%grids)
       if (g%grids(ch)%num_calls >= 2) then
          call vamp_sample_grid0 &
               (rng, g%grids(ch), func, data, &
                ch, weights, g%grids, exc = exc)
       end if
    end do
  end subroutine vamp_warmup_grids
  subroutine vamp_integrate_grid &
       (rng, g, func, data, calls, integral, std_dev, avg_chi2, num_div, &
        stratified, quadrupole, accuracy, exc, history)
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grid), intent(inout) :: g
    class(vamp_data_t), intent(in) :: data
    integer, dimension(:,:), intent(in) :: calls
    real(kind=default), intent(out), optional :: integral, std_dev, avg_chi2
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    real(kind=default), intent(in), optional :: accuracy
    type(exception), intent(inout), optional :: exc
    type(vamp_history), dimension(:), intent(inout), optional :: history
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    character(len=*), parameter :: FN = "vamp_integrate_grid"
    integer :: step, last_step, it
    last_step = size (calls, dim = 2)
    it = 1
    do step = 1, last_step - 1
       call vamp_discard_integral (g, calls(2,step), num_div, &
                                  stratified, quadrupole, exc = exc)
       call vamp_sample_grid (rng, g, func, data, calls(1,step), &
                              exc = exc, history = history(it:))
       if (present (exc)) then
          if (exc%level > EXC_WARN) then
             return
          end if
       end if
       it = it + calls(1,step)
    end do
    call vamp_discard_integral (g, calls(2,last_step), exc = exc)
    call vamp_sample_grid (rng, g, func, data, calls(1,last_step), &
                           integral, std_dev, avg_chi2, accuracy, exc = exc, &
                           history = history(it:))
  end subroutine vamp_integrate_grid
  subroutine vamp_integrate_region &
       (rng, region, func, data, calls, &
        integral, std_dev, avg_chi2, num_div, &
        stratified, quadrupole, accuracy, map, covariance, exc, history)
    type(tao_random_state), intent(inout) :: rng
    real(kind=default), dimension(:,:), intent(in) :: region
    class(vamp_data_t), intent(in) :: data
    integer, dimension(:,:), intent(in) :: calls
    real(kind=default), intent(out), optional :: integral, std_dev, avg_chi2
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    real(kind=default), intent(in), optional :: accuracy
    real(kind=default), dimension(:,:), intent(in), optional :: map
    real(kind=default), dimension(:,:), intent(out), optional :: covariance
    type(exception), intent(inout), optional :: exc
    type(vamp_history), dimension(:), intent(inout), optional :: history
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    character(len=*), parameter :: FN = "vamp_integrate_region"
    type(vamp_grid) :: g
    call vamp_create_grid &
         (g, region, calls(2,1), num_div, &
          stratified, quadrupole, present (covariance), map, exc)
    call vamp_integrate_grid &
         (rng, g, func, data, calls, &
          integral, std_dev, avg_chi2, num_div, &
          accuracy = accuracy, exc = exc, history = history)
    if (present (covariance)) then
       covariance = vamp_get_covariance (g)
    end if
    call vamp_delete_grid (g)
  end subroutine vamp_integrate_region
  subroutine vamp_integratex_region &
       (rng, region, func, data, calls, integral, std_dev, avg_chi2, &
        num_div, stratified, quadrupole, accuracy, pancake, cigar, &
        exc, history)
    type(tao_random_state), intent(inout) :: rng
    real(kind=default), dimension(:,:), intent(in) :: region
    class(vamp_data_t), intent(in) :: data
    integer, dimension(:,:,:), intent(in) :: calls
    real(kind=default), intent(out), optional :: integral, std_dev, avg_chi2
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    real(kind=default), intent(in), optional :: accuracy
    integer, intent(in), optional :: pancake, cigar
    type(exception), intent(inout), optional :: exc
    type(vamp_history), dimension(:), intent(inout), optional :: history
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    real(kind=default), dimension(size(region,dim=2)) :: eval
    real(kind=default), dimension(size(region,dim=2),size(region,dim=2)) :: evec
    type(vamp_grid) :: g
    integer :: step, last_step, it
    it = 1
    call vamp_create_grid &
         (g, region, calls(2,1,1), num_div, &
          stratified, quadrupole, covariance = .true., exc = exc)
    call vamp_integrate_grid &
         (rng, g, func, data, calls(:,:,1), num_div = num_div, &
          exc = exc, history = history(it:))
    if (present (exc)) then
       if (exc%level > EXC_WARN) then
          return
       end if
    end if
    it = it + sum (calls(1,:,1))
    last_step = size (calls, dim = 3)
    do step = 2, last_step - 1
       call diagonalize_real_symmetric (vamp_get_covariance(g), eval, evec)
       call sort (eval, evec)
       call select_rotation_axis (vamp_get_covariance(g), evec, pancake, cigar)
       call vamp_delete_grid (g)
       call vamp_create_grid &
            (g, region, calls(2,1,step), num_div, stratified, quadrupole, &
             covariance = .true., map = evec, exc = exc)
       call vamp_integrate_grid &
            (rng, g, func, data, calls(:,:,step), num_div = num_div, &
             exc = exc, history = history(it:))
       if (present (exc)) then
          if (exc%level > EXC_WARN) then
             return
          end if
       end if
       it = it + sum (calls(1,:,step))
    end do
    call diagonalize_real_symmetric (vamp_get_covariance(g), eval, evec)
    call sort (eval, evec)
    call select_rotation_axis (vamp_get_covariance(g), evec, pancake, cigar)
    call vamp_delete_grid (g)
    call vamp_create_grid &
         (g, region, calls(2,1,last_step), num_div, stratified, quadrupole, &
          covariance = .true., map = evec, exc = exc)
    call vamp_integrate_grid &
         (rng, g, func, data, calls(:,:,last_step), &
          integral, std_dev, avg_chi2, &
          num_div = num_div, exc = exc, history = history(it:))
    call vamp_delete_grid (g)
  end subroutine vamp_integratex_region
  subroutine write_grid_unit (g, unit, write_integrals)
    type(vamp_grid), intent(in) :: g
    integer, intent(in) :: unit
    logical, intent(in), optional :: write_integrals
    integer :: i, j
    write (unit = unit, fmt = descr_fmt) "begin type(vamp_grid) :: g"
    write (unit = unit, fmt = integer_fmt) "size (g%div) = ", size (g%div)
    write (unit = unit, fmt = integer_fmt) "num_calls = ", g%num_calls
    write (unit = unit, fmt = integer_fmt) "calls_per_cell = ", g%calls_per_cell
    write (unit = unit, fmt = logical_fmt) "stratified = ", g%stratified
    write (unit = unit, fmt = logical_fmt) "all_stratified = ", g%all_stratified
    write (unit = unit, fmt = logical_fmt) "quadrupole = ", g%quadrupole
    write (unit = unit, fmt = double_fmt) "mu(1) = ", g%mu(1)
    write (unit = unit, fmt = double_fmt) "mu(2) = ", g%mu(2)
    write (unit = unit, fmt = double_fmt) "mu_plus(1) = ", g%mu_plus(1)
    write (unit = unit, fmt = double_fmt) "mu_plus(2) = ", g%mu_plus(2)
    write (unit = unit, fmt = double_fmt) "mu_minus(1) = ", g%mu_minus(1)
    write (unit = unit, fmt = double_fmt) "mu_minus(2) = ", g%mu_minus(2)
    write (unit = unit, fmt = double_fmt) "sum_integral = ", g%sum_integral
    write (unit = unit, fmt = double_fmt) "sum_weights = ", g%sum_weights
    write (unit = unit, fmt = double_fmt) "sum_chi2 = ", g%sum_chi2
    write (unit = unit, fmt = double_fmt) "calls = ", g%calls
    write (unit = unit, fmt = double_fmt) "dv2g = ", g%dv2g
    write (unit = unit, fmt = double_fmt) "jacobi = ", g%jacobi
    write (unit = unit, fmt = double_fmt) "f_min = ", g%f_min
    write (unit = unit, fmt = double_fmt) "f_max = ", g%f_max
    write (unit = unit, fmt = double_fmt) "mu_gi = ", g%mu_gi
    write (unit = unit, fmt = double_fmt) "sum_mu_gi = ", g%sum_mu_gi
    write (unit = unit, fmt = descr_fmt) "begin g%num_div"
    do i = 1, size (g%div)
       write (unit = unit, fmt = integer_array_fmt) i, g%num_div(i)
    end do
    write (unit = unit, fmt = descr_fmt) "end g%num_div"
    write (unit = unit, fmt = descr_fmt) "begin g%div"
    do i = 1, size (g%div)
       call write_division (g%div(i), unit, write_integrals)
    end do
    write (unit = unit, fmt = descr_fmt) "end g%div"
    if (associated (g%map)) then
       write (unit = unit, fmt = descr_fmt) "begin g%map"
       do i = 1, size (g%div)
          do j = 1, size (g%div)
             write (unit = unit, fmt = double_array2_fmt) i, j, g%map(i,j)
          end do
       end do
       write (unit = unit, fmt = descr_fmt) "end g%map"
    else
       write (unit = unit, fmt = descr_fmt) "empty g%map"
    end if
    if (associated (g%mu_x)) then
       write (unit = unit, fmt = descr_fmt) "begin g%mu_x"
       do i = 1, size (g%div)
          write (unit = unit, fmt = double_array_fmt) i, g%mu_x(i)
          write (unit = unit, fmt = double_array_fmt) i, g%sum_mu_x(i)
          do j = 1, size (g%div)
             write (unit = unit, fmt = double_array2_fmt) i, j, g%mu_xx(i,j)
             write (unit = unit, fmt = double_array2_fmt) i, j, g%sum_mu_xx(i,j)
          end do
       end do
       write (unit = unit, fmt = descr_fmt) "end g%mu_x"
    else
       write (unit = unit, fmt = descr_fmt) "empty g%mu_x"
    end if
    write (unit = unit, fmt = descr_fmt) "end type(vamp_grid)"
  end subroutine write_grid_unit
  subroutine read_grid_unit (g, unit, read_integrals)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in) :: unit
    logical, intent(in), optional :: read_integrals
    character(len=*), parameter :: FN = "vamp_read_grid"
    character(len=80) :: chdum
    integer :: ndim, i, j, idum, jdum
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = integer_fmt) chdum, ndim
    if (associated (g%div)) then
       if (size (g%div) /= ndim) then
          call delete_division (g%div)
          deallocate (g%div)
          allocate (g%div(ndim))
          call create_empty_division (g%div)
       end if
    else
       allocate (g%div(ndim))
       call create_empty_division (g%div)
    end if
    call create_array_pointer (g%num_div, ndim)
    read (unit = unit, fmt = integer_fmt) chdum, g%num_calls
    read (unit = unit, fmt = integer_fmt) chdum, g%calls_per_cell
    read (unit = unit, fmt = logical_fmt) chdum, g%stratified
    read (unit = unit, fmt = logical_fmt) chdum, g%all_stratified
    read (unit = unit, fmt = logical_fmt) chdum, g%quadrupole
    read (unit = unit, fmt = double_fmt) chdum, g%mu(1)
    read (unit = unit, fmt = double_fmt) chdum, g%mu(2)
    read (unit = unit, fmt = double_fmt) chdum, g%mu_plus(1)
    read (unit = unit, fmt = double_fmt) chdum, g%mu_plus(2)
    read (unit = unit, fmt = double_fmt) chdum, g%mu_minus(1)
    read (unit = unit, fmt = double_fmt) chdum, g%mu_minus(2)
    read (unit = unit, fmt = double_fmt) chdum, g%sum_integral
    read (unit = unit, fmt = double_fmt) chdum, g%sum_weights
    read (unit = unit, fmt = double_fmt) chdum, g%sum_chi2
    read (unit = unit, fmt = double_fmt) chdum, g%calls
    read (unit = unit, fmt = double_fmt) chdum, g%dv2g
    read (unit = unit, fmt = double_fmt) chdum, g%jacobi
    read (unit = unit, fmt = double_fmt) chdum, g%f_min
    read (unit = unit, fmt = double_fmt) chdum, g%f_max
    read (unit = unit, fmt = double_fmt) chdum, g%mu_gi
    read (unit = unit, fmt = double_fmt) chdum, g%sum_mu_gi
    read (unit = unit, fmt = descr_fmt) chdum
    do i = 1, size (g%div)
       read (unit = unit, fmt = integer_array_fmt) idum, g%num_div(i)
    end do
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = descr_fmt) chdum
    do i = 1, size (g%div)
       call read_division (g%div(i), unit, read_integrals)
    end do
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = descr_fmt) chdum
    if (chdum == "begin g%map") then
       call create_array_pointer (g%map, (/ ndim, ndim /))
       do i = 1, size (g%div)
          do j = 1, size (g%div)
             read (unit = unit, fmt = double_array2_fmt) idum, jdum, g%map(i,j)
          end do
       end do
       read (unit = unit, fmt = descr_fmt) chdum
    else
       if (associated (g%map)) then
          deallocate (g%map)
       end if
    end if
    read (unit = unit, fmt = descr_fmt) chdum
    if (chdum == "begin g%mu_x") then
       call create_array_pointer (g%mu_x, ndim )
       call create_array_pointer (g%sum_mu_x, ndim)
       call create_array_pointer (g%mu_xx, (/ ndim, ndim /))
       call create_array_pointer (g%sum_mu_xx, (/ ndim, ndim /))
       do i = 1, size (g%div)
          read (unit = unit, fmt = double_array_fmt) idum, jdum, g%mu_x(i)
          read (unit = unit, fmt = double_array_fmt) idum, jdum, g%sum_mu_x(i)
          do j = 1, size (g%div)
             read (unit = unit, fmt = double_array2_fmt) &
                  idum, jdum, g%mu_xx(i,j)
             read (unit = unit, fmt = double_array2_fmt) &
                  idum, jdum, g%sum_mu_xx(i,j)
          end do
       end do
       read (unit = unit, fmt = descr_fmt) chdum
    else
       if (associated (g%mu_x)) then
          deallocate (g%mu_x)
       end if
       if (associated (g%mu_xx)) then
          deallocate (g%mu_xx)
       end if
       if (associated (g%sum_mu_x)) then
          deallocate (g%sum_mu_x)
       end if
       if (associated (g%sum_mu_xx)) then
          deallocate (g%sum_mu_xx)
       end if
    end if
    read (unit = unit, fmt = descr_fmt) chdum
  end subroutine read_grid_unit
  subroutine write_grid_name (g, name, write_integrals)
    type(vamp_grid), intent(inout) :: g
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: write_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "write", status = "replace", file = name)
    call write_grid_unit (g, unit, write_integrals)
    close (unit = unit)
  end subroutine write_grid_name
  subroutine read_grid_name (g, name, read_integrals)
    type(vamp_grid), intent(inout) :: g
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: read_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "read", status = "old", file = name)
    call read_grid_unit (g, unit, read_integrals)
    close (unit = unit)
  end subroutine read_grid_name
  subroutine write_grids_unit (g, unit, write_integrals)
    type(vamp_grids), intent(in) :: g
    integer, intent(in) :: unit
    logical, intent(in), optional :: write_integrals
    integer :: i
    write (unit = unit, fmt = descr_fmt) "begin type(vamp_grids) :: g"
    write (unit = unit, fmt = integer_fmt) "size (g%grids) = ", size (g%grids)
    write (unit = unit, fmt = double_fmt) "sum_integral = ", g%sum_integral
    write (unit = unit, fmt = double_fmt) "sum_weights = ", g%sum_weights
    write (unit = unit, fmt = double_fmt) "sum_chi2 = ", g%sum_chi2
    write (unit = unit, fmt = descr_fmt) "begin g%weights"
    do i = 1, size (g%grids)
       write (unit = unit, fmt = double_array_fmt) i, g%weights(i)
    end do
    write (unit = unit, fmt = descr_fmt) "end g%weights"
    write (unit = unit, fmt = descr_fmt) "begin g%num_calls"
    do i = 1, size (g%grids)
       write (unit = unit, fmt = integer_array_fmt) i, g%num_calls(i)
    end do
    write (unit = unit, fmt = descr_fmt) "end g%num_calls"
    write (unit = unit, fmt = descr_fmt) "begin g%grids"
    do i = 1, size (g%grids)
       call write_grid_unit (g%grids(i), unit, write_integrals)
    end do
    write (unit = unit, fmt = descr_fmt) "end g%grids"
    write (unit = unit, fmt = descr_fmt) "end type(vamp_grids)"
  end subroutine write_grids_unit
  subroutine read_grids_unit (g, unit, read_integrals)
    type(vamp_grids), intent(inout) :: g
    integer, intent(in) :: unit
    logical, intent(in), optional :: read_integrals
    character(len=*), parameter :: FN = "vamp_read_grids"
    character(len=80) :: chdum
    integer :: i, nch, idum
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = integer_fmt) chdum, nch
    if (associated (g%grids)) then
       if (size (g%grids) /= nch) then
          call vamp_delete_grid (g%grids)
          deallocate (g%grids, g%weights, g%num_calls)
          allocate (g%grids(nch), g%weights(nch), g%num_calls(nch))
          call vamp_create_empty_grid (g%grids)
       end if
    else
       allocate (g%grids(nch), g%weights(nch), g%num_calls(nch))
       call vamp_create_empty_grid (g%grids)
    end if
    read (unit = unit, fmt = double_fmt) chdum, g%sum_integral
    read (unit = unit, fmt = double_fmt) chdum, g%sum_weights
    read (unit = unit, fmt = double_fmt) chdum, g%sum_chi2
    read (unit = unit, fmt = descr_fmt) chdum
    do i = 1, nch
       read (unit = unit, fmt = double_array_fmt) idum, g%weights(i)
    end do
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = descr_fmt) chdum
    do i = 1, nch
       read (unit = unit, fmt = integer_array_fmt) idum, g%num_calls(i)
    end do
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = descr_fmt) chdum
    do i = 1, nch
       call read_grid_unit (g%grids(i), unit, read_integrals)
    end do
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = descr_fmt) chdum
  end subroutine read_grids_unit
  subroutine write_grids_name (g, name, write_integrals)
    type(vamp_grids), intent(inout) :: g
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: write_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "write", status = "replace", file = name)
    call write_grids_unit (g, unit, write_integrals)
    close (unit = unit)
  end subroutine write_grids_name
  subroutine read_grids_name (g, name, read_integrals)
    type(vamp_grids), intent(inout) :: g
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: read_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "read", status = "old", file = name)
    call read_grids_unit (g, unit, read_integrals)
    close (unit = unit)
  end subroutine read_grids_name
  subroutine write_grid_raw_unit (g, unit, write_integrals)
    type(vamp_grid), intent(in) :: g
    integer, intent(in) :: unit
    logical, intent(in), optional :: write_integrals
    integer :: i, j
    write (unit = unit) MAGIC_GRID_BEGIN
    write (unit = unit) size (g%div)
    write (unit = unit) g%num_calls
    write (unit = unit) g%calls_per_cell
    write (unit = unit) g%stratified
    write (unit = unit) g%all_stratified
    write (unit = unit) g%quadrupole
    write (unit = unit) g%mu(1)
    write (unit = unit) g%mu(2)
    write (unit = unit) g%mu_plus(1)
    write (unit = unit) g%mu_plus(2)
    write (unit = unit) g%mu_minus(1)
    write (unit = unit) g%mu_minus(2)
    write (unit = unit) g%sum_integral
    write (unit = unit) g%sum_weights
    write (unit = unit) g%sum_chi2
    write (unit = unit) g%calls
    write (unit = unit) g%dv2g
    write (unit = unit) g%jacobi
    write (unit = unit) g%f_min
    write (unit = unit) g%f_max
    write (unit = unit) g%mu_gi
    write (unit = unit) g%sum_mu_gi
    do i = 1, size (g%div)
       write (unit = unit) g%num_div(i)
    end do
    do i = 1, size (g%div)
       call write_division_raw (g%div(i), unit, write_integrals)
    end do
    if (associated (g%map)) then
       write (unit = unit) MAGIC_GRID_MAP
       do i = 1, size (g%div)
          do j = 1, size (g%div)
             write (unit = unit) g%map(i,j)
          end do
       end do
    else
       write (unit = unit) MAGIC_GRID_EMPTY
    end if
    if (associated (g%mu_x)) then
       write (unit = unit) MAGIC_GRID_MU_X
       do i = 1, size (g%div)
          write (unit = unit) g%mu_x(i)
          write (unit = unit) g%sum_mu_x(i)
          do j = 1, size (g%div)
             write (unit = unit) g%mu_xx(i,j)
             write (unit = unit) g%sum_mu_xx(i,j)
          end do
       end do
    else
       write (unit = unit) MAGIC_GRID_EMPTY
    end if
    write (unit = unit) MAGIC_GRID_END
  end subroutine write_grid_raw_unit
  subroutine read_grid_raw_unit (g, unit, read_integrals)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in) :: unit
    logical, intent(in), optional :: read_integrals
    character(len=*), parameter :: FN = "vamp_read_raw_grid"
    integer :: ndim, i, j, magic
    read (unit = unit) magic
    if (magic /= MAGIC_GRID_BEGIN) then
       print *, FN, " fatal: expecting magic ", MAGIC_GRID_BEGIN, &
                    ", found ", magic
       stop
    end if
    read (unit = unit) ndim
    if (associated (g%div)) then
       if (size (g%div) /= ndim) then
          call delete_division (g%div)
          deallocate (g%div)
          allocate (g%div(ndim))
          call create_empty_division (g%div)
       end if
    else
       allocate (g%div(ndim))
       call create_empty_division (g%div)
    end if
    call create_array_pointer (g%num_div, ndim)
    read (unit = unit) g%num_calls
    read (unit = unit) g%calls_per_cell
    read (unit = unit) g%stratified
    read (unit = unit) g%all_stratified
    read (unit = unit) g%quadrupole
    read (unit = unit) g%mu(1)
    read (unit = unit) g%mu(2)
    read (unit = unit) g%mu_plus(1)
    read (unit = unit) g%mu_plus(2)
    read (unit = unit) g%mu_minus(1)
    read (unit = unit) g%mu_minus(2)
    read (unit = unit) g%sum_integral
    read (unit = unit) g%sum_weights
    read (unit = unit) g%sum_chi2
    read (unit = unit) g%calls
    read (unit = unit) g%dv2g
    read (unit = unit) g%jacobi
    read (unit = unit) g%f_min
    read (unit = unit) g%f_max
    read (unit = unit) g%mu_gi
    read (unit = unit) g%sum_mu_gi
    do i = 1, size (g%div)
       read (unit = unit) g%num_div(i)
    end do
    do i = 1, size (g%div)
       call read_division_raw (g%div(i), unit, read_integrals)
    end do
    read (unit = unit) magic
    if (magic == MAGIC_GRID_MAP) then
       call create_array_pointer (g%map, (/ ndim, ndim /))
       do i = 1, size (g%div)
          do j = 1, size (g%div)
             read (unit = unit) g%map(i,j)
          end do
       end do
    else if (magic == MAGIC_GRID_EMPTY) then
       if (associated (g%map)) then
          deallocate (g%map)
       end if
    else
       print *, FN, " fatal: expecting magic ", MAGIC_GRID_EMPTY, &
                    " or ", MAGIC_GRID_MAP, ", found ", magic
       stop
    end if
    read (unit = unit) magic
    if (magic == MAGIC_GRID_MU_X) then
       call create_array_pointer (g%mu_x, ndim )
       call create_array_pointer (g%sum_mu_x, ndim)
       call create_array_pointer (g%mu_xx, (/ ndim, ndim /))
       call create_array_pointer (g%sum_mu_xx, (/ ndim, ndim /))
       do i = 1, size (g%div)
          read (unit = unit) g%mu_x(i)
          read (unit = unit) g%sum_mu_x(i)
          do j = 1, size (g%div)
             read (unit = unit) g%mu_xx(i,j)
             read (unit = unit) g%sum_mu_xx(i,j)
          end do
       end do
    else if (magic == MAGIC_GRID_EMPTY) then
       if (associated (g%mu_x)) then
          deallocate (g%mu_x)
       end if
       if (associated (g%mu_xx)) then
          deallocate (g%mu_xx)
       end if
       if (associated (g%sum_mu_x)) then
          deallocate (g%sum_mu_x)
       end if
       if (associated (g%sum_mu_xx)) then
          deallocate (g%sum_mu_xx)
       end if
    else
       print *, FN, " fatal: expecting magic ", MAGIC_GRID_EMPTY, &
                    " or ", MAGIC_GRID_MU_X, ", found ", magic
       stop
    end if
    read (unit = unit) magic
    if (magic /= MAGIC_GRID_END) then
       print *, FN, " fatal: expecting magic ", MAGIC_GRID_END, &
                    " found ", magic
       stop
    end if
  end subroutine read_grid_raw_unit
  subroutine write_grid_raw_name (g, name, write_integrals)
    type(vamp_grid), intent(inout) :: g
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: write_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "write", status = "replace", &
          form = "unformatted", file = name)
    call write_grid_raw_unit (g, unit, write_integrals)
    close (unit = unit)
  end subroutine write_grid_raw_name
  subroutine read_grid_raw_name (g, name, read_integrals)
    type(vamp_grid), intent(inout) :: g
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: read_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "read", status = "old", &
          form = "unformatted",  file = name)
    call read_grid_raw_unit (g, unit, read_integrals)
    close (unit = unit)
  end subroutine read_grid_raw_name
  subroutine write_grids_raw_unit (g, unit, write_integrals)
    type(vamp_grids), intent(in) :: g
    integer, intent(in) :: unit
    logical, intent(in), optional :: write_integrals
    integer :: i
    write (unit = unit) MAGIC_GRIDS_BEGIN
    write (unit = unit) size (g%grids)
    write (unit = unit) g%sum_integral
    write (unit = unit) g%sum_weights
    write (unit = unit) g%sum_chi2
    do i = 1, size (g%grids)
       write (unit = unit) g%weights(i)
    end do
    do i = 1, size (g%grids)
       write (unit = unit) g%num_calls(i)
    end do
    do i = 1, size (g%grids)
       call write_grid_raw_unit (g%grids(i), unit, write_integrals)
    end do
    write (unit = unit) MAGIC_GRIDS_END
  end subroutine write_grids_raw_unit
  subroutine read_grids_raw_unit (g, unit, read_integrals)
    type(vamp_grids), intent(inout) :: g
    integer, intent(in) :: unit
    logical, intent(in), optional :: read_integrals
    character(len=*), parameter :: FN = "vamp_read_grids_raw"
    integer :: i, nch, magic
    read (unit = unit) magic
    if (magic /= MAGIC_GRIDS_BEGIN) then
       print *, FN, " fatal: expecting magic ", MAGIC_GRIDS_BEGIN, &
                    " found ", magic
       stop
    end if
    read (unit = unit) nch
    if (associated (g%grids)) then
       if (size (g%grids) /= nch) then
          call vamp_delete_grid (g%grids)
          deallocate (g%grids, g%weights, g%num_calls)
          allocate (g%grids(nch), g%weights(nch), g%num_calls(nch))
          call vamp_create_empty_grid (g%grids)
       end if
    else
       allocate (g%grids(nch), g%weights(nch), g%num_calls(nch))
       call vamp_create_empty_grid (g%grids)
    end if
    read (unit = unit) g%sum_integral
    read (unit = unit) g%sum_weights
    read (unit = unit) g%sum_chi2
    do i = 1, nch
       read (unit = unit) g%weights(i)
    end do
    do i = 1, nch
       read (unit = unit) g%num_calls(i)
    end do
    do i = 1, nch
       call read_grid_raw_unit (g%grids(i), unit, read_integrals)
    end do
    read (unit = unit) magic
    if (magic /= MAGIC_GRIDS_END) then
       print *, FN, " fatal: expecting magic ", MAGIC_GRIDS_END, &
                    " found ", magic
       stop
    end if
  end subroutine read_grids_raw_unit
  subroutine write_grids_raw_name (g, name, write_integrals)
    type(vamp_grids), intent(inout) :: g
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: write_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "write", status = "replace", &
          form = "unformatted", file = name)
    call write_grids_raw_unit (g, unit, write_integrals)
    close (unit = unit)
  end subroutine write_grids_raw_name
  subroutine read_grids_raw_name (g, name, read_integrals)
    type(vamp_grids), intent(inout) :: g
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: read_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "read", status = "old", &
          form = "unformatted", file = name)
    call read_grids_raw_unit (g, unit, read_integrals)
    close (unit = unit)
  end subroutine read_grids_raw_name
  pure subroutine vamp_marshal_grid (g, ibuf, dbuf)
    type(vamp_grid), intent(in) :: g
    integer, dimension(:), intent(inout) :: ibuf
    real(kind=default), dimension(:), intent(inout) :: dbuf
    integer :: i, iwords, dwords, iidx, didx, ndim
    ndim = size (g%div)
    ibuf(1) = g%num_calls
    ibuf(2) = g%calls_per_cell
    ibuf(3) = ndim
    if (g%stratified) then
       ibuf(4) = 1
    else
       ibuf(4) = 0
    end if
    if (g%all_stratified) then
       ibuf(5) = 1
    else
       ibuf(5) = 0
    end if
    if (g%quadrupole) then
       ibuf(6) = 1
    else
       ibuf(6) = 0
    end if
    dbuf(1:2) = g%mu
    dbuf(3) = g%sum_integral
    dbuf(4) = g%sum_weights
    dbuf(5) = g%sum_chi2
    dbuf(6) = g%calls
    dbuf(7) = g%dv2g
    dbuf(8) = g%jacobi
    dbuf(9) = g%f_min
    dbuf(10) = g%f_max
    dbuf(11) = g%mu_gi
    dbuf(12) = g%sum_mu_gi
    ibuf(7:6+ndim) = g%num_div
    iidx = 7 + ndim
    didx = 13
    do i = 1, ndim
       call marshal_division_size (g%div(i), iwords, dwords)
       ibuf(iidx) = iwords
       ibuf(iidx+1) = dwords
       iidx = iidx + 2
       call marshal_division (g%div(i), ibuf(iidx:iidx-1+iwords), &
                                        dbuf(didx:didx-1+dwords))
       iidx = iidx + iwords
       didx = didx + dwords
    end do
    if (associated (g%map)) then
       ibuf(iidx) = 1
       dbuf(didx:didx-1+ndim**2) = reshape (g%map, (/ ndim**2 /))
       didx = didx + ndim**2
    else
       ibuf(iidx) = 0
    end if
    iidx = iidx + 1
    if (associated (g%mu_x)) then
       ibuf(iidx) = 1
       dbuf(didx:didx-1+ndim) = g%mu_x
       didx = didx + ndim
       dbuf(didx:didx-1+ndim) = g%sum_mu_x
       didx = didx + ndim
       dbuf(didx:didx-1+ndim**2) = reshape (g%mu_xx, (/ ndim**2 /))
       didx = didx + ndim**2
       dbuf(didx:didx-1+ndim**2) = reshape (g%sum_mu_xx, (/ ndim**2 /))
       didx = didx + ndim**2
    else
       ibuf(iidx) = 0
    end if
    iidx = iidx + 1
  end subroutine vamp_marshal_grid
  pure subroutine vamp_marshal_grid_size (g, iwords, dwords)
    type(vamp_grid), intent(in) :: g
    integer, intent(out) :: iwords, dwords
    integer :: i, ndim, iw, dw
    ndim = size (g%div)
    iwords = 6 + ndim
    dwords = 12
    do i = 1, ndim
       call marshal_division_size (g%div(i), iw, dw)
       iwords = iwords + 2 + iw
       dwords = dwords + dw
    end do
    iwords = iwords + 1
    if (associated (g%map)) then
       dwords = dwords + ndim**2
    end if
    iwords = iwords + 1
    if (associated (g%mu_x)) then
       dwords = dwords + 2 * (ndim + ndim**2)
    end if
  end subroutine vamp_marshal_grid_size
  pure subroutine vamp_unmarshal_grid (g, ibuf, dbuf)
    type(vamp_grid), intent(inout) :: g
    integer, dimension(:), intent(in) :: ibuf
    real(kind=default), dimension(:), intent(in) :: dbuf
    integer :: i, iwords, dwords, iidx, didx, ndim
    g%num_calls = ibuf(1)
    g%calls_per_cell = ibuf(2)
    ndim = ibuf(3)
    g%stratified = ibuf(4) /= 0
    g%all_stratified = ibuf(5) /= 0
    g%quadrupole = ibuf(6) /= 0
    g%mu = dbuf(1:2)
    g%sum_integral = dbuf(3)
    g%sum_weights = dbuf(4)
    g%sum_chi2 = dbuf(5)
    g%calls = dbuf(6)
    g%dv2g = dbuf(7)
    g%jacobi = dbuf(8)
    g%f_min = dbuf(9)
    g%f_max = dbuf(10)
    g%mu_gi = dbuf(11)
    g%sum_mu_gi = dbuf(12)
    call copy_array_pointer (g%num_div, ibuf(7:6+ndim))
    if (associated (g%div)) then
       if (size (g%div) /= ndim) then
          call delete_division (g%div)
          deallocate (g%div)
          allocate (g%div(ndim))
          call create_empty_division (g%div)
       end if
    else
       allocate (g%div(ndim))
       call create_empty_division (g%div)
    end if
    iidx = 7 + ndim
    didx = 13
    do i = 1, ndim
       iwords = ibuf(iidx)
       dwords = ibuf(iidx+1)
       iidx = iidx + 2
       call unmarshal_division (g%div(i), ibuf(iidx:iidx-1+iwords), &
                                          dbuf(didx:didx-1+dwords))
       iidx = iidx + iwords
       didx = didx + dwords
    end do
    if (ibuf(iidx) > 0) then
       call copy_array_pointer &
            (g%map, reshape (dbuf(didx:didx-1+ibuf(iidx)), (/ ndim, ndim /)))
       didx = didx + ibuf(iidx)
    else
       if (associated (g%map)) then
          deallocate (g%map)
       end if
    end if
    iidx = iidx + 1
    if (ibuf(iidx) > 0) then
       call copy_array_pointer (g%mu_x, dbuf(didx:didx-1+ndim))
       didx = didx + ndim
       call copy_array_pointer (g%sum_mu_x, dbuf(didx:didx-1+ndim))
       didx = didx + ndim
       call copy_array_pointer &
            (g%mu_xx, reshape (dbuf(didx:didx-1+ndim**2), (/ ndim, ndim /)))
       didx = didx + ndim**2
       call copy_array_pointer &
            (g%sum_mu_xx, reshape (dbuf(didx:didx-1+ndim**2), (/ ndim, ndim /)))
       didx = didx + ndim**2
    else
       if (associated (g%mu_x)) then
          deallocate (g%mu_x)
       end if
       if (associated (g%mu_xx)) then
          deallocate (g%mu_xx)
       end if
       if (associated (g%sum_mu_x)) then
          deallocate (g%sum_mu_x)
       end if
       if (associated (g%sum_mu_xx)) then
          deallocate (g%sum_mu_xx)
       end if
    end if
    iidx = iidx + 1
  end subroutine vamp_unmarshal_grid
  pure subroutine vamp_marshal_history (h, ibuf, dbuf)
    type(vamp_history), intent(in) :: h
    integer, dimension(:), intent(inout) :: ibuf
    real(kind=default), dimension(:), intent(inout) :: dbuf
    integer :: j, ndim, iidx, didx, iwords, dwords
    if (h%verbose .and. (associated (h%div))) then
       ndim = size (h%div)
    else
       ndim = 0
    end if
    ibuf(1) = ndim
    ibuf(2) = h%calls
    if (h%stratified) then
       ibuf(3) = 1
    else
       ibuf(3) = 0
    end if
    dbuf(1) = h%integral
    dbuf(2) = h%std_dev
    dbuf(3) = h%avg_integral
    dbuf(4) = h%avg_std_dev
    dbuf(5) = h%avg_chi2
    dbuf(6) = h%f_min
    dbuf(7) = h%f_max
    iidx = 4
    didx = 8
    do j = 1, ndim
       call marshal_div_history_size (h%div(j), iwords, dwords)
       ibuf(iidx) = iwords
       ibuf(iidx+1) = dwords
       iidx = iidx + 2
       call marshal_div_history (h%div(j), ibuf(iidx:iidx-1+iwords), &
                                           dbuf(didx:didx-1+dwords))
       iidx = iidx + iwords
       didx = didx + dwords
    end do
  end subroutine vamp_marshal_history
  pure subroutine vamp_marshal_history_size (h, iwords, dwords)
    type(vamp_history), intent(in) :: h
    integer, intent(out) :: iwords, dwords
    integer :: i, ndim, iw, dw
    if (h%verbose .and. (associated (h%div))) then
       ndim = size (h%div)
    else
       ndim = 0
    end if
    iwords = 3
    dwords = 7
    do i = 1, ndim
       call marshal_div_history_size (h%div(i), iw, dw)
       iwords = iwords + 2 + iw
       dwords = dwords + dw
    end do
  end subroutine vamp_marshal_history_size
  pure subroutine vamp_unmarshal_history (h, ibuf, dbuf)
    type(vamp_history), intent(inout) :: h
    integer, dimension(:), intent(in) :: ibuf
    real(kind=default), dimension(:), intent(in) :: dbuf
    integer :: j, ndim, iidx, didx, iwords, dwords
    ndim = ibuf(1)
    h%calls = ibuf(2)
    h%stratified = ibuf(3) /= 0
    h%integral = dbuf(1)
    h%std_dev = dbuf(2)
    h%avg_integral = dbuf(3)
    h%avg_std_dev = dbuf(4)
    h%avg_chi2 = dbuf(5)
    h%f_min = dbuf(6)
    h%f_max = dbuf(7)
    if (ndim > 0) then
       if (associated (h%div)) then
          if (size (h%div) /= ndim) then
             deallocate (h%div)
             allocate (h%div(ndim))
          end if
       else
          allocate (h%div(ndim))
       end if
       iidx = 4
       didx = 8
       do j = 1, ndim
          iwords = ibuf(iidx)
          dwords = ibuf(iidx+1)
          iidx = iidx + 2
          call unmarshal_div_history (h%div(j), ibuf(iidx:iidx-1+iwords), &
                                                dbuf(didx:didx-1+dwords))
          iidx = iidx + iwords
          didx = didx + dwords
       end do
    end if
  end subroutine vamp_unmarshal_history
  elemental subroutine vamp_copy_grid (lhs, rhs)
    type(vamp_grid), intent(inout) :: lhs
    type(vamp_grid), intent(in) :: rhs
    integer :: ndim
    ndim = size (rhs%div)
    lhs%mu = rhs%mu
    lhs%mu_plus = rhs%mu_plus
    lhs%mu_minus = rhs%mu_minus
    lhs%sum_integral = rhs%sum_integral
    lhs%sum_weights = rhs%sum_weights
    lhs%sum_chi2 = rhs%sum_chi2
    lhs%calls = rhs%calls
    lhs%num_calls = rhs%num_calls
    call copy_array_pointer (lhs%num_div, rhs%num_div)
    lhs%dv2g = rhs%dv2g
    lhs%jacobi = rhs%jacobi
    lhs%f_min = rhs%f_min
    lhs%f_max = rhs%f_max
    lhs%mu_gi = rhs%mu_gi
    lhs%sum_mu_gi = rhs%sum_mu_gi
    lhs%calls_per_cell = rhs%calls_per_cell
    lhs%stratified = rhs%stratified
    lhs%all_stratified = rhs%all_stratified
    lhs%quadrupole = rhs%quadrupole
    if (associated (lhs%div)) then
       if (size (lhs%div) /= ndim) then
          call delete_division (lhs%div)
          deallocate (lhs%div)
          allocate (lhs%div(ndim))
       end if
    else
       allocate (lhs%div(ndim))
    end if
    call copy_division (lhs%div, rhs%div)
    if (associated (rhs%map)) then
       call copy_array_pointer (lhs%map, rhs%map)
    else if (associated (lhs%map)) then
       deallocate (lhs%map)
    end if
    if (associated (rhs%mu_x)) then
       call copy_array_pointer (lhs%mu_x, rhs%mu_x)
       call copy_array_pointer (lhs%mu_xx, rhs%mu_xx)
       call copy_array_pointer (lhs%sum_mu_x, rhs%sum_mu_x)
       call copy_array_pointer (lhs%sum_mu_xx, rhs%sum_mu_xx)
    else if (associated (lhs%mu_x)) then
       deallocate (lhs%mu_x, lhs%mu_xx, lhs%sum_mu_x, lhs%sum_mu_xx)
    end if
  end subroutine vamp_copy_grid
  elemental subroutine vamp_delete_grid (g)
    type(vamp_grid), intent(inout) :: g
    if (associated (g%div)) then
       call delete_division (g%div)
       deallocate (g%div, g%num_div)
    end if
    if (associated (g%map)) then
       deallocate (g%map)
    end if
    if (associated (g%mu_x)) then
       deallocate (g%mu_x, g%mu_xx, g%sum_mu_x, g%sum_mu_xx)
    end if
  end subroutine vamp_delete_grid
  elemental subroutine vamp_copy_grids (lhs, rhs)
    type(vamp_grids), intent(inout) :: lhs
    type(vamp_grids), intent(in) :: rhs
    integer :: nch
    nch = size (rhs%grids)
    lhs%sum_integral = rhs%sum_integral
    lhs%sum_chi2 = rhs%sum_chi2
    lhs%sum_weights = rhs%sum_weights
    if (associated (lhs%grids)) then
       if (size (lhs%grids) /= nch) then
          deallocate (lhs%grids)
          allocate (lhs%grids(nch))
          call vamp_create_empty_grid (lhs%grids(nch))
       end if
    else
       allocate (lhs%grids(nch))
       call vamp_create_empty_grid (lhs%grids(nch))
    end if
    call vamp_copy_grid (lhs%grids, rhs%grids)
    call copy_array_pointer (lhs%weights, rhs%weights)
    call copy_array_pointer (lhs%num_calls, rhs%num_calls)
  end subroutine vamp_copy_grids
  elemental subroutine vamp_delete_grids (g)
    type(vamp_grids), intent(inout) :: g
    if (associated (g%grids)) then
       call vamp_delete_grid (g%grids)
       deallocate (g%weights, g%grids, g%num_calls)
    end if
  end subroutine vamp_delete_grids
    elemental subroutine vamp_copy_history (lhs, rhs)
      type(vamp_history), intent(inout) :: lhs
      type(vamp_history), intent(in) :: rhs
      lhs%calls = rhs%calls
      lhs%stratified = rhs%stratified
      lhs%verbose = rhs%verbose
      lhs%integral = rhs%integral
      lhs%std_dev = rhs%std_dev
      lhs%avg_integral = rhs%avg_integral
      lhs%avg_std_dev = rhs%avg_std_dev
      lhs%avg_chi2 = rhs%avg_chi2
      lhs%f_min = rhs%f_min
      lhs%f_max = rhs%f_max
      if (rhs%verbose) then
         if (associated (lhs%div)) then
            if (size (lhs%div) /= size (rhs%div)) then
               deallocate (lhs%div)
               allocate (lhs%div(size(rhs%div)))
            end if
         else
            allocate (lhs%div(size(rhs%div)))
         end if
         call copy_history (lhs%div, rhs%div)
      end if
    end subroutine vamp_copy_history

    elemental subroutine vamp_delete_history (h)
      type(vamp_history), intent(inout) :: h
      if (associated (h%div)) then
         deallocate (h%div)
      end if
    end subroutine vamp_delete_history
end module vamp_rest
module vamp
  use vamp_grid_type    !NODEP!
  use vamp_rest         !NODEP!
  use vamp_equivalences !NODEP!
  public
end module vamp
