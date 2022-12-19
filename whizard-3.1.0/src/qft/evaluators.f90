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

module evaluators

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use state_matrices
  use interactions

  implicit none
  private

  public :: pairing_array_write
  public :: evaluator_t
  public :: assignment(=)
  public :: evaluator_reassign_links
  public :: evaluator_get_int_in_ptr

  integer, parameter :: &
       EVAL_UNDEFINED = 0, &
       EVAL_PRODUCT = 1, &
       EVAL_SQUARED_FLOWS = 2, &
       EVAL_SQUARE_WITH_COLOR_FACTORS = 3, &
       EVAL_COLOR_CONTRACTION = 4, &
       EVAL_IDENTITY = 5, &
       EVAL_QN_SUM = 6

  type :: pairing_array_t
     integer, dimension(:), allocatable :: i1, i2
     complex(default), dimension(:), allocatable :: factor
  end type pairing_array_t

  type, extends (interaction_t) :: evaluator_t
     private
     integer :: type = EVAL_UNDEFINED
     class(interaction_t), pointer :: int_in1 => null ()
     class(interaction_t), pointer :: int_in2 => null ()
     type(pairing_array_t), dimension(:), allocatable :: pairing_array
   contains
     procedure :: write => evaluator_write
     procedure :: write_pairing_array => evaluator_write_pairing_array
     procedure :: init_product => evaluator_init_product
     procedure :: init_square => evaluator_init_square
     procedure :: init_square_diag => evaluator_init_square_diag
     procedure :: init_square_nondiag => evaluator_init_square_nondiag
     procedure :: init_color_contractions => evaluator_init_color_contractions
     procedure :: init_identity => evaluator_init_identity
     procedure :: init_qn_sum => evaluator_init_qn_sum
     procedure :: evaluate => evaluator_evaluate
  end type evaluator_t

  type :: index_map_t
     integer, dimension(:), allocatable :: entry
  end type index_map_t

  type :: index_map2_t
     integer :: s = 0
     integer, dimension(:,:), allocatable :: entry
  end type index_map2_t

  type :: prt_mask_t
     logical, dimension(:), allocatable :: entry
  end type prt_mask_t

  type :: qn_list_t
     type(quantum_numbers_t), dimension(:,:), allocatable :: qn
  end type qn_list_t

  type :: qn_mask_array_t
     type(quantum_numbers_mask_t), dimension(:), allocatable :: mask
  end type qn_mask_array_t

  type :: connection_entry_t
     type(quantum_numbers_t), dimension(:), allocatable :: qn_conn
     integer, dimension(:), allocatable :: n_index
     integer, dimension(:), allocatable :: count
     type(index_map_t), dimension(:), allocatable :: index_in
     type(qn_list_t), dimension(:), allocatable :: qn_in_list
  end type connection_entry_t

  type :: connection_table_t
     integer :: n_conn = 0
     integer, dimension(2) :: n_rest = 0
     integer :: n_tot = 0
     integer :: n_me_conn = 0
     type(state_matrix_t) :: state
     type(index_map_t), dimension(:), allocatable :: index_conn
     type(connection_entry_t), dimension(:), allocatable :: entry
     type(index_map_t) :: index_result
  end type connection_table_t

  type :: connection_table_diag_t
    integer :: n_tot = 0
    integer :: n_me_conn = 0
    type(state_matrix_t) :: state
    type(index_map_t) :: index_conn
    type(connection_entry_t), dimension(:), allocatable :: entry
    type(index_map_t) :: index_result
  end type connection_table_diag_t

  type :: connection_table_nondiag_t
    integer :: n_tot = 0
    integer :: n_me_conn = 0
    type(state_matrix_t) :: state
    type(index_map2_t) :: index_conn
    type(connection_entry_t), dimension(:), allocatable :: entry
    type(index_map_t) :: index_result
  end type connection_table_nondiag_t

  type :: color_table_t
     integer, dimension(:), allocatable :: index
     type(color_t), dimension(:,:), allocatable :: col
     logical, dimension(:,:), allocatable :: factor_is_known
     complex(default), dimension(:,:), allocatable :: factor
  end type color_table_t


  interface assignment(=)
     module procedure evaluator_assign
  end interface

  interface size
     module procedure index_map_size
  end interface

  interface assignment(=)
     module procedure index_map_assign_int
     module procedure index_map_assign_array
  end interface

  interface size
     module procedure index_map2_size
  end interface

  interface assignment(=)
     module procedure index_map2_assign_int
  end interface

  interface size
     module procedure prt_mask_size
  end interface

  interface evaluator_reassign_links
     module procedure evaluator_reassign_links_eval
     module procedure evaluator_reassign_links_int
  end interface


  interface
    elemental module subroutine pairing_array_init (pa, n, has_i2, has_factor)
      type(pairing_array_t), intent(out) :: pa
      integer, intent(in) :: n
      logical, intent(in) :: has_i2, has_factor
    end subroutine pairing_array_init
    module subroutine pairing_array_write (pa, unit)
      type(pairing_array_t), intent(in) :: pa
      integer, intent(in), optional :: unit
    end subroutine pairing_array_write
    module subroutine evaluator_write (eval, unit, &
         verbose, show_momentum_sum, show_mass, show_state, show_table, &
         col_verbose, testflag)
      class(evaluator_t), intent(in) :: eval
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, show_momentum_sum, show_mass
      logical, intent(in), optional :: show_state, show_table, col_verbose
      logical, intent(in), optional :: testflag
    end subroutine evaluator_write
    module subroutine evaluator_write_pairing_array (eval, conjugate, square, unit)
      class(evaluator_t), intent(in) :: eval
      logical, intent(in) :: conjugate, square
      integer, intent(in), optional :: unit
    end subroutine evaluator_write_pairing_array
    module subroutine evaluator_assign (eval_out, eval_in)
      type(evaluator_t), intent(out) :: eval_out
      type(evaluator_t), intent(in) :: eval_in
    end subroutine evaluator_assign
    elemental module subroutine index_map_init (map, n)
      type(index_map_t), intent(out) :: map
      integer, intent(in) :: n
    end subroutine index_map_init
    module function index_map_size (map) result (s)
      integer :: s
      type(index_map_t), intent(in) :: map
    end function index_map_size
    elemental module subroutine index_map_assign_int (map, ival)
      type(index_map_t), intent(inout) :: map
      integer, intent(in) :: ival
    end subroutine index_map_assign_int
    module subroutine index_map_assign_array (map, array)
      type(index_map_t), intent(inout) :: map
      integer, dimension(:), intent(in) :: array
    end subroutine index_map_assign_array
    elemental module subroutine index_map_set_entry (map, i, ival)
      type(index_map_t), intent(inout) :: map
      integer, intent(in) :: i
      integer, intent(in) :: ival
    end subroutine index_map_set_entry
    elemental module function index_map_get_entry (map, i) result (ival)
      integer :: ival
      type(index_map_t), intent(in) :: map
      integer, intent(in) :: i
    end function index_map_get_entry
    module function index_map2_size (map) result (s)
      integer :: s
      type(index_map2_t), intent(in) :: map
    end function index_map2_size
    elemental module subroutine index_map2_assign_int (map, ival)
      type(index_map2_t), intent(inout) :: map
      integer, intent(in) :: ival
    end subroutine index_map2_assign_int
    module function prt_mask_size (mask) result (s)
      integer :: s
      type(prt_mask_t), intent(in) :: mask
    end function prt_mask_size
    module subroutine evaluator_init_product &
         (eval, int_in1, int_in2, qn_mask_conn, qn_filter_conn, qn_mask_rest, &
          connections_are_resonant, ignore_sub_for_qn)
      class(evaluator_t), intent(out), target :: eval
      class(interaction_t), intent(in), target :: int_in1, int_in2
      type(quantum_numbers_mask_t), intent(in) :: qn_mask_conn
      type(quantum_numbers_t), intent(in), optional :: qn_filter_conn
      type(quantum_numbers_mask_t), intent(in), optional :: qn_mask_rest
      logical, intent(in), optional :: connections_are_resonant
      logical, intent(in), optional :: ignore_sub_for_qn
    end subroutine evaluator_init_product
    module subroutine evaluator_init_square (eval, int_in, qn_mask, &
         col_flow_index, col_factor, col_index_hi, expand_color_flows, nc)
      class(evaluator_t), intent(out), target :: eval
      class(interaction_t), intent(in), target :: int_in
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      integer, dimension(:,:), intent(in), optional :: col_flow_index
      complex(default), dimension(:), intent(in), optional :: col_factor
      integer, dimension(:), intent(in), optional :: col_index_hi
      logical, intent(in), optional :: expand_color_flows
      integer, intent(in), optional :: nc
    end subroutine evaluator_init_square
    module subroutine evaluator_init_square_diag (eval, int_in, qn_mask, &
         col_flow_index, col_factor, col_index_hi, expand_color_flows, nc)
      class(evaluator_t), intent(out), target :: eval
      class(interaction_t), intent(in), target :: int_in
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      integer, dimension(:,:), intent(in), optional :: col_flow_index
      complex(default), dimension(:), intent(in), optional :: col_factor
      integer, dimension(:), intent(in), optional :: col_index_hi
      logical, intent(in), optional :: expand_color_flows
      integer, intent(in), optional :: nc
    end subroutine evaluator_init_square_diag
    module subroutine evaluator_init_square_nondiag (eval, int_in, qn_mask, &
         col_flow_index, col_factor, col_index_hi, expand_color_flows, nc)
      class(evaluator_t), intent(out), target :: eval
      class(interaction_t), intent(in), target :: int_in
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      integer, dimension(:,:), intent(in), optional :: col_flow_index
      complex(default), dimension(:), intent(in), optional :: col_factor
      integer, dimension(:), intent(in), optional :: col_index_hi
      logical, intent(in), optional :: expand_color_flows
      integer, intent(in), optional :: nc
    end subroutine evaluator_init_square_nondiag
    module subroutine evaluator_init_color_contractions (eval, int_in)
      class(evaluator_t), intent(out), target :: eval
      type(interaction_t), intent(in), target :: int_in
    end subroutine evaluator_init_color_contractions
    module subroutine evaluator_reassign_links_eval (eval, eval_src, eval_target)
      type(evaluator_t), intent(inout) :: eval
      type(evaluator_t), intent(in) :: eval_src
      type(evaluator_t), intent(in), target :: eval_target
    end subroutine evaluator_reassign_links_eval
    module subroutine evaluator_reassign_links_int (eval, int_src, int_target)
      type(evaluator_t), intent(inout) :: eval
      type(interaction_t), intent(in) :: int_src
      type(interaction_t), intent(in), target :: int_target
    end subroutine evaluator_reassign_links_int
    module function evaluator_get_int_in_ptr (eval, i) result (int_in)
      class(interaction_t), pointer :: int_in
      type(evaluator_t), intent(in), target :: eval
      integer, intent(in) :: i
    end function evaluator_get_int_in_ptr
    module subroutine evaluator_init_identity (eval, int)
      class(evaluator_t), intent(out), target :: eval
      class(interaction_t), intent(in), target :: int
    end subroutine evaluator_init_identity
    module subroutine evaluator_init_qn_sum (eval, int, qn_mask, drop)
      class(evaluator_t), intent(out), target :: eval
      class(interaction_t), target, intent(in) :: int
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      logical, intent(in), optional, dimension(:) :: drop
    end subroutine evaluator_init_qn_sum
    module subroutine evaluator_evaluate (eval)
      class(evaluator_t), intent(inout), target :: eval
    end subroutine evaluator_evaluate
  end interface

end module evaluators
