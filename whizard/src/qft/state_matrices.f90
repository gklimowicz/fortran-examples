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

module state_matrices

  use kinds, only: default
  use model_data
  use flavors
  use colors
  use helicities
  use quantum_numbers

  implicit none
  private

  public :: state_matrix_t
  public :: state_iterator_t
  public :: assignment(=)
  public :: merge_state_matrices
  public :: outer_multiply
  public :: state_flv_content_t

  integer, parameter, public :: FM_IGNORE_HELICITY = 1
  integer, parameter, public :: FM_SELECT_HELICITY = 2
  integer, parameter, public :: FM_FACTOR_HELICITY = 3
  integer, parameter, public :: FM_CORRELATED_HELICITY = 4


  type :: node_t
     private
     type(quantum_numbers_t) :: qn
     type(node_t), pointer :: parent => null ()
     type(node_t), pointer :: child_first => null ()
     type(node_t), pointer :: child_last => null ()
     type(node_t), pointer :: next => null ()
     type(node_t), pointer :: previous => null ()
     integer :: me_index = 0
     integer, dimension(:), allocatable :: me_count
     complex(default) :: me = 0
  end type node_t

  type :: state_matrix_t
     private
     type(node_t), pointer :: root => null ()
     integer :: depth = 0
     integer :: n_matrix_elements = 0
     logical :: leaf_nodes_store_values = .false.
     integer :: n_counters = 0
     complex(default), dimension(:), allocatable :: me
     real(default) :: norm = 1
     integer :: n_sub = -1
   contains
     procedure :: init => state_matrix_init
     procedure :: final => state_matrix_final
     procedure :: write => state_matrix_write
     procedure :: write_raw => state_matrix_write_raw
     procedure :: read_raw => state_matrix_read_raw
     procedure :: set_model => state_matrix_set_model
     procedure :: tag_hard_process => state_matrix_tag_hard_process
     procedure :: is_defined => state_matrix_is_defined
     procedure :: is_empty => state_matrix_is_empty
     generic :: get_n_matrix_elements => get_n_matrix_elements_all, get_n_matrix_elements_mask
     procedure :: get_n_matrix_elements_all => state_matrix_get_n_matrix_elements_all
     procedure :: get_n_matrix_elements_mask => state_matrix_get_n_matrix_elements_mask
     procedure :: get_me_size => state_matrix_get_me_size
     procedure :: compute_n_sub => state_matrix_compute_n_sub
     procedure :: set_n_sub => state_matrix_set_n_sub
     procedure :: get_n_sub => state_matrix_get_n_sub
     procedure :: get_n_leaves => state_matrix_get_n_leaves
     procedure :: get_depth => state_matrix_get_depth
     procedure :: get_norm => state_matrix_get_norm
     procedure :: get_quantum_number => &
        state_matrix_get_quantum_number
     generic :: get_quantum_numbers => get_quantum_numbers_all, get_quantum_numbers_mask
     procedure :: get_quantum_numbers_all => state_matrix_get_quantum_numbers_all
     procedure :: get_quantum_numbers_mask => state_matrix_get_quantum_numbers_mask
     procedure :: get_flavors => state_matrix_get_flavors
     generic :: get_matrix_element => get_matrix_element_single
     generic :: get_matrix_element => get_matrix_element_array
     procedure :: get_matrix_element_single => &
       state_matrix_get_matrix_element_single
     procedure :: get_matrix_element_array => &
       state_matrix_get_matrix_element_array
     procedure :: get_max_color_value => state_matrix_get_max_color_value
     procedure :: add_state => state_matrix_add_state
     procedure :: collapse => state_matrix_collapse
     procedure :: reduce => state_matrix_reduce
     procedure :: reorder_me => state_matrix_reorder_me
     procedure :: set_duplicate_flv_zero => state_matrix_set_duplicate_flv_zero
     procedure :: freeze => state_matrix_freeze
     generic :: set_matrix_element => set_matrix_element_qn
     generic :: set_matrix_element => set_matrix_element_all
     generic :: set_matrix_element => set_matrix_element_array
     generic :: set_matrix_element => set_matrix_element_single
     generic :: set_matrix_element => set_matrix_element_clone
     procedure :: set_matrix_element_qn => state_matrix_set_matrix_element_qn
     procedure :: set_matrix_element_all => state_matrix_set_matrix_element_all
     procedure :: set_matrix_element_array => &
          state_matrix_set_matrix_element_array
     procedure :: set_matrix_element_single => &
          state_matrix_set_matrix_element_single
     procedure :: set_matrix_element_clone => &
        state_matrix_set_matrix_element_clone
     procedure :: add_to_matrix_element => state_matrix_add_to_matrix_element
     procedure :: get_diagonal_entries => state_matrix_get_diagonal_entries
     procedure :: renormalize => state_matrix_renormalize
     procedure :: normalize_by_trace => state_matrix_normalize_by_trace
     procedure :: normalize_by_max => state_matrix_normalize_by_max
     procedure :: set_norm => state_matrix_set_norm
     procedure :: sum => state_matrix_sum
     procedure :: trace => state_matrix_trace
     procedure :: add_color_contractions => state_matrix_add_color_contractions
     procedure :: evaluate_product => state_matrix_evaluate_product
     procedure :: evaluate_product_cf => state_matrix_evaluate_product_cf
     procedure :: evaluate_square_c => state_matrix_evaluate_square_c
     procedure :: evaluate_sum => state_matrix_evaluate_sum
     procedure :: evaluate_me_sum => state_matrix_evaluate_me_sum
     procedure :: factorize => state_matrix_factorize
     procedure :: get_polarization_density_matrix &
        => state_matrix_get_polarization_density_matrix
  end type state_matrix_t

  type :: state_iterator_t
     private
     integer :: depth = 0
     type(state_matrix_t), pointer :: state => null ()
     type(node_t), pointer :: node => null ()
   contains
     procedure :: init => state_iterator_init
     procedure :: advance => state_iterator_advance
     procedure :: is_valid => state_iterator_is_valid
     procedure :: get_me_index => state_iterator_get_me_index
     procedure :: get_me_count => state_iterator_get_me_count
     procedure :: get_depth => state_iterator_get_depth
     procedure :: go_to_qn => state_iterator_go_to_qn
     generic :: get_quantum_numbers => get_qn_multi, get_qn_slice, &
          get_qn_range, get_qn_single
     generic :: get_flavor => get_flv_multi, get_flv_slice, &
          get_flv_range, get_flv_single
     generic :: get_color => get_col_multi, get_col_slice, &
          get_col_range, get_col_single
     generic :: get_helicity => get_hel_multi, get_hel_slice, &
          get_hel_range, get_hel_single
     procedure :: get_qn_multi => state_iterator_get_qn_multi
     procedure :: get_qn_slice => state_iterator_get_qn_slice
     procedure :: get_qn_range => state_iterator_get_qn_range
     procedure :: get_qn_single => state_iterator_get_qn_single
     procedure :: get_flv_multi => state_iterator_get_flv_multi
     procedure :: get_flv_slice => state_iterator_get_flv_slice
     procedure :: get_flv_range => state_iterator_get_flv_range
     procedure :: get_flv_single => state_iterator_get_flv_single
     procedure :: get_col_multi => state_iterator_get_col_multi
     procedure :: get_col_slice => state_iterator_get_col_slice
     procedure :: get_col_range => state_iterator_get_col_range
     procedure :: get_col_single => state_iterator_get_col_single
     procedure :: get_hel_multi => state_iterator_get_hel_multi
     procedure :: get_hel_slice => state_iterator_get_hel_slice
     procedure :: get_hel_range => state_iterator_get_hel_range
     procedure :: get_hel_single => state_iterator_get_hel_single
     procedure :: set_model => state_iterator_set_model
     procedure :: retag_hard_process => state_iterator_retag_hard_process
     procedure :: get_matrix_element => state_iterator_get_matrix_element
     procedure :: set_matrix_element => state_iterator_set_matrix_element
     procedure :: add_to_matrix_element => state_iterator_add_to_matrix_element
  end type state_iterator_t

  type :: state_flv_content_t
     private
     integer, dimension(:,:), allocatable :: pdg
     integer, dimension(:,:), allocatable :: map
     logical, dimension(:), allocatable :: mask
   contains
     procedure :: write => state_flv_content_write
     procedure :: init => state_flv_content_init
     procedure :: set_entry => state_flv_content_set_entry
     procedure :: fill => state_flv_content_fill
     procedure :: match => state_flv_content_match
     procedure :: contains => state_flv_content_contains
     procedure :: find_duplicates => state_flv_content_find_duplicates
  end type state_flv_content_t


  interface assignment(=)
     module procedure state_matrix_assign
  end interface

  interface outer_multiply
     module procedure outer_multiply_pair
     module procedure outer_multiply_array
  end interface


  interface
    module subroutine state_matrix_init (state, store_values, n_counters)
      class(state_matrix_t), intent(out) :: state
      logical, intent(in), optional :: store_values
      integer, intent(in), optional :: n_counters
    end subroutine state_matrix_init
    module subroutine state_matrix_final (state)
      class(state_matrix_t), intent(inout) :: state
    end subroutine state_matrix_final
    module subroutine state_matrix_write (state, unit, write_value_list, &
           verbose, col_verbose, testflag)
      class(state_matrix_t), intent(in) :: state
      logical, intent(in), optional :: write_value_list, verbose, col_verbose
      logical, intent(in), optional :: testflag
      integer, intent(in), optional :: unit
    end subroutine state_matrix_write
    module subroutine state_matrix_write_raw (state, u)
      class(state_matrix_t), intent(in), target :: state
      integer, intent(in) :: u
    end subroutine state_matrix_write_raw
    module subroutine state_matrix_read_raw (state, u, iostat)
      class(state_matrix_t), intent(out) :: state
      integer, intent(in) :: u
      integer, intent(out) :: iostat
    end subroutine state_matrix_read_raw
    module subroutine state_matrix_set_model (state, model)
      class(state_matrix_t), intent(inout), target :: state
      class(model_data_t), intent(in), target :: model
    end subroutine state_matrix_set_model
    module subroutine state_matrix_tag_hard_process (state, tagged_state, tag)
      class(state_matrix_t), intent(in), target :: state
      type(state_matrix_t), intent(out) :: tagged_state
      integer, dimension(:), intent(in), optional :: tag
    end subroutine state_matrix_tag_hard_process
    elemental module function state_matrix_is_defined (state) result (defined)
      logical :: defined
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_is_defined
    elemental module function state_matrix_is_empty (state) result (flag)
      logical :: flag
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_is_empty
    pure module function state_matrix_get_n_matrix_elements_all (state) result (n)
      integer :: n
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_get_n_matrix_elements_all
    module function state_matrix_get_n_matrix_elements_mask (state, qn_mask) result (n)
      integer :: n
      class(state_matrix_t), intent(in) :: state
      type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask
    end function state_matrix_get_n_matrix_elements_mask
    pure module function state_matrix_get_me_size (state) result (n)
      integer :: n
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_get_me_size
    module function state_matrix_compute_n_sub (state) result (n_sub)
      integer :: n_sub
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_compute_n_sub
    module subroutine state_matrix_set_n_sub (state)
      class(state_matrix_t), intent(inout) :: state
    end subroutine state_matrix_set_n_sub
    module function state_matrix_get_n_sub (state) result (n_sub)
      integer :: n_sub
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_get_n_sub
    module function state_matrix_get_n_leaves (state) result (n)
      integer :: n
      class(state_matrix_t), intent(in) :: state
      type(state_iterator_t) :: it
    end function state_matrix_get_n_leaves
    pure module function state_matrix_get_depth (state) result (depth)
      integer :: depth
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_get_depth
    pure module function state_matrix_get_norm (state) result (norm)
      real(default) :: norm
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_get_norm
    module function state_matrix_get_quantum_number (state, i, by_me_index) result (qn)
      class(state_matrix_t), intent(in), target :: state
      integer, intent(in) :: i
      logical, intent(in), optional :: by_me_index
      type(quantum_numbers_t), dimension(state%depth) :: qn
    end function state_matrix_get_quantum_number
    module subroutine state_matrix_get_quantum_numbers_all (state, qn)
      class(state_matrix_t), intent(in), target :: state
      type(quantum_numbers_t), intent(out), dimension(:,:), allocatable :: qn
    end subroutine state_matrix_get_quantum_numbers_all
    module subroutine state_matrix_get_quantum_numbers_mask (state, qn_mask, qn)
      class(state_matrix_t), intent(in), target :: state
      type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask
      type(quantum_numbers_t), intent(out), dimension(:,:), allocatable :: qn
    end subroutine state_matrix_get_quantum_numbers_mask
    module subroutine state_matrix_get_flavors (state, only_elementary, qn_mask, flv)
      class(state_matrix_t), intent(in), target :: state
      logical, intent(in) :: only_elementary
      type(quantum_numbers_mask_t), intent(in), dimension(:), optional :: qn_mask
      integer, intent(out), dimension(:,:), allocatable :: flv
    end subroutine state_matrix_get_flavors
    elemental module function state_matrix_get_matrix_element_single (state, i) result (me)
      complex(default) :: me
      class(state_matrix_t), intent(in) :: state
      integer, intent(in) :: i
    end function state_matrix_get_matrix_element_single
    module function state_matrix_get_matrix_element_array (state) result (me)
       complex(default), dimension(:), allocatable :: me
       class(state_matrix_t), intent(in) :: state
    end function state_matrix_get_matrix_element_array
    module function state_matrix_get_max_color_value (state) result (cmax)
      integer :: cmax
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_get_max_color_value
    module subroutine state_matrix_add_state (state, qn, index, value, &
           sum_values, counter_index, ignore_sub_for_qn, me_index)
      class(state_matrix_t), intent(inout) :: state
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      integer, intent(in), optional :: index
      complex(default), intent(in), optional :: value
      logical, intent(in), optional :: sum_values
      integer, intent(in), optional :: counter_index
      logical, intent(in), optional :: ignore_sub_for_qn
      integer, intent(out), optional :: me_index
    end subroutine state_matrix_add_state
    module subroutine state_matrix_collapse (state, mask)
      class(state_matrix_t), intent(inout) :: state
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
    end subroutine state_matrix_collapse
    module subroutine state_matrix_reduce (state, mask, red_state, keep_me_index)
      class(state_matrix_t), intent(in), target :: state
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
      type(state_matrix_t), intent(out) :: red_state
      logical, optional, intent(in)  :: keep_me_index
    end subroutine state_matrix_reduce
    module subroutine state_matrix_reorder_me (state, ordered_state)
      class(state_matrix_t), intent(in), target :: state
      type(state_matrix_t), intent(out) :: ordered_state
    end subroutine state_matrix_reorder_me
    module subroutine state_matrix_set_duplicate_flv_zero (state)
      class(state_matrix_t), intent(inout), target :: state
    end subroutine state_matrix_set_duplicate_flv_zero
    module subroutine state_matrix_freeze (state)
      class(state_matrix_t), intent(inout), target :: state
    end subroutine state_matrix_freeze
    module subroutine state_matrix_set_matrix_element_qn (state, qn, value)
      class(state_matrix_t), intent(inout), target :: state
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      complex(default), intent(in) :: value
    end subroutine state_matrix_set_matrix_element_qn
    module subroutine state_matrix_set_matrix_element_all (state, value)
      class(state_matrix_t), intent(inout) :: state
      complex(default), intent(in) :: value
    end subroutine state_matrix_set_matrix_element_all
    module subroutine state_matrix_set_matrix_element_array (state, value, range)
      class(state_matrix_t), intent(inout) :: state
      complex(default), intent(in), dimension(:) :: value
      integer, intent(in), dimension(:), optional :: range
    end subroutine state_matrix_set_matrix_element_array
    pure module subroutine state_matrix_set_matrix_element_single (state, i, value)
      class(state_matrix_t), intent(inout) :: state
      integer, intent(in) :: i
      complex(default), intent(in) :: value
    end subroutine state_matrix_set_matrix_element_single
    module subroutine state_matrix_set_matrix_element_clone (state, state1)
      class(state_matrix_t), intent(inout) :: state
      type(state_matrix_t), intent(in) :: state1
    end subroutine state_matrix_set_matrix_element_clone
    module subroutine state_matrix_add_to_matrix_element (state, qn, value, match_only_flavor)
      class(state_matrix_t), intent(inout), target :: state
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      complex(default), intent(in) :: value
      logical, intent(in), optional :: match_only_flavor
    end subroutine state_matrix_add_to_matrix_element
    module subroutine state_iterator_init (it, state)
      class(state_iterator_t), intent(out) :: it
      type(state_matrix_t), intent(in), target :: state
    end subroutine state_iterator_init
    module subroutine state_iterator_advance (it)
      class(state_iterator_t), intent(inout) :: it
    end subroutine state_iterator_advance
    module function state_iterator_is_valid (it) result (defined)
      logical :: defined
      class(state_iterator_t), intent(in) :: it
    end function state_iterator_is_valid
    module function state_iterator_get_me_index (it) result (n)
      integer :: n
      class(state_iterator_t), intent(in) :: it
    end function state_iterator_get_me_index
    module function state_iterator_get_me_count (it) result (n)
      integer, dimension(:), allocatable :: n
      class(state_iterator_t), intent(in) :: it
    end function state_iterator_get_me_count
    pure module function state_iterator_get_depth (state_iterator) result (depth)
      integer :: depth
      class(state_iterator_t), intent(in) :: state_iterator
    end function state_iterator_get_depth
    module subroutine state_iterator_go_to_qn (it, qn, match_only_flavor)
      class(state_iterator_t), intent(inout) :: it
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      logical, intent(in), optional :: match_only_flavor
    end subroutine state_iterator_go_to_qn
    module function state_iterator_get_qn_multi (it) result (qn)
      class(state_iterator_t), intent(in) :: it
      type(quantum_numbers_t), dimension(it%depth) :: qn
    end function state_iterator_get_qn_multi
    module function state_iterator_get_flv_multi (it) result (flv)
      class(state_iterator_t), intent(in) :: it
      type(flavor_t), dimension(it%depth) :: flv
    end function state_iterator_get_flv_multi
    module function state_iterator_get_col_multi (it) result (col)
      class(state_iterator_t), intent(in) :: it
      type(color_t), dimension(it%depth) :: col
    end function state_iterator_get_col_multi
    module function state_iterator_get_hel_multi (it) result (hel)
      class(state_iterator_t), intent(in) :: it
      type(helicity_t), dimension(it%depth) :: hel
    end function state_iterator_get_hel_multi
    module function state_iterator_get_qn_slice (it, index) result (qn)
      class(state_iterator_t), intent(in) :: it
      integer, dimension(:), intent(in) :: index
      type(quantum_numbers_t), dimension(size(index)) :: qn
    end function state_iterator_get_qn_slice
    module function state_iterator_get_flv_slice (it, index) result (flv)
      class(state_iterator_t), intent(in) :: it
      integer, dimension(:), intent(in) :: index
      type(flavor_t), dimension(size(index)) :: flv
    end function state_iterator_get_flv_slice
    module function state_iterator_get_col_slice (it, index) result (col)
      class(state_iterator_t), intent(in) :: it
      integer, dimension(:), intent(in) :: index
      type(color_t), dimension(size(index)) :: col
    end function state_iterator_get_col_slice
    module function state_iterator_get_hel_slice (it, index) result (hel)
      class(state_iterator_t), intent(in) :: it
      integer, dimension(:), intent(in) :: index
      type(helicity_t), dimension(size(index)) :: hel
    end function state_iterator_get_hel_slice
    module function state_iterator_get_qn_range (it, k1, k2) result (qn)
      class(state_iterator_t), intent(in) :: it
      integer, intent(in) :: k1, k2
      type(quantum_numbers_t), dimension(k2-k1+1) :: qn
    end function state_iterator_get_qn_range
    module function state_iterator_get_flv_range (it, k1, k2) result (flv)
      class(state_iterator_t), intent(in) :: it
      integer, intent(in) :: k1, k2
      type(flavor_t), dimension(k2-k1+1) :: flv
    end function state_iterator_get_flv_range
    module function state_iterator_get_col_range (it, k1, k2) result (col)
      class(state_iterator_t), intent(in) :: it
      integer, intent(in) :: k1, k2
      type(color_t), dimension(k2-k1+1) :: col
    end function state_iterator_get_col_range
    module function state_iterator_get_hel_range (it, k1, k2) result (hel)
      class(state_iterator_t), intent(in) :: it
      integer, intent(in) :: k1, k2
      type(helicity_t), dimension(k2-k1+1) :: hel
    end function state_iterator_get_hel_range
    module function state_iterator_get_qn_single (it, k) result (qn)
      class(state_iterator_t), intent(in) :: it
      integer, intent(in) :: k
      type(quantum_numbers_t) :: qn
    end function state_iterator_get_qn_single
    module function state_iterator_get_flv_single (it, k) result (flv)
      class(state_iterator_t), intent(in) :: it
      integer, intent(in) :: k
      type(flavor_t) :: flv
    end function state_iterator_get_flv_single
    module function state_iterator_get_col_single (it, k) result (col)
      class(state_iterator_t), intent(in) :: it
      integer, intent(in) :: k
      type(color_t) :: col
    end function state_iterator_get_col_single
    module function state_iterator_get_hel_single (it, k) result (hel)
      class(state_iterator_t), intent(in) :: it
      integer, intent(in) :: k
      type(helicity_t) :: hel
    end function state_iterator_get_hel_single
    module subroutine state_iterator_set_model (it, model)
      class(state_iterator_t), intent(inout) :: it
      class(model_data_t), intent(in), target :: model
    end subroutine state_iterator_set_model
    module subroutine state_iterator_retag_hard_process (it, i, hard)
      class(state_iterator_t), intent(inout) :: it
      integer, intent(in) :: i
      logical, intent(in) :: hard
    end subroutine state_iterator_retag_hard_process
    module function state_iterator_get_matrix_element (it) result (me)
      complex(default) :: me
      class(state_iterator_t), intent(in) :: it
    end function state_iterator_get_matrix_element
    module subroutine state_iterator_set_matrix_element (it, value)
      class(state_iterator_t), intent(inout) :: it
      complex(default), intent(in) :: value
    end subroutine state_iterator_set_matrix_element
    module subroutine state_iterator_add_to_matrix_element (it, value)
      class(state_iterator_t), intent(inout) :: it
      complex(default), intent(in) :: value
    end subroutine state_iterator_add_to_matrix_element
    module subroutine state_matrix_assign (state_out, state_in)
      type(state_matrix_t), intent(out) :: state_out
      type(state_matrix_t), intent(in), target :: state_in
    end subroutine state_matrix_assign
    module subroutine state_matrix_get_diagonal_entries (state, i)
      class(state_matrix_t), intent(in) :: state
      integer, dimension(:), allocatable, intent(out) :: i
    end subroutine state_matrix_get_diagonal_entries
    module subroutine state_matrix_renormalize (state, factor)
      class(state_matrix_t), intent(inout) :: state
      complex(default), intent(in) :: factor
    end subroutine state_matrix_renormalize
    module subroutine state_matrix_normalize_by_trace (state)
      class(state_matrix_t), intent(inout) :: state
    end subroutine state_matrix_normalize_by_trace
    module subroutine state_matrix_normalize_by_max (state)
      class(state_matrix_t), intent(inout) :: state
    end subroutine state_matrix_normalize_by_max
    module subroutine state_matrix_set_norm (state, norm)
      class(state_matrix_t), intent(inout) :: state
      real(default), intent(in) :: norm
    end subroutine state_matrix_set_norm
    pure module function state_matrix_sum (state) result (value)
      complex(default) :: value
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_sum
    module function state_matrix_trace (state, qn_in) result (trace)
      complex(default) :: trace
      class(state_matrix_t), intent(in), target :: state
      type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_in
    end function state_matrix_trace
    module subroutine state_matrix_add_color_contractions (state)
      class(state_matrix_t), intent(inout), target :: state
    end subroutine state_matrix_add_color_contractions
    module subroutine merge_state_matrices (state1, state2, state3)
      type(state_matrix_t), intent(in), target :: state1, state2
      type(state_matrix_t), intent(out) :: state3
    end subroutine merge_state_matrices
    pure module subroutine state_matrix_evaluate_product &
         (state, i, state1, state2, index1, index2)
      class(state_matrix_t), intent(inout) :: state
      integer, intent(in) :: i
      type(state_matrix_t), intent(in) :: state1, state2
      integer, dimension(:), intent(in) :: index1, index2
    end subroutine state_matrix_evaluate_product
    pure module subroutine state_matrix_evaluate_product_cf &
         (state, i, state1, state2, index1, index2, factor)
      class(state_matrix_t), intent(inout) :: state
      integer, intent(in) :: i
      type(state_matrix_t), intent(in) :: state1, state2
      integer, dimension(:), intent(in) :: index1, index2
      complex(default), dimension(:), intent(in) :: factor
    end subroutine state_matrix_evaluate_product_cf
    pure module subroutine state_matrix_evaluate_square_c (state, i, state1, index1)
      class(state_matrix_t), intent(inout) :: state
      integer, intent(in) :: i
      type(state_matrix_t), intent(in) :: state1
      integer, dimension(:), intent(in) :: index1
    end subroutine state_matrix_evaluate_square_c
    pure module subroutine state_matrix_evaluate_sum (state, i, state1, index1)
      class(state_matrix_t), intent(inout) :: state
      integer, intent(in) :: i
      type(state_matrix_t), intent(in) :: state1
      integer, dimension(:), intent(in) :: index1
    end subroutine state_matrix_evaluate_sum
    pure module subroutine state_matrix_evaluate_me_sum (state, i, state1, index1)
      class(state_matrix_t), intent(inout) :: state
      integer, intent(in) :: i
      type(state_matrix_t), intent(in) :: state1
      integer, dimension(:), intent(in) :: index1
    end subroutine state_matrix_evaluate_me_sum
    module subroutine outer_multiply_pair (state1, state2, state3)
      type(state_matrix_t), intent(in), target :: state1, state2
      type(state_matrix_t), intent(out) :: state3
    end subroutine outer_multiply_pair
    module subroutine outer_multiply_array (state_in, state_out)
      type(state_matrix_t), dimension(:), intent(in), target :: state_in
      type(state_matrix_t), intent(out) :: state_out
    end subroutine outer_multiply_array
    module subroutine state_matrix_factorize &
         (state, mode, x, ok, single_state, correlated_state, qn_in)
      class(state_matrix_t), intent(in), target :: state
      integer, intent(in) :: mode
      real(default), intent(in) :: x
      logical, intent(out) :: ok
      type(state_matrix_t), &
           dimension(:), allocatable, intent(out) :: single_state
      type(state_matrix_t), intent(out), optional :: correlated_state
      type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_in
    end subroutine state_matrix_factorize
    module function state_matrix_get_polarization_density_matrix &
         (state) result (pol_matrix)
      real(default), dimension(:,:), allocatable :: pol_matrix
      class(state_matrix_t), intent(in) :: state
    end function state_matrix_get_polarization_density_matrix
    module subroutine state_flv_content_write (state_flv, unit)
      class(state_flv_content_t), intent(in), target :: state_flv
      integer, intent(in), optional :: unit
    end subroutine state_flv_content_write
    module subroutine state_flv_content_init (state_flv, n, mask)
      class(state_flv_content_t), intent(out) :: state_flv
      integer, intent(in) :: n
      logical, dimension(:), intent(in) :: mask
    end subroutine state_flv_content_init
    module subroutine state_flv_content_set_entry (state_flv, i, pdg, map)
      class(state_flv_content_t), intent(inout) :: state_flv
      integer, intent(in) :: i
      integer, dimension(:), intent(in) :: pdg, map
    end subroutine state_flv_content_set_entry
    module subroutine state_flv_content_fill &
         (state_flv, state_full, mask)
      class(state_flv_content_t), intent(out) :: state_flv
      type(state_matrix_t), intent(in), target :: state_full
      logical, dimension(:), intent(in) :: mask
    end subroutine state_flv_content_fill
    module subroutine state_flv_content_match (state_flv, pdg, success, map)
      class(state_flv_content_t), intent(in) :: state_flv
      integer, dimension(:), intent(in) :: pdg
      logical, intent(out) :: success
      integer, dimension(:), intent(out) :: map
    end subroutine state_flv_content_match
    module function state_flv_content_contains (state_flv, pdg) result (success)
      class(state_flv_content_t), intent(in) :: state_flv
      integer, intent(in) :: pdg
      logical :: success
    end function state_flv_content_contains
    module subroutine state_flv_content_find_duplicates (state_flv, duplicate_mask)
      class(state_flv_content_t), intent(in) :: state_flv
      logical, dimension(:), allocatable, intent(out) :: duplicate_mask
    end subroutine state_flv_content_find_duplicates
  end interface

end module state_matrices
