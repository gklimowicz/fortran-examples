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

module interactions

  use kinds, only: default
  use lorentz
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use state_matrices

  implicit none
  private

  public :: qn_index_map_t
  public :: external_link_get_ptr
  public :: external_link_get_index
  public :: interaction_t
  public :: reset_interaction_counter
  public :: assignment(=)
  public :: interaction_reassign_links
  public :: interaction_find_link
  public :: find_connections

  type :: qn_index_map_t
     private
     type(quantum_numbers_t), dimension(:, :), allocatable :: qn_flv
     type(quantum_numbers_t), dimension(:, :), allocatable :: qn_hel
     logical :: flip_hel = .false.
     integer :: n_flv = 0, n_hel = 0, n_sub = 0
     integer, dimension(:, :, :), allocatable :: index
     integer, dimension(:,:), allocatable :: sf_index_born, sf_index_real
   contains
     generic :: init => init_trivial, &
                        init_involved, &
                        init_sf
     procedure, private :: init_trivial => qn_index_map_init_trivial
     procedure, private :: init_involved => qn_index_map_init_involved
     procedure, private :: init_sf => qn_index_map_init_sf
     procedure :: write => qn_index_map_write
     procedure :: set_helicity_flip => qn_index_map_set_helicity_flip
     procedure :: get_index => qn_index_map_get_index
     procedure :: get_n_flv => qn_index_map_get_n_flv
     procedure :: get_n_hel => qn_index_map_get_n_hel
     procedure :: get_n_sub => qn_index_map_get_n_sub
     procedure :: get_index_by_qn => qn_index_map_get_index_by_qn
     procedure :: get_sf_index_born => qn_index_map_get_sf_index_born
     procedure :: get_sf_index_real => qn_index_map_get_sf_index_real
  end type qn_index_map_t

  type :: external_link_t
     private
     type(interaction_t), pointer :: int => null ()
     integer :: i
  end type external_link_t

  type :: internal_link_list_t
     private
     integer :: length = 0
     integer, dimension(:), allocatable :: link
   contains
     procedure :: write => internal_link_list_write
     procedure :: append => internal_link_list_append
     procedure :: has_entries => internal_link_list_has_entries
     procedure :: get_length => internal_link_list_get_length
     procedure :: get_link => internal_link_list_get_link
  end type internal_link_list_t

  type :: interaction_t
     private
     integer :: tag = 0
     type(state_matrix_t) :: state_matrix
     integer :: n_in = 0
     integer :: n_vir = 0
     integer :: n_out = 0
     integer :: n_tot = 0
     logical, dimension(:), allocatable :: p_is_known
     type(vector4_t), dimension(:), allocatable :: p
     type(external_link_t), dimension(:), allocatable :: source
     type(internal_link_list_t), dimension(:), allocatable :: parents
     type(internal_link_list_t), dimension(:), allocatable :: children
     logical, dimension(:), allocatable :: resonant
     type(quantum_numbers_mask_t), dimension(:), allocatable :: mask
     integer, dimension(:), allocatable :: hel_lock
     logical :: update_state_matrix = .false.
     logical :: update_values = .false.
     type(qn_index_map_t) :: qn_index
   contains
     procedure :: basic_init => interaction_init
     generic :: init_qn_index => init_qn_index_trivial, &
                                 init_qn_index_involved, &
                                 init_qn_index_sf
     procedure :: init_qn_index_trivial => interaction_init_qn_index_trivial
     procedure :: init_qn_index_involved => interaction_init_qn_index_involved
     procedure :: init_qn_index_sf => interaction_init_qn_index_sf
     procedure :: set_qn_index_helicity_flip => interaction_set_qn_index_helicity_flip
     procedure :: get_qn_index => interaction_get_qn_index
     procedure :: get_sf_qn_index_born => interaction_get_sf_qn_index_born
     procedure :: get_sf_qn_index_real => interaction_get_sf_qn_index_real
     procedure :: get_qn_index_n_flv => interaction_get_qn_index_n_flv
     procedure :: get_qn_index_n_hel => interaction_get_qn_index_n_hel
     procedure :: get_qn_index_n_sub => interaction_get_qn_index_n_sub
     procedure :: final => interaction_final
     procedure :: basic_write => interaction_write
     procedure :: write_state_matrix => interaction_write_state_matrix
     procedure :: reduce_state_matrix => interaction_reduce_state_matrix
     procedure :: add_state => interaction_add_state
     procedure :: set_duplicate_flv_zero => interaction_set_duplicate_flv_zero
     procedure :: freeze => interaction_freeze
     procedure :: is_empty => interaction_is_empty
     procedure :: get_n_matrix_elements => &
          interaction_get_n_matrix_elements
     procedure :: get_state_depth => interaction_get_state_depth
     procedure :: get_n_in_helicities => interaction_get_n_in_helicities
     procedure :: get_me_size => interaction_get_me_size
     procedure :: get_norm => interaction_get_norm
     procedure :: get_n_sub => interaction_get_n_sub
     generic :: get_quantum_numbers => get_quantum_numbers_single, &
                                       get_quantum_numbers_all, &
                                       get_quantum_numbers_all_qn_mask
     procedure :: get_quantum_numbers_single => &
        interaction_get_quantum_numbers_single
     procedure :: get_quantum_numbers_all => &
        interaction_get_quantum_numbers_all
     procedure :: get_quantum_numbers_all_qn_mask => &
        interaction_get_quantum_numbers_all_qn_mask
     procedure :: get_quantum_numbers_all_sub => interaction_get_quantum_numbers_all_sub
     procedure :: get_flavors => interaction_get_flavors
     procedure :: get_quantum_numbers_mask => interaction_get_quantum_numbers_mask
     generic :: get_matrix_element => get_matrix_element_single
     generic :: get_matrix_element => get_matrix_element_array
     procedure :: get_matrix_element_single => &
        interaction_get_matrix_element_single
     procedure :: get_matrix_element_array => &
        interaction_get_matrix_element_array
     generic :: set_matrix_element => interaction_set_matrix_element_qn, &
          interaction_set_matrix_element_all, &
          interaction_set_matrix_element_array, &
          interaction_set_matrix_element_single, &
          interaction_set_matrix_element_clone
     procedure :: interaction_set_matrix_element_qn
     procedure :: interaction_set_matrix_element_all
     procedure :: interaction_set_matrix_element_array
     procedure :: interaction_set_matrix_element_single
     procedure :: interaction_set_matrix_element_clone
     procedure :: set_only_matrix_element => interaction_set_only_matrix_element
     procedure :: add_to_matrix_element => interaction_add_to_matrix_element
     procedure :: get_diagonal_entries => interaction_get_diagonal_entries
     procedure :: normalize_by_trace => interaction_normalize_by_trace
     procedure :: normalize_by_max => interaction_normalize_by_max
     procedure :: set_norm => interaction_set_norm
     procedure :: set_state_matrix => interaction_set_state_matrix
     procedure :: get_max_color_value => &
          interaction_get_max_color_value
     procedure :: factorize => interaction_factorize
     procedure :: sum => interaction_sum
     procedure :: add_color_contractions => &
          interaction_add_color_contractions
     procedure :: evaluate_product => interaction_evaluate_product
     procedure :: evaluate_product_cf => interaction_evaluate_product_cf
     procedure :: evaluate_square_c => interaction_evaluate_square_c
     procedure :: evaluate_sum => interaction_evaluate_sum
     procedure :: evaluate_me_sum => interaction_evaluate_me_sum
     procedure :: tag_hard_process => interaction_tag_hard_process
     procedure :: retag_hard_process => interaction_retag_hard_process
     procedure :: get_tag => interaction_get_tag
     procedure :: get_n_tot => interaction_get_n_tot
     procedure :: get_n_in => interaction_get_n_in
     procedure :: get_n_vir => interaction_get_n_vir
     procedure :: get_n_out => interaction_get_n_out
     generic :: get_momenta => get_momenta_all, get_momenta_idx
     procedure :: get_momentum => interaction_get_momentum
     procedure :: get_momenta_all => interaction_get_momenta_all
     procedure :: get_momenta_idx => interaction_get_momenta_idx
     procedure :: get_state_matrix_ptr => &
          interaction_get_state_matrix_ptr
     procedure :: get_resonance_flags => interaction_get_resonance_flags
     generic :: get_mask => get_mask_all, get_mask_slice
     procedure :: get_mask_all => interaction_get_mask_all
     procedure :: get_mask_slice => interaction_get_mask_slice
     procedure :: get_s => interaction_get_s
     procedure :: get_cm_transformation => interaction_get_cm_transformation
     procedure :: get_unstable_particle => interaction_get_unstable_particle
     procedure :: get_flv_out => interaction_get_flv_out
     procedure :: get_flv_content => interaction_get_flv_content
     procedure :: set_mask => interaction_set_mask
     procedure :: reset_momenta => interaction_reset_momenta
     procedure :: set_momenta => interaction_set_momenta
     procedure :: set_momentum => interaction_set_momentum
     procedure :: set_flavored_values => interaction_set_flavored_values
     procedure :: relate => interaction_relate
     procedure :: transfer_relations => interaction_transfer_relations
     procedure :: relate_connections => interaction_relate_connections
     procedure :: get_n_children => interaction_get_n_children
     procedure :: get_n_parents => interaction_get_n_parents
     procedure :: get_children => interaction_get_children
     procedure :: get_parents => interaction_get_parents
     procedure :: set_source_link => interaction_set_source_link
     procedure :: find_source => interaction_find_source
     procedure :: exchange_mask => interaction_exchange_mask
     procedure :: receive_momenta => interaction_receive_momenta
     procedure :: send_momenta => interaction_send_momenta
     procedure :: pacify_momenta => interaction_pacify_momenta
     procedure :: declare_subtraction => interaction_declare_subtraction
  end type interaction_t


  interface assignment(=)
     module procedure interaction_assign
  end interface


  interface
    module subroutine qn_index_map_init_trivial (self, int)
      class(qn_index_map_t), intent(out) :: self
      class(interaction_t), intent(in) :: int
    end subroutine qn_index_map_init_trivial
    module subroutine qn_index_map_init_involved (self, int, qn_flv, n_sub, qn_hel)
      class(qn_index_map_t), intent(out) :: self
      type(interaction_t), intent(in) :: int
      type(quantum_numbers_t), dimension(:, :), intent(in) :: qn_flv
      integer, intent(in) :: n_sub
      type(quantum_numbers_t), dimension(:, :), intent(in), optional :: qn_hel
    end subroutine qn_index_map_init_involved
    module subroutine qn_index_map_init_sf (self, int, qn_flv, n_flv_born, n_flv_real)
      class(qn_index_map_t), intent(out) :: self
      type(interaction_t), intent(in) :: int
      integer, intent(in) :: n_flv_born, n_flv_real
      type(quantum_numbers_t), dimension(:,:), intent(in) :: qn_flv
    end subroutine qn_index_map_init_sf
    module subroutine qn_index_map_write (self, unit)
      class(qn_index_map_t), intent(in) :: self
      integer, intent(in), optional :: unit
    end subroutine qn_index_map_write
    module subroutine qn_index_map_set_helicity_flip (self, yorn)
      class(qn_index_map_t), intent(inout) :: self
      logical, intent(in) :: yorn
    end subroutine qn_index_map_set_helicity_flip
    module function qn_index_map_get_index (self, i_flv, i_hel, i_sub) result (index)
      class(qn_index_map_t), intent(in) :: self
      integer :: index
      integer, intent(in) :: i_flv
      integer, intent(in), optional :: i_hel
      integer, intent(in), optional :: i_sub
    end function qn_index_map_get_index
    module function qn_index_map_get_n_flv (self) result (n_flv)
      class(qn_index_map_t), intent(in) :: self
      integer :: n_flv
    end function qn_index_map_get_n_flv
    module function qn_index_map_get_n_hel (self) result (n_hel)
      class(qn_index_map_t), intent(in) :: self
      integer :: n_hel
    end function qn_index_map_get_n_hel
    module function qn_index_map_get_n_sub (self) result (n_sub)
      class(qn_index_map_t), intent(in) :: self
      integer :: n_sub
    end function qn_index_map_get_n_sub
    module function qn_index_map_get_index_by_qn (self, qn, i_sub) result (index)
      class(qn_index_map_t), intent(in) :: self
      integer :: index
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      integer, intent(in), optional :: i_sub
    end function qn_index_map_get_index_by_qn
    module function qn_index_map_get_sf_index_born (self, i_born, i_sub) result (index)
      class(qn_index_map_t), intent(in) :: self
      integer, intent(in) :: i_born, i_sub
      integer :: index
    end function qn_index_map_get_sf_index_born
    module function qn_index_map_get_sf_index_real (self, i_real, i_sub) result (index)
      class(qn_index_map_t), intent(in) :: self
      integer, intent(in) :: i_real, i_sub
      integer :: index
    end function qn_index_map_get_sf_index_real
    module subroutine external_link_set (link, int, i)
      type(external_link_t), intent(out) :: link
      type(interaction_t), target, intent(in) :: int
      integer, intent(in) :: i
    end subroutine external_link_set
    module subroutine external_link_reassign (link, int_src, int_target)
      type(external_link_t), intent(inout) :: link
      type(interaction_t), intent(in) :: int_src
      type(interaction_t), intent(in), target :: int_target
    end subroutine external_link_reassign
    module function external_link_is_set (link) result (flag)
      logical :: flag
      type(external_link_t), intent(in) :: link
    end function external_link_is_set
    module function external_link_get_ptr (link) result (int)
      type(interaction_t), pointer :: int
      type(external_link_t), intent(in) :: link
    end function external_link_get_ptr
    module function external_link_get_index (link) result (i)
      integer :: i
      type(external_link_t), intent(in) :: link
    end function external_link_get_index
    module function external_link_get_momentum_ptr (link) result (p)
      type(vector4_t), pointer :: p
      type(external_link_t), intent(in) :: link
    end function external_link_get_momentum_ptr
    module subroutine internal_link_list_write (object, unit)
      class(internal_link_list_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine internal_link_list_write
    module subroutine internal_link_list_append (link_list, link)
      class(internal_link_list_t), intent(inout) :: link_list
      integer, intent(in) :: link
    end subroutine internal_link_list_append
    module function internal_link_list_has_entries (link_list) result (flag)
      class(internal_link_list_t), intent(in) :: link_list
      logical :: flag
    end function internal_link_list_has_entries
    module function internal_link_list_get_length (link_list) result (length)
      class(internal_link_list_t), intent(in) :: link_list
      integer :: length
    end function internal_link_list_get_length
    module function internal_link_list_get_link (link_list, i) result (link)
      class(internal_link_list_t), intent(in) :: link_list
      integer, intent(in) :: i
      integer :: link
    end function internal_link_list_get_link
    module subroutine interaction_init &
         (int, n_in, n_vir, n_out, &
          tag, resonant, mask, hel_lock, set_relations, store_values)
      class(interaction_t), intent(out) :: int
      integer, intent(in) :: n_in, n_vir, n_out
      integer, intent(in), optional :: tag
      logical, dimension(:), intent(in), optional :: resonant
      type(quantum_numbers_mask_t), dimension(:), intent(in), optional :: mask
      integer, dimension(:), intent(in), optional :: hel_lock
      logical, intent(in), optional :: set_relations, store_values
    end subroutine interaction_init
    module subroutine interaction_init_qn_index_trivial (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_init_qn_index_trivial
    module subroutine interaction_init_qn_index_involved (int, qn_flv, n_sub, qn_hel)
      class(interaction_t), intent(inout) :: int
      type(quantum_numbers_t), dimension(:, :), intent(in) :: qn_flv
      integer, intent(in) :: n_sub
      type(quantum_numbers_t), dimension(:, :), intent(in), optional :: qn_hel
    end subroutine interaction_init_qn_index_involved
    module subroutine interaction_init_qn_index_sf (int, qn_flv, n_flv_born, n_flv_real)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: n_flv_born, n_flv_real
      type(quantum_numbers_t), dimension(:,:), intent(in) :: qn_flv
    end subroutine interaction_init_qn_index_sf
    module subroutine interaction_set_qn_index_helicity_flip (int, yorn)
      class(interaction_t), intent(inout) :: int
      logical, intent(in) :: yorn
    end subroutine interaction_set_qn_index_helicity_flip
    module function interaction_get_qn_index (int, i_flv, i_hel, i_sub) result (index)
      class(interaction_t), intent(in) :: int
      integer :: index
      integer, intent(in) :: i_flv
      integer, intent(in), optional :: i_hel
      integer, intent(in), optional :: i_sub
    end function interaction_get_qn_index
    module function interaction_get_sf_qn_index_born (int, i_born, i_sub) result (index)
      class(interaction_t), intent(in) :: int
      integer :: index
      integer, intent(in) :: i_born, i_sub
    end function interaction_get_sf_qn_index_born
    module function interaction_get_sf_qn_index_real (int, i_real, i_sub) result (index)
      class(interaction_t), intent(in) :: int
      integer :: index
      integer, intent(in) :: i_real, i_sub
    end function interaction_get_sf_qn_index_real
    module function interaction_get_qn_index_n_flv (int) result (index)
      class(interaction_t), intent(in) :: int
      integer :: index
    end function interaction_get_qn_index_n_flv
    module function interaction_get_qn_index_n_hel (int) result (index)
      class(interaction_t), intent(in) :: int
      integer :: index
    end function interaction_get_qn_index_n_hel
    module function interaction_get_qn_index_n_sub (int) result (index)
      class(interaction_t), intent(in) :: int
      integer :: index
    end function interaction_get_qn_index_n_sub
    module subroutine interaction_set_tag (int, tag)
      type(interaction_t), intent(inout), optional :: int
      integer, intent(in), optional :: tag
    end subroutine interaction_set_tag
    module subroutine reset_interaction_counter (tag)
      integer, intent(in), optional :: tag
    end subroutine reset_interaction_counter
    module subroutine interaction_final (object)
      class(interaction_t), intent(inout) :: object
    end subroutine interaction_final
    module subroutine interaction_write &
         (int, unit, verbose, show_momentum_sum, show_mass, show_state, &
         col_verbose, testflag)
      class(interaction_t), intent(in) :: int
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, show_momentum_sum, show_mass
      logical, intent(in), optional :: show_state, col_verbose, testflag
    end subroutine interaction_write
    module subroutine interaction_write_state_matrix (int, unit, write_value_list, &
       verbose, col_verbose, testflag)
      class(interaction_t), intent(in) :: int
      logical, intent(in), optional :: write_value_list, verbose, col_verbose
      logical, intent(in), optional :: testflag
      integer, intent(in), optional :: unit
    end subroutine interaction_write_state_matrix
    module subroutine interaction_reduce_state_matrix (int, qn_mask, keep_order)
      class(interaction_t), intent(inout) :: int
      type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask
      logical, optional, intent(in) :: keep_order
    end subroutine interaction_reduce_state_matrix
    module subroutine interaction_assign (int_out, int_in)
      type(interaction_t), intent(out) :: int_out
      type(interaction_t), intent(in), target :: int_in
    end subroutine interaction_assign
    module subroutine interaction_add_state &
         (int, qn, index, value, sum_values, counter_index, ignore_sub_for_qn, me_index)
      class(interaction_t), intent(inout) :: int
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      integer, intent(in), optional :: index
      complex(default), intent(in), optional :: value
      logical, intent(in), optional :: sum_values
      integer, intent(in), optional :: counter_index
      logical, intent(in), optional :: ignore_sub_for_qn
      integer, intent(out), optional :: me_index
    end subroutine interaction_add_state
    module subroutine interaction_set_duplicate_flv_zero (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_set_duplicate_flv_zero
    module subroutine interaction_freeze (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_freeze
    pure module function interaction_is_empty (int) result (flag)
      logical :: flag
      class(interaction_t), intent(in) :: int
    end function interaction_is_empty
    pure module function interaction_get_n_matrix_elements (int) result (n)
      integer :: n
      class(interaction_t), intent(in) :: int
    end function interaction_get_n_matrix_elements
    module function interaction_get_state_depth (int) result (n)
      integer :: n
      class(interaction_t), intent(in) :: int
    end function interaction_get_state_depth
    module function interaction_get_n_in_helicities (int) result (n_hel)
      integer :: n_hel
      class(interaction_t), intent(in) :: int
    end function interaction_get_n_in_helicities
    pure module function interaction_get_me_size (int) result (n)
      integer :: n
      class(interaction_t), intent(in) :: int
    end function interaction_get_me_size
    pure module function interaction_get_norm (int) result (norm)
      real(default) :: norm
      class(interaction_t), intent(in) :: int
    end function interaction_get_norm
    module function interaction_get_n_sub (int) result (n_sub)
      integer :: n_sub
      class(interaction_t), intent(in) :: int
    end function interaction_get_n_sub
    module function interaction_get_quantum_numbers_single (int, i, by_me_index) result (qn)
      type(quantum_numbers_t), dimension(:), allocatable :: qn
      class(interaction_t), intent(in), target :: int
      integer, intent(in) :: i
      logical, intent(in), optional :: by_me_index
    end function interaction_get_quantum_numbers_single
    module function interaction_get_quantum_numbers_all (int) result (qn)
      type(quantum_numbers_t), dimension(:,:), allocatable :: qn
      class(interaction_t), intent(in), target :: int
    end function interaction_get_quantum_numbers_all
    module function interaction_get_quantum_numbers_all_qn_mask (int, qn_mask) &
       result (qn)
      type(quantum_numbers_t), dimension(:,:), allocatable :: qn
      class(interaction_t), intent(in) :: int
      type(quantum_numbers_mask_t), intent(in) :: qn_mask
    end function interaction_get_quantum_numbers_all_qn_mask
    module subroutine interaction_get_quantum_numbers_all_sub (int, qn)
      class(interaction_t), intent(in) :: int
      type(quantum_numbers_t), dimension(:,:), allocatable, intent(out) :: qn
    end subroutine interaction_get_quantum_numbers_all_sub
    module subroutine interaction_get_flavors (int, only_elementary, qn_mask, flv)
      class(interaction_t), intent(in), target :: int
      logical, intent(in) :: only_elementary
      type(quantum_numbers_mask_t), intent(in), dimension(:), optional :: qn_mask
      integer, intent(out), dimension(:,:), allocatable :: flv
    end subroutine interaction_get_flavors
    module subroutine interaction_get_quantum_numbers_mask (int, qn_mask, qn)
      class(interaction_t), intent(in) :: int
      type(quantum_numbers_mask_t), intent(in) :: qn_mask
      type(quantum_numbers_t), dimension(:,:), allocatable, intent(out) :: qn
    end subroutine interaction_get_quantum_numbers_mask
    elemental module function interaction_get_matrix_element_single (int, i) result (me)
      complex(default) :: me
      class(interaction_t), intent(in) :: int
      integer, intent(in) :: i
    end function interaction_get_matrix_element_single
    module function interaction_get_matrix_element_array (int) result (me)
      complex(default), dimension(:), allocatable :: me
      class(interaction_t), intent(in) :: int
    end function interaction_get_matrix_element_array
    module subroutine interaction_set_matrix_element_qn (int, qn, val)
      class(interaction_t), intent(inout) :: int
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      complex(default), intent(in) :: val
    end subroutine interaction_set_matrix_element_qn
    module subroutine interaction_set_matrix_element_all (int, value)
      class(interaction_t), intent(inout) :: int
      complex(default), intent(in) :: value
    end subroutine interaction_set_matrix_element_all
    module subroutine interaction_set_matrix_element_array (int, value, range)
      class(interaction_t), intent(inout) :: int
      complex(default), intent(in), dimension(:) :: value
      integer, intent(in), dimension(:), optional :: range
    end subroutine interaction_set_matrix_element_array
    pure module subroutine interaction_set_matrix_element_single (int, i, value)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: i
      complex(default), intent(in) :: value
    end subroutine interaction_set_matrix_element_single
    module subroutine interaction_set_matrix_element_clone (int, int1)
      class(interaction_t), intent(inout) :: int
      class(interaction_t), intent(in) :: int1
    end subroutine interaction_set_matrix_element_clone
    module subroutine interaction_set_only_matrix_element (int, i, value)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: i
      complex(default), intent(in) :: value
    end subroutine interaction_set_only_matrix_element
    module subroutine interaction_add_to_matrix_element (int, qn, value, match_only_flavor)
      class(interaction_t), intent(inout) :: int
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      complex(default), intent(in) :: value
      logical, intent(in), optional :: match_only_flavor
    end subroutine interaction_add_to_matrix_element
    module subroutine interaction_get_diagonal_entries (int, i)
      class(interaction_t), intent(in) :: int
      integer, dimension(:), allocatable, intent(out) :: i
    end subroutine interaction_get_diagonal_entries
    module subroutine interaction_normalize_by_trace (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_normalize_by_trace
    module subroutine interaction_normalize_by_max (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_normalize_by_max
    module subroutine interaction_set_norm (int, norm)
      class(interaction_t), intent(inout) :: int
      real(default), intent(in) :: norm
    end subroutine interaction_set_norm
    module subroutine interaction_set_state_matrix (int, state)
      class(interaction_t), intent(inout) :: int
      type(state_matrix_t), intent(in) :: state
    end subroutine interaction_set_state_matrix
    module function interaction_get_max_color_value (int) result (cmax)
      class(interaction_t), intent(in) :: int
      integer :: cmax
    end function interaction_get_max_color_value
    module subroutine interaction_factorize &
         (int, mode, x, ok, single_state, correlated_state, qn_in)
      class(interaction_t), intent(in), target :: int
      integer, intent(in) :: mode
      real(default), intent(in) :: x
      logical, intent(out) :: ok
      type(state_matrix_t), &
           dimension(:), allocatable, intent(out) :: single_state
      type(state_matrix_t), intent(out), optional :: correlated_state
      type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_in
    end subroutine interaction_factorize
    module function interaction_sum (int) result (value)
      class(interaction_t), intent(in) :: int
      complex(default) :: value
    end function interaction_sum
    module subroutine interaction_add_color_contractions (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_add_color_contractions
    pure module subroutine interaction_evaluate_product &
         (int, i, int1, int2, index1, index2)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: i
      type(interaction_t), intent(in) :: int1, int2
      integer, dimension(:), intent(in) :: index1, index2
    end subroutine interaction_evaluate_product
    pure module subroutine interaction_evaluate_product_cf &
         (int, i, int1, int2, index1, index2, factor)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: i
      type(interaction_t), intent(in) :: int1, int2
      integer, dimension(:), intent(in) :: index1, index2
      complex(default), dimension(:), intent(in) :: factor
    end subroutine interaction_evaluate_product_cf
    pure module subroutine interaction_evaluate_square_c (int, i, int1, index1)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: i
      type(interaction_t), intent(in) :: int1
      integer, dimension(:), intent(in) :: index1
    end subroutine interaction_evaluate_square_c
    pure module subroutine interaction_evaluate_sum (int, i, int1, index1)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: i
      type(interaction_t), intent(in) :: int1
      integer, dimension(:), intent(in) :: index1
    end subroutine interaction_evaluate_sum
    pure module subroutine interaction_evaluate_me_sum (int, i, int1, index1)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: i
      type(interaction_t), intent(in) :: int1
      integer, dimension(:), intent(in) :: index1
    end subroutine interaction_evaluate_me_sum
    module subroutine interaction_tag_hard_process (int, tag)
      class(interaction_t), intent(inout) :: int
      integer, dimension(:), intent(in), optional :: tag
    end subroutine interaction_tag_hard_process
    module subroutine interaction_retag_hard_process (int, i, hard)
      class(interaction_t), intent(inout), target :: int
      integer, intent(in) :: i
      logical, intent(in) :: hard
    end subroutine interaction_retag_hard_process
    module function interaction_get_tag (int) result (tag)
      class(interaction_t), intent(in) :: int
      integer :: tag
    end function interaction_get_tag
    pure module function interaction_get_n_tot (object) result (n_tot)
      class(interaction_t), intent(in) :: object
      integer :: n_tot
    end function interaction_get_n_tot
    pure module function interaction_get_n_in (object) result (n_in)
      class(interaction_t), intent(in) :: object
      integer :: n_in
    end function interaction_get_n_in
    pure module function interaction_get_n_vir (object) result (n_vir)
      class(interaction_t), intent(in) :: object
      integer :: n_vir
    end function interaction_get_n_vir
    pure module function interaction_get_n_out (object) result (n_out)
      class(interaction_t), intent(in) :: object
      integer :: n_out
    end function interaction_get_n_out
    module function idx (int, i, outgoing)
      integer :: idx
      type(interaction_t), intent(in) :: int
      integer, intent(in) :: i
      logical, intent(in), optional :: outgoing
    end function idx
    module function interaction_get_momenta_all (int, outgoing) result (p)
      class(interaction_t), intent(in) :: int
      type(vector4_t), dimension(:), allocatable :: p
      logical, intent(in), optional :: outgoing
    end function interaction_get_momenta_all
    module function interaction_get_momenta_idx (int, jj) result (p)
      class(interaction_t), intent(in) :: int
      type(vector4_t), dimension(:), allocatable :: p
      integer, dimension(:), intent(in) :: jj
    end function interaction_get_momenta_idx
    module function interaction_get_momentum (int, i, outgoing) result (p)
      class(interaction_t), intent(in) :: int
      type(vector4_t) :: p
      integer, intent(in) :: i
      logical, intent(in), optional :: outgoing
    end function interaction_get_momentum
    module function interaction_get_state_matrix_ptr (int) result (state)
      class(interaction_t), intent(in), target :: int
      type(state_matrix_t), pointer :: state
    end function interaction_get_state_matrix_ptr
    module function interaction_get_resonance_flags (int) result (resonant)
      class(interaction_t), intent(in) :: int
      logical, dimension(size(int%resonant)) :: resonant
    end function interaction_get_resonance_flags
    module function interaction_get_mask_all (int) result (mask)
      class(interaction_t), intent(in) :: int
      type(quantum_numbers_mask_t), dimension(size(int%mask)) :: mask
    end function interaction_get_mask_all
    module function interaction_get_mask_slice (int, index) result (mask)
      class(interaction_t), intent(in) :: int
      integer, dimension(:), intent(in) :: index
      type(quantum_numbers_mask_t), dimension(size(index)) :: mask
    end function interaction_get_mask_slice
    module function interaction_get_s (int) result (s)
      real(default) :: s
      class(interaction_t), intent(in) :: int
    end function interaction_get_s
    module function interaction_get_cm_transformation (int) result (lt)
      type(lorentz_transformation_t) :: lt
      class(interaction_t), intent(in) :: int
    end function interaction_get_cm_transformation
    module subroutine interaction_get_unstable_particle (int, flv, p, i)
      class(interaction_t), intent(in), target :: int
      type(flavor_t), intent(out) :: flv
      type(vector4_t), intent(out) :: p
      integer, intent(out) :: i
    end subroutine interaction_get_unstable_particle
    module subroutine interaction_get_flv_out (int, flv)
      class(interaction_t), intent(in), target :: int
      type(flavor_t), dimension(:,:), allocatable, intent(out) :: flv
    end subroutine interaction_get_flv_out
    module subroutine interaction_get_flv_content (int, state_flv, n_out_hard)
      class(interaction_t), intent(in), target :: int
      type(state_flv_content_t), intent(out) :: state_flv
      integer, intent(in) :: n_out_hard
    end subroutine interaction_get_flv_content
    module subroutine interaction_set_mask (int, mask)
      class(interaction_t), intent(inout) :: int
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
    end subroutine interaction_set_mask
    module subroutine interaction_reset_momenta (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_reset_momenta
    module subroutine interaction_set_momenta (int, p, outgoing)
      class(interaction_t), intent(inout) :: int
      type(vector4_t), dimension(:), intent(in) :: p
      logical, intent(in), optional :: outgoing
    end subroutine interaction_set_momenta
    module subroutine interaction_set_momentum (int, p, i, outgoing)
      class(interaction_t), intent(inout) :: int
      type(vector4_t), intent(in) :: p
      integer, intent(in) :: i
      logical, intent(in), optional :: outgoing
    end subroutine interaction_set_momentum
    module subroutine interaction_set_flavored_values (int, value, flv_in, pos)
      class(interaction_t), intent(inout) :: int
      complex(default), dimension(:), intent(in) :: value
      type(flavor_t), dimension(:), intent(in) :: flv_in
      integer, intent(in) :: pos
    end subroutine interaction_set_flavored_values
    module subroutine interaction_relate (int, i1, i2)
      class(interaction_t), intent(inout), target :: int
      integer, intent(in) :: i1, i2
    end subroutine interaction_relate
    module subroutine interaction_transfer_relations (int1, int2, map)
      class(interaction_t), intent(in) :: int1
      class(interaction_t), intent(inout), target :: int2
      integer, dimension(:), intent(in) :: map
    end subroutine interaction_transfer_relations
    module subroutine interaction_relate_connections &
         (int, int_in, connection_index, &
          map, map_connections, resonant)
      class(interaction_t), intent(inout), target :: int
      class(interaction_t), intent(in) :: int_in
      integer, dimension(:), intent(in) :: connection_index
      integer, dimension(:), intent(in) :: map, map_connections
      logical, intent(in), optional :: resonant
    end subroutine interaction_relate_connections
    module function interaction_get_n_children (int, i) result (n)
      integer :: n
      class(interaction_t), intent(in) :: int
      integer, intent(in) :: i
    end function interaction_get_n_children
    module function interaction_get_n_parents (int, i) result (n)
      integer :: n
      class(interaction_t), intent(in) :: int
      integer, intent(in) :: i
    end function interaction_get_n_parents
    module function interaction_get_children (int, i) result (idx)
      integer, dimension(:), allocatable :: idx
      class(interaction_t), intent(in) :: int
      integer, intent(in) :: i
    end function interaction_get_children
    module function interaction_get_parents (int, i) result (idx)
      integer, dimension(:), allocatable :: idx
      class(interaction_t), intent(in) :: int
      integer, intent(in) :: i
    end function interaction_get_parents
    module subroutine interaction_set_source_link (int, i, int1, i1)
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: i
      class(interaction_t), intent(in), target :: int1
      integer, intent(in) :: i1
    end subroutine interaction_set_source_link
    module subroutine interaction_reassign_links (int, int_src, int_target)
      type(interaction_t), intent(inout) :: int
      type(interaction_t), intent(in) :: int_src
      type(interaction_t), intent(in), target :: int_target
    end subroutine interaction_reassign_links
    module function interaction_find_link (int, int1, i1) result (i)
      integer :: i
      type(interaction_t), intent(in) :: int, int1
      integer, intent(in) :: i1
    end function interaction_find_link
    module subroutine interaction_find_source (int, i, int1, i1)
      class(interaction_t), intent(in) :: int
      integer, intent(in) :: i
      type(interaction_t), intent(out), pointer :: int1
      integer, intent(out) :: i1
    end subroutine interaction_find_source
    module function interaction_get_ultimate_source (int, i) result (link)
      type(external_link_t) :: link
      type(interaction_t), intent(in) :: int
      integer, intent(in) :: i
    end function interaction_get_ultimate_source
    module subroutine interaction_exchange_mask (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_exchange_mask
    module subroutine interaction_receive_momenta (int)
      class(interaction_t), intent(inout) :: int
    end subroutine interaction_receive_momenta
    module subroutine interaction_send_momenta (int)
      class(interaction_t), intent(in) :: int
    end subroutine interaction_send_momenta
    module subroutine interaction_pacify_momenta (int, acc)
      class(interaction_t), intent(inout) :: int
      real(default), intent(in) :: acc
    end subroutine interaction_pacify_momenta
    module subroutine interaction_declare_subtraction (int, n_sub)
      class(interaction_t), intent(inout), target :: int
      integer, intent(in) :: n_sub
    end subroutine interaction_declare_subtraction
    module subroutine find_connections (int1, int2, n, connection_index)
      class(interaction_t), intent(in) :: int1, int2
      integer, intent(out) :: n
      integer, dimension(:,:), intent(out), allocatable :: connection_index
      integer, dimension(:,:), allocatable :: conn_index_tmp
      integer, dimension(:), allocatable :: ordering
    end subroutine find_connections
  end interface

end module interactions
