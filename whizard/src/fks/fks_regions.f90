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

module fks_regions

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use diagnostics
  use os_interface
  use constants
  use process_constants
  use lorentz
  use models
  use resonances, only: resonance_contributors_t, resonance_history_t
  use phs_fks, only: phs_identifier_t, check_for_phs_identifier

  use nlo_data

  implicit none
  private

  public :: ftuple_t
  public :: flv_structure_t
  public :: singular_region_t
  public :: fks_mapping_default_t
  public :: fks_mapping_resonances_t
  public :: operator(.equiv.)
  public :: operator(.equivtag.)
  public :: region_data_t
  public :: assignment(=)
  public :: create_resonance_histories_for_threshold
  public :: setup_region_data_for_test

  integer, parameter :: UNDEFINED_SPLITTING = 0
  integer, parameter :: F_TO_FV = 1
  integer, parameter :: V_TO_VV = 2
  integer, parameter :: V_TO_FF = 3
  integer, parameter :: F_TO_VF = 4


  type :: ftuple_t
    integer, dimension(2) :: ireg = [-1,-1]
    integer :: i_res = 0
    integer :: splitting_type
    logical :: pseudo_isr = .false.
    logical :: qcd_split = .false.
  contains
    procedure :: write => ftuple_write
    procedure :: get => ftuple_get
    procedure :: set => ftuple_set
    procedure :: determine_splitting_type_fsr => &
         ftuple_determine_splitting_type_fsr
    procedure :: determine_splitting_type_isr => &
         ftuple_determine_splitting_type_isr
    procedure :: determine_sub_correction_type => &
         ftuple_determine_sub_correction_type
    procedure :: has_negative_elements => ftuple_has_negative_elements
    procedure :: has_identical_elements => ftuple_has_identical_elements
  end type ftuple_t

  type :: ftuple_list_t
    integer :: index = 0
    type(ftuple_t) :: ftuple
    type(ftuple_list_t), pointer :: next => null ()
    type(ftuple_list_t), pointer :: prev => null ()
    type(ftuple_list_t), pointer :: equiv => null ()
  contains
     procedure :: write => ftuple_list_write
     procedure :: append => ftuple_list_append
     procedure :: get_n_tuples => ftuple_list_get_n_tuples
     procedure :: get_entry => ftuple_list_get_entry
     procedure :: get_ftuple => ftuple_list_get_ftuple
     procedure :: set_equiv => ftuple_list_set_equiv
     procedure :: check_equiv => ftuple_list_check_equiv
     procedure :: to_array => ftuple_list_to_array
  end type ftuple_list_t

  type :: flv_structure_t
    integer, dimension(:), allocatable :: flst
    integer, dimension(:), allocatable :: tag
    integer :: nlegs = 0
    integer :: n_in = 0
    logical, dimension(:), allocatable :: massive
    logical, dimension(:), allocatable :: colored
    real(default), dimension(:), allocatable :: charge
    real(default) :: prt_symm_fs = 1._default
    integer :: eqv_index = 0
  contains
    procedure :: valid_pair => flv_structure_valid_pair
    procedure :: remove_particle => flv_structure_remove_particle
    procedure :: insert_particle_fsr => flv_structure_insert_particle_fsr
    procedure :: insert_particle_isr => flv_structure_insert_particle_isr
    procedure :: insert_particle => flv_structure_insert_particle
    procedure :: count_particle => flv_structure_count_particle
    procedure :: init => flv_structure_init
    procedure :: compute_prt_symm_fs => flv_structure_compute_prt_symm_fs
    procedure :: write => flv_structure_write
    procedure :: to_string => flv_structure_to_string
    procedure :: create_uborn => flv_structure_create_uborn
    procedure :: init_mass_color_and_charge => &
         flv_structure_init_mass_color_and_charge
    procedure :: get_last_two => flv_structure_get_last_two
    procedure :: final => flv_structure_final
  end type flv_structure_t

  type :: flv_perm_t
    integer, dimension(:,:), allocatable :: perms
  contains
    procedure :: init => flv_perm_init
    procedure :: write => flv_perm_write
    procedure :: reset => flv_perm_final
    procedure :: final => flv_perm_final
    generic :: apply => &
           apply_flv_structure, apply_index, apply_ftuple
    procedure :: apply_flv_structure => flv_perm_apply_flv_structure
    procedure :: apply_index => flv_perm_apply_index
    procedure :: apply_ftuple => flv_perm_apply_ftuple
    procedure :: eqv => flv_perm_eqv
  end type flv_perm_t

  type :: singular_region_t
    integer :: alr
    integer :: i_res
    type(flv_structure_t) :: flst_real
    type(flv_structure_t) :: flst_uborn
    integer :: mult
    integer :: emitter
    integer :: nregions
    integer :: real_index
    type(ftuple_t), dimension(:), allocatable :: ftuples
    integer :: uborn_index
    logical :: double_fsr = .false.
    logical :: soft_divergence = .false.
    logical :: coll_divergence = .false.
    type(string_t) :: nlo_correction_type
    integer, dimension(:), allocatable :: i_reg_to_i_con
    logical :: pseudo_isr = .false.
    logical :: sc_required = .false.
    integer :: eqv_index = 0
  contains
    procedure :: init => singular_region_init
    procedure :: write => singular_region_write
    procedure :: write_latex => singular_region_write_latex
    procedure :: set_splitting_info => singular_region_set_splitting_info
    procedure :: double_fsr_factor => singular_region_double_fsr_factor
    procedure :: has_soft_divergence => singular_region_has_soft_divergence
    procedure :: has_collinear_divergence => &
           singular_region_has_collinear_divergence
    procedure :: has_identical_ftuples => singular_region_has_identical_ftuples
  end type singular_region_t

  type :: resonance_mapping_t
    type(resonance_history_t), dimension(:), allocatable :: res_histories
    integer, dimension(:), allocatable :: alr_to_i_res
    integer, dimension(:,:), allocatable :: i_res_to_alr
    type(vector4_t), dimension(:), allocatable :: p_res
  contains
    procedure :: init => resonance_mapping_init
    procedure :: set_alr_to_i_res => resonance_mapping_set_alr_to_i_res
    procedure :: get_resonance_history => resonance_mapping_get_resonance_history
    procedure :: write => resonance_mapping_write
    procedure :: get_resonance_value => resonance_mapping_get_resonance_value
    procedure :: get_resonance_all => resonance_mapping_get_resonance_all
    procedure :: get_weight => resonance_mapping_get_weight
    procedure :: get_resonance_alr => resonance_mapping_get_resonance_alr
  end type resonance_mapping_t

  type, abstract :: fks_mapping_t
     real(default) :: sumdij
     real(default) :: sumdij_soft
     logical :: pseudo_isr = .false.
     real(default) :: normalization_factor = one
  contains
    procedure (fks_mapping_dij), deferred :: dij
    procedure (fks_mapping_compute_sumdij), deferred :: compute_sumdij
    procedure (fks_mapping_svalue), deferred :: svalue
    procedure (fks_mapping_dij_soft), deferred :: dij_soft
    procedure (fks_mapping_compute_sumdij_soft), deferred :: compute_sumdij_soft
    procedure (fks_mapping_svalue_soft), deferred :: svalue_soft
  end type fks_mapping_t

  type, extends (fks_mapping_t) :: fks_mapping_default_t
    real(default) :: exp_1, exp_2
    integer :: n_in
  contains
    procedure :: set_parameter => fks_mapping_default_set_parameter
    procedure :: dij => fks_mapping_default_dij
    procedure :: compute_sumdij => fks_mapping_default_compute_sumdij
    procedure :: svalue => fks_mapping_default_svalue
    procedure :: dij_soft => fks_mapping_default_dij_soft
    procedure :: compute_sumdij_soft => fks_mapping_default_compute_sumdij_soft
    procedure :: svalue_soft => fks_mapping_default_svalue_soft
  end type fks_mapping_default_t

  type, extends (fks_mapping_t) :: fks_mapping_resonances_t
    real(default) :: exp_1, exp_2
    type(resonance_mapping_t) :: res_map
    integer :: i_con = 0
  contains
    procedure :: dij => fks_mapping_resonances_dij
    procedure :: compute_sumdij => fks_mapping_resonances_compute_sumdij
    procedure :: svalue => fks_mapping_resonances_svalue
    procedure :: get_resonance_weight => &
         fks_mapping_resonances_get_resonance_weight
    procedure :: dij_soft => fks_mapping_resonances_dij_soft
    procedure :: compute_sumdij_soft => fks_mapping_resonances_compute_sumdij_soft
    procedure :: svalue_soft => fks_mapping_resonances_svalue_soft
    procedure :: set_resonance_momentum => &
         fks_mapping_resonances_set_resonance_momentum
    procedure :: set_resonance_momenta => &
         fks_mapping_resonances_set_resonance_momenta
  end type fks_mapping_resonances_t

  type :: region_data_t
    type(singular_region_t), dimension(:), allocatable :: regions
    type(flv_structure_t), dimension(:), allocatable :: flv_born
    type(flv_structure_t), dimension(:), allocatable :: flv_real
    integer, dimension(:), allocatable :: eqv_flv_index_born
    integer, dimension(:), allocatable :: eqv_flv_index_real
    integer, dimension(:), allocatable :: emitters
    integer :: n_regions = 0
    integer :: n_emitters = 0
    integer :: n_flv_born = 0
    integer :: n_flv_real = 0
    integer :: n_in = 0
    integer :: n_legs_born = 0
    integer :: n_legs_real = 0
    integer :: n_phs = 0
    integer :: alpha_power = 0
    integer :: alphas_power = 0
    type(string_t) :: nlo_correction_type
    class(fks_mapping_t), allocatable :: fks_mapping
    integer, dimension(:), allocatable :: resonances
    type(resonance_contributors_t), dimension(:), allocatable :: alr_contributors
    integer, dimension(:), allocatable :: alr_to_i_contributor
    integer, dimension(:), allocatable :: i_phs_to_i_con
  contains
    procedure :: allocate_fks_mappings => region_data_allocate_fks_mappings
    procedure :: init => region_data_init
    procedure :: init_resonance_information => &
         region_data_init_resonance_information
    procedure :: set_resonance_mappings => region_data_set_resonance_mappings
    procedure :: setup_fks_mappings => region_data_setup_fks_mappings
    procedure :: enlarge_singular_regions_with_resonances &
       => region_data_enlarge_singular_regions_with_resonances
    procedure :: set_isr_pseudo_regions => region_data_set_isr_pseudo_regions
    procedure :: split_up_interference_regions_for_threshold => &
         region_data_split_up_interference_regions_for_threshold
    procedure :: set_mass_color_and_charge => &
         region_data_set_mass_color_and_charge
    procedure :: uses_resonances => region_data_uses_resonances
    procedure :: get_emitter_list => region_data_get_emitter_list
    procedure :: get_n_emitters_sc => region_data_get_n_emitters_sc
    procedure :: get_associated_resonances => &
         region_data_get_associated_resonances
    procedure :: emitter_is_compatible_with_resonance => &
       region_data_emitter_is_compatible_with_resonance
    procedure :: emitter_is_in_resonance => region_data_emitter_is_in_resonance
    procedure :: get_contributors => region_data_get_contributors
    procedure :: get_emitter => region_data_get_emitter
    procedure :: map_real_to_born_index => region_data_map_real_to_born_index
    generic :: get_flv_states_born => get_flv_states_born_single, &
         get_flv_states_born_array
    procedure :: get_flv_states_born_single => &
         region_data_get_flv_states_born_single
    procedure :: get_flv_states_born_array => &
         region_data_get_flv_states_born_array
    generic :: get_flv_states_real => &
         get_flv_states_real_single, get_flv_states_real_array
    procedure :: get_flv_states_real_single => &
         region_data_get_flv_states_real_single
    procedure :: get_flv_states_real_array => &
         region_data_get_flv_states_real_array
    procedure :: get_all_flv_states => region_data_get_all_flv_states
    procedure :: get_n_in => region_data_get_n_in
    procedure :: get_n_legs_real => region_data_get_n_legs_real
    procedure :: get_n_legs_born => region_data_get_n_legs_born
    procedure :: get_n_flv_real => region_data_get_n_flv_real
    procedure :: get_n_flv_born => region_data_get_n_flv_born
    generic :: get_svalue => get_svalue_last_pos, get_svalue_ij
    procedure :: get_svalue_last_pos => region_data_get_svalue_last_pos
    procedure :: get_svalue_ij => region_data_get_svalue_ij
    procedure :: get_svalue_soft => region_data_get_svalue_soft
    procedure :: find_regions => region_data_find_regions
    procedure :: find_eqv_regions => region_data_find_eqv_regions
    procedure :: init_singular_regions => region_data_init_singular_regions
    procedure :: find_emitters => region_data_find_emitters
    procedure :: find_resonances => region_data_find_resonances
    procedure :: set_i_phs_to_i_con => region_data_set_i_phs_to_i_con
    procedure :: set_alr_to_i_phs => region_data_set_alr_to_i_phs
    procedure :: set_contributors => region_data_set_contributors
    procedure :: extend_ftuples => region_data_extend_ftuples
    procedure :: get_flavor_indices => region_data_get_flavor_indices
    procedure :: get_matrix_element_index => region_data_get_matrix_element_index
    procedure :: compute_number_of_phase_spaces &
       => region_data_compute_number_of_phase_spaces
    procedure :: get_n_phs => region_data_get_n_phs
    procedure :: set_splitting_info => region_data_set_splitting_info
    procedure :: init_phs_identifiers => region_data_init_phs_identifiers
    procedure :: get_all_ftuples => region_data_get_all_ftuples
    procedure :: write_to_file => region_data_write_to_file
    procedure :: write_latex => region_data_write_latex
    procedure :: write => region_data_write
    procedure :: has_pseudo_isr => region_data_has_pseudo_isr
    procedure :: check_consistency => region_data_check_consistency
    procedure :: requires_spin_correlations => &
         region_data_requires_spin_correlations
    procedure :: born_to_real_symm_factor_fs => &
         region_data_born_to_real_symm_factor_fs
    procedure :: final => region_data_final
  end type region_data_t


  interface assignment(=)
     module procedure ftuple_assign
  end interface

  interface operator(==)
     module procedure ftuple_equal
  end interface

  interface operator(>)
     module procedure ftuple_greater
  end interface

  interface operator(<)
     module procedure ftuple_less
  end interface

  interface assignment(=)
     module procedure singular_region_assign
  end interface

  interface operator(.match.)
     module procedure singular_region_match
  end interface

  interface assignment(=)
     module procedure resonance_mapping_assign
  end interface

  interface operator(.equiv.)
    module procedure flv_structure_equivalent_no_tag
  end interface

  interface operator(.equivtag.)
    module procedure flv_structure_equivalent_with_tag
  end interface

  interface assignment(=)
    module procedure flv_structure_assign_flv
    module procedure flv_structure_assign_integer
  end interface

  interface assignment(=)
     module procedure region_data_assign
  end interface

  abstract interface
    function fks_mapping_dij (map, p, i, j, i_con) result (d)
      import
      real(default) :: d
      class(fks_mapping_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: i, j
      integer, intent(in), optional :: i_con
    end function fks_mapping_dij
  end interface

  abstract interface
    subroutine fks_mapping_compute_sumdij (map, sregion, p_real)
      import
      class(fks_mapping_t), intent(inout) :: map
      type(singular_region_t), intent(in) :: sregion
      type(vector4_t), intent(in), dimension(:) :: p_real
    end subroutine fks_mapping_compute_sumdij
  end interface

  abstract interface
    function fks_mapping_svalue (map, p, i, j, i_res) result (value)
      import
      real(default) :: value
      class(fks_mapping_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: i, j
      integer, intent(in), optional :: i_res
    end function fks_mapping_svalue
  end interface

  abstract interface
    function fks_mapping_dij_soft (map, p_born, p_soft, em, i_con) result (d)
      import
      real(default) :: d
      class(fks_mapping_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
      integer, intent(in) :: em
      integer, intent(in), optional :: i_con
    end function fks_mapping_dij_soft
  end interface

  abstract interface
    subroutine fks_mapping_compute_sumdij_soft (map, sregion, p_born, p_soft)
      import
      class(fks_mapping_t), intent(inout) :: map
      type(singular_region_t), intent(in) :: sregion
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
    end subroutine fks_mapping_compute_sumdij_soft
  end interface
  abstract interface
    function fks_mapping_svalue_soft (map, p_born, p_soft, em, i_res) result (value)
      import
      real(default) :: value
      class(fks_mapping_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
      integer, intent(in) :: em
      integer, intent(in), optional :: i_res
    end function fks_mapping_svalue_soft
  end interface

  interface assignment(=)
     module procedure fks_mapping_default_assign
  end interface

  interface assignment(=)
     module procedure fks_mapping_resonances_assign
  end interface


  interface
    pure module subroutine ftuple_assign (ftuple_out, ftuple_in)
      type(ftuple_t), intent(out) :: ftuple_out
      type(ftuple_t), intent(in) :: ftuple_in
    end subroutine ftuple_assign
    elemental module function ftuple_equal (f1, f2) result (value)
      logical :: value
      type(ftuple_t), intent(in) :: f1, f2
    end function ftuple_equal
    elemental module function ftuple_greater (f1, f2) result (greater)
      logical :: greater
      type(ftuple_t), intent(in) :: f1, f2
    end function ftuple_greater
    elemental module function ftuple_less (f1, f2) result (less)
      logical :: less
      type(ftuple_t), intent(in) :: f1, f2
    end function ftuple_less
    module subroutine ftuple_write (ftuple, unit, newline)
      class(ftuple_t), intent(in) :: ftuple
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: newline
    end subroutine ftuple_write
    module subroutine ftuple_get (ftuple, pos1, pos2)
      class(ftuple_t), intent(in) :: ftuple
      integer, intent(out) :: pos1, pos2
    end subroutine ftuple_get
    module subroutine ftuple_set (ftuple, pos1, pos2)
      class(ftuple_t), intent(inout) :: ftuple
      integer, intent(in) ::  pos1, pos2
    end subroutine ftuple_set
    module subroutine ftuple_determine_splitting_type_fsr (ftuple, flv, i, j)
      class(ftuple_t), intent(inout) :: ftuple
      type(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: i, j
    end subroutine ftuple_determine_splitting_type_fsr
    module subroutine ftuple_determine_splitting_type_isr (ftuple, flv, i, j)
      class(ftuple_t), intent(inout) :: ftuple
      type(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: i, j
    end subroutine ftuple_determine_splitting_type_isr
    module subroutine ftuple_determine_sub_correction_type &
         (ftuple, flv_born, flv_real, i, j)
      class(ftuple_t), intent(inout) :: ftuple
      type(flv_structure_t), intent(in) :: flv_born, flv_real
      integer, intent(in) :: i, j
    end subroutine ftuple_determine_sub_correction_type
    elemental module function ftuple_has_negative_elements &
         (ftuple) result (value)
      logical :: value
      class(ftuple_t), intent(in) :: ftuple
    end function ftuple_has_negative_elements
    elemental module function ftuple_has_identical_elements &
         (ftuple) result (value)
      logical :: value
      class(ftuple_t), intent(in) :: ftuple
    end function ftuple_has_identical_elements
    module subroutine ftuple_list_write (list, unit, verbose)
      class(ftuple_list_t), intent(in), target :: list
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine ftuple_list_write
    module subroutine ftuple_list_append (list, ftuple)
      class(ftuple_list_t), intent(inout), target :: list
      type(ftuple_t), intent(in) :: ftuple
    end subroutine ftuple_list_append
    impure elemental module function ftuple_list_get_n_tuples &
         (list) result(n_tuples)
      integer :: n_tuples
      class(ftuple_list_t), intent(in), target :: list
    end function ftuple_list_get_n_tuples
    module function ftuple_list_get_entry (list, index) result (entry)
      type(ftuple_list_t), pointer :: entry
      class(ftuple_list_t), intent(in), target :: list
      integer, intent(in) :: index
    end function ftuple_list_get_entry
    module function ftuple_list_get_ftuple (list, index)  result (ftuple)
      type(ftuple_t) :: ftuple
      class(ftuple_list_t), intent(in), target :: list
      integer, intent(in) :: index
    end function ftuple_list_get_ftuple
    module subroutine ftuple_list_set_equiv (list, i1, i2)
      class(ftuple_list_t), intent(in) :: list
      integer, intent(in) :: i1, i2
    end subroutine ftuple_list_set_equiv
    module function ftuple_list_check_equiv(list, i1, i2) result (eq)
      class(ftuple_list_t), intent(in) :: list
      integer, intent(in) :: i1, i2
      logical :: eq
    end function ftuple_list_check_equiv
    module subroutine ftuple_list_to_array &
         (ftuple_list, ftuple_array, equivalences, ordered)
      class(ftuple_list_t), intent(in), target :: ftuple_list
      type(ftuple_t), intent(out), dimension(:), allocatable :: ftuple_array
      logical, intent(out), dimension(:,:), allocatable :: equivalences
      logical, intent(in) :: ordered
    end subroutine ftuple_list_to_array
    module function flv_structure_valid_pair &
       (flv, i, j, flv_ref, model) result (valid)
      logical :: valid
      class(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: i,j
      type(flv_structure_t), intent(in) :: flv_ref
      type(model_t), intent(in) :: model
    end function flv_structure_valid_pair
    module function flv_structure_equivalent_no_tag (flv1, flv2) result(equiv)
      logical :: equiv
      type(flv_structure_t), intent(in) :: flv1, flv2
    end function flv_structure_equivalent_no_tag
    module function flv_structure_equivalent_with_tag (flv1, flv2) result(equiv)
      logical :: equiv
      type(flv_structure_t), intent(in) :: flv1, flv2
    endfunction flv_structure_equivalent_with_tag
    pure module subroutine flv_structure_assign_flv (flv_out, flv_in)
      type(flv_structure_t), intent(out) :: flv_out
      type(flv_structure_t), intent(in) :: flv_in
    end subroutine flv_structure_assign_flv
    pure module subroutine flv_structure_assign_integer (flv_out, iarray)
      type(flv_structure_t), intent(out) :: flv_out
      integer, intent(in), dimension(:) :: iarray
    end subroutine flv_structure_assign_integer
    module function flv_structure_remove_particle (flv, index) result (flv_new)
      type(flv_structure_t) :: flv_new
      class(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: index
    end function flv_structure_remove_particle
    module function flv_structure_insert_particle_fsr &
         (flv, i1, i2, flv_add) result (flv_new)
      type(flv_structure_t) :: flv_new
      class(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: i1, i2, flv_add
    end function flv_structure_insert_particle_fsr
    module function flv_structure_insert_particle_isr &
         (flv, i_in, i_out, flv_add) result (flv_new)
      type(flv_structure_t) :: flv_new
      class(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: i_in, i_out, flv_add
    end function flv_structure_insert_particle_isr
    module function flv_structure_insert_particle &
         (flv, i1, i2, particle) result (flv_new)
      type(flv_structure_t) :: flv_new
      class(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: i1, i2, particle
    end function flv_structure_insert_particle
    module function flv_structure_count_particle (flv, part) result (n)
      class(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: part
      integer :: n
    end function flv_structure_count_particle
    module subroutine flv_structure_init (flv, aval, n_in, tags)
      class(flv_structure_t), intent(inout) :: flv
      integer, intent(in), dimension(:) :: aval
      integer, intent(in) :: n_in
      integer, intent(in), dimension(:), optional :: tags
    end subroutine flv_structure_init
    module subroutine flv_structure_compute_prt_symm_fs (flv, n_in)
      class(flv_structure_t), intent(inout) :: flv
      integer, intent(in) :: n_in
    end subroutine flv_structure_compute_prt_symm_fs
    module subroutine flv_structure_write (flv, unit)
      class(flv_structure_t), intent(in) :: flv
      integer, intent(in), optional :: unit
    end subroutine flv_structure_write
    module function flv_structure_to_string (flv) result (flv_string)
      type(string_t) :: flv_string
      class(flv_structure_t), intent(in) :: flv
    end function flv_structure_to_string
    module function flv_structure_create_uborn &
         (flv, emitter, nlo_correction_type) result(flv_uborn)
      type(flv_structure_t) :: flv_uborn
      class(flv_structure_t), intent(in) :: flv
      type(string_t), intent(in) :: nlo_correction_type
      integer, intent(in) :: emitter
    end function flv_structure_create_uborn
    module subroutine flv_structure_init_mass_color_and_charge (flv, model)
      class(flv_structure_t), intent(inout) :: flv
      type(model_t), intent(in) :: model
    end subroutine flv_structure_init_mass_color_and_charge
    module function flv_structure_get_last_two (flv, n) result (flst_last)
      integer, dimension(2) :: flst_last
      class(flv_structure_t), intent(in) :: flv
      integer, intent(in) :: n
    end function flv_structure_get_last_two
    module subroutine flv_structure_final (flv)
      class(flv_structure_t), intent(inout) :: flv
    end subroutine flv_structure_final
    module subroutine flv_perm_init &
         (perm, flv_in, flv_ref, n_first, n_last, with_tag)
      class(flv_perm_t), intent(out) :: perm
      type(flv_structure_t), intent(in) :: flv_in, flv_ref
      integer, intent(in) :: n_first, n_last
      logical, intent(in) :: with_tag
    end subroutine flv_perm_init
    module subroutine flv_perm_write (perm, unit)
      class(flv_perm_t), intent(in) :: perm
      integer, intent(in), optional :: unit
    end subroutine flv_perm_write
    module subroutine flv_perm_final (perm)
      class(flv_perm_t), intent(inout) :: perm
    end subroutine flv_perm_final
    elemental module function flv_perm_apply_flv_structure &
         (perm, flv_in, invert) result (flv_out)
      type(flv_structure_t) :: flv_out
      class(flv_perm_t), intent(in) :: perm
      type(flv_structure_t), intent(in) :: flv_in
      logical, intent(in), optional :: invert
    end function flv_perm_apply_flv_structure
    elemental module function flv_perm_apply_index &
         (perm, i_in) result (i_out)
      integer :: i_out
      class(flv_perm_t), intent(in) :: perm
      integer, intent(in) :: i_in
    end function flv_perm_apply_index
    elemental module function flv_perm_apply_ftuple &
         (perm, f_in) result (f_out)
      type(ftuple_t) :: f_out
      class(flv_perm_t), intent(in) :: perm
      type(ftuple_t), intent(in) :: f_in
    end function flv_perm_apply_ftuple
    module function flv_perm_eqv &
         (perm, flv1, flv2, with_tag) result (valid)
      logical :: valid
      class(flv_perm_t), intent(in) :: perm
      type(flv_structure_t), intent(in) :: flv1, flv2
      logical, intent(in) :: with_tag
    end function flv_perm_eqv
    module subroutine singular_region_init (sregion, alr, mult, i_res, &
         flst_real, flst_uborn, flv_born, emitter, ftuples, equivalences, &
         nlo_correction_type)
      class(singular_region_t), intent(out) :: sregion
      integer, intent(in) :: alr, mult, i_res
      type(flv_structure_t), intent(in) :: flst_real
      type(flv_structure_t), intent(in) :: flst_uborn
      type(flv_structure_t), dimension(:), intent(in) :: flv_born
      integer, intent(in) :: emitter
      type(ftuple_t), intent(inout), dimension(:) :: ftuples
      logical, intent(inout), dimension(:,:) :: equivalences
      type(string_t), intent(in) :: nlo_correction_type
    end subroutine singular_region_init
    module subroutine singular_region_write (sregion, unit, maxnregions)
      class(singular_region_t), intent(in) :: sregion
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: maxnregions
    end subroutine singular_region_write
    module subroutine singular_region_write_latex (region, unit)
      class(singular_region_t), intent(in) :: region
      integer, intent(in), optional :: unit
    end subroutine singular_region_write_latex
    module subroutine singular_region_set_splitting_info (region, n_in)
      class(singular_region_t), intent(inout) :: region
      integer, intent(in) :: n_in
    end subroutine singular_region_set_splitting_info
    module function singular_region_double_fsr_factor (region, p) result (val)
      class(singular_region_t), intent(in) :: region
      type(vector4_t), intent(in), dimension(:) :: p
      real(default) :: val
    end function singular_region_double_fsr_factor
    module function singular_region_has_soft_divergence (region) result (div)
      logical :: div
      class(singular_region_t), intent(in) :: region
    end function singular_region_has_soft_divergence
    module function singular_region_has_collinear_divergence &
         (region) result (div)
      logical :: div
      class(singular_region_t), intent(in) :: region
    end function singular_region_has_collinear_divergence
    elemental module function singular_region_has_identical_ftuples &
         (sregion) result (value)
      logical :: value
      class(singular_region_t), intent(in) :: sregion
    end function singular_region_has_identical_ftuples
    module subroutine singular_region_assign (reg_out, reg_in)
      type(singular_region_t), intent(out) :: reg_out
      type(singular_region_t), intent(in) :: reg_in
    end subroutine singular_region_assign
    module function singular_region_match (reg1, reg2) result (match)
      type(singular_region_t), intent(in) :: reg1, reg2
      logical :: match
    end function singular_region_match
    module subroutine resonance_mapping_init (res_map, res_hist)
      class(resonance_mapping_t), intent(inout) :: res_map
      type(resonance_history_t), intent(in), dimension(:) :: res_hist
    end subroutine resonance_mapping_init
    module subroutine resonance_mapping_set_alr_to_i_res &
         (res_map, regions, alr_new_to_old)
      class(resonance_mapping_t), intent(inout) :: res_map
      type(singular_region_t), intent(in), dimension(:) :: regions
      integer, intent(out), dimension(:), allocatable :: alr_new_to_old
    end subroutine resonance_mapping_set_alr_to_i_res
    module function resonance_mapping_get_resonance_history &
         (res_map, alr) result (res_hist)
      type(resonance_history_t) :: res_hist
      class(resonance_mapping_t), intent(in) :: res_map
      integer, intent(in) :: alr
    end function resonance_mapping_get_resonance_history
    module subroutine resonance_mapping_write (res_map)
      class(resonance_mapping_t), intent(in) :: res_map
    end subroutine resonance_mapping_write
    module function resonance_mapping_get_resonance_value &
         (res_map, i_res, p, i_gluon) result (p_map)
      real(default) :: p_map
      class(resonance_mapping_t), intent(in) :: res_map
      integer, intent(in) :: i_res
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in), optional :: i_gluon
    end function resonance_mapping_get_resonance_value
    module function resonance_mapping_get_resonance_all &
         (res_map, alr, p, i_gluon) result (p_map)
      real(default) :: p_map
      class(resonance_mapping_t), intent(in) :: res_map
      integer, intent(in) :: alr
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in), optional :: i_gluon
    end function resonance_mapping_get_resonance_all
    module function resonance_mapping_get_weight (res_map, alr, p) result (pfr)
      real(default) :: pfr
      class(resonance_mapping_t), intent(in) :: res_map
      integer, intent(in) :: alr
      type(vector4_t), intent(in), dimension(:) :: p
    end function resonance_mapping_get_weight
    module function resonance_mapping_get_resonance_alr &
         (res_map, alr, p, i_gluon) result (p_map)
      real(default) :: p_map
      class(resonance_mapping_t), intent(in) :: res_map
      integer, intent(in) :: alr
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in), optional :: i_gluon
    end function resonance_mapping_get_resonance_alr
    module subroutine resonance_mapping_assign (res_map_out, res_map_in)
      type(resonance_mapping_t), intent(out) :: res_map_out
      type(resonance_mapping_t), intent(in) :: res_map_in
    end subroutine resonance_mapping_assign
    module subroutine region_data_init (reg_data, n_in, model, flavor_born, &
         flavor_real, nlo_correction_type, alpha_pow, alphas_pow)
      class(region_data_t), intent(out) :: reg_data
      integer, intent(in) :: n_in, alpha_pow, alphas_pow
      type(model_t), intent(in) :: model
      integer, intent(in), dimension(:,:) :: flavor_born, flavor_real
      type(string_t), intent(in) :: nlo_correction_type
    end subroutine region_data_init
    module subroutine region_data_init_resonance_information (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_init_resonance_information
    module subroutine region_data_set_resonance_mappings &
         (reg_data, resonance_histories)
      class(region_data_t), intent(inout) :: reg_data
      type(resonance_history_t), intent(in), dimension(:) :: resonance_histories
    end subroutine region_data_set_resonance_mappings
    module subroutine region_data_setup_fks_mappings (reg_data, template, n_in)
      class(region_data_t), intent(inout) :: reg_data
      type(fks_template_t), intent(in) :: template
      integer, intent(in) :: n_in
    end subroutine region_data_setup_fks_mappings
    module subroutine region_data_enlarge_singular_regions_with_resonances &
         (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_enlarge_singular_regions_with_resonances
    module subroutine region_data_set_isr_pseudo_regions (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_set_isr_pseudo_regions
    module subroutine region_data_split_up_interference_regions_for_threshold &
         (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_split_up_interference_regions_for_threshold
    module subroutine region_data_set_mass_color_and_charge (reg_data, model)
      class(region_data_t), intent(inout) :: reg_data
      type(model_t), intent(in) :: model
    end subroutine region_data_set_mass_color_and_charge
    module function region_data_uses_resonances (reg_data) result (val)
      logical :: val
      class(region_data_t), intent(in) :: reg_data
    end function region_data_uses_resonances
    pure module function region_data_get_emitter_list &
         (reg_data) result (emitters)
      class(region_data_t), intent(in) :: reg_data
      integer, dimension(:), allocatable :: emitters
    end function region_data_get_emitter_list
    module function region_data_get_n_emitters_sc &
         (reg_data) result (n_emitters_sc)
      class(region_data_t), intent(in) :: reg_data
      integer :: n_emitters_sc
    end function region_data_get_n_emitters_sc
    module function region_data_get_associated_resonances &
         (reg_data, emitter) result (res)
      integer, dimension(:), allocatable :: res
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: emitter
    end function region_data_get_associated_resonances
    module function region_data_emitter_is_compatible_with_resonance &
         (reg_data, i_res, emitter) result (compatible)
      logical :: compatible
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: i_res, emitter
    end function region_data_emitter_is_compatible_with_resonance
    module function region_data_emitter_is_in_resonance &
         (reg_data, i_res, emitter) result (exist)
      logical :: exist
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: i_res, emitter
    end function region_data_emitter_is_in_resonance
    module subroutine region_data_get_contributors &
         (reg_data, i_res, emitter, c, success)
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: i_res, emitter
      integer, intent(inout), dimension(:), allocatable :: c
      logical, intent(out) :: success
    end subroutine region_data_get_contributors
    pure module function region_data_get_emitter &
         (reg_data, alr) result (emitter)
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: alr
      integer :: emitter
    end function region_data_get_emitter
    module function region_data_map_real_to_born_index &
         (reg_data, real_index) result (uborn_index)
      integer :: uborn_index
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: real_index
    end function region_data_map_real_to_born_index
    module function region_data_get_flv_states_born_array &
         (reg_data) result (flv_states)
      integer, dimension(:,:), allocatable :: flv_states
      class(region_data_t), intent(in) :: reg_data
    end function region_data_get_flv_states_born_array
    module function region_data_get_flv_states_born_single &
         (reg_data, i_flv) result (flv_states)
      integer, dimension(:), allocatable :: flv_states
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: i_flv
    end function region_data_get_flv_states_born_single
    module function region_data_get_flv_states_real_single &
         (reg_data, i_flv) result (flv_states)
      integer, dimension(:), allocatable :: flv_states
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: i_flv
    end function region_data_get_flv_states_real_single
    module function region_data_get_flv_states_real_array &
         (reg_data) result (flv_states)
      integer, dimension(:,:), allocatable :: flv_states
      class(region_data_t), intent(in) :: reg_data
    end function region_data_get_flv_states_real_array
    module subroutine region_data_get_all_flv_states &
         (reg_data, flv_born, flv_real)
      class(region_data_t), intent(in) :: reg_data
      integer, dimension(:,:), allocatable, intent(out) :: flv_born, flv_real
    end subroutine region_data_get_all_flv_states
    module function region_data_get_n_in (reg_data) result (n_in)
      integer :: n_in
      class(region_data_t), intent(in) :: reg_data
    end function region_data_get_n_in
    module function region_data_get_n_legs_real (reg_data) result (n_legs)
      integer :: n_legs
      class(region_data_t), intent(in) :: reg_data
    end function region_data_get_n_legs_real
    module function region_data_get_n_legs_born (reg_data) result (n_legs)
      integer :: n_legs
      class(region_data_t), intent(in) :: reg_data
    end function region_data_get_n_legs_born
    module function region_data_get_n_flv_real (reg_data) result (n_flv)
      integer :: n_flv
      class(region_data_t), intent(in) :: reg_data
    end function region_data_get_n_flv_real
    module function region_data_get_n_flv_born (reg_data) result (n_flv)
      integer :: n_flv
      class(region_data_t), intent(in) :: reg_data
    end function region_data_get_n_flv_born
    module function region_data_get_svalue_ij &
         (reg_data, p_real, alr, i, j, i_res) result (sval)
      class(region_data_t), intent(inout) :: reg_data
      type(vector4_t), intent(in), dimension(:) :: p_real
      integer, intent(in) :: alr, i, j
      integer, intent(in) :: i_res
      real(default) :: sval
    end function region_data_get_svalue_ij
    module function region_data_get_svalue_last_pos &
         (reg_data, p, alr, emitter, i_res) result (sval)
      class(region_data_t), intent(inout) :: reg_data
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: alr, emitter
      integer, intent(in) :: i_res
      real(default) :: sval
    end function region_data_get_svalue_last_pos
    module function region_data_get_svalue_soft &
         (reg_data, p_born, p_soft, alr, emitter, i_res) result (sval)
      class(region_data_t), intent(inout) :: reg_data
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
      integer, intent(in) :: alr, emitter, i_res
      real(default) :: sval
    end function region_data_get_svalue_soft
    module subroutine region_data_find_regions &
         (reg_data, model, ftuples, emitters, flst_alr)
      class(region_data_t), intent(in) :: reg_data
      type(model_t), intent(in) :: model
      type(ftuple_list_t), intent(out), dimension(:), allocatable :: ftuples
      integer, intent(out), dimension(:), allocatable :: emitters
      type(flv_structure_t), intent(out), dimension(:), allocatable :: flst_alr
    end subroutine region_data_find_regions
    module subroutine region_data_find_eqv_regions (reg_data, optimize)
      class(region_data_t), intent(inout) :: reg_data
      logical, intent(in) :: optimize
    end subroutine region_data_find_eqv_regions
    module subroutine region_data_init_singular_regions &
         (reg_data, ftuples, emitter, flv_alr, nlo_correction_type)
      class(region_data_t), intent(inout) :: reg_data
      type(ftuple_list_t), intent(inout), dimension(:), allocatable :: ftuples
      type(string_t), intent(in) :: nlo_correction_type
      integer, intent(in), dimension(:) :: emitter
      type(flv_structure_t), intent(in), dimension(:) :: flv_alr
    end subroutine region_data_init_singular_regions
    module subroutine region_data_find_emitters (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_find_emitters
    module subroutine region_data_find_resonances (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_find_resonances
    module subroutine region_data_set_i_phs_to_i_con (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_set_i_phs_to_i_con
    module subroutine region_data_set_alr_to_i_phs &
         (reg_data, phs_identifiers, alr_to_i_phs)
      class(region_data_t), intent(inout) :: reg_data
      type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
      integer, intent(out), dimension(:) :: alr_to_i_phs
    end subroutine region_data_set_alr_to_i_phs
    module subroutine region_data_set_contributors (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_set_contributors
    module subroutine region_data_extend_ftuples (reg_data, n_res)
      class(region_data_t), intent(inout) :: reg_data
      integer, intent(in) :: n_res
    end subroutine region_data_extend_ftuples
    module function region_data_get_flavor_indices &
         (reg_data, born) result (i_flv)
      integer, dimension(:), allocatable :: i_flv
      class(region_data_t), intent(in) :: reg_data
      logical, intent(in) :: born
    end function region_data_get_flavor_indices
    module function region_data_get_matrix_element_index &
         (reg_data, i_reg) result (i_me)
      integer :: i_me
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: i_reg
    end function region_data_get_matrix_element_index
    module subroutine region_data_compute_number_of_phase_spaces (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_compute_number_of_phase_spaces
    module function region_data_get_n_phs (reg_data) result (n_phs)
      integer :: n_phs
      class(region_data_t), intent(in) :: reg_data
    end function region_data_get_n_phs
    module subroutine region_data_set_splitting_info (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_set_splitting_info
    module subroutine region_data_init_phs_identifiers (reg_data, phs_id)
      class(region_data_t), intent(in) :: reg_data
      type(phs_identifier_t), intent(out), dimension(:), allocatable :: phs_id
    end subroutine region_data_init_phs_identifiers
    module subroutine region_data_get_all_ftuples (reg_data, ftuples)
      class(region_data_t), intent(in) :: reg_data
      type(ftuple_t), intent(inout), dimension(:), allocatable :: ftuples
    end subroutine region_data_get_all_ftuples
    module subroutine region_data_write_to_file &
         (reg_data, proc_id, latex, os_data)
      class(region_data_t), intent(inout) :: reg_data
      type(string_t), intent(in) :: proc_id
      logical, intent(in) :: latex
      type(os_data_t), intent(in) :: os_data
    end subroutine region_data_write_to_file
    module subroutine region_data_write_latex (reg_data, unit)
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in), optional :: unit
    end subroutine region_data_write_latex
    module subroutine region_data_write (reg_data, unit)
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in), optional :: unit
    end subroutine region_data_write
    module subroutine region_data_assign (reg_data_out, reg_data_in)
      type(region_data_t), intent(out) :: reg_data_out
      type(region_data_t), intent(in) :: reg_data_in
    end subroutine region_data_assign
    module function region_data_has_pseudo_isr (reg_data) result (flag)
      logical :: flag
      class(region_data_t), intent(in) :: reg_data
    end function region_data_has_pseudo_isr
    module subroutine region_data_check_consistency (reg_data, fail_fatal, unit)
      class(region_data_t), intent(in) :: reg_data
      logical, intent(in) :: fail_fatal
      integer, intent(in), optional :: unit
    end subroutine region_data_check_consistency
    module function region_data_requires_spin_correlations &
         (reg_data) result (flag)
      class(region_data_t), intent(in) :: reg_data
      logical :: flag
    end function region_data_requires_spin_correlations
    module function region_data_born_to_real_symm_factor_fs &
         (reg_data, alr) result (factor)
      class(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: alr
      real(default) :: factor
    end function region_data_born_to_real_symm_factor_fs
    module subroutine region_data_final (reg_data)
      class(region_data_t), intent(inout) :: reg_data
    end subroutine region_data_final
    module subroutine fks_mapping_default_set_parameter &
         (map, n_in, dij_exp1, dij_exp2)
      class(fks_mapping_default_t), intent(inout) :: map
      integer, intent(in) :: n_in
      real(default), intent(in) :: dij_exp1, dij_exp2
    end subroutine fks_mapping_default_set_parameter
    module function fks_mapping_default_dij (map, p, i, j, i_con) result (d)
      real(default) :: d
      class(fks_mapping_default_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: i, j
      integer, intent(in), optional :: i_con
    end function fks_mapping_default_dij
    module subroutine fks_mapping_default_compute_sumdij (map, sregion, p_real)
      class(fks_mapping_default_t), intent(inout) :: map
      type(singular_region_t), intent(in) :: sregion
      type(vector4_t), intent(in), dimension(:) :: p_real
    end subroutine fks_mapping_default_compute_sumdij
    module function fks_mapping_default_svalue &
         (map, p, i, j, i_res) result (value)
      real(default) :: value
      class(fks_mapping_default_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: i, j
      integer, intent(in), optional :: i_res
    end function fks_mapping_default_svalue
    module function fks_mapping_default_dij_soft &
         (map, p_born, p_soft, em, i_con) result (d)
      real(default) :: d
      class(fks_mapping_default_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
      integer, intent(in) :: em
      integer, intent(in), optional :: i_con
    end function fks_mapping_default_dij_soft
    module subroutine fks_mapping_default_compute_sumdij_soft &
         (map, sregion, p_born, p_soft)
      class(fks_mapping_default_t), intent(inout) :: map
      type(singular_region_t), intent(in) :: sregion
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
    end subroutine fks_mapping_default_compute_sumdij_soft
    module function fks_mapping_default_svalue_soft &
         (map, p_born, p_soft, em, i_res) result (value)
      real(default) :: value
      class(fks_mapping_default_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
      integer, intent(in) :: em
      integer, intent(in), optional :: i_res
    end function fks_mapping_default_svalue_soft
    module subroutine fks_mapping_default_assign (fks_map_out, fks_map_in)
      type(fks_mapping_default_t), intent(out) :: fks_map_out
      type(fks_mapping_default_t), intent(in) :: fks_map_in
    end subroutine fks_mapping_default_assign
    module function fks_mapping_resonances_dij (map, p, i, j, i_con) result (d)
      real(default) :: d
      class(fks_mapping_resonances_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: i, j
      integer, intent(in), optional :: i_con
    end function fks_mapping_resonances_dij
    module subroutine fks_mapping_resonances_compute_sumdij &
         (map, sregion, p_real)
      class(fks_mapping_resonances_t), intent(inout) :: map
      type(singular_region_t), intent(in) :: sregion
      type(vector4_t), intent(in), dimension(:) :: p_real
    end subroutine fks_mapping_resonances_compute_sumdij
    module function fks_mapping_resonances_svalue &
         (map, p, i, j, i_res) result (value)
      real(default) :: value
      class(fks_mapping_resonances_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: i, j
      integer, intent(in), optional :: i_res
    end function fks_mapping_resonances_svalue
    module function fks_mapping_resonances_get_resonance_weight &
         (map, alr, p) result (pfr)
      real(default) :: pfr
      class(fks_mapping_resonances_t), intent(in) :: map
      integer, intent(in) :: alr
      type(vector4_t), intent(in), dimension(:) :: p
    end function fks_mapping_resonances_get_resonance_weight
    module function fks_mapping_resonances_dij_soft &
         (map, p_born, p_soft, em, i_con) result (d)
      real(default) :: d
      class(fks_mapping_resonances_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
      integer, intent(in) :: em
      integer, intent(in), optional :: i_con
    end function fks_mapping_resonances_dij_soft
    module subroutine fks_mapping_resonances_compute_sumdij_soft &
         (map, sregion, p_born, p_soft)
      class(fks_mapping_resonances_t), intent(inout) :: map
      type(singular_region_t), intent(in) :: sregion
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
    end subroutine fks_mapping_resonances_compute_sumdij_soft
    module function fks_mapping_resonances_svalue_soft &
         (map, p_born, p_soft, em, i_res) result (value)
      real(default) :: value
      class(fks_mapping_resonances_t), intent(in) :: map
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
      integer, intent(in) :: em
      integer, intent(in), optional :: i_res
    end function fks_mapping_resonances_svalue_soft
    module subroutine fks_mapping_resonances_set_resonance_momentum (map, p)
      class(fks_mapping_resonances_t), intent(inout) :: map
      type(vector4_t), intent(in) :: p
    end subroutine fks_mapping_resonances_set_resonance_momentum
    module subroutine fks_mapping_resonances_set_resonance_momenta (map, p)
      class(fks_mapping_resonances_t), intent(inout) :: map
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine fks_mapping_resonances_set_resonance_momenta
    module subroutine fks_mapping_resonances_assign (fks_map_out, fks_map_in)
      type(fks_mapping_resonances_t), intent(out) :: fks_map_out
      type(fks_mapping_resonances_t), intent(in) :: fks_map_in
    end subroutine fks_mapping_resonances_assign
    module function create_resonance_histories_for_threshold &
         () result (res_history)
      type(resonance_history_t) :: res_history
    end function create_resonance_histories_for_threshold
    module subroutine setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, &
         nlo_corr_type, alpha_pow, alphas_pow)
      integer, intent(in) :: n_in, alpha_pow, alphas_pow
      integer, intent(in), dimension(:,:) :: flv_born, flv_real
      type(string_t), intent(in) :: nlo_corr_type
      type(region_data_t), intent(out) :: reg_data
    end subroutine setup_region_data_for_test
  end interface

contains

  subroutine region_data_allocate_fks_mappings (reg_data, mapping_type)
    class(region_data_t), intent(inout) :: reg_data
    integer, intent(in) :: mapping_type
    select case (mapping_type)
    case (FKS_DEFAULT)
       allocate (fks_mapping_default_t :: reg_data%fks_mapping)
    case (FKS_RESONANCES)
       allocate (fks_mapping_resonances_t :: reg_data%fks_mapping)
    case default
       call msg_fatal ("Init region_data: FKS mapping not implemented!")
    end select
  end subroutine region_data_allocate_fks_mappings


end module fks_regions

