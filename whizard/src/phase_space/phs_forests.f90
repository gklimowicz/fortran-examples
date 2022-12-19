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

module phs_forests

  use kinds, only: default
  use kinds, only: TC
  use iso_varying_string, string_t => varying_string
  use lorentz
  use permutations
  use syntax_rules
  use parser
  use model_data
  use flavors
  use interactions

  use phs_base
  use resonances, only: resonance_history_t
  use resonances, only: resonance_history_set_t
  use mappings
  use phs_trees

  implicit none
  private

  public :: phs_parameters_t
  public :: phs_forest_t
  public :: assignment(=)
  public :: syntax_phs_forest_init
  public :: syntax_phs_forest_final
  public :: syntax_phs_forest_write

  type :: phs_parameters_t
     real(default) :: sqrts = 0
     real(default) :: m_threshold_s = 50._default
     real(default) :: m_threshold_t = 100._default
     integer :: off_shell = 1
     integer :: t_channel = 2
     logical :: keep_nonresonant = .true.
   contains
     procedure :: write => phs_parameters_write
     procedure :: read => phs_parameters_read
  end type phs_parameters_t

  type :: equivalence_t
     private
     integer :: left, right
     type(permutation_t) :: perm
     type(permutation_t) :: msq_perm, angle_perm
     logical, dimension(:), allocatable :: angle_sig
     type(equivalence_t), pointer :: next => null ()
  end type equivalence_t

  type :: equivalence_list_t
     private
     integer :: length = 0
     type(equivalence_t), pointer :: first => null ()
     type(equivalence_t), pointer :: last => null ()
  end type equivalence_list_t

  type :: phs_grove_t
     private
     integer :: tree_count_offset
     type(phs_tree_t), dimension(:), allocatable :: tree
     type(equivalence_list_t) :: equivalence_list
  end type phs_grove_t

  type :: phs_forest_t
     private
     integer :: n_in, n_out, n_tot
     integer :: n_masses, n_angles, n_dimensions
     integer :: n_trees, n_equivalences
     type(flavor_t), dimension(:), allocatable :: flv
     type(phs_grove_t), dimension(:), allocatable :: grove
     integer, dimension(:), allocatable :: grove_lookup
     type(phs_prt_t), dimension(:), allocatable :: prt_in
     type(phs_prt_t), dimension(:), allocatable :: prt_out
     type(phs_prt_t), dimension(:), allocatable :: prt
     integer(TC), dimension(:,:), allocatable :: prt_combination
     type(mapping_t), dimension(:), allocatable :: s_mapping
   contains
     procedure :: init => phs_forest_init
     procedure :: set_s_mappings => phs_forest_set_s_mappings
     procedure :: final => phs_forest_final
     procedure :: write => phs_forest_write
     procedure :: get_n_parameters => phs_forest_get_n_parameters
     procedure :: get_n_channels => phs_forest_get_n_channels
     procedure :: get_n_groves => phs_forest_get_n_groves
     procedure :: get_grove_bounds => phs_forest_get_grove_bounds
     procedure :: get_n_equivalences => phs_forest_get_n_equivalences
     procedure :: get_s_mapping => phs_forest_get_s_mapping
     procedure :: get_on_shell => phs_forest_get_on_shell
     procedure :: extract_resonance_history_set &
          => phs_forest_extract_resonance_history_set
     generic :: read => read_file, read_unit, read_parse_tree
     procedure :: read_file => phs_forest_read_file
     procedure :: read_unit => phs_forest_read_unit
     procedure :: read_parse_tree => phs_forest_read_parse_tree
     procedure :: set_flavors => phs_forest_set_flavors
     procedure :: set_momentum_links => phs_forest_set_momentum_links
     procedure :: set_parameters => phs_forest_set_parameters
     procedure :: setup_prt_combinations => phs_forest_setup_prt_combinations
     generic :: set_prt_in => set_prt_in_int, set_prt_in_mom
     procedure :: set_prt_in_int => phs_forest_set_prt_in_int
     procedure :: set_prt_in_mom => phs_forest_set_prt_in_mom
     generic :: set_prt_out => set_prt_out_int, set_prt_out_mom
     procedure :: set_prt_out_int => phs_forest_set_prt_out_int
     procedure :: set_prt_out_mom => phs_forest_set_prt_out_mom
     procedure :: combine_particles => phs_forest_combine_particles
     procedure :: get_prt_out => phs_forest_get_prt_out
     procedure :: get_momenta_out => phs_forest_get_momenta_out
     procedure :: set_equivalences => phs_forest_set_equivalences
     procedure :: get_equivalences => phs_forest_get_equivalences
     procedure :: evaluate_selected_channel => phs_forest_evaluate_selected_channel
     procedure :: evaluate_other_channels => phs_forest_evaluate_other_channels
     procedure :: recover_channel => phs_forest_recover_channel
  end type phs_forest_t


  interface operator(==)
     module procedure phs_parameters_eq
  end interface
  interface operator(/=)
     module procedure phs_parameters_ne
  end interface
  interface assignment(=)
     module procedure equivalence_list_assign
  end interface

  interface assignment(=)
     module procedure phs_grove_assign0
     module procedure phs_grove_assign1
  end interface

  interface assignment(=)
     module procedure phs_forest_assign
  end interface


  type(syntax_t), target, save :: syntax_phs_forest


  interface
    module subroutine phs_parameters_write (phs_par, unit)
      class(phs_parameters_t), intent(in) :: phs_par
      integer, intent(in), optional :: unit
    end subroutine phs_parameters_write
    module subroutine phs_parameters_read (phs_par, unit)
      class(phs_parameters_t), intent(out) :: phs_par
      integer, intent(in) :: unit
    end subroutine phs_parameters_read
    module function phs_parameters_eq (phs_par1, phs_par2) result (equal)
      logical :: equal
      type(phs_parameters_t), intent(in) :: phs_par1, phs_par2
    end function phs_parameters_eq
    module function phs_parameters_ne (phs_par1, phs_par2) result (ne)
      logical :: ne
      type(phs_parameters_t), intent(in) :: phs_par1, phs_par2
    end function phs_parameters_ne
    module subroutine phs_forest_init (forest, n_tree, n_in, n_out)
      class(phs_forest_t), intent(inout) :: forest
      integer, dimension(:), intent(in) :: n_tree
      integer, intent(in) :: n_in, n_out
    end subroutine phs_forest_init
    module subroutine phs_forest_set_s_mappings (forest)
      class(phs_forest_t), intent(inout) :: forest
    end subroutine phs_forest_set_s_mappings
    module subroutine phs_forest_final (forest)
      class(phs_forest_t), intent(inout) :: forest
    end subroutine phs_forest_final
    module subroutine phs_forest_write (forest, unit)
      class(phs_forest_t), intent(in) :: forest
      integer, intent(in), optional :: unit
    end subroutine phs_forest_write
    module subroutine phs_forest_assign (forest_out, forest_in)
      type(phs_forest_t), intent(out) :: forest_out
      type(phs_forest_t), intent(in) :: forest_in
    end subroutine phs_forest_assign
    module function phs_forest_get_n_parameters (forest) result (n)
      integer :: n
      class(phs_forest_t), intent(in) :: forest
    end function phs_forest_get_n_parameters
    module function phs_forest_get_n_channels (forest) result (n)
      integer :: n
      class(phs_forest_t), intent(in) :: forest
    end function phs_forest_get_n_channels
    module function phs_forest_get_n_groves (forest) result (n)
      integer :: n
      class(phs_forest_t), intent(in) :: forest
    end function phs_forest_get_n_groves
    module subroutine phs_forest_get_grove_bounds (forest, g, i0, i1, n)
      class(phs_forest_t), intent(in) :: forest
      integer, intent(in) :: g
      integer, intent(out) :: i0, i1, n
    end subroutine phs_forest_get_grove_bounds
    module function phs_forest_get_n_equivalences (forest) result (n)
      integer :: n
      class(phs_forest_t), intent(in) :: forest
    end function phs_forest_get_n_equivalences
    module subroutine phs_forest_get_s_mapping &
         (forest, channel, flag, mass, width)
      class(phs_forest_t), intent(in) :: forest
      integer, intent(in) :: channel
      logical, intent(out) :: flag
      real(default), intent(out) :: mass, width
    end subroutine phs_forest_get_s_mapping
    module subroutine phs_forest_get_on_shell (forest, channel, flag, mass)
      class(phs_forest_t), intent(in) :: forest
      integer, intent(in) :: channel
      logical, intent(out) :: flag
      real(default), intent(out) :: mass
    end subroutine phs_forest_get_on_shell
    module subroutine phs_forest_extract_resonance_history_set &
         (forest, res_set, include_trivial)
      class(phs_forest_t), intent(in) :: forest
      type(resonance_history_set_t), intent(out) :: res_set
      logical, intent(in), optional :: include_trivial
    end subroutine phs_forest_extract_resonance_history_set
    module subroutine syntax_phs_forest_init ()
    end subroutine syntax_phs_forest_init
    module subroutine syntax_phs_forest_final ()
    end subroutine syntax_phs_forest_final
    module subroutine syntax_phs_forest_write (unit)
      integer, intent(in), optional :: unit
    end subroutine syntax_phs_forest_write
    module subroutine phs_forest_read_file &
         (forest, filename, process_id, n_in, n_out, model, found, &
          md5sum_process, md5sum_model_par, &
          md5sum_phs_config, phs_par, match)
      class(phs_forest_t), intent(out) :: forest
      type(string_t), intent(in) :: filename
      type(string_t), intent(in) :: process_id
      integer, intent(in) :: n_in, n_out
      class(model_data_t), intent(in), target :: model
      logical, intent(out) :: found
      character(32), intent(in), optional :: &
           md5sum_process, md5sum_model_par, md5sum_phs_config
      type(phs_parameters_t), intent(in), optional :: phs_par
      logical, intent(out), optional :: match
    end subroutine phs_forest_read_file
    module subroutine phs_forest_read_unit &
         (forest, unit, process_id, n_in, n_out, model, found, &
          md5sum_process, md5sum_model_par, md5sum_phs_config, &
          phs_par, match)
      class(phs_forest_t), intent(out) :: forest
      integer, intent(in) :: unit
      type(string_t), intent(in) :: process_id
      integer, intent(in) :: n_in, n_out
      class(model_data_t), intent(in), target :: model
      logical, intent(out) :: found
      character(32), intent(in), optional :: &
           md5sum_process, md5sum_model_par, md5sum_phs_config
      type(phs_parameters_t), intent(in), optional :: phs_par
      logical, intent(out), optional :: match
    end subroutine phs_forest_read_unit
    module subroutine phs_forest_read_parse_tree &
         (forest, parse_tree, process_id, n_in, n_out, model, found, &
          md5sum_process, md5sum_model_par, md5sum_phs_config, &
          phs_par, match)
      class(phs_forest_t), intent(out) :: forest
      type(parse_tree_t), intent(in), target :: parse_tree
      type(string_t), intent(in) :: process_id
      integer, intent(in) :: n_in, n_out
      class(model_data_t), intent(in), target :: model
      logical, intent(out) :: found
      character(32), intent(in), optional :: &
           md5sum_process, md5sum_model_par, md5sum_phs_config
      type(phs_parameters_t), intent(in), optional :: phs_par
      logical, intent(out), optional :: match
    end subroutine phs_forest_read_parse_tree
    module subroutine phs_forest_set_flavors (forest, flv, reshuffle, flv_extra)
      class(phs_forest_t), intent(inout) :: forest
      type(flavor_t), dimension(:), intent(in) :: flv
      integer, intent(in), dimension(:), allocatable, optional :: reshuffle
      type(flavor_t), intent(in), optional :: flv_extra
    end subroutine phs_forest_set_flavors
    module subroutine phs_forest_set_momentum_links (forest, list)
      class(phs_forest_t), intent(inout) :: forest
      integer, intent(in), dimension(:), allocatable :: list
    end subroutine phs_forest_set_momentum_links
    module subroutine phs_forest_set_parameters &
         (forest, mapping_defaults, variable_limits)
      class(phs_forest_t), intent(inout) :: forest
      type(mapping_defaults_t), intent(in) :: mapping_defaults
      logical, intent(in) :: variable_limits
    end subroutine phs_forest_set_parameters
    module subroutine phs_forest_setup_prt_combinations (forest)
      class(phs_forest_t), intent(inout) :: forest
    end subroutine phs_forest_setup_prt_combinations
    module subroutine phs_forest_set_prt_in_int (forest, int, lt_cm_to_lab)
      class(phs_forest_t), intent(inout) :: forest
      type(interaction_t), intent(in) :: int
      type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    end subroutine phs_forest_set_prt_in_int
    module subroutine phs_forest_set_prt_in_mom (forest, mom, lt_cm_to_lab)
      class(phs_forest_t), intent(inout) :: forest
      type(vector4_t), dimension(size (forest%prt_in)), intent(in) :: mom
      type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    end subroutine phs_forest_set_prt_in_mom
    module subroutine phs_forest_set_prt_out_int (forest, int, lt_cm_to_lab)
      class(phs_forest_t), intent(inout) :: forest
      type(interaction_t), intent(in) :: int
      type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    end subroutine phs_forest_set_prt_out_int
    module subroutine phs_forest_set_prt_out_mom (forest, mom, lt_cm_to_lab)
      class(phs_forest_t), intent(inout) :: forest
      type(vector4_t), dimension(size (forest%prt_out)), intent(in) :: mom
      type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    end subroutine phs_forest_set_prt_out_mom
    module subroutine phs_forest_combine_particles (forest)
      class(phs_forest_t), intent(inout) :: forest
    end subroutine phs_forest_combine_particles
    module subroutine phs_forest_get_prt_out (forest, int, lt_cm_to_lab)
      class(phs_forest_t), intent(in) :: forest
      type(interaction_t), intent(inout) :: int
      type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    end subroutine phs_forest_get_prt_out
    module function phs_forest_get_momenta_out (forest, lt_cm_to_lab) result (p)
      class(phs_forest_t), intent(in) :: forest
      type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
      type(vector4_t), dimension(size (forest%prt_out)) :: p
    end function phs_forest_get_momenta_out
    module subroutine phs_forest_set_equivalences (forest)
      class(phs_forest_t), intent(inout) :: forest
    end subroutine phs_forest_set_equivalences
    module subroutine phs_forest_get_equivalences &
         (forest, channel, azimuthal_dependence)
      class(phs_forest_t), intent(in) :: forest
      type(phs_channel_t), dimension(:), intent(out) :: channel
      logical, intent(in) :: azimuthal_dependence
    end subroutine phs_forest_get_equivalences
    module subroutine phs_forest_evaluate_selected_channel &
         (forest, channel, active, sqrts, x, phs_factor, volume, ok)
      class(phs_forest_t), intent(inout) :: forest
      integer, intent(in) :: channel
      logical, dimension(:), intent(in) :: active
      real(default), intent(in) :: sqrts
      real(default), dimension(:,:), intent(inout) :: x
      real(default), dimension(:), intent(out) :: phs_factor
      real(default), intent(out) :: volume
      logical, intent(out) :: ok
    end subroutine phs_forest_evaluate_selected_channel
    module subroutine phs_forest_evaluate_other_channels &
         (forest, channel, active, sqrts, x, phs_factor, combine)
      class(phs_forest_t), intent(inout) :: forest
      integer, intent(in) :: channel
      logical, dimension(:), intent(in) :: active
      real(default), intent(in) :: sqrts
      real(default), dimension(:,:), intent(inout) :: x
      real(default), dimension(:), intent(inout) :: phs_factor
      logical, intent(in) :: combine
    end subroutine phs_forest_evaluate_other_channels
    module subroutine phs_forest_recover_channel &
         (forest, channel, sqrts, x, phs_factor, volume)
      class(phs_forest_t), intent(inout) :: forest
      integer, intent(in) :: channel
      real(default), intent(in) :: sqrts
      real(default), dimension(:,:), intent(inout) :: x
      real(default), dimension(:), intent(inout) :: phs_factor
      real(default), intent(out) :: volume
    end subroutine phs_forest_recover_channel
  end interface

contains

  subroutine equivalence_list_assign (eql_out, eql_in)
    type(equivalence_list_t), intent(out) :: eql_out
    type(equivalence_list_t), intent(in) :: eql_in
    type(equivalence_t), pointer :: eq, eq_copy
    eq => eql_in%first
    do while (associated (eq))
       allocate (eq_copy)
       eq_copy = eq
       eq_copy%next => null ()
       if (associated (eql_out%first)) then
          eql_out%last%next => eq_copy
       else
          eql_out%first => eq_copy
       end if
       eql_out%last => eq_copy
       eq => eq%next
    end do
  end subroutine equivalence_list_assign

  subroutine phs_grove_assign0 (grove_out, grove_in)
    type(phs_grove_t), intent(out) :: grove_out
    type(phs_grove_t), intent(in) :: grove_in
    grove_out%tree_count_offset = grove_in%tree_count_offset
    if (allocated (grove_in%tree)) then
       allocate (grove_out%tree (size (grove_in%tree)))
       grove_out%tree = grove_in%tree
    end if
    grove_out%equivalence_list = grove_in%equivalence_list
  end subroutine phs_grove_assign0

  subroutine phs_grove_assign1 (grove_out, grove_in)
    type(phs_grove_t), dimension(:), intent(out) :: grove_out
    type(phs_grove_t), dimension(:), intent(in) :: grove_in
    integer :: i
    do i = 1, size (grove_in)
       call phs_grove_assign0 (grove_out(i), grove_in(i))
    end do
  end subroutine phs_grove_assign1


end module phs_forests
