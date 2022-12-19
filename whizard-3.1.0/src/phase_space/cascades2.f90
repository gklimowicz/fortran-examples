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

module cascades2

  use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
  use kinds, only: default
  use kinds, only: TC, i8
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use diagnostics
  use flavors
  use model_data
  use phs_forests, only: phs_parameters_t
  use resonances, only: resonance_info_t
  use resonances, only: resonance_history_t
  use resonances, only: resonance_history_set_t
  use cascades2_lexer

  implicit none
  private

  public :: DECAY, SCATTERING
  public :: feyngraph_set_t
  public :: init_sm_full_test
  public :: feyngraph_set_write_file_format
  public :: feyngraph_set_generate_single
  public :: feyngraph_set_write_process_bincode_format
  public :: feyngraph_set_write_graph_format
  public :: feyngraph_set_generate
  public :: feyngraph_set_is_valid
  public :: grove_list_get_n_trees
  public :: feyngraph_set_get_resonance_histories

  integer, parameter :: LABEL_LEN=30
  integer, parameter :: &
       & NONRESONANT = -2, EXTERNAL_PRT = -1, &
       & NO_MAPPING = 0, S_CHANNEL = 1, T_CHANNEL =  2, U_CHANNEL = 3, &
       & RADIATION = 4, COLLINEAR = 5, INFRARED = 6, &
       & STEP_MAPPING_E = 11, STEP_MAPPING_H = 12, &
       & ON_SHELL = 99
  integer, parameter :: FEYNGRAPH_LEN=300
  integer, parameter :: DECAY=1, SCATTERING=2

  integer, parameter :: BUFFER_LEN = 1000
  integer, parameter :: STACK_SIZE = 100
  integer, parameter :: PRT_ARRAY_SIZE = 200
  integer, parameter :: DAG_STACK_SIZE = 1000
  integer, parameter :: EMPTY = -999

  type :: part_prop_t
     character(len=LABEL_LEN) :: particle_label
     integer :: pdg = 0
     real(default) :: mass = 0.
     real :: width = 0.
     integer :: spin_type = 0
     logical :: is_vector = .false.
     logical :: empty = .true.
     type(part_prop_t), pointer :: anti => null ()
     type(string_t) :: tex_name
   contains
     procedure :: final => part_prop_final
     procedure :: init => part_prop_init
  end type part_prop_t

  type :: grove_prop_t
     integer :: multiplicity = 0
     integer :: n_resonances = 0
     integer :: n_log_enhanced = 0
     integer :: n_off_shell = 0
     integer :: n_t_channel = 0
     integer :: res_hash = 0
  end type grove_prop_t

  type :: tree_t
     integer(TC), dimension(:), allocatable :: bc
     integer, dimension(:), allocatable :: pdg
     integer, dimension(:), allocatable :: mapping
     integer :: n_entries = 0
     logical :: keep = .true.
     logical :: empty = .true.
   contains
     procedure :: final => tree_final
     procedure :: add_entry_from_numbers => tree_add_entry_from_numbers
     procedure :: add_entry_from_node => tree_add_entry_from_node
     generic :: add_entry =>  add_entry_from_numbers, add_entry_from_node
     procedure :: sort => tree_sort
  end type tree_t

  type, abstract :: graph_t
     integer :: index = 0
     integer :: n_nodes = 0
     logical :: keep = .true.
  end type graph_t

  type, extends (graph_t) :: feyngraph_t
     type(string_t) :: omega_feyngraph_output
     type(f_node_t), pointer :: root => null ()
     type(feyngraph_t), pointer :: next => null()
     type(kingraph_t), pointer :: kin_first => null ()
     type(kingraph_t), pointer :: kin_last => null ()
   contains
     procedure :: final => feyngraph_final
     procedure :: make_kingraphs => feyngraph_make_kingraphs
     procedure :: make_inverse_kingraphs => feyngraph_make_inverse_kingraphs
     procedure :: compute_mappings => feyngraph_compute_mappings
     procedure :: make_invertible => feyngraph_make_invertible
  end type feyngraph_t

  type :: feyngraph_ptr_t
     type(feyngraph_t), pointer :: graph => null ()
  end type feyngraph_ptr_t

  type, extends (graph_t) :: kingraph_t
     type(k_node_t), pointer :: root => null ()
     type(kingraph_t), pointer :: next => null()
     type(kingraph_t), pointer :: grove_next => null ()
     type(tree_t) :: tree
     type(grove_prop_t) :: grove_prop
     logical :: inverse = .false.
     integer :: prc_component = 0
     contains
     procedure :: final => kingraph_final
     procedure :: write_file_format => kingraph_write_file_format
     procedure :: make_inverse_copy => kingraph_make_inverse_copy
     procedure :: assign_resonance_hash => kingraph_assign_resonance_hash
     procedure :: extract_resonance_history => kingraph_extract_resonance_history
  end type kingraph_t

  type :: kingraph_ptr_t
     type(kingraph_t), pointer :: graph => null ()
  end type kingraph_ptr_t

  type, abstract :: node_t
     type(part_prop_t), pointer :: particle => null ()
     logical :: incoming = .false.
     logical :: t_line = .false.
     integer :: index = 0
     logical :: keep = .true.
     integer :: n_subtree_nodes = 1
  end type node_t

  type, abstract :: list_t
     integer :: n_entries = 0
  end type list_t

  type :: k_node_entry_t
     type(k_node_t), pointer :: node => null ()
     type(k_node_entry_t), pointer :: next => null ()
     logical :: recycle = .false.
   contains
     procedure :: final => k_node_entry_final
     procedure :: write => k_node_entry_write
  end type k_node_entry_t

  type, extends (list_t) :: k_node_list_t
     type(k_node_entry_t), pointer :: first => null ()
     type(k_node_entry_t), pointer :: last => null ()
     integer :: n_recycle
     logical :: observer = .false.
   contains
     procedure :: final => k_node_list_final
     procedure :: add_entry => k_node_list_add_entry
     procedure :: add_pointer => k_node_list_add_pointer
     procedure :: check_subtree_equivalences => &
          k_node_list_check_subtree_equivalences
     procedure :: get_nodes => k_node_list_get_nodes
  end type k_node_list_t

  type, extends (node_t) :: f_node_t
     type(f_node_t), pointer :: daughter1 => null ()
     type(f_node_t), pointer :: daughter2 => null ()
     character(len=LABEL_LEN) :: particle_label
     type(k_node_list_t) :: k_node_list
   contains
     procedure :: final => f_node_final
     procedure :: set_index => f_node_set_index
     procedure :: assign_particle_properties => f_node_assign_particle_properties
  end type f_node_t

  type :: f_node_ptr_t
     type(f_node_t), pointer :: node => null ()
     contains
         procedure :: final => f_node_ptr_final
  end type f_node_ptr_t

  type :: k_node_ptr_t
     type(k_node_t), pointer :: node => null ()
  end type k_node_ptr_t

  type, extends (node_t) :: k_node_t
     type(k_node_t), pointer :: daughter1 => null ()
     type(k_node_t), pointer :: daughter2 => null ()
     type(k_node_t), pointer :: inverse_daughter1 => null ()
     type(k_node_t), pointer :: inverse_daughter2 => null ()
     type(f_node_t), pointer :: f_node => null ()
     type(tree_t) :: subtree
     real (default) :: ext_mass_sum = 0.
     real (default) :: effective_mass = 0.
     logical :: resonant = .false.
     logical :: on_shell = .false.
     logical :: log_enhanced = .false.
     integer :: mapping = NO_MAPPING
     integer(TC) :: bincode = 0
     logical :: mapping_assigned = .false.
     logical :: is_nonresonant_copy = .false.
     logical :: subtree_checked = .false.
     integer :: n_off_shell = 0
     integer :: n_log_enhanced = 0
     integer :: n_resonances = 0
     integer :: multiplicity = 0
     integer :: n_t_channel = 0
     integer :: f_node_index = 0
   contains
     procedure :: final => k_node_final
     procedure :: set_index => k_node_set_index
  end type k_node_t

  type :: f_node_entry_t
     character(len=FEYNGRAPH_LEN) :: subtree_string
     integer :: string_len = 0
     type(f_node_t), pointer :: node => null ()
     type(f_node_entry_t), pointer :: next => null ()
     integer :: subtree_size = 0
   contains
     procedure :: final => f_node_entry_final
     procedure :: write => f_node_entry_write
  end type f_node_entry_t

  type, extends (list_t) :: f_node_list_t
     type(f_node_entry_t), pointer :: first => null ()
     type(f_node_entry_t), pointer :: last => null ()
     type(k_node_list_t), pointer :: k_node_list => null ()
     integer :: max_tree_size = 0
   contains
     procedure :: add_entry => f_node_list_add_entry
     procedure :: write => f_node_list_write
     procedure :: final => f_node_list_final
  end type f_node_list_t

  type :: compare_tree_entry_t
     type(compare_tree_entry_t), dimension(:), pointer :: entry => null ()
     type(kingraph_ptr_t), dimension(:), allocatable :: graph_entry
   contains
       procedure :: final => compare_tree_entry_final
       procedure :: check_kingraph => compare_tree_entry_check_kingraph
  end type compare_tree_entry_t

  type :: compare_tree_t
    integer :: depth = 3
    type(compare_tree_entry_t), dimension(:), pointer :: entry => null ()
  contains
      procedure :: final => compare_tree_final
      procedure :: check_kingraph => compare_tree_check_kingraph
  end type compare_tree_t

  type :: grove_t
     type(grove_prop_t) :: grove_prop
     type(grove_t), pointer :: next => null ()
     type(kingraph_t), pointer :: first => null ()
     type(kingraph_t), pointer :: last => null ()
     type(compare_tree_t) :: compare_tree
   contains
     procedure :: final => grove_final
     procedure :: write_file_format => grove_write_file_format
  end type grove_t

  type :: grove_ptr_t
     type(grove_t), pointer :: grove => null ()
  end type grove_ptr_t

  type :: grove_list_t
     type(grove_t), pointer :: first => null ()
   contains
     procedure :: final => grove_list_final
     procedure :: get_grove => grove_list_get_grove
     procedure :: add_kingraph => grove_list_add_kingraph
     procedure :: add_feyngraph => grove_list_add_feyngraph
     procedure :: merge => grove_list_merge
     procedure :: rebuild => grove_list_rebuild
  end type grove_list_t

  type :: feyngraph_set_t
     type(model_data_t), pointer :: model => null ()
     type(flavor_t), dimension(:,:), allocatable :: flv
     integer :: n_in = 0
     integer :: n_out = 0
     integer :: process_type = DECAY
     type(phs_parameters_t) :: phs_par
     logical :: fatal_beam_decay = .true.
     type(part_prop_t), dimension (:), pointer :: particle => null ()
     type(f_node_list_t) :: f_node_list
     type(feyngraph_t), pointer :: first => null ()
     type(feyngraph_t), pointer :: last => null ()
     integer :: n_graphs = 0
     type(grove_list_t), pointer :: grove_list => null ()
     logical :: use_dag = .true.
     type(dag_t), pointer :: dag => null ()
     type(feyngraph_set_t), dimension (:), pointer :: fset => null ()
   contains
     procedure :: final => feyngraph_set_final
     procedure :: build => feyngraph_set_build
  end type feyngraph_set_t

  type :: dag_node_t
     integer :: string_len
     type(dag_string_t) :: string
     logical :: leaf = .false.
     type(f_node_ptr_t), dimension (:), allocatable :: f_node
     integer :: subtree_size = 0
   contains
       procedure :: final => dag_node_final
       procedure :: make_f_nodes => dag_node_make_f_nodes
  end type dag_node_t

  type :: dag_options_t
     integer :: string_len
     type(dag_string_t) :: string
     type(f_node_ptr_t), dimension (:), allocatable :: f_node_ptr1
     type(f_node_ptr_t), dimension (:), allocatable :: f_node_ptr2
   contains
       procedure :: final => dag_options_final
       procedure :: make_f_nodes => dag_options_make_f_nodes
  end type dag_options_t

  type :: dag_combination_t
     integer :: string_len
     type(dag_string_t) :: string
     integer, dimension (2) :: combination
     type(f_node_ptr_t), dimension (:), allocatable :: f_node_ptr1
     type(f_node_ptr_t), dimension (:), allocatable :: f_node_ptr2
   contains
       procedure :: final => dag_combination_final
       procedure :: make_f_nodes => dag_combination_make_f_nodes
  end type dag_combination_t

  type :: dag_t
     type(dag_string_t) :: string
     type(dag_node_t), dimension (:), allocatable :: node
     type(dag_options_t), dimension (:), allocatable :: options
     type(dag_combination_t), dimension (:), allocatable :: combination
     integer :: n_nodes = 0
     integer :: n_options = 0
     integer :: n_combinations = 0
   contains
       procedure :: read_string => dag_read_string
       procedure :: final => dag_final
       procedure :: construct => dag_construct
       procedure :: get_nodes_and_combinations => dag_get_nodes_and_combinations
       procedure :: get_options => dag_get_options
       procedure :: add_node => dag_add_node
       procedure :: add_options => dag_add_options
       procedure :: add_combination => dag_add_combination
       procedure :: make_feyngraphs => dag_make_feyngraphs
       procedure :: write => dag_write
  end type dag_t


  interface assignment (=)
     module procedure tree_assign
  end interface assignment (=)

  interface assignment (=)
     module procedure f_node_ptr_assign
  end interface assignment (=)
  interface assignment (=)
     module procedure k_node_assign
  end interface assignment (=)
  interface assignment (=)
     module procedure f_node_entry_assign
  end interface assignment (=)
  interface assignment (=)
     module procedure k_node_entry_assign
  end interface assignment (=)
  interface operator (.match.)
     module procedure grove_prop_match
  end interface operator (.match.)
  interface operator (==)
     module procedure grove_prop_equal
  end interface operator (==)
  interface operator (==)
     module procedure tree_equal
  end interface operator (==)
  interface operator (.eqv.)
     module procedure subtree_eqv
  end interface operator (.eqv.)

  interface
    module subroutine part_prop_final (part)
      class(part_prop_t), intent(inout) :: part
    end subroutine part_prop_final
    module subroutine tree_final (tree)
      class(tree_t), intent(inout) :: tree
    end subroutine tree_final
    module subroutine tree_assign (tree1, tree2)
      type(tree_t), intent(inout) :: tree1
      type(tree_t), intent(in) :: tree2
    end subroutine tree_assign
    module subroutine tree_add_entry_from_numbers (tree, bincode, pdg, mapping)
      class(tree_t), intent(inout) :: tree
      integer(TC), intent(in) :: bincode
      integer, intent(in) :: pdg
      integer, intent(in) :: mapping
    end subroutine tree_add_entry_from_numbers
    module subroutine tree_add_entry_from_node (tree, node)
      class(tree_t), intent(inout) :: tree
      type(k_node_t), intent(in) :: node
    end subroutine tree_add_entry_from_node
    module subroutine tree_sort (tree)
      class(tree_t), intent(inout) :: tree
    end subroutine tree_sort
    module subroutine feyngraph_final (graph)
      class(feyngraph_t), intent(inout) :: graph
    end subroutine feyngraph_final
    module subroutine kingraph_final (graph)
      class(kingraph_t), intent(inout) :: graph
    end subroutine kingraph_final
    module subroutine k_node_entry_final (entry)
      class(k_node_entry_t), intent(inout) :: entry
    end subroutine k_node_entry_final
    module subroutine k_node_entry_write (k_node_entry, u)
      class(k_node_entry_t), intent(in) :: k_node_entry
      integer, intent(in) :: u
    end subroutine k_node_entry_write
    module subroutine k_node_list_final (list)
      class(k_node_list_t), intent(inout) :: list
    end subroutine k_node_list_final
    recursive module subroutine f_node_final (node)
      class(f_node_t), intent(inout) :: node
    end subroutine f_node_final
    module subroutine f_node_entry_final (entry)
      class(f_node_entry_t), intent(inout) :: entry
    end subroutine f_node_entry_final
    module subroutine f_node_set_index (f_node)
      class(f_node_t), intent(inout) :: f_node
    end subroutine f_node_set_index
    module subroutine f_node_ptr_final (f_node_ptr)
      class(f_node_ptr_t), intent(inout) :: f_node_ptr
    end subroutine f_node_ptr_final
    module subroutine f_node_ptr_assign (ptr1, ptr2)
      type(f_node_ptr_t), intent(out) :: ptr1
      type(f_node_ptr_t), intent(in) :: ptr2
    end subroutine f_node_ptr_assign
    module subroutine k_node_assign (k_node1, k_node2)
      type(k_node_t), intent(inout) :: k_node1
      type(k_node_t), intent(in) :: k_node2
    end subroutine k_node_assign
    recursive module subroutine k_node_final (k_node)
      class(k_node_t), intent(inout) :: k_node
    end subroutine k_node_final
    module subroutine k_node_set_index (k_node)
      class(k_node_t), intent(inout) :: k_node
    end subroutine k_node_set_index
    module subroutine f_node_entry_write (f_node_entry, u)
      class(f_node_entry_t), intent(in) :: f_node_entry
      integer, intent(in) :: u
    end subroutine f_node_entry_write
    module subroutine f_node_entry_assign (entry1, entry2)
      type(f_node_entry_t), intent(out) :: entry1
      type(f_node_entry_t), intent(in) :: entry2
    end subroutine f_node_entry_assign
    module subroutine f_node_list_add_entry (list, subtree_string, &
         ptr_to_node, recycle, subtree_size)
      class(f_node_list_t), intent(inout) :: list
      character(len=*), intent(in) :: subtree_string
      type(f_node_t), pointer, intent(out) :: ptr_to_node
      logical, intent(in) :: recycle
      integer, intent(in), optional :: subtree_size
    end subroutine f_node_list_add_entry
    module subroutine f_node_list_write (f_node_list, u)
      class(f_node_list_t), intent(in) :: f_node_list
      integer, intent(in) :: u
    end subroutine f_node_list_write
    module subroutine k_node_entry_assign (entry1, entry2)
      type(k_node_entry_t), intent(out) :: entry1
      type(k_node_entry_t), intent(in) :: entry2
    end subroutine k_node_entry_assign
    recursive module subroutine k_node_list_add_entry &
         (list, ptr_to_node, recycle)
      class(k_node_list_t), intent(inout) :: list
      type(k_node_t), pointer, intent(out) :: ptr_to_node
      logical, intent(in) :: recycle
    end subroutine k_node_list_add_entry
    module subroutine k_node_list_add_pointer (list, ptr_to_node, recycle)
      class(k_node_list_t), intent(inout) :: list
      type(k_node_t), pointer, intent(in) :: ptr_to_node
      logical, optional, intent(in) :: recycle
    end subroutine k_node_list_add_pointer
    module subroutine k_node_list_check_subtree_equivalences (list, model)
      class(k_node_list_t), intent(inout) :: list
      type(model_data_t), intent(in) :: model
    end subroutine k_node_list_check_subtree_equivalences
    module subroutine k_node_list_get_nodes (list, nodes)
      class(k_node_list_t), intent(inout) :: list
      type(k_node_ptr_t), dimension(:), allocatable, intent(out) :: nodes
    end subroutine k_node_list_get_nodes
    module subroutine compare_tree_final (ctree)
      class(compare_tree_t), intent(inout) :: ctree
    end subroutine compare_tree_final
    recursive module subroutine compare_tree_entry_final (ct_entry)
      class(compare_tree_entry_t), intent(inout) :: ct_entry
    end subroutine compare_tree_entry_final
    module subroutine compare_tree_check_kingraph &
         (ctree, kingraph, model, preliminary)
      class(compare_tree_t), intent(inout) :: ctree
      type(kingraph_t), intent(inout), pointer :: kingraph
      type(model_data_t), intent(in) :: model
      logical, intent(in) :: preliminary
    end subroutine compare_tree_check_kingraph
    recursive module subroutine compare_tree_entry_check_kingraph (ct_entry, &
         kingraph, model, preliminary, subtree_size, identifier)
      class(compare_tree_entry_t), intent(inout) :: ct_entry
      type(kingraph_t), pointer, intent(inout) :: kingraph
      type(model_data_t), intent(in) :: model
      logical, intent(in) :: preliminary
      integer, intent(in), optional :: subtree_size
      integer, dimension (:), intent(in), optional :: identifier
    end subroutine compare_tree_entry_check_kingraph
    module subroutine grove_final (grove)
      class(grove_t), intent(inout) :: grove
    end subroutine grove_final
    module subroutine feyngraph_set_build (feyngraph_set, u_in)
      class(feyngraph_set_t), intent(inout) :: feyngraph_set
      integer, intent(in) :: u_in
    end subroutine feyngraph_set_build
    module subroutine dag_read_string (dag, u_in, flv)
      class(dag_t), intent(inout) :: dag
      integer, intent(in) :: u_in
      type(flavor_t), dimension(:), intent(in) :: flv
    end subroutine dag_read_string
    module subroutine init_sm_full_test (model)
      class(model_data_t), intent(out) :: model
    end subroutine init_sm_full_test
    recursive module subroutine part_prop_init &
         (part_prop, feyngraph_set, particle_label)
      class(part_prop_t), intent(out), target :: part_prop
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
      character(len=*), intent(in) :: particle_label
    end subroutine part_prop_init
    module subroutine f_node_assign_particle_properties (node, feyngraph_set)
      class(f_node_t), intent(inout ) :: node
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
    end subroutine f_node_assign_particle_properties
    module subroutine dag_node_final (dag_node)
      class(dag_node_t), intent(inout) :: dag_node
    end subroutine dag_node_final
    module subroutine dag_options_final (dag_options)
      class(dag_options_t), intent(inout) :: dag_options
    end subroutine dag_options_final
    module subroutine dag_combination_final (dag_combination)
      class(dag_combination_t), intent(inout) :: dag_combination
    end subroutine dag_combination_final
    module subroutine dag_final (dag)
      class(dag_t), intent(inout) :: dag
    end subroutine dag_final
    module subroutine dag_construct (dag, feyngraph_set)
      class(dag_t), intent(inout) :: dag
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
    end subroutine dag_construct
      module subroutine dag_get_nodes_and_combinations (dag, leaves)
        class(dag_t), intent(inout) :: dag
        logical, intent(in) :: leaves
      end subroutine dag_get_nodes_and_combinations
      module subroutine dag_get_options (dag)
        class(dag_t), intent(inout) :: dag
      end subroutine dag_get_options
    module subroutine dag_add_node (dag, string, leaf, i_node)
      class(dag_t), intent(inout) :: dag
      type(dag_token_t), dimension (:), intent(in) :: string
      logical, intent(in) :: leaf
      integer, intent(out), optional :: i_node
    end subroutine dag_add_node
    module subroutine dag_add_options (dag, string, i_options)
      class(dag_t), intent(inout) :: dag
      type(dag_token_t), dimension (:), intent(in) :: string
      integer, intent(out), optional :: i_options
    end subroutine dag_add_options
    module subroutine dag_add_combination (dag, string, i_combination)
      class(dag_t), intent(inout) :: dag
      type(dag_token_t), dimension (:), intent(in) :: string
      integer, intent(out), optional :: i_combination
    end subroutine dag_add_combination
    module subroutine dag_node_make_f_nodes (dag_node, feyngraph_set, dag)
      class(dag_node_t), intent(inout) :: dag_node
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
      type(dag_t), intent(inout) :: dag
    end subroutine dag_node_make_f_nodes
    module subroutine dag_options_make_f_nodes (dag_options, &
         feyngraph_set, dag)
      class(dag_options_t), intent(inout) :: dag_options
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
      type(dag_t), intent(inout) :: dag
    end subroutine dag_options_make_f_nodes
    module subroutine dag_combination_make_f_nodes (dag_combination, &
         feyngraph_set, dag)
      class(dag_combination_t), intent(inout) :: dag_combination
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
      type(dag_t), intent(inout) :: dag
    end subroutine dag_combination_make_f_nodes
    module subroutine dag_make_feyngraphs (dag, feyngraph_set)
      class(dag_t), intent(inout) :: dag
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
    end subroutine dag_make_feyngraphs
    module subroutine dag_write (dag, u)
      class(dag_t), intent(in) :: dag
      integer, intent(in) :: u
    end subroutine dag_write
    module subroutine feyngraph_make_kingraphs (feyngraph, feyngraph_set)
      class(feyngraph_t), intent(inout) :: feyngraph
      type(feyngraph_set_t), intent(in) :: feyngraph_set
    end subroutine feyngraph_make_kingraphs
    module subroutine feyngraph_make_inverse_kingraphs (feyngraph)
      class(feyngraph_t), intent(inout) :: feyngraph
    end subroutine feyngraph_make_inverse_kingraphs
    module subroutine feyngraph_compute_mappings (feyngraph, feyngraph_set)
      class(feyngraph_t), intent(inout) :: feyngraph
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
    end subroutine feyngraph_compute_mappings
    module subroutine grove_list_get_grove (grove_list, kingraph, &
         return_grove, preliminary)
      class(grove_list_t), intent(inout) :: grove_list
      type(kingraph_t), intent(in), pointer :: kingraph
      type(grove_t), intent(inout), pointer :: return_grove
      logical, intent(in) :: preliminary
    end subroutine grove_list_get_grove
    module subroutine grove_list_add_kingraph (grove_list, kingraph, &
         preliminary, check, model)
      class(grove_list_t), intent(inout) :: grove_list
      type(kingraph_t), pointer, intent(inout) :: kingraph
      logical, intent(in) :: preliminary
      logical, intent(in) :: check
      type(model_data_t), optional, intent(in) :: model
    end subroutine grove_list_add_kingraph
    module subroutine grove_list_add_feyngraph (grove_list, feyngraph, model)
      class(grove_list_t), intent(inout) :: grove_list
      type(feyngraph_t), intent(inout) :: feyngraph
      type(model_data_t), intent(in) :: model
    end subroutine grove_list_add_feyngraph
    module function grove_prop_match (grove_prop1, grove_prop2) &
         result (gp_match)
      type(grove_prop_t), intent(in) :: grove_prop1
      type(grove_prop_t), intent(in) :: grove_prop2
      logical :: gp_match
    end function grove_prop_match
    module function grove_prop_equal (grove_prop1, grove_prop2) &
         result (gp_equal)
      type(grove_prop_t), intent(in) :: grove_prop1
      type(grove_prop_t), intent(in) :: grove_prop2
      logical :: gp_equal
    end function grove_prop_equal
    module subroutine grove_list_merge (target_list, grove_list, model, &
         prc_component)
      class(grove_list_t), intent(inout) :: target_list
      type(grove_list_t), intent(inout) :: grove_list
      type(model_data_t), intent(in) :: model
      integer, intent(in) :: prc_component
    end subroutine grove_list_merge
    module subroutine grove_list_rebuild (grove_list)
      class(grove_list_t), intent(inout) :: grove_list
    end subroutine grove_list_rebuild
    module subroutine feyngraph_set_write_file_format (feyngraph_set, u)
      type(feyngraph_set_t), intent(in) :: feyngraph_set
      integer, intent(in) :: u
    end subroutine feyngraph_set_write_file_format
    recursive module subroutine grove_write_file_format &
         (grove, feyngraph_set, gr_number, ch_number, u)
      class(grove_t), intent(in) :: grove
      type(feyngraph_set_t), intent(in) :: feyngraph_set
      integer, intent(in) :: u
      integer, intent(inout) :: gr_number
      integer, intent(inout) :: ch_number
    end subroutine grove_write_file_format
    module subroutine kingraph_write_file_format &
         (kingraph, feyngraph_set, ch_number, u)
      class(kingraph_t), intent(in) :: kingraph
      type(feyngraph_set_t), intent(in) :: feyngraph_set
      integer, intent(in) :: ch_number
      integer, intent(in) :: u
    end subroutine kingraph_write_file_format
    module subroutine feyngraph_make_invertible (feyngraph)
      class(feyngraph_t), intent(inout) :: feyngraph
    end subroutine feyngraph_make_invertible
    module subroutine kingraph_make_inverse_copy (original_kingraph, feyngraph)
      class(kingraph_t), intent(inout) :: original_kingraph
      type(feyngraph_t), intent(inout) :: feyngraph
    end subroutine kingraph_make_inverse_copy
    module subroutine feyngraph_set_generate_single (feyngraph_set, model, &
         n_in, n_out, phs_par, fatal_beam_decay, u_in)
      type(feyngraph_set_t), intent(inout) :: feyngraph_set
      type(model_data_t), target, intent(in) :: model
      integer, intent(in) :: n_in, n_out
      type(phs_parameters_t), intent(in) :: phs_par
      logical, intent(in) :: fatal_beam_decay
      integer, intent(in) :: u_in
    end subroutine feyngraph_set_generate_single
    elemental module function tree_equal (tree1, tree2) result (flag)
      type(tree_t), intent(in) :: tree1, tree2
      logical :: flag
    end function tree_equal
    pure module function subtree_eqv (subtree1, subtree2) result (eqv)
      type(tree_t), intent(in) :: subtree1, subtree2
      logical :: eqv
    end function subtree_eqv
    module subroutine kingraph_assign_resonance_hash (kingraph)
      class(kingraph_t), intent(inout) :: kingraph
    end subroutine kingraph_assign_resonance_hash
    module subroutine feyngraph_set_write_process_bincode_format &
         (feyngraph_set, unit)
      type(feyngraph_set_t), intent(in), target :: feyngraph_set
      integer, intent(in), optional :: unit
    end subroutine feyngraph_set_write_process_bincode_format
    module subroutine feyngraph_set_write_graph_format &
         (feyngraph_set, filename, process_id, unit)
      type(feyngraph_set_t), intent(in), target :: feyngraph_set
      type(string_t), intent(in) :: filename, process_id
      integer, intent(in), optional :: unit
    end subroutine feyngraph_set_write_graph_format
    module subroutine feyngraph_set_generate &
      (feyngraph_set, model, n_in, n_out, flv, phs_par, fatal_beam_decay, &
      u_in, vis_channels, use_dag)
      type(feyngraph_set_t), intent(out) :: feyngraph_set
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in, n_out
      type(flavor_t), dimension(:,:), intent(in) :: flv
      type(phs_parameters_t), intent(in) :: phs_par
      logical, intent(in) :: fatal_beam_decay
      integer, intent(in) :: u_in
      logical, intent(in) :: vis_channels
      logical, optional, intent(in) :: use_dag
    end subroutine feyngraph_set_generate
    module function feyngraph_set_is_valid (feyngraph_set) result (flag)
      class(feyngraph_set_t), intent(in) :: feyngraph_set
      logical :: flag
    end function feyngraph_set_is_valid
    module subroutine kingraph_extract_resonance_history &
         (kingraph, res_hist, model, n_out)
      class(kingraph_t), intent(in), target :: kingraph
      type(resonance_history_t), intent(out) :: res_hist
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_out
    end subroutine kingraph_extract_resonance_history
    module function grove_list_get_n_trees (grove_list) result (n)
      class(grove_list_t), intent(in) :: grove_list
      integer :: n
    end function grove_list_get_n_trees
    module subroutine feyngraph_set_get_resonance_histories &
         (feyngraph_set, n_filter, res_hists)
      type(feyngraph_set_t), intent(in), target :: feyngraph_set
      integer, intent(in), optional :: n_filter
      type(resonance_history_t), dimension(:), allocatable, intent(out) :: &
           res_hists
    end subroutine feyngraph_set_get_resonance_histories
  end interface

contains

  subroutine f_node_list_final (list)
    class(f_node_list_t) :: list
    type(f_node_entry_t), pointer :: current
    list%k_node_list => null ()
    do while (associated (list%first))
       current => list%first
       list%first => list%first%next
       call current%final ()
       deallocate (current)
    end do
  end subroutine f_node_list_final

  subroutine grove_list_final (list)
    class(grove_list_t), intent(inout) :: list
    class(grove_t), pointer :: current
    do while (associated (list%first))
       current => list%first
       list%first => list%first%next
       call current%final ()
       deallocate (current)
    end do
  end subroutine grove_list_final

  recursive subroutine feyngraph_set_final (set)
    class(feyngraph_set_t), intent(inout) :: set
    class(feyngraph_t), pointer :: current
    integer :: i
    if (associated (set%fset)) then
       do i=1, size (set%fset)
          call set%fset(i)%final ()
       end do
       deallocate (set%fset)
    else
       set%particle => null ()
       set%grove_list => null ()
    end if
    set%model => null ()
    if (allocated (set%flv)) deallocate (set%flv)
    set%last => null ()
    do while (associated (set%first))
       current => set%first
       set%first => set%first%next
       call current%final ()
       deallocate (current)
    end do
    if (associated (set%particle)) then
       do i = 1, size (set%particle)
          call set%particle(i)%final ()
       end do
       deallocate (set%particle)
    end if
    if (associated (set%grove_list)) then
       if (debug_on) call msg_debug (D_PHASESPACE, "grove_list: final")
       call set%grove_list%final ()
       deallocate (set%grove_list)
    end if
    if (debug_on) call msg_debug (D_PHASESPACE, "f_node_list: final")
    call set%f_node_list%final ()
    if (associated (set%dag)) then
       if (debug_on) call msg_debug (D_PHASESPACE, "dag: final")
       if (associated (set%dag)) then
          call set%dag%final ()
          deallocate (set%dag)
       end if
    end if
  end subroutine feyngraph_set_final


end module cascades2

