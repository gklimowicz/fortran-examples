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

module decays

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use diagnostics
  use flavors
  use interactions
  use evaluators
  use variables, only: var_list_t
  use model_data
  use rng_base
  use selectors
  use parton_states
  use process, only: process_t
  use instances, only: process_instance_t, pacify
  use process_stacks
  use event_transforms

  implicit none
  private

  public :: decay_root_config_t
  public :: decay_root_t
  public :: unstable_config_t
  public :: unstable_t
  public :: decay_chain_t
  public :: evt_decay_t
  public :: pacify

  type, abstract :: any_config_t
     private
   contains
     procedure (any_config_final), deferred :: final
     procedure (any_config_write), deferred :: write
  end type any_config_t

  type :: particle_config_t
     private
     class(any_config_t), allocatable :: c
  end type particle_config_t

  type, abstract :: any_t
     private
   contains
     procedure (any_final), deferred :: final
     procedure (any_write), deferred :: write
  end type any_t

  type :: particle_out_t
     private
     class(any_t), allocatable :: c
  end type particle_out_t

  type :: decay_term_config_t
     private
     type(particle_config_t), dimension(:), allocatable :: prt
   contains
     procedure :: final => decay_term_config_final
     procedure :: write => decay_term_config_write
     procedure :: init => decay_term_config_init
     procedure :: compute => decay_term_config_compute
  end type decay_term_config_t

  type :: decay_term_t
     private
     type(decay_term_config_t), pointer :: config => null ()
     type(particle_out_t), dimension(:), allocatable :: particle_out
   contains
     procedure :: final => decay_term_final
     procedure :: write => decay_term_write
     procedure :: write_process_instances => decay_term_write_process_instances
     procedure :: init => decay_term_init
     procedure :: make_rng => decay_term_make_rng
     procedure :: link_interactions => decay_term_link_interactions
     procedure :: select_chain => decay_term_select_chain
     procedure :: generate => decay_term_generate
  end type decay_term_t

  type :: decay_root_config_t
     private
     type(string_t) :: process_id
     type(process_t), pointer :: process => null ()
     class(model_data_t), pointer :: model => null ()
     type(decay_term_config_t), dimension(:), allocatable :: term_config
   contains
     procedure :: final => decay_root_config_final
     procedure :: write => decay_root_config_write
     procedure :: write_header => decay_root_config_write_header
     procedure :: write_terms => decay_root_config_write_terms
     procedure :: init => decay_root_config_init
     procedure :: init_term => decay_root_config_init_term
     procedure :: connect => decay_root_config_connect
     procedure :: compute => decay_root_config_compute
  end type decay_root_config_t

  type, abstract :: decay_gen_t
     private
     type(decay_term_t), dimension(:), allocatable :: term
     type(process_instance_t), pointer :: process_instance => null ()
     integer :: selected_mci = 0
     integer :: selected_term = 0
   contains
     procedure :: base_final => decay_gen_final
     procedure :: write_process_instances => decay_gen_write_process_instances
     procedure :: base_init => decay_gen_init
     procedure :: set_mci => decay_gen_set_mci
     procedure :: set_term => decay_gen_set_term
     procedure :: get_mci => decay_gen_get_mci
     procedure :: get_term => decay_gen_get_term
     procedure :: make_term_rng => decay_gen_make_term_rng
     procedure :: link_term_interactions => decay_gen_link_term_interactions
  end type decay_gen_t

  type, extends (decay_gen_t) :: decay_root_t
     private
     type(decay_root_config_t), pointer :: config => null ()
   contains
     procedure :: final => decay_root_final
     procedure :: write => decay_root_write
     procedure :: init => decay_root_init
     procedure :: select_chain => decay_root_select_chain
     procedure :: generate => decay_root_generate
  end type decay_root_t

  type, extends (decay_root_config_t) :: decay_config_t
     private
     type(flavor_t) :: flv
     real(default) :: weight = 0
     real(default) :: integral = 0
     real(default) :: abs_error = 0
     real(default) :: rel_error = 0
     type(selector_t) :: mci_selector
   contains
     procedure :: write => decay_config_write
     procedure :: connect => decay_config_connect
     procedure :: set_flv => decay_config_set_flv
     procedure :: compute => decay_config_compute
  end type decay_config_t

  type, extends (decay_gen_t) :: decay_t
     private
     type(decay_config_t), pointer :: config => null ()
     class(rng_t), allocatable :: rng
   contains
     procedure :: final => decay_final
     procedure :: write => decay_write
     procedure :: init => decay_init
     procedure :: link_interactions => decay_link_interactions
     procedure :: select_chain => decay_select_chain
     procedure :: generate => decay_generate
  end type decay_t

  type, extends (any_config_t) :: stable_config_t
     private
     type(flavor_t), dimension(:), allocatable :: flv
   contains
     procedure :: final => stable_config_final
     procedure :: write => stable_config_write
     procedure :: init => stable_config_init
  end type stable_config_t

  type, extends (any_t) :: stable_t
     private
     type(stable_config_t), pointer :: config => null ()
   contains
     procedure :: final => stable_final
     procedure :: write => stable_write
     procedure :: init => stable_init
  end type stable_t

  type, extends (any_config_t) :: unstable_config_t
     private
     type(flavor_t) :: flv
     real(default) :: integral = 0
     real(default) :: abs_error = 0
     real(default) :: rel_error = 0
     type(selector_t) :: selector
     type(decay_config_t), dimension(:), allocatable :: decay_config
   contains
     procedure :: final => unstable_config_final
     procedure :: write => unstable_config_write
     procedure :: init => unstable_config_init
     procedure :: init_decays => unstable_config_init_decays
     procedure :: connect_decay => unstable_config_connect_decay
     procedure :: compute => unstable_config_compute
     procedure :: init_test_case1
     procedure :: init_test_case2
  end type unstable_config_t

  type, extends (any_t) :: unstable_t
     private
     type(unstable_config_t), pointer :: config => null ()
     class(rng_t), allocatable :: rng
     integer :: selected_decay = 0
     type(decay_t), dimension(:), allocatable :: decay
   contains
     procedure :: final => unstable_final
     procedure :: write => unstable_write
     procedure :: write_process_instances => unstable_write_process_instances
     procedure :: init => unstable_init
     procedure :: link_interactions => unstable_link_interactions
     procedure :: import_rng => unstable_import_rng
     procedure :: select_chain => unstable_select_chain
     procedure :: generate => unstable_generate
  end type unstable_t

  type, extends (connected_state_t) :: decay_chain_entry_t
     private
     integer :: index = 0
     type(decay_config_t), pointer :: config => null ()
     integer :: selected_mci = 0
     integer :: selected_term = 0
     type(decay_chain_entry_t), pointer :: previous => null ()
  end type decay_chain_entry_t

  type :: decay_chain_t
     private
     type(process_instance_t), pointer :: process_instance => null ()
     integer :: selected_term = 0
     type(evaluator_t) :: correlated_trace
     type(decay_chain_entry_t), pointer :: last => null ()
   contains
     procedure :: final => decay_chain_final
     procedure :: write => decay_chain_write
     procedure :: build => decay_chain_build
     procedure :: build_term_entries => decay_chain_build_term_entries
     procedure :: build_decay_entries => decay_chain_build_decay_entries
     procedure :: evaluate => decay_chain_evaluate
     procedure :: get_probability => decay_chain_get_probability
  end type decay_chain_t

  type, extends (evt_t) :: evt_decay_t
     private
     type(decay_root_config_t) :: decay_root_config
     type(decay_root_t) :: decay_root
     type(decay_chain_t) :: decay_chain
     type(var_list_t), pointer :: var_list => null ()
   contains
     procedure :: write_name => evt_decay_write_name
     procedure :: write => evt_decay_write
     procedure :: set_var_list => evt_decay_set_var_list
     procedure :: connect => evt_decay_connect
     procedure :: prepare_new_event => evt_decay_prepare_new_event
     procedure :: generate_weighted => evt_decay_generate_weighted
     procedure :: make_particle_set => evt_decay_make_particle_set
  end type evt_decay_t


  interface
     subroutine any_config_final (object)
       import
       class(any_config_t), intent(inout) :: object
     end subroutine any_config_final
  end interface

  interface
     subroutine any_config_write (object, unit, indent, verbose)
       import
       class(any_config_t), intent(in) :: object
       integer, intent(in), optional :: unit, indent
       logical, intent(in), optional :: verbose
     end subroutine any_config_write
  end interface

  interface
     subroutine any_final (object)
       import
       class(any_t), intent(inout) :: object
     end subroutine any_final
  end interface

  interface
     subroutine any_write (object, unit, indent)
       import
       class(any_t), intent(in) :: object
       integer, intent(in), optional :: unit, indent
     end subroutine any_write
  end interface

  interface pacify
     module procedure pacify_decay
     module procedure pacify_decay_gen
     module procedure pacify_term
     module procedure pacify_unstable
  end interface pacify

  interface
    recursive module subroutine decay_term_config_final (object)
      class(decay_term_config_t), intent(inout) :: object
    end subroutine decay_term_config_final
    recursive module subroutine decay_term_config_write &
         (object, unit, indent, verbose)
      class(decay_term_config_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
      logical, intent(in), optional :: verbose
    end subroutine decay_term_config_write
    recursive module subroutine decay_term_config_compute (term)
      class(decay_term_config_t), intent(inout) :: term
    end subroutine decay_term_config_compute
    recursive module subroutine decay_term_final (object)
      class(decay_term_t), intent(inout) :: object
    end subroutine decay_term_final
    recursive module subroutine decay_term_write (object, unit, indent)
      class(decay_term_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
    end subroutine decay_term_write
    recursive module subroutine decay_term_write_process_instances &
         (term, unit, verbose)
      class(decay_term_t), intent(in) :: term
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine decay_term_write_process_instances
    module subroutine decay_term_make_rng (term, process)
      class(decay_term_t), intent(inout) :: term
      type(process_t), intent(inout) :: process
      class(rng_t), allocatable :: rng
    end subroutine decay_term_make_rng
    recursive module subroutine decay_term_link_interactions (term, trace)
      class(decay_term_t), intent(inout) :: term
      type(interaction_t), intent(in), target :: trace
    end subroutine decay_term_link_interactions
    recursive module subroutine decay_term_select_chain (term)
      class(decay_term_t), intent(inout) :: term
    end subroutine decay_term_select_chain
    recursive module subroutine decay_term_generate (term)
      class(decay_term_t), intent(inout) :: term
    end subroutine decay_term_generate
    recursive module subroutine decay_root_config_final (object)
      class(decay_root_config_t), intent(inout) :: object
    end subroutine decay_root_config_final
    recursive module subroutine decay_root_config_write &
         (object, unit, indent, verbose)
      class(decay_root_config_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
      logical, intent(in), optional :: verbose
    end subroutine decay_root_config_write
    module subroutine decay_root_config_write_header (object, unit, indent)
      class(decay_root_config_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
    end subroutine decay_root_config_write_header
    module recursive subroutine decay_root_config_write_terms &
         (object, unit, indent, verbose)
      class(decay_root_config_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
      logical, intent(in), optional :: verbose
    end subroutine decay_root_config_write_terms
    module subroutine decay_root_config_init (decay, model, process_id, n_terms)
      class(decay_root_config_t), intent(out) :: decay
      class(model_data_t), intent(in), target :: model
      type(string_t), intent(in) :: process_id
      integer, intent(in), optional :: n_terms
    end subroutine decay_root_config_init
    recursive module subroutine decay_root_config_init_term &
         (decay, i, flv, stable, model, process_stack, var_list)
      class(decay_root_config_t), intent(inout) :: decay
      integer, intent(in) :: i
      type(flavor_t), dimension(:,:), intent(in) :: flv
      logical, dimension(:), intent(in) :: stable
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
      type(var_list_t), intent(in), optional, target :: var_list
    end subroutine decay_root_config_init_term
    recursive module subroutine decay_root_config_connect &
         (decay, process, model, process_stack, process_instance, var_list)
      class(decay_root_config_t), intent(out) :: decay
      type(process_t), intent(in), target :: process
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
      type(process_instance_t), intent(in), optional, target :: process_instance
      type(var_list_t), intent(in), optional, target :: var_list
    end subroutine decay_root_config_connect
    recursive module subroutine decay_root_config_compute (decay)
      class(decay_root_config_t), intent(inout) :: decay
    end subroutine decay_root_config_compute
    recursive module subroutine decay_gen_final (object)
      class(decay_gen_t), intent(inout) :: object
    end subroutine decay_gen_final
    module subroutine decay_root_final (object)
      class(decay_root_t), intent(inout) :: object
    end subroutine decay_root_final
    module subroutine decay_root_write (object, unit)
      class(decay_root_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine decay_root_write
    recursive module subroutine decay_gen_write_process_instances &
         (decay, unit, verbose)
      class(decay_gen_t), intent(in) :: decay
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine decay_gen_write_process_instances
    recursive module subroutine decay_gen_init (decay, term_config)
      class(decay_gen_t), intent(out) :: decay
      type(decay_term_config_t), dimension(:), intent(in), target :: term_config
    end subroutine decay_gen_init
    module subroutine decay_root_init (decay_root, config, process_instance)
      class(decay_root_t), intent(out) :: decay_root
      type(decay_root_config_t), intent(in), target :: config
      type(process_instance_t), intent(in), target :: process_instance
    end subroutine decay_root_init
    module subroutine decay_gen_set_mci (decay, i)
      class(decay_gen_t), intent(inout) :: decay
      integer, intent(in) :: i
    end subroutine decay_gen_set_mci
    module subroutine decay_gen_set_term (decay, i)
      class(decay_gen_t), intent(inout) :: decay
      integer, intent(in) :: i
    end subroutine decay_gen_set_term
    module function decay_gen_get_mci (decay) result (i)
      class(decay_gen_t), intent(inout) :: decay
      integer :: i
    end function decay_gen_get_mci
    module function decay_gen_get_term (decay) result (i)
      class(decay_gen_t), intent(inout) :: decay
      integer :: i
    end function decay_gen_get_term
    module subroutine decay_gen_make_term_rng (decay, process)
      class(decay_gen_t), intent(inout) :: decay
      type(process_t), intent(in), pointer :: process
    end subroutine decay_gen_make_term_rng
    recursive module subroutine decay_gen_link_term_interactions (decay)
      class(decay_gen_t), intent(inout) :: decay
    end subroutine decay_gen_link_term_interactions
    module subroutine decay_root_select_chain (decay_root)
      class(decay_root_t), intent(inout) :: decay_root
    end subroutine decay_root_select_chain
    module subroutine decay_root_generate (decay_root)
      class(decay_root_t), intent(inout) :: decay_root
    end subroutine decay_root_generate
    recursive module subroutine decay_config_write &
         (object, unit, indent, verbose)
      class(decay_config_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
      logical, intent(in), optional :: verbose
    end subroutine decay_config_write
    recursive module subroutine decay_config_connect &
         (decay, process, model, process_stack, process_instance, var_list)
      class(decay_config_t), intent(out) :: decay
      type(process_t), intent(in), target :: process
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
      type(process_instance_t), intent(in), optional, target :: process_instance
      type(var_list_t), intent(in), optional, target :: var_list
    end subroutine decay_config_connect
    module subroutine decay_config_set_flv (decay, flv)
      class(decay_config_t), intent(inout) :: decay
      type(flavor_t), intent(in) :: flv
    end subroutine decay_config_set_flv
    recursive module subroutine decay_config_compute (decay)
      class(decay_config_t), intent(inout) :: decay
    end subroutine decay_config_compute
    recursive module subroutine decay_final (object)
      class(decay_t), intent(inout) :: object
    end subroutine decay_final
    recursive module subroutine decay_write (object, unit, indent, recursive)
      class(decay_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent, recursive
    end subroutine decay_write
    recursive module subroutine decay_init (decay, config)
      class(decay_t), intent(out) :: decay
      type(decay_config_t), intent(in), target :: config
    end subroutine decay_init
    recursive module subroutine decay_link_interactions (decay, i_prt, trace)
      class(decay_t), intent(inout) :: decay
      integer, intent(in) :: i_prt
      type(interaction_t), intent(in), target :: trace
    end subroutine decay_link_interactions
    recursive module subroutine decay_select_chain (decay)
      class(decay_t), intent(inout) :: decay
    end subroutine decay_select_chain
    recursive module subroutine decay_generate (decay)
      class(decay_t), intent(inout) :: decay
    end subroutine decay_generate
    module subroutine stable_config_final (object)
      class(stable_config_t), intent(inout) :: object
    end subroutine stable_config_final
    recursive module subroutine stable_config_write &
         (object, unit, indent, verbose)
      class(stable_config_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
      logical, intent(in), optional :: verbose
    end subroutine stable_config_write
    module subroutine stable_config_init (config, flv)
      class(stable_config_t), intent(out) :: config
      type(flavor_t), dimension(:), intent(in) :: flv
    end subroutine stable_config_init
    module subroutine stable_final (object)
      class(stable_t), intent(inout) :: object
    end subroutine stable_final
    module subroutine stable_write (object, unit, indent)
      class(stable_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
    end subroutine stable_write
    module subroutine stable_init (stable, config)
      class(stable_t), intent(out) :: stable
      type(stable_config_t), intent(in), target :: config
    end subroutine stable_init
    recursive module subroutine unstable_config_final (object)
      class(unstable_config_t), intent(inout) :: object
    end subroutine unstable_config_final
    recursive module subroutine unstable_config_write &
         (object, unit, indent, verbose)
      class(unstable_config_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
      logical, intent(in), optional :: verbose
    end subroutine unstable_config_write
    module subroutine unstable_config_init (unstable, flv, set_decays, model)
      class(unstable_config_t), intent(out) :: unstable
      type(flavor_t), intent(in) :: flv
      logical, intent(in), optional :: set_decays
      class(model_data_t), intent(in), optional, target :: model
    end subroutine unstable_config_init
    recursive module subroutine unstable_config_init_decays &
         (unstable, decay_id, model, process_stack, var_list)
      class(unstable_config_t), intent(inout) :: unstable
      type(string_t), dimension(:), intent(in) :: decay_id
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
      type(var_list_t), intent(in), optional :: var_list
    end subroutine unstable_config_init_decays
    module subroutine unstable_config_connect_decay &
         (unstable, i, process, model)
      class(unstable_config_t), intent(inout) :: unstable
      integer, intent(in) :: i
      type(process_t), intent(in), target :: process
      class(model_data_t), intent(in), target :: model
    end subroutine unstable_config_connect_decay
    recursive module subroutine unstable_config_compute (unstable)
      class(unstable_config_t), intent(inout) :: unstable
    end subroutine unstable_config_compute
    recursive module subroutine unstable_final (object)
      class(unstable_t), intent(inout) :: object
    end subroutine unstable_final
    recursive module subroutine unstable_write (object, unit, indent)
      class(unstable_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
    end subroutine unstable_write
    recursive module subroutine unstable_write_process_instances &
         (unstable, unit, verbose)
      class(unstable_t), intent(in) :: unstable
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine unstable_write_process_instances
    recursive module subroutine unstable_init (unstable, config)
      class(unstable_t), intent(out) :: unstable
      type(unstable_config_t), intent(in), target :: config
    end subroutine unstable_init
    recursive module subroutine unstable_link_interactions &
         (unstable, i_prt, trace)
      class(unstable_t), intent(inout) :: unstable
      integer, intent(in) :: i_prt
      type(interaction_t), intent(in), target :: trace
    end subroutine unstable_link_interactions
    module subroutine unstable_import_rng (unstable, rng)
      class(unstable_t), intent(inout) :: unstable
      class(rng_t), intent(inout), allocatable :: rng
    end subroutine unstable_import_rng
    recursive module subroutine unstable_select_chain (unstable)
      class(unstable_t), intent(inout) :: unstable
    end subroutine unstable_select_chain
    recursive module subroutine unstable_generate (unstable)
      class(unstable_t), intent(inout) :: unstable
    end subroutine unstable_generate
    module subroutine decay_chain_final (object)
      class(decay_chain_t), intent(inout) :: object
    end subroutine decay_chain_final
    module subroutine decay_chain_write (object, unit)
      class(decay_chain_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine decay_chain_write
    module subroutine decay_chain_build (chain, decay_root)
      class(decay_chain_t), intent(inout), target :: chain
      type(decay_root_t), intent(in) :: decay_root
    end subroutine decay_chain_build
    recursive module subroutine decay_chain_build_term_entries (chain, term)
      class(decay_chain_t), intent(inout) :: chain
      type(decay_term_t), intent(in) :: term
    end subroutine decay_chain_build_term_entries
    recursive module subroutine decay_chain_build_decay_entries (chain, decay)
      class(decay_chain_t), intent(inout) :: chain
      type(decay_t), intent(in) :: decay
    end subroutine decay_chain_build_decay_entries
    module subroutine decay_chain_evaluate (chain)
      class(decay_chain_t), intent(inout) :: chain
    end subroutine decay_chain_evaluate
    module function decay_chain_get_probability (chain) result (x)
      class(decay_chain_t), intent(in) :: chain
      real(default) :: x
    end function decay_chain_get_probability
    module subroutine evt_decay_write_name (evt, unit)
      class(evt_decay_t), intent(in) :: evt
      integer, intent(in), optional :: unit
    end subroutine evt_decay_write_name
    module subroutine evt_decay_write &
         (evt, unit, verbose, more_verbose, testflag)
      class(evt_decay_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, more_verbose, testflag
    end subroutine evt_decay_write
    module subroutine evt_decay_set_var_list (evt, var_list)
      class(evt_decay_t), intent(inout) :: evt
      type(var_list_t), intent(in), target :: var_list
    end subroutine evt_decay_set_var_list
    module subroutine evt_decay_connect &
         (evt, process_instance, model, process_stack)
      class(evt_decay_t), intent(inout), target :: evt
      type(process_instance_t), intent(in), target :: process_instance
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
    end subroutine evt_decay_connect
    module subroutine evt_decay_prepare_new_event (evt, i_mci, i_term)
      class(evt_decay_t), intent(inout) :: evt
      integer, intent(in) :: i_mci, i_term
    end subroutine evt_decay_prepare_new_event
    module subroutine evt_decay_generate_weighted (evt, probability)
      class(evt_decay_t), intent(inout) :: evt
      real(default), intent(inout) :: probability
    end subroutine evt_decay_generate_weighted
    module subroutine evt_decay_make_particle_set &
         (evt, factorization_mode, keep_correlations, r)
      class(evt_decay_t), intent(inout) :: evt
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
    end subroutine evt_decay_make_particle_set
    module subroutine pacify_decay (evt)
      class(evt_decay_t), intent(inout) :: evt
    end subroutine pacify_decay
    recursive module subroutine pacify_decay_gen (decay)
      class(decay_gen_t), intent(inout) :: decay
    end subroutine pacify_decay_gen
    recursive module subroutine pacify_term (term)
      class(decay_term_t), intent(inout) :: term
    end subroutine pacify_term
    recursive module subroutine pacify_unstable (unstable)
      class(unstable_t), intent(inout) :: unstable
    end subroutine pacify_unstable
    module subroutine init_test_case1 &
         (unstable, i, flv, integral, relerr, model)
      class(unstable_config_t), intent(inout) :: unstable
      integer, intent(in) :: i
      type(flavor_t), dimension(:,:), intent(in) :: flv
      real(default), intent(in) :: integral
      real(default), intent(in) :: relerr
      class(model_data_t), intent(in), target :: model
    end subroutine init_test_case1
    module subroutine init_test_case2 (unstable, flv1, flv21, flv22, model)
      class(unstable_config_t), intent(inout) :: unstable
      type(flavor_t), dimension(:,:), intent(in) :: flv1, flv21, flv22
      class(model_data_t), intent(in), target :: model
    end subroutine init_test_case2
  end interface

contains

  recursive subroutine decay_term_config_init &
       (term, flv, stable, model, process_stack, var_list)
    class(decay_term_config_t), intent(out) :: term
    type(flavor_t), dimension(:,:), intent(in) :: flv
    logical, dimension(:), intent(in) :: stable
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    type(var_list_t), intent(in), optional :: var_list
    type(string_t), dimension(:), allocatable :: decay
    integer :: i
    allocate (term%prt (size (flv, 1)))
    do i = 1, size (flv, 1)
       associate (prt => term%prt(i))
         if (stable(i)) then
            allocate (stable_config_t :: prt%c)
         else
            allocate (unstable_config_t :: prt%c)
         end if
         select type (prt_config => prt%c)
         type is (stable_config_t)
            call prt_config%init (flv(i,:))
         type is (unstable_config_t)
            if (all (flv(i,:) == flv(i,1))) then
               call prt_config%init (flv(i,1))
               call flv(i,1)%get_decays (decay)
               call prt_config%init_decays &
                    (decay, model, process_stack, var_list)
            else
               call prt_config%write ()
               call msg_fatal ("Decay configuration: &
                    &unstable product must be unique")
            end if
         end select
       end associate
    end do
  end subroutine decay_term_config_init

  recursive subroutine decay_term_init (term, config)
    class(decay_term_t), intent(out) :: term
    type(decay_term_config_t), intent(in), target :: config
    integer :: i
    term%config => config
    allocate (term%particle_out (size (config%prt)))
    do i = 1, size (config%prt)
       select type (prt_config => config%prt(i)%c)
       type is (stable_config_t)
          allocate (stable_t :: term%particle_out(i)%c)
          select type (stable => term%particle_out(i)%c)
          type is (stable_t)
             call stable%init (prt_config)
          end select
       type is (unstable_config_t)
          allocate (unstable_t :: term%particle_out(i)%c)
          select type (unstable => term%particle_out(i)%c)
          type is (unstable_t)
             call unstable%init (prt_config)
          end select
       end select
    end do
  end subroutine decay_term_init


end module decays

