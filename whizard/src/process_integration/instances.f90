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

module instances

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use mci_base
  use particles
  use sm_qcd, only: qcd_t
  use quantum_numbers
  use interactions
  use model_data
  use variables
  use sf_base
  use pdf, only: pdf_data_t
  use physics_defs
  use process_constants
  use state_matrices
  use phs_base
  use prc_core, only: prc_core_t, prc_core_state_t
  !!! local modules
  use parton_states
  use process_counter
  use pcm_base
  use pcm
  use process_config
  use process_mci
  use process
  use kinematics

  implicit none
  private

  public :: setup_interaction_qn_index
  public :: setup_interaction_qn_hel
  public :: process_instance_t
  public :: process_instance_ptr_t
  public :: process_instance_hook_t
  public :: process_instance_hook_final, process_instance_hook_evaluate
  public :: pacify

  type :: term_instance_t
     type(process_term_t), pointer :: config => null ()
     class(pcm_t), pointer :: pcm => null ()
     class(pcm_workspace_t), pointer :: pcm_work => null ()
     logical :: active = .false.
     complex(default), dimension(:), allocatable :: amp
     type(interaction_t) :: int_hard
     type(isolated_state_t) :: isolated
     type(connected_state_t) :: connected
     class(prc_core_state_t), allocatable :: core_state
     logical :: checked = .false.
     logical :: passed = .false.
     logical, dimension(:), allocatable :: passed_array
     integer, dimension(:), allocatable :: i_flv_to_i_flv_rep
     real(default) :: scale = 0
     real(default), allocatable :: fac_scale
     real(default), allocatable :: ren_scale
     real(default), allocatable :: es_scale
     real(default), allocatable :: alpha_qcd_forced
     real(default) :: weight = 1
     type(vector4_t), dimension(:), allocatable :: p_seed
     type(vector4_t), dimension(:), allocatable :: p_hard
     integer :: nlo_type = BORN
     integer, dimension(:), allocatable :: same_kinematics
     logical :: negative_sf = .false.
     logical :: flv_dep_cut_eval = .false.
   contains
     procedure :: write => term_instance_write
     procedure :: final => term_instance_final
     procedure :: configure => term_instance_configure
     procedure :: init => term_instance_init
     procedure :: setup_dynamics => term_instance_setup_dynamics
     procedure :: init_eqv_expr_classes => term_instance_init_eqv_expr_classes
     procedure :: init_interaction_qn_index => &
          term_instance_init_interaction_qn_index
     procedure :: setup_fks_kinematics => term_instance_setup_fks_kinematics
     procedure :: compute_seed_kinematics => term_instance_compute_seed_kinematics
     procedure :: evaluate_projections => term_instance_evaluate_projections
     procedure :: compute_hard_kinematics => &
          term_instance_compute_hard_kinematics
     procedure :: recover_seed_kinematics => &
          term_instance_recover_seed_kinematics
     procedure :: apply_real_partition => term_instance_apply_real_partition
     procedure :: get_p_hard => term_instance_get_p_hard
     procedure :: set_emitter => term_instance_set_emitter
     procedure :: setup_expressions => term_instance_setup_expressions
     procedure :: setup_event_data => term_instance_setup_event_data
     procedure :: evaluate_color_correlations => &
        term_instance_evaluate_color_correlations
     procedure :: evaluate_charge_correlations => &
        term_instance_evaluate_charge_correlations
     procedure :: evaluate_spin_correlations => &
          term_instance_evaluate_spin_correlations
     procedure :: apply_fks => term_instance_apply_fks
     procedure :: evaluate_sqme_virt => term_instance_evaluate_sqme_virt
     procedure :: evaluate_sqme_mismatch => term_instance_evaluate_sqme_mismatch
     procedure :: evaluate_sqme_dglap => term_instance_evaluate_sqme_dglap
     procedure :: reset => term_instance_reset
     procedure :: set_alpha_qcd_forced => term_instance_set_alpha_qcd_forced
     procedure :: compute_eff_kinematics => &
          term_instance_compute_eff_kinematics
     procedure :: recover_hard_kinematics => &
          term_instance_recover_hard_kinematics
     procedure :: evaluate_expressions => &
          term_instance_evaluate_expressions
     procedure :: evaluate_interaction => term_instance_evaluate_interaction
     procedure :: evaluate_interaction_default &
        => term_instance_evaluate_interaction_default
     procedure :: evaluate_interaction_external &
        => term_instance_evaluate_interaction_external
     procedure :: evaluate_interaction_external_tree &
        => term_instance_evaluate_interaction_external_tree
     procedure :: evaluate_interaction_external_loop &
        => term_instance_evaluate_interaction_external_loop
     procedure :: evaluate_trace => term_instance_evaluate_trace
     procedure :: evaluate_scaled_sf_chains => &
          term_instance_evaluate_scaled_sf_chains
     procedure :: evaluate_event_data => term_instance_evaluate_event_data
     procedure :: set_fac_scale => term_instance_set_fac_scale
     procedure :: get_fac_scale => term_instance_get_fac_scale
     procedure :: get_ren_scale => term_instance_get_ren_scale
     procedure :: get_alpha_s => term_instance_get_alpha_s
     procedure :: get_helicities_for_openloops => &
          term_instance_get_helicities_for_openloops
     procedure :: get_i_term_global => term_instance_get_i_term_global
     procedure :: is_subtraction => term_instance_is_subtraction
     procedure :: get_n_sub => term_instance_get_n_sub
     procedure :: get_n_sub_color => term_instance_get_n_sub_color
     procedure :: get_n_sub_spin => term_instance_get_n_sub_spin
     procedure :: set_born_sqmes => term_instance_set_born_sqmes
     procedure :: set_sf_factors => term_instance_set_sf_factors
  end type term_instance_t

  type, extends (mci_sampler_t) :: process_instance_t
     type(process_t), pointer :: process => null ()
     class(pcm_t), pointer :: pcm => null ()
     class(pcm_workspace_t), allocatable :: pcm_work
     integer :: evaluation_status = STAT_UNDEFINED
     real(default) :: sqme = 0
     real(default) :: weight = 0
     real(default) :: excess = 0
     integer :: n_dropped = 0
     integer :: i_mci = 0
     integer :: selected_channel = 0
     type(sf_chain_t) :: sf_chain
     type(kinematics_t), dimension(:), allocatable :: kin
     type(term_instance_t), dimension(:), allocatable :: term
     type(mci_work_t), dimension(:), allocatable :: mci_work
     class(process_instance_hook_t), pointer :: hook => null ()
   contains
     procedure :: write_header => process_instance_write_header
     procedure :: write => process_instance_write
     procedure :: init => process_instance_init
     procedure :: final => process_instance_final
     procedure :: reset => process_instance_reset
     procedure :: sampler_test => process_instance_sampler_test
     procedure :: generate_weighted_event => &
          process_instance_generate_weighted_event
     procedure :: generate_unweighted_event => &
          process_instance_generate_unweighted_event
     procedure :: recover_event => process_instance_recover_event
     procedure :: activate => process_instance_activate
     procedure :: find_same_kinematics => process_instance_find_same_kinematics
     procedure :: transfer_same_kinematics => &
          process_instance_transfer_same_kinematics
     procedure :: redo_sf_chains => process_instance_redo_sf_chains
     procedure :: integrate => process_instance_integrate
     procedure :: setup_sf_chain => process_instance_setup_sf_chain
     procedure :: setup_event_data => process_instance_setup_event_data
     procedure :: choose_mci => process_instance_choose_mci
     procedure :: set_mcpar => process_instance_set_mcpar
     procedure :: receive_beam_momenta => process_instance_receive_beam_momenta
     procedure :: set_beam_momenta => process_instance_set_beam_momenta
     procedure :: recover_beam_momenta => process_instance_recover_beam_momenta
     procedure :: select_channel => process_instance_select_channel
     procedure :: compute_seed_kinematics => &
          process_instance_compute_seed_kinematics
     procedure :: get_x_process => process_instance_get_x_process
     procedure :: get_active_component_type => &
          process_instance_get_active_component_type
     procedure :: recover_mcpar => process_instance_recover_mcpar
     procedure :: recover_sfchain => process_instance_recover_sfchain
     procedure :: compute_hard_kinematics => &
          process_instance_compute_hard_kinematics
     procedure :: recover_seed_kinematics => &
          process_instance_recover_seed_kinematics
     procedure :: compute_eff_kinematics => &
          process_instance_compute_eff_kinematics
     procedure :: recover_hard_kinematics => &
          process_instance_recover_hard_kinematics
     procedure :: evaluate_expressions => &
          process_instance_evaluate_expressions
     procedure :: compute_other_channels => &
          process_instance_compute_other_channels
     procedure :: reset_core_kinematics => process_instance_reset_core_kinematics
     procedure :: evaluate_trace => process_instance_evaluate_trace
     procedure :: apply_real_partition => process_instance_apply_real_partition
     procedure :: set_i_mci_to_real_component => &
          process_instance_set_i_mci_to_real_component
     procedure :: evaluate_event_data => process_instance_evaluate_event_data
     procedure :: compute_sqme_rad => process_instance_compute_sqme_rad
     procedure :: normalize_weight => process_instance_normalize_weight
     procedure :: evaluate_sqme => process_instance_evaluate_sqme
     procedure :: recover => process_instance_recover
     procedure :: evaluate => process_instance_evaluate
     procedure :: is_valid => process_instance_is_valid
     procedure :: append_after_hook => process_instance_append_after_hook
     procedure :: evaluate_after_hook => process_instance_evaluate_after_hook
     procedure :: rebuild => process_instance_rebuild
     procedure :: fetch => process_instance_fetch
     procedure :: init_simulation => process_instance_init_simulation
     procedure :: final_simulation => process_instance_final_simulation
     procedure :: get_mcpar => process_instance_get_mcpar
     procedure :: has_evaluated_trace => process_instance_has_evaluated_trace
     procedure :: is_complete_event => process_instance_is_complete_event
     procedure :: select_i_term => process_instance_select_i_term
     procedure :: get_beam_int_ptr => process_instance_get_beam_int_ptr
     procedure :: get_trace_int_ptr => process_instance_get_trace_int_ptr
     procedure :: get_matrix_int_ptr => process_instance_get_matrix_int_ptr
     procedure :: get_flows_int_ptr => process_instance_get_flows_int_ptr
     procedure :: get_state_flv => process_instance_get_state_flv
     procedure :: get_isolated_state_ptr => &
          process_instance_get_isolated_state_ptr
     procedure :: get_connected_state_ptr => &
          process_instance_get_connected_state_ptr
     procedure :: get_beam_index => process_instance_get_beam_index
     procedure :: get_in_index => process_instance_get_in_index
     procedure :: get_sqme => process_instance_get_sqme
     procedure :: get_weight => process_instance_get_weight
     procedure :: get_excess => process_instance_get_excess
     procedure :: get_n_dropped => process_instance_get_n_dropped
     procedure :: get_channel => process_instance_get_channel
     procedure :: set_fac_scale => process_instance_set_fac_scale
     procedure :: get_fac_scale => process_instance_get_fac_scale
     procedure :: get_alpha_s => process_instance_get_alpha_s
     procedure :: get_qcd => process_instance_get_qcd
     procedure :: reset_counter => process_instance_reset_counter
     procedure :: record_call => process_instance_record_call
     procedure :: get_counter => process_instance_get_counter
     procedure :: get_actual_calls_total => process_instance_get_actual_calls_total
     procedure :: reset_matrix_elements => process_instance_reset_matrix_elements
     procedure :: get_test_phase_space_point &
        => process_instance_get_test_phase_space_point
     procedure :: get_p_hard => process_instance_get_p_hard
     procedure :: get_first_active_i_term => &
          process_instance_get_first_active_i_term
     procedure :: get_real_of_mci => process_instance_get_real_of_mci
     procedure :: get_connected_states => process_instance_get_connected_states
     procedure :: get_sqrts => process_instance_get_sqrts
     procedure :: get_polarization => process_instance_get_polarization
     procedure :: get_beam_file => process_instance_get_beam_file
     procedure :: get_process_name => process_instance_get_process_name
     procedure :: get_trace => process_instance_get_trace
     procedure :: set_trace => process_instance_set_trace
     procedure :: set_alpha_qcd_forced => process_instance_set_alpha_qcd_forced
     procedure :: has_nlo_component => process_instance_has_nlo_component
     procedure :: keep_failed_events => process_instance_keep_failed_events
     procedure :: get_term_indices => process_instance_get_term_indices
     procedure :: get_boost_to_lab => process_instance_get_boost_to_lab
     procedure :: get_boost_to_cms => process_instance_get_boost_to_cms
     procedure :: lab_is_cm => process_instance_lab_is_cm
  end type process_instance_t

  type :: process_instance_ptr_t
     type(process_instance_t), pointer :: p => null ()
  end type process_instance_ptr_t

  type, abstract :: process_instance_hook_t
     class(process_instance_hook_t), pointer :: next => null ()
   contains
     procedure(process_instance_hook_init), deferred :: init
     procedure(process_instance_hook_final), deferred :: final
     procedure(process_instance_hook_evaluate), deferred :: evaluate
  end type process_instance_hook_t


  abstract interface
     subroutine process_instance_hook_init (hook, var_list, instance, pdf_data)
       import :: process_instance_hook_t, var_list_t, process_instance_t, pdf_data_t
       class(process_instance_hook_t), intent(inout), target :: hook
       type(var_list_t), intent(in) :: var_list
       class(process_instance_t), intent(in), target :: instance
       type(pdf_data_t), intent(in), optional :: pdf_data
     end subroutine process_instance_hook_init

     subroutine process_instance_hook_final (hook)
       import :: process_instance_hook_t
       class(process_instance_hook_t), intent(inout) :: hook
     end subroutine process_instance_hook_final

     subroutine process_instance_hook_evaluate (hook, instance)
       import :: process_instance_hook_t, process_instance_t
       class(process_instance_hook_t), intent(inout) :: hook
       class(process_instance_t), intent(in), target :: instance
     end subroutine process_instance_hook_evaluate
  end interface

  interface pacify
     module procedure pacify_process_instance
  end interface pacify


  interface
    module subroutine term_instance_write &
         (term, unit, kin, show_eff_state, testflag)
      class(term_instance_t), intent(in) :: term
      integer, intent(in), optional :: unit
      type(kinematics_t), intent(in), optional :: kin
      logical, intent(in), optional :: show_eff_state
      logical, intent(in), optional :: testflag
    end subroutine term_instance_write
    module subroutine term_instance_final (term)
      class(term_instance_t), intent(inout) :: term
    end subroutine term_instance_final
    module subroutine term_instance_configure &
         (term_instance, process, i, pcm_work, sf_chain, kin)
      class(term_instance_t), intent(out), target :: term_instance
      type(process_t), intent(in), target :: process
      integer, intent(in) :: i
      class(pcm_workspace_t), intent(in), target :: pcm_work
      type(sf_chain_t), intent(in), target :: sf_chain
      type(kinematics_t), intent(inout), target :: kin
    end subroutine term_instance_configure
    module subroutine term_instance_init &
         (term_instance, pcm, pcm_work, nlo_type)
      class(term_instance_t), intent(out) :: term_instance
      class(pcm_t), intent(in), target :: pcm
      class(pcm_workspace_t), intent(in), target :: pcm_work
      integer, intent(in) :: nlo_type
    end subroutine term_instance_init
    module subroutine term_instance_setup_dynamics &
         (term, process, i_term, kin, real_finite)
      class(term_instance_t), intent(inout), target :: term
      type(process_t), intent(in), target:: process
      integer, intent(in) :: i_term
      type(kinematics_t), intent(in) :: kin
      logical, intent(in), optional :: real_finite
    end subroutine term_instance_setup_dynamics
    module subroutine setup_interaction_qn_index &
         (int, data, qn_config, n_sub, is_polarized)
      class(interaction_t), intent(inout) :: int
      class(process_constants_t), intent(in) :: data
      type(quantum_numbers_t), dimension(:, :), intent(in) :: qn_config
      integer, intent(in) :: n_sub
      logical, intent(in) :: is_polarized
    end subroutine setup_interaction_qn_index
    module subroutine setup_interaction_qn_hel (int, data, qn_hel)
      class(interaction_t), intent(in) :: int
      class(process_constants_t), intent(in) :: data
      type(quantum_numbers_t), dimension(:, :), allocatable, intent(out) :: &
           qn_hel
    end subroutine setup_interaction_qn_hel
    module subroutine term_instance_init_eqv_expr_classes (term)
      class(term_instance_t), intent(inout), target :: term
    end subroutine term_instance_init_eqv_expr_classes
    module subroutine term_instance_init_interaction_qn_index (term, core, &
         int, n_sub, model, is_polarized)
      class(term_instance_t), intent(inout), target :: term
      class(prc_core_t), intent(in) :: core
      class(interaction_t), intent(inout) :: int
      integer, intent(in) :: n_sub
      class(model_data_t), intent(in) :: model
      logical, intent(in), optional :: is_polarized
    end subroutine term_instance_init_interaction_qn_index
    module subroutine term_instance_setup_fks_kinematics &
         (term, kin, var_list, beam_config)
      class(term_instance_t), intent(inout), target :: term
      type(kinematics_t), intent(inout) :: kin
      type(var_list_t), intent(in) :: var_list
      type(process_beam_config_t), intent(in) :: beam_config
    end subroutine term_instance_setup_fks_kinematics
    module subroutine term_instance_compute_seed_kinematics &
         (term, kin, mci_work, phs_channel, success)
      class(term_instance_t), intent(inout), target :: term
      type(kinematics_t), intent(inout) :: kin
      type(mci_work_t), intent(in) :: mci_work
      integer, intent(in) :: phs_channel
      logical, intent(out) :: success
    end subroutine term_instance_compute_seed_kinematics
    module subroutine term_instance_evaluate_projections (term, kin)
      class(term_instance_t), intent(inout) :: term
      type(kinematics_t), intent(inout) :: kin
    end subroutine term_instance_evaluate_projections
    module subroutine term_instance_compute_hard_kinematics &
         (term, kin, recover, skip_term, success)
      class(term_instance_t), intent(inout) :: term
      type(kinematics_t), intent(inout) :: kin
      integer, intent(in), optional :: skip_term
      logical, intent(in), optional :: recover
      logical, intent(out) :: success
    end subroutine term_instance_compute_hard_kinematics
    module subroutine term_instance_recover_seed_kinematics &
         (term, kin, p_seed_ref)
      class(term_instance_t), intent(inout) :: term
      type(kinematics_t), intent(in) :: kin
      type(vector4_t), dimension(:), intent(in), optional :: p_seed_ref
    end subroutine term_instance_recover_seed_kinematics
    module subroutine term_instance_apply_real_partition (term, kin)
      class(term_instance_t), intent(inout) :: term
      type(kinematics_t), intent(in) :: kin
    end subroutine term_instance_apply_real_partition
    pure module function term_instance_get_p_hard &
         (term_instance) result (p_hard)
      type(vector4_t), dimension(:), allocatable :: p_hard
      class(term_instance_t), intent(in) :: term_instance
    end function term_instance_get_p_hard
    module subroutine term_instance_set_emitter (term, kin)
      class(term_instance_t), intent(inout) :: term
      type(kinematics_t), intent(inout) :: kin
    end subroutine term_instance_set_emitter
    module subroutine term_instance_setup_expressions (term, meta, config)
      class(term_instance_t), intent(inout), target :: term
      type(process_metadata_t), intent(in), target :: meta
      type(process_config_data_t), intent(in) :: config
    end subroutine term_instance_setup_expressions
    module subroutine term_instance_setup_event_data (term, kin, core, model)
      class(term_instance_t), intent(inout), target :: term
      type(kinematics_t), intent(in) :: kin
      class(prc_core_t), intent(in) :: core
      class(model_data_t), intent(in), target :: model
    end subroutine term_instance_setup_event_data
    module subroutine term_instance_evaluate_color_correlations (term, core)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(inout) :: core
    end subroutine term_instance_evaluate_color_correlations
    module subroutine term_instance_evaluate_charge_correlations (term, core)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(inout) :: core
    end subroutine term_instance_evaluate_charge_correlations
    module subroutine term_instance_evaluate_spin_correlations (term, core)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(inout) :: core
    end subroutine term_instance_evaluate_spin_correlations
    module subroutine term_instance_apply_fks &
         (term, kin, alpha_s_sub, alpha_qed_sub)
      class(term_instance_t), intent(inout) :: term
      class(kinematics_t), intent(inout) :: kin
      real(default), intent(in) :: alpha_s_sub, alpha_qed_sub
    end subroutine term_instance_apply_fks
    module subroutine term_instance_evaluate_sqme_virt &
         (term, alpha_s, alpha_qed)
      class(term_instance_t), intent(inout) :: term
      real(default), intent(in) :: alpha_s, alpha_qed
    end subroutine term_instance_evaluate_sqme_virt
    module subroutine term_instance_evaluate_sqme_mismatch (term, alpha_s)
      class(term_instance_t), intent(inout) :: term
      real(default), intent(in) :: alpha_s
    end subroutine term_instance_evaluate_sqme_mismatch
    module subroutine term_instance_evaluate_sqme_dglap &
         (term, alpha_s, alpha_qed)
      class(term_instance_t), intent(inout) :: term
      real(default), intent(in) :: alpha_s, alpha_qed
    end subroutine term_instance_evaluate_sqme_dglap
    module subroutine term_instance_reset (term)
      class(term_instance_t), intent(inout) :: term
    end subroutine term_instance_reset
    module subroutine term_instance_set_alpha_qcd_forced (term, alpha_qcd)
      class(term_instance_t), intent(inout) :: term
      real(default), intent(in) :: alpha_qcd
    end subroutine term_instance_set_alpha_qcd_forced
    module subroutine term_instance_compute_eff_kinematics (term)
      class(term_instance_t), intent(inout) :: term
    end subroutine term_instance_compute_eff_kinematics
    module subroutine term_instance_recover_hard_kinematics (term)
      class(term_instance_t), intent(inout) :: term
    end subroutine term_instance_recover_hard_kinematics
    module subroutine term_instance_evaluate_expressions &
         (term, config, scale_forced)
      class(term_instance_t), intent(inout) :: term
      type(process_beam_config_t), intent(in) :: config
      real(default), intent(in), allocatable, optional :: scale_forced
    end subroutine term_instance_evaluate_expressions
    module subroutine term_instance_evaluate_interaction (term, core, kin)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(in), pointer :: core
      type(kinematics_t), intent(inout) :: kin
    end subroutine term_instance_evaluate_interaction
    module subroutine term_instance_evaluate_interaction_default (term, core)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(in) :: core
    end subroutine term_instance_evaluate_interaction_default
    module subroutine term_instance_evaluate_interaction_external &
         (term, core, kin)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(inout) :: core
      type(kinematics_t), intent(inout) :: kin
    end subroutine term_instance_evaluate_interaction_external
    module subroutine term_instance_evaluate_interaction_external_tree &
         (term, core)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(inout) :: core
    end subroutine term_instance_evaluate_interaction_external_tree
    module subroutine term_instance_evaluate_interaction_external_loop &
         (term, core)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(in) :: core
    end subroutine term_instance_evaluate_interaction_external_loop
    module subroutine term_instance_evaluate_trace (term, kin)
      class(term_instance_t), intent(inout) :: term
      type(kinematics_t), intent(inout) :: kin
    end subroutine term_instance_evaluate_trace
    module subroutine term_instance_evaluate_event_data (term)
      class(term_instance_t), intent(inout) :: term
    end subroutine term_instance_evaluate_event_data
    module subroutine term_instance_set_fac_scale (term, fac_scale)
      class(term_instance_t), intent(inout) :: term
      real(default), intent(in) :: fac_scale
    end subroutine term_instance_set_fac_scale
    module function term_instance_get_fac_scale (term) result (fac_scale)
      class(term_instance_t), intent(in) :: term
      real(default) :: fac_scale
    end function term_instance_get_fac_scale
    module function term_instance_get_ren_scale (term) result (ren_scale)
      class(term_instance_t), intent(in) :: term
      real(default) :: ren_scale
    end function term_instance_get_ren_scale
    module function term_instance_get_alpha_s (term, core) result (alpha_s)
      class(term_instance_t), intent(in) :: term
      class(prc_core_t), intent(in) :: core
      real(default) :: alpha_s
    end function term_instance_get_alpha_s
    module subroutine term_instance_get_helicities_for_openloops &
         (term, helicities)
      class(term_instance_t), intent(in) :: term
      integer, dimension(:,:), allocatable, intent(out) :: helicities
    end subroutine term_instance_get_helicities_for_openloops
    elemental module function term_instance_get_i_term_global &
         (term) result (i_term)
      integer :: i_term
      class(term_instance_t), intent(in) :: term
    end function term_instance_get_i_term_global
    elemental module function term_instance_is_subtraction (term) result (sub)
      logical :: sub
      class(term_instance_t), intent(in) :: term
    end function term_instance_is_subtraction
    module function term_instance_get_n_sub (term) result (n_sub)
      integer :: n_sub
      class(term_instance_t), intent(in) :: term
    end function term_instance_get_n_sub
    module function term_instance_get_n_sub_color (term) result (n_sub_color)
      integer :: n_sub_color
      class(term_instance_t), intent(in) :: term
    end function term_instance_get_n_sub_color
    module function term_instance_get_n_sub_spin (term) result (n_sub_spin)
      integer :: n_sub_spin
      class(term_instance_t), intent(in) :: term
    end function term_instance_get_n_sub_spin
    module subroutine process_instance_write_header (object, unit, testflag)
      class(process_instance_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine process_instance_write_header
    module subroutine process_instance_write (object, unit, testflag)
      class(process_instance_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine process_instance_write
    module subroutine process_instance_init (instance, process)
      class(process_instance_t), intent(out), target :: instance
      type(process_t), intent(inout), target :: process
    end subroutine process_instance_init
    module subroutine process_instance_final (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_final
    module subroutine process_instance_reset (instance, reset_mci)
      class(process_instance_t), intent(inout), target :: instance
      logical, intent(in), optional :: reset_mci
    end subroutine process_instance_reset
    module subroutine process_instance_sampler_test (instance, i_mci, n_calls)
      class(process_instance_t), intent(inout), target :: instance
      integer, intent(in) :: i_mci
      integer, intent(in) :: n_calls
    end subroutine process_instance_sampler_test
    module subroutine process_instance_generate_weighted_event (instance, i_mci)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_mci
    end subroutine process_instance_generate_weighted_event
    module subroutine process_instance_generate_unweighted_event &
         (instance, i_mci)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_mci
    end subroutine process_instance_generate_unweighted_event
    module subroutine process_instance_recover_event (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_recover_event
    module subroutine process_instance_activate (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_activate
    module subroutine process_instance_find_same_kinematics (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_find_same_kinematics
    module subroutine process_instance_transfer_same_kinematics &
         (instance, i_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_term
    end subroutine process_instance_transfer_same_kinematics
    module subroutine process_instance_redo_sf_chains &
         (instance, i_term, phs_channel)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in), dimension(:) :: i_term
      integer, intent(in) :: phs_channel
    end subroutine process_instance_redo_sf_chains
    module subroutine process_instance_integrate (instance, i_mci, &
         n_it, n_calls, adapt_grids, adapt_weights, final, pacify)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_mci
      integer, intent(in) :: n_it
      integer, intent(in) :: n_calls
      logical, intent(in), optional :: adapt_grids
      logical, intent(in), optional :: adapt_weights
      logical, intent(in), optional :: final, pacify
    end subroutine process_instance_integrate
    module subroutine process_instance_setup_sf_chain (instance, config)
      class(process_instance_t), intent(inout) :: instance
      type(process_beam_config_t), intent(in), target :: config
    end subroutine process_instance_setup_sf_chain
    module subroutine process_instance_setup_event_data &
         (instance, model, i_core)
      class(process_instance_t), intent(inout), target :: instance
      class(model_data_t), intent(in), optional, target :: model
      integer, intent(in), optional :: i_core
    end subroutine process_instance_setup_event_data
    module subroutine process_instance_choose_mci (instance, i_mci)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_mci
    end subroutine process_instance_choose_mci
    module subroutine process_instance_set_mcpar (instance, x, warmup_flag)
      class(process_instance_t), intent(inout) :: instance
      real(default), dimension(:), intent(in) :: x
      logical, intent(in), optional :: warmup_flag
    end subroutine process_instance_set_mcpar
    module subroutine process_instance_receive_beam_momenta (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_receive_beam_momenta
    module subroutine process_instance_set_beam_momenta (instance, p)
      class(process_instance_t), intent(inout) :: instance
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine process_instance_set_beam_momenta
    module subroutine process_instance_recover_beam_momenta (instance, i_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_term
    end subroutine process_instance_recover_beam_momenta
    module subroutine process_instance_select_channel (instance, channel)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: channel
    end subroutine process_instance_select_channel
    module subroutine process_instance_compute_seed_kinematics &
         (instance, recover, skip_term)
      class(process_instance_t), intent(inout) :: instance
      logical, intent(in), optional :: recover
      integer, intent(in), optional :: skip_term
    end subroutine process_instance_compute_seed_kinematics
    pure module function process_instance_get_x_process (instance) result (x)
      real(default), dimension(:), allocatable :: x
      class(process_instance_t), intent(in) :: instance
    end function process_instance_get_x_process
    pure module function process_instance_get_active_component_type &
         (instance) result (nlo_type)
      integer :: nlo_type
      class(process_instance_t), intent(in) :: instance
    end function process_instance_get_active_component_type
    module subroutine process_instance_recover_mcpar (instance, i_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_term
    end subroutine process_instance_recover_mcpar
    module subroutine process_instance_recover_sfchain (instance, i_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_term
    end subroutine process_instance_recover_sfchain
    module subroutine process_instance_compute_hard_kinematics &
         (instance, recover, skip_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in), optional :: skip_term
      logical, intent(in), optional :: recover
    end subroutine process_instance_compute_hard_kinematics
    module subroutine process_instance_recover_seed_kinematics &
         (instance, i_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_term
    end subroutine process_instance_recover_seed_kinematics
    module subroutine process_instance_compute_eff_kinematics &
         (instance, skip_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in), optional :: skip_term
    end subroutine process_instance_compute_eff_kinematics
    module subroutine process_instance_recover_hard_kinematics &
         (instance, i_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_term
    end subroutine process_instance_recover_hard_kinematics
    module subroutine process_instance_evaluate_expressions &
         (instance, scale_forced)
      class(process_instance_t), intent(inout) :: instance
      real(default), intent(in), allocatable, optional :: scale_forced
    end subroutine process_instance_evaluate_expressions
    module subroutine process_instance_compute_other_channels &
         (instance, skip_term)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in), optional :: skip_term
    end subroutine process_instance_compute_other_channels
    module subroutine process_instance_reset_core_kinematics (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_reset_core_kinematics
    module subroutine process_instance_evaluate_trace (instance, recover)
      class(process_instance_t), intent(inout) :: instance
      logical, intent(in), optional :: recover
    end subroutine process_instance_evaluate_trace
    module subroutine term_instance_set_born_sqmes (term, core)
      class(term_instance_t), intent(inout) :: term
      class(prc_core_t), intent(in) :: core
    end subroutine term_instance_set_born_sqmes
    module subroutine term_instance_set_sf_factors (term, kin, has_pdfs)
      class(term_instance_t), intent(inout) :: term
      type(kinematics_t), intent(inout) :: kin
      logical, intent(in) :: has_pdfs
    end subroutine term_instance_set_sf_factors
    module subroutine process_instance_apply_real_partition (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_apply_real_partition
    module subroutine process_instance_set_i_mci_to_real_component (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_set_i_mci_to_real_component
    module subroutine process_instance_evaluate_event_data (instance, weight)
      class(process_instance_t), intent(inout) :: instance
      real(default), intent(in), optional :: weight
    end subroutine process_instance_evaluate_event_data
    module subroutine process_instance_compute_sqme_rad (instance, &
         i_term, i_phs, is_subtraction, alpha_s_external, scale_forced)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_term, i_phs
      logical, intent(in) :: is_subtraction
      real(default), intent(in), optional :: alpha_s_external
      real(default), intent(in), allocatable, optional :: scale_forced
    end subroutine process_instance_compute_sqme_rad
    module subroutine process_instance_normalize_weight (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_normalize_weight
    module subroutine process_instance_evaluate_sqme (instance, channel, x)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: channel
      real(default), dimension(:), intent(in) :: x
    end subroutine process_instance_evaluate_sqme
    module subroutine process_instance_recover &
         (instance, channel, i_term, update_sqme, recover_phs, scale_forced)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: channel
      integer, intent(in) :: i_term
      logical, intent(in) :: update_sqme
      logical, intent(in) :: recover_phs
      real(default), intent(in), allocatable, optional :: scale_forced
    end subroutine process_instance_recover
    module subroutine process_instance_evaluate (sampler, c, x_in, val, x, f)
      class(process_instance_t), intent(inout) :: sampler
      integer, intent(in) :: c
      real(default), dimension(:), intent(in) :: x_in
      real(default), intent(out) :: val
      real(default), dimension(:,:), intent(out) :: x
      real(default), dimension(:), intent(out) :: f
    end subroutine process_instance_evaluate
    module function process_instance_is_valid (sampler) result (valid)
      class(process_instance_t), intent(in) :: sampler
      logical :: valid
    end function process_instance_is_valid
    module subroutine process_instance_append_after_hook (sampler, new_hook)
      class(process_instance_t), intent(inout), target :: sampler
      class(process_instance_hook_t), intent(inout), target :: new_hook
    end subroutine process_instance_append_after_hook
    module subroutine process_instance_evaluate_after_hook (sampler)
      class(process_instance_t), intent(in) :: sampler
    end subroutine process_instance_evaluate_after_hook
    module subroutine process_instance_rebuild (sampler, c, x_in, val, x, f)
      class(process_instance_t), intent(inout) :: sampler
      integer, intent(in) :: c
      real(default), dimension(:), intent(in) :: x_in
      real(default), intent(in) :: val
      real(default), dimension(:,:), intent(out) :: x
      real(default), dimension(:), intent(out) :: f
    end subroutine process_instance_rebuild
    module subroutine process_instance_fetch (sampler, val, x, f)
      class(process_instance_t), intent(in) :: sampler
      real(default), intent(out) :: val
      real(default), dimension(:,:), intent(out) :: x
      real(default), dimension(:), intent(out) :: f
    end subroutine process_instance_fetch
    module subroutine process_instance_init_simulation (instance, i_mci, &
       safety_factor, keep_failed_events)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_mci
      real(default), intent(in), optional :: safety_factor
      logical, intent(in), optional :: keep_failed_events
    end subroutine process_instance_init_simulation
    module subroutine process_instance_final_simulation (instance, i_mci)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_mci
    end subroutine process_instance_final_simulation
    module subroutine process_instance_get_mcpar (instance, channel, x)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: channel
      real(default), dimension(:), intent(out) :: x
    end subroutine process_instance_get_mcpar
    module function process_instance_has_evaluated_trace &
         (instance) result (flag)
      class(process_instance_t), intent(in) :: instance
      logical :: flag
    end function process_instance_has_evaluated_trace
    module function process_instance_is_complete_event (instance) result (flag)
      class(process_instance_t), intent(in) :: instance
      logical :: flag
    end function process_instance_is_complete_event
    module function process_instance_select_i_term (instance) result (i_term)
      integer :: i_term
      class(process_instance_t), intent(in) :: instance
    end function process_instance_select_i_term
    module function process_instance_get_beam_int_ptr (instance) result (ptr)
      class(process_instance_t), intent(in), target :: instance
      type(interaction_t), pointer :: ptr
    end function process_instance_get_beam_int_ptr
    module function process_instance_get_trace_int_ptr &
         (instance, i_term) result (ptr)
      class(process_instance_t), intent(in), target :: instance
      integer, intent(in) :: i_term
      type(interaction_t), pointer :: ptr
    end function process_instance_get_trace_int_ptr
    module function process_instance_get_matrix_int_ptr &
         (instance, i_term) result (ptr)
      class(process_instance_t), intent(in), target :: instance
      integer, intent(in) :: i_term
      type(interaction_t), pointer :: ptr
    end function process_instance_get_matrix_int_ptr
    module function process_instance_get_flows_int_ptr &
         (instance, i_term) result (ptr)
      class(process_instance_t), intent(in), target :: instance
      integer, intent(in) :: i_term
      type(interaction_t), pointer :: ptr
    end function process_instance_get_flows_int_ptr
    module function process_instance_get_state_flv &
         (instance, i_term) result (state_flv)
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_term
      type(state_flv_content_t) :: state_flv
    end function process_instance_get_state_flv
    module function process_instance_get_isolated_state_ptr &
         (instance, i_term) result (ptr)
      class(process_instance_t), intent(in), target :: instance
      integer, intent(in) :: i_term
      type(isolated_state_t), pointer :: ptr
    end function process_instance_get_isolated_state_ptr
    module function process_instance_get_connected_state_ptr &
         (instance, i_term) result (ptr)
      class(process_instance_t), intent(in), target :: instance
      integer, intent(in) :: i_term
      type(connected_state_t), pointer :: ptr
    end function process_instance_get_connected_state_ptr
    module subroutine process_instance_get_beam_index (instance, i_term, i_beam)
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_term
      integer, dimension(:), intent(out) :: i_beam
    end subroutine process_instance_get_beam_index
    module subroutine process_instance_get_in_index (instance, i_term, i_in)
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_term
      integer, dimension(:), intent(out) :: i_in
    end subroutine process_instance_get_in_index
    module function process_instance_get_sqme (instance, i_term) result (sqme)
      real(default) :: sqme
      class(process_instance_t), intent(in) :: instance
      integer, intent(in), optional :: i_term
    end function process_instance_get_sqme
    module function process_instance_get_weight (instance) result (weight)
      real(default) :: weight
      class(process_instance_t), intent(in) :: instance
    end function process_instance_get_weight
    module function process_instance_get_excess (instance) result (excess)
      real(default) :: excess
      class(process_instance_t), intent(in) :: instance
    end function process_instance_get_excess
    module function process_instance_get_n_dropped (instance) result (n_dropped)
      integer :: n_dropped
      class(process_instance_t), intent(in) :: instance
    end function process_instance_get_n_dropped
    module function process_instance_get_channel (instance) result (channel)
      integer :: channel
      class(process_instance_t), intent(in) :: instance
    end function process_instance_get_channel
    module subroutine process_instance_set_fac_scale (instance, fac_scale)
      class(process_instance_t), intent(inout) :: instance
      real(default), intent(in) :: fac_scale
    end subroutine process_instance_set_fac_scale
    module function process_instance_get_fac_scale &
         (instance, i_term) result (fac_scale)
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_term
      real(default) :: fac_scale
    end function process_instance_get_fac_scale
    module function process_instance_get_alpha_s &
         (instance, i_term) result (alpha_s)
      real(default) :: alpha_s
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_term
    end function process_instance_get_alpha_s
    module function process_instance_get_qcd (process_instance) result (qcd)
      type(qcd_t) :: qcd
      class(process_instance_t), intent(in) :: process_instance
    end function process_instance_get_qcd
    module subroutine process_instance_reset_counter (process_instance)
      class(process_instance_t), intent(inout) :: process_instance
    end subroutine process_instance_reset_counter
    module subroutine process_instance_record_call (process_instance)
      class(process_instance_t), intent(inout) :: process_instance
    end subroutine process_instance_record_call
    pure module function process_instance_get_counter &
         (process_instance) result (counter)
      class(process_instance_t), intent(in) :: process_instance
      type(process_counter_t) :: counter
    end function process_instance_get_counter
    pure module function process_instance_get_actual_calls_total &
         (process_instance) result (n)
      class(process_instance_t), intent(in) :: process_instance
      integer :: n
    end function process_instance_get_actual_calls_total
    module subroutine process_instance_reset_matrix_elements (instance)
      class(process_instance_t), intent(inout) :: instance
    end subroutine process_instance_reset_matrix_elements
    module subroutine process_instance_get_test_phase_space_point (instance, &
           i_component, i_core, p)
      type(vector4_t), dimension(:), allocatable, intent(out) :: p
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_component, i_core
    end subroutine process_instance_get_test_phase_space_point
    pure module function process_instance_get_p_hard &
         (process_instance, i_term) result (p_hard)
      type(vector4_t), dimension(:), allocatable :: p_hard
      class(process_instance_t), intent(in) :: process_instance
      integer, intent(in) :: i_term
    end function process_instance_get_p_hard
    module function process_instance_get_first_active_i_term &
         (instance) result (i_term)
      integer :: i_term
      class(process_instance_t), intent(in) :: instance
    end function process_instance_get_first_active_i_term
    module function process_instance_get_real_of_mci (instance) result (i_real)
      integer :: i_real
      class(process_instance_t), intent(in) :: instance
    end function process_instance_get_real_of_mci
    module function process_instance_get_connected_states &
         (instance, i_component) result (connected)
      type(connected_state_t), dimension(:), allocatable :: connected
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_component
    end function process_instance_get_connected_states
    module function process_instance_get_sqrts (instance) result (sqrts)
      class(process_instance_t), intent(in) :: instance
      real(default) :: sqrts
    end function process_instance_get_sqrts
    module function process_instance_get_polarization (instance) result (pol)
      class(process_instance_t), intent(in) :: instance
      real(default), dimension(:), allocatable :: pol
    end function process_instance_get_polarization
    module function process_instance_get_beam_file (instance) result (file)
      class(process_instance_t), intent(in) :: instance
      type(string_t) :: file
    end function process_instance_get_beam_file
    module function process_instance_get_process_name (instance) result (name)
      class(process_instance_t), intent(in) :: instance
      type(string_t) :: name
    end function process_instance_get_process_name
    module subroutine process_instance_get_trace &
         (instance, pset, i_term, n_incoming)
      class(process_instance_t), intent(in), target :: instance
      type(particle_set_t), intent(out) :: pset
      integer, intent(in) :: i_term
      integer, intent(in), optional :: n_incoming
    end subroutine process_instance_get_trace
    module subroutine process_instance_set_trace &
         (instance, pset, i_term, recover_beams, check_match, success)
      class(process_instance_t), intent(inout), target :: instance
      type(particle_set_t), intent(in) :: pset
      integer, intent(in) :: i_term
      logical, intent(in), optional :: recover_beams, check_match
      logical, intent(out), optional :: success
    end subroutine process_instance_set_trace
    module subroutine process_instance_set_alpha_qcd_forced &
         (instance, i_term, alpha_qcd)
      class(process_instance_t), intent(inout) :: instance
      integer, intent(in) :: i_term
      real(default), intent(in) :: alpha_qcd
    end subroutine process_instance_set_alpha_qcd_forced
    module function process_instance_has_nlo_component (instance) result (nlo)
      class(process_instance_t), intent(in) :: instance
      logical :: nlo
    end function process_instance_has_nlo_component
    module function process_instance_keep_failed_events (instance) result (keep)
      logical :: keep
      class(process_instance_t), intent(in) :: instance
    end function process_instance_keep_failed_events
    module function process_instance_get_term_indices &
         (instance, nlo_type) result (i_term)
      integer, dimension(:), allocatable :: i_term
      class(process_instance_t), intent(in) :: instance
      integer :: nlo_type
    end function process_instance_get_term_indices
    module function process_instance_get_boost_to_lab &
         (instance, i_term) result (lt)
      type(lorentz_transformation_t) :: lt
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_term
    end function process_instance_get_boost_to_lab
    module function process_instance_get_boost_to_cms &
         (instance, i_term) result (lt)
      type(lorentz_transformation_t) :: lt
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_term
    end function process_instance_get_boost_to_cms
    module function process_instance_lab_is_cm &
         (instance, i_term) result (lab_is_cm)
      logical :: lab_is_cm
      class(process_instance_t), intent(in) :: instance
      integer, intent(in) :: i_term
    end function process_instance_lab_is_cm
    module subroutine pacify_process_instance (instance)
      type(process_instance_t), intent(inout) :: instance
    end subroutine pacify_process_instance
  end interface

contains

  subroutine term_instance_evaluate_scaled_sf_chains (term, kin)
    class(term_instance_t), intent(inout) :: term
    type(kinematics_t), intent(inout) :: kin
    class(sf_rescale_t), allocatable :: sf_rescale
    if (.not. term%pcm%has_pdfs) return
    if (term%nlo_type == NLO_REAL) then
       if (term%is_subtraction ()) then
          allocate (sf_rescale_collinear_t :: sf_rescale)
          select type (pcm_work => term%pcm_work)
          type is (pcm_nlo_workspace_t)
             select type (sf_rescale)
             type is (sf_rescale_collinear_t)
                call sf_rescale%set (pcm_work%real_kinematics%xi_tilde)
             end select
          end select
          call kin%sf_chain%evaluate (term%get_fac_scale (), &
               term%negative_sf, sf_rescale)
          deallocate (sf_rescale)
       else if (kin%emitter >= 0 .and. kin%emitter <= kin%n_in) then
          allocate (sf_rescale_real_t :: sf_rescale)
          select type (pcm_work => term%pcm_work)
          type is (pcm_nlo_workspace_t)
             select type (sf_rescale)
             type is (sf_rescale_real_t)
                call sf_rescale%set (pcm_work%real_kinematics%xi_tilde * &
                     pcm_work%real_kinematics%xi_max (kin%i_phs), &
                     pcm_work%real_kinematics%y (kin%i_phs))
             end select
          end select
          call kin%sf_chain%evaluate (term%get_fac_scale (), &
               term%negative_sf, sf_rescale)
          deallocate (sf_rescale)
       else
          call kin%sf_chain%evaluate (term%get_fac_scale (), term%negative_sf)
       end if
    else if (term%nlo_type == NLO_DGLAP) then
       allocate (sf_rescale_dglap_t :: sf_rescale)
       select type (pcm_work => term%pcm_work)
       type is (pcm_nlo_workspace_t)
          select type (sf_rescale)
          type is (sf_rescale_dglap_t)
             call sf_rescale%set (pcm_work%isr_kinematics%z)
          end select
       end select
       call kin%sf_chain%evaluate (term%get_fac_scale (), &
            term%negative_sf, sf_rescale)
       deallocate (sf_rescale)
    end if
  end subroutine term_instance_evaluate_scaled_sf_chains


end module instances
