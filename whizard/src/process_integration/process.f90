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

module process

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use diagnostics
  use lorentz
  use rng_base
  use dispatch_rng, only: dispatch_rng_factory
  use dispatch_rng, only: update_rng_seed_in_var_list
  use os_interface
  use sm_qcd
  use mci_base
  use flavors
  use model_data
  use models
  use process_libraries
  use process_constants
  use variables
  use beam_structures
  use beams
  use pdg_arrays
  use expr_base
  use sf_base
  use sf_mappings
  use resonances, only: resonance_history_t, resonance_history_set_t

  use prc_test_core, only: test_t
  use prc_core_def, only: prc_core_def_t
  use prc_core, only: prc_core_t, helicity_selection_t
  use phs_base

  use parton_states, only: connected_state_t
  use pcm_base
  use pcm
  use process_counter
  use process_config
  use process_mci

  implicit none
  private

  public :: process_t
  public :: process_ptr_t

  type :: process_status_t
     private
  end type process_status_t

  type :: process_results_t
     private
  end type process_results_t

  type :: process_t
     private
     type(process_metadata_t) :: &
          meta
     type(process_environment_t) :: &
          env
     type(process_config_data_t) :: &
          config
     class(pcm_t), allocatable :: &
          pcm
     type(process_component_t), dimension(:), allocatable :: &
          component
     type(process_phs_config_t), dimension(:), allocatable :: &
          phs_entry
     type(core_entry_t), dimension(:), allocatable :: &
          core_entry
     type(process_mci_entry_t), dimension(:), allocatable :: &
          mci_entry
     class(rng_factory_t), allocatable :: &
          rng_factory
     type(process_beam_config_t) :: &
          beam_config
     type(process_term_t), dimension(:), allocatable :: &
          term
     type(process_status_t) :: &
          status
     type(process_results_t) :: &
          result
   contains
     procedure :: write => process_write
     ! generic :: write (formatted) => write_formatted
     procedure :: write_formatted => process_write_formatted
     procedure :: write_meta => process_write_meta
     procedure :: show => process_show
     procedure :: final => process_final
     procedure :: init => process_init
     procedure :: complete_pcm_setup => process_complete_pcm_setup
     procedure :: setup_cores => process_setup_cores
     procedure :: prepare_blha_cores => process_prepare_blha_cores
     procedure :: create_blha_interface => process_create_blha_interface
     procedure :: init_components => process_init_components
     procedure :: record_inactive_components => process_record_inactive_components
     procedure :: setup_terms => process_setup_terms
     procedure :: setup_beams_sqrts => process_setup_beams_sqrts
     procedure :: setup_beams_decay => process_setup_beams_decay
     procedure :: check_masses => process_check_masses
     procedure :: optimize_nlo_singular_regions => &
          process_optimize_nlo_singular_regions
     procedure :: get_pdg_in => process_get_pdg_in
     procedure :: get_phs_config => process_get_phs_config
     procedure :: extract_resonance_history_set &
          => process_extract_resonance_history_set
     procedure :: setup_beams_beam_structure => process_setup_beams_beam_structure
     procedure :: beams_startup_message => process_beams_startup_message
     procedure :: init_phs_config => process_init_phs_config
     procedure :: configure_phs => process_configure_phs
     procedure :: print_phs_startup_message => process_print_phs_startup_message
     procedure :: init_sf_chain => process_init_sf_chain
     generic :: set_sf_channel => set_sf_channel_single
     procedure :: set_sf_channel_single => process_set_sf_channel
     generic :: set_sf_channel => set_sf_channel_array
     procedure :: set_sf_channel_array => process_set_sf_channel_array
     procedure :: sf_startup_message => process_sf_startup_message
     procedure :: collect_channels => process_collect_channels
     procedure :: contains_trivial_component => process_contains_trivial_component
     procedure :: get_master_component => process_get_master_component
     procedure :: setup_mci => process_setup_mci
     procedure :: set_cuts => process_set_cuts
     procedure :: set_scale => process_set_scale
     procedure :: set_fac_scale => process_set_fac_scale
     procedure :: set_ren_scale => process_set_ren_scale
     procedure :: set_weight => process_set_weight
     procedure :: compute_md5sum => process_compute_md5sum
     procedure :: sampler_test => process_sampler_test
     procedure :: final_integration => process_final_integration
     procedure :: integrate_dummy => process_integrate_dummy
     procedure :: integrate => process_integrate
     procedure :: generate_weighted_event => process_generate_weighted_event
     procedure :: generate_unweighted_event => process_generate_unweighted_event
     procedure :: display_summed_results => process_display_summed_results
     procedure :: display_integration_history => &
          process_display_integration_history
     procedure :: write_logfile => process_write_logfile
     procedure :: write_state_summary => process_write_state_summary
     procedure :: prepare_simulation => process_prepare_simulation
     generic :: has_integral => has_integral_tot, has_integral_mci
     procedure :: has_integral_tot => process_has_integral_tot
     procedure :: has_integral_mci => process_has_integral_mci
     generic :: get_integral => get_integral_tot, get_integral_mci
     generic :: get_error => get_error_tot, get_error_mci
     generic :: get_efficiency => get_efficiency_tot, get_efficiency_mci
     procedure :: get_integral_tot => process_get_integral_tot
     procedure :: get_integral_mci => process_get_integral_mci
     procedure :: get_error_tot => process_get_error_tot
     procedure :: get_error_mci => process_get_error_mci
     procedure :: get_efficiency_tot => process_get_efficiency_tot
     procedure :: get_efficiency_mci => process_get_efficiency_mci
     procedure :: get_correction => process_get_correction
     procedure :: get_correction_error => process_get_correction_error
     procedure :: lab_is_cm => process_lab_is_cm
     procedure :: get_component_ptr => process_get_component_ptr
     procedure :: get_qcd => process_get_qcd
     generic :: get_component_type => get_component_type_single
     procedure :: get_component_type_single => process_get_component_type_single
     generic :: get_component_type => get_component_type_all
     procedure :: get_component_type_all => process_get_component_type_all
     procedure :: get_component_i_terms => process_get_component_i_terms
     procedure :: get_n_allowed_born => process_get_n_allowed_born
     procedure :: get_pcm_ptr => process_get_pcm_ptr
     generic :: component_can_be_integrated => component_can_be_integrated_single
     generic :: component_can_be_integrated => component_can_be_integrated_all
     procedure :: component_can_be_integrated_single => &
          process_component_can_be_integrated_single
     procedure :: component_can_be_integrated_all => &
          process_component_can_be_integrated_all
     procedure :: reset_selected_cores => process_reset_selected_cores
     procedure :: select_components => process_select_components
     procedure :: component_is_selected => process_component_is_selected
     procedure :: get_coupling_powers => process_get_coupling_powers
     procedure :: get_real_component => process_get_real_component
     procedure :: extract_active_component_mci => &
          process_extract_active_component_mci
     procedure :: uses_real_partition => process_uses_real_partition
     procedure :: get_md5sum_prc => process_get_md5sum_prc
     procedure :: get_md5sum_mci => process_get_md5sum_mci
     procedure :: get_md5sum_cfg => process_get_md5sum_cfg
     procedure :: get_n_cores => process_get_n_cores
     procedure :: get_base_i_term => process_get_base_i_term
     procedure :: get_core_term => process_get_core_term
     procedure :: get_core_ptr => process_get_core_ptr
     procedure :: get_term_ptr => process_get_term_ptr
     procedure :: get_i_term => process_get_i_term
     procedure :: get_i_core => process_get_i_core
     procedure :: set_i_mci_work => process_set_i_mci_work
     procedure :: get_i_mci_work => process_get_i_mci_work
     procedure :: get_i_sub => process_get_i_sub
     procedure :: get_i_term_virtual => process_get_i_term_virtual
     generic :: component_is_active => component_is_active_single
     procedure :: component_is_active_single => process_component_is_active_single
     generic :: component_is_active => component_is_active_all
     procedure :: component_is_active_all => process_component_is_active_all
     procedure :: get_n_pass_default => process_get_n_pass_default
     procedure :: adapt_grids_default => process_adapt_grids_default
     procedure :: adapt_weights_default => process_adapt_weights_default
     procedure :: get_n_it_default => process_get_n_it_default
     procedure :: get_n_calls_default => process_get_n_calls_default
     procedure :: set_run_id => process_set_run_id
     procedure :: get_id => process_get_id
     procedure :: get_num_id => process_get_num_id
     procedure :: get_run_id => process_get_run_id
     procedure :: get_library_name => process_get_library_name
     procedure :: get_n_in => process_get_n_in
     procedure :: get_n_mci => process_get_n_mci
     procedure :: get_n_components => process_get_n_components
     procedure :: get_n_terms => process_get_n_terms
     procedure :: get_i_component => process_get_i_component
     procedure :: get_component_id => process_get_component_id
     procedure :: get_component_def_ptr => process_get_component_def_ptr
     procedure :: extract_core => process_extract_core
     procedure :: restore_core => process_restore_core
     procedure :: get_constants => process_get_constants
     procedure :: get_config => process_get_config
     procedure :: get_md5sum_constants => process_get_md5sum_constants
     procedure :: get_term_flv_out => process_get_term_flv_out
     procedure :: contains_unstable => process_contains_unstable
     procedure :: get_sqrts => process_get_sqrts
     procedure :: get_energy => process_get_energy
     procedure :: get_polarization => process_get_polarization
     procedure :: get_meta => process_get_meta
     procedure :: has_matrix_element => process_has_matrix_element
     procedure :: get_beam_data_ptr => process_get_beam_data_ptr
     procedure :: get_beam_config => process_get_beam_config
     procedure :: get_beam_config_ptr => process_get_beam_config_ptr
     procedure :: get_pdf_set => process_get_pdf_set
     procedure :: pcm_contains_pdfs => process_pcm_contains_pdfs
     procedure :: get_beam_file => process_get_beam_file
     procedure :: get_var_list_ptr => process_get_var_list_ptr
     procedure :: get_model_ptr => process_get_model_ptr
     procedure :: make_rng => process_make_rng
     procedure :: compute_amplitude => process_compute_amplitude
     procedure :: check_library_sanity => process_check_library_sanity
     procedure :: reset_library_ptr => process_reset_library_ptr
     procedure :: set_counter_mci_entry => process_set_counter_mci_entry
     procedure :: pacify => process_pacify
     procedure :: test_allocate_sf_channels
     procedure :: test_set_component_sf_channel
     procedure :: test_get_mci_ptr
     procedure :: init_mci_work => process_init_mci_work
     procedure :: setup_test_cores => process_setup_test_cores
     procedure :: get_connected_states => process_get_connected_states
     procedure :: init_nlo_settings => process_init_nlo_settings
     generic :: get_nlo_type_component => get_nlo_type_component_single
     procedure :: get_nlo_type_component_single => &
          process_get_nlo_type_component_single
     generic :: get_nlo_type_component => get_nlo_type_component_all
     procedure :: get_nlo_type_component_all => process_get_nlo_type_component_all
     procedure :: is_nlo_calculation => process_is_nlo_calculation
     procedure :: get_negative_sf => process_get_negative_sf
     procedure :: is_combined_nlo_integration &
          => process_is_combined_nlo_integration
     procedure :: component_is_real_finite => process_component_is_real_finite
     procedure :: get_component_nlo_type => process_get_component_nlo_type
     procedure :: get_component_core_ptr => process_get_component_core_ptr
     procedure :: get_component_associated_born &
               => process_get_component_associated_born
     procedure :: get_first_real_component => process_get_first_real_component
     procedure :: get_first_real_term => process_get_first_real_term
     procedure :: get_associated_real_fin => process_get_associated_real_fin
     procedure :: select_i_term => process_select_i_term
     procedure :: prepare_any_external_code &
        => process_prepare_any_external_code
  end type process_t

  type :: process_ptr_t
     type(process_t), pointer :: p => null ()
  end type process_ptr_t


  abstract interface
     subroutine dispatch_core_proc (core, core_def, model, &
          helicity_selection, qcd, use_color_factors, has_beam_pol)
       import
       class(prc_core_t), allocatable, intent(inout) :: core
       class(prc_core_def_t), intent(in) :: core_def
       class(model_data_t), intent(in), target, optional :: model
       type(helicity_selection_t), intent(in), optional :: helicity_selection
       type(qcd_t), intent(in), optional :: qcd
       logical, intent(in), optional :: use_color_factors
       logical, intent(in), optional :: has_beam_pol
     end subroutine dispatch_core_proc
  end interface


  interface
    module subroutine process_write (process, screen, unit, &
         show_os_data, show_var_list, show_rng, show_expressions, pacify)
      class(process_t), intent(in) :: process
      logical, intent(in) :: screen
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_os_data
      logical, intent(in), optional :: show_var_list
      logical, intent(in), optional :: show_rng
      logical, intent(in), optional :: show_expressions
      logical, intent(in), optional :: pacify
    end subroutine process_write
    module subroutine process_write_formatted (dtv, unit, iotype, &
         v_list, iostat, iomsg)
      class(process_t), intent(in) :: dtv
      integer, intent(in) :: unit
      character(*), intent(in) :: iotype
      integer, dimension(:), intent(in) :: v_list
      integer, intent(out) :: iostat
      character(*), intent(inout) :: iomsg
    end subroutine process_write_formatted
    module subroutine process_write_meta (process, unit, testflag)
      class(process_t), intent(in) :: process
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine process_write_meta
    module subroutine process_show (object, unit, verbose)
      class(process_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine process_show
    module subroutine process_final (process)
      class(process_t), intent(inout) :: process
    end subroutine process_final
    module subroutine process_complete_pcm_setup (process)
      class(process_t), intent(inout) :: process
    end subroutine process_complete_pcm_setup
    module subroutine process_setup_cores (process, dispatch_core, &
         helicity_selection, use_color_factors, has_beam_pol)
      class(process_t), intent(inout) :: process
      procedure(dispatch_core_proc) :: dispatch_core
      type(helicity_selection_t), intent(in), optional :: helicity_selection
      logical, intent(in), optional :: use_color_factors
      logical, intent(in), optional :: has_beam_pol
    end subroutine process_setup_cores
    module subroutine process_prepare_blha_cores (process)
      class(process_t), intent(inout), target :: process
    end subroutine process_prepare_blha_cores
    module subroutine process_create_blha_interface (process)
      class(process_t), intent(inout) :: process
    end subroutine process_create_blha_interface
    module subroutine process_init_components (process, phs_config)
      class(process_t), intent(inout), target :: process
      class(phs_config_t), allocatable, intent(in), optional :: phs_config
    end subroutine process_init_components
    module subroutine process_record_inactive_components (process)
      class(process_t), intent(inout) :: process
    end subroutine process_record_inactive_components
    module subroutine process_setup_terms (process, with_beams)
      class(process_t), intent(inout), target :: process
      logical, intent(in), optional :: with_beams
    end subroutine process_setup_terms
    module subroutine process_setup_beams_sqrts &
         (process, sqrts, beam_structure, i_core)
      class(process_t), intent(inout) :: process
      real(default), intent(in) :: sqrts
      type(beam_structure_t), intent(in), optional :: beam_structure
      integer, intent(in), optional :: i_core
    end subroutine process_setup_beams_sqrts
    module subroutine process_setup_beams_decay &
         (process, rest_frame, beam_structure, i_core)
      class(process_t), intent(inout), target :: process
      logical, intent(in), optional :: rest_frame
      type(beam_structure_t), intent(in), optional :: beam_structure
      integer, intent(in), optional :: i_core
    end subroutine process_setup_beams_decay
    module subroutine process_check_masses (process)
      class(process_t), intent(in) :: process
    end subroutine process_check_masses
    module subroutine process_optimize_nlo_singular_regions (process)
      class(process_t), intent(inout) :: process
    end subroutine process_optimize_nlo_singular_regions
    module subroutine process_get_pdg_in (process, pdg_in)
      class(process_t), intent(in), target :: process
      type(pdg_array_t), dimension(:,:), allocatable, intent(out) :: pdg_in
    end subroutine process_get_pdg_in
    module function process_get_phs_config &
         (process, i_component) result (phs_config)
      class(phs_config_t), pointer :: phs_config
      class(process_t), intent(in), target :: process
      integer, intent(in) :: i_component
    end function process_get_phs_config
    module subroutine process_extract_resonance_history_set &
         (process, res_set, include_trivial, i_component)
      class(process_t), intent(in), target :: process
      type(resonance_history_set_t), intent(out) :: res_set
      logical, intent(in), optional :: include_trivial
      integer, intent(in), optional :: i_component
    end subroutine process_extract_resonance_history_set
    module subroutine process_setup_beams_beam_structure &
         (process, beam_structure, sqrts, decay_rest_frame)
      class(process_t), intent(inout) :: process
      type(beam_structure_t), intent(in) :: beam_structure
      real(default), intent(in) :: sqrts
      logical, intent(in), optional :: decay_rest_frame
    end subroutine process_setup_beams_beam_structure
    module subroutine process_beams_startup_message &
         (process, unit, beam_structure)
      class(process_t), intent(in) :: process
      integer, intent(in), optional :: unit
      type(beam_structure_t), intent(in), optional :: beam_structure
    end subroutine process_beams_startup_message
    module subroutine process_init_phs_config (process)
      class(process_t), intent(inout) :: process
    end subroutine process_init_phs_config
    module subroutine process_configure_phs (process, rebuild, &
         ignore_mismatch, combined_integration, subdir)
      class(process_t), intent(inout) :: process
      logical, intent(in), optional :: rebuild
      logical, intent(in), optional :: ignore_mismatch
      logical, intent(in), optional :: combined_integration
      type(string_t), intent(in), optional :: subdir
    end subroutine process_configure_phs
    module subroutine process_print_phs_startup_message (process)
      class(process_t), intent(in) :: process
    end subroutine process_print_phs_startup_message
    module subroutine process_init_sf_chain (process, sf_config, sf_trace_file)
      class(process_t), intent(inout) :: process
      type(sf_config_t), dimension(:), intent(in) :: sf_config
      type(string_t), intent(in), optional :: sf_trace_file
    end subroutine process_init_sf_chain
    module subroutine process_set_sf_channel (process, c, sf_channel)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: c
      type(sf_channel_t), intent(in) :: sf_channel
    end subroutine process_set_sf_channel
    module subroutine process_set_sf_channel_array (process, sf_channel)
      class(process_t), intent(inout) :: process
      type(sf_channel_t), dimension(:), intent(in) :: sf_channel
    end subroutine process_set_sf_channel_array
    module subroutine process_sf_startup_message (process, sf_string, unit)
      class(process_t), intent(in) :: process
      type(string_t), intent(in) :: sf_string
      integer, intent(in), optional :: unit
    end subroutine process_sf_startup_message
    module subroutine process_collect_channels (process, coll)
      class(process_t), intent(inout) :: process
      type(phs_channel_collection_t), intent(inout) :: coll
    end subroutine process_collect_channels
    module function process_contains_trivial_component (process) result (flag)
      class(process_t), intent(in) :: process
      logical :: flag
    end function process_contains_trivial_component
    module function process_get_master_component &
         (process, i_mci) result (i_component)
      integer :: i_component
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
    end function process_get_master_component
    module subroutine process_setup_mci (process, dispatch_mci)
      class(process_t), intent(inout) :: process
      procedure(dispatch_mci_proc) :: dispatch_mci
    end subroutine process_setup_mci
    module subroutine process_set_cuts (process, ef_cuts)
      class(process_t), intent(inout) :: process
      class(expr_factory_t), intent(in) :: ef_cuts
    end subroutine process_set_cuts
    module subroutine process_set_scale (process, ef_scale)
      class(process_t), intent(inout) :: process
      class(expr_factory_t), intent(in) :: ef_scale
    end subroutine process_set_scale
    module subroutine process_set_weight (process, ef_weight)
      class(process_t), intent(inout) :: process
      class(expr_factory_t), intent(in) :: ef_weight
    end subroutine process_set_weight
    module subroutine process_set_fac_scale (process, ef_fac_scale)
      class(process_t), intent(inout) :: process
      class(expr_factory_t), intent(in) :: ef_fac_scale
    end subroutine process_set_fac_scale
    module subroutine process_set_ren_scale (process, ef_ren_scale)
      class(process_t), intent(inout) :: process
      class(expr_factory_t), intent(in) :: ef_ren_scale
    end subroutine process_set_ren_scale
    module subroutine process_compute_md5sum (process)
      class(process_t), intent(inout) :: process
    end subroutine process_compute_md5sum
    module subroutine process_sampler_test (process, sampler, n_calls, i_mci)
      class(process_t), intent(inout) :: process
      class(mci_sampler_t), intent(inout) :: sampler
      integer, intent(in) :: n_calls, i_mci
    end subroutine process_sampler_test
    module subroutine process_final_integration (process, i_mci)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
    end subroutine process_final_integration
    module subroutine process_integrate_dummy (process)
      class(process_t), intent(inout) :: process
    end subroutine process_integrate_dummy
    module subroutine process_integrate (process, i_mci, mci_work, &
       mci_sampler, n_it, n_calls, adapt_grids, adapt_weights, final, &
       pacify, nlo_type)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
      type(mci_work_t), intent(inout) :: mci_work
      class(mci_sampler_t), intent(inout) :: mci_sampler
      integer, intent(in) :: n_it, n_calls
      logical, intent(in), optional :: adapt_grids, adapt_weights
      logical, intent(in), optional :: final
      logical, intent(in), optional :: pacify
      integer, intent(in), optional :: nlo_type
    end subroutine process_integrate
    module subroutine process_generate_weighted_event (process, i_mci, &
         mci_work, mci_sampler, keep_failed_events)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
      type(mci_work_t), intent(inout) :: mci_work
      class(mci_sampler_t), intent(inout) :: mci_sampler
      logical, intent(in) :: keep_failed_events
    end subroutine process_generate_weighted_event
    module subroutine process_generate_unweighted_event (process, i_mci, &
       mci_work, mci_sampler)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
      type(mci_work_t), intent(inout) :: mci_work
      class(mci_sampler_t), intent(inout) :: mci_sampler
    end subroutine process_generate_unweighted_event
    module subroutine process_display_summed_results (process, pacify)
      class(process_t), intent(inout) :: process
      logical, intent(in) :: pacify
    end subroutine process_display_summed_results
    module subroutine process_display_integration_history &
         (process, i_mci, filename, os_data, eff_reset)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
      type(string_t), intent(in) :: filename
      type(os_data_t), intent(in) :: os_data
      logical, intent(in), optional :: eff_reset
    end subroutine process_display_integration_history
    module subroutine process_write_logfile (process, i_mci, filename)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
      type(string_t), intent(in) :: filename
    end subroutine process_write_logfile
    module subroutine process_write_state_summary (process, unit)
      class(process_t), intent(in) :: process
      integer, intent(in), optional :: unit
    end subroutine process_write_state_summary
    module subroutine process_prepare_simulation (process, i_mci)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
    end subroutine process_prepare_simulation
    module function process_has_integral_tot (process) result (flag)
      logical :: flag
      class(process_t), intent(in) :: process
    end function process_has_integral_tot
    module function process_has_integral_mci (process, i_mci) result (flag)
      logical :: flag
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
    end function process_has_integral_mci
    module function process_get_integral_mci (process, i_mci) result (integral)
      real(default) :: integral
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
    end function process_get_integral_mci
    module function process_get_error_mci (process, i_mci) result (error)
      real(default) :: error
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
    end function process_get_error_mci
    module function process_get_efficiency_mci &
         (process, i_mci) result (efficiency)
      real(default) :: efficiency
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
    end function process_get_efficiency_mci
    module function process_get_integral_tot (process) result (integral)
      real(default) :: integral
      class(process_t), intent(in) :: process
    end function process_get_integral_tot
    module function process_get_error_tot (process) result (error)
      real(default) :: variance
      class(process_t), intent(in) :: process
      real(default) :: error
    end function process_get_error_tot
    module function process_get_efficiency_tot (process) result (efficiency)
      real(default) :: efficiency
      class(process_t), intent(in) :: process
    end function process_get_efficiency_tot
    module function process_get_correction (process) result (ratio)
      real(default) :: ratio
      class(process_t), intent(in) :: process
    end function process_get_correction
    module function process_get_correction_error (process) result (error)
      real(default) :: error
      class(process_t), intent(in) :: process
    end function process_get_correction_error
    pure module function process_lab_is_cm (process) result (lab_is_cm)
      logical :: lab_is_cm
      class(process_t), intent(in) :: process
    end function process_lab_is_cm
    module function process_get_component_ptr (process, i) result (component)
      type(process_component_t), pointer :: component
      class(process_t), intent(in), target :: process
      integer, intent(in) :: i
    end function process_get_component_ptr
    module function process_get_qcd (process) result (qcd)
      type(qcd_t) :: qcd
      class(process_t), intent(in) :: process
    end function process_get_qcd
    elemental module function process_get_component_type_single &
       (process, i_component) result (comp_type)
      integer :: comp_type
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
    end function process_get_component_type_single
    module function process_get_component_type_all &
       (process) result (comp_type)
      integer, dimension(:), allocatable :: comp_type
      class(process_t), intent(in) :: process
    end function process_get_component_type_all
    module function process_get_component_i_terms &
         (process, i_component) result (i_term)
       integer, dimension(:), allocatable :: i_term
       class(process_t), intent(in) :: process
       integer, intent(in) :: i_component
    end function process_get_component_i_terms
    module function process_get_n_allowed_born (process, i_born) result (n_born)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_born
      integer :: n_born
    end function process_get_n_allowed_born
    module function process_get_pcm_ptr (process) result (pcm)
      class(pcm_t), pointer :: pcm
      class(process_t), intent(in), target :: process
    end function process_get_pcm_ptr
    module function process_component_can_be_integrated_single &
         (process, i_component) result (active)
      logical :: active
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
    end function process_component_can_be_integrated_single
    module function process_component_can_be_integrated_all &
         (process) result (val)
      logical, dimension(:), allocatable :: val
      class(process_t), intent(in) :: process
    end function process_component_can_be_integrated_all
    pure module subroutine process_reset_selected_cores (process)
      class(process_t), intent(inout) :: process
    end subroutine process_reset_selected_cores
    pure module subroutine process_select_components (process, indices)
      class(process_t), intent(inout) :: process
      integer, dimension(:), intent(in) :: indices
    end subroutine process_select_components
    pure module function process_component_is_selected &
         (process, index) result (val)
      logical :: val
      class(process_t), intent(in) :: process
      integer, intent(in) :: index
    end function process_component_is_selected
    pure module subroutine process_get_coupling_powers &
         (process, alpha_power, alphas_power)
      class(process_t), intent(in) :: process
      integer, intent(out) :: alpha_power, alphas_power
    end subroutine process_get_coupling_powers
    module function process_get_real_component (process) result (i_real)
      integer :: i_real
      class(process_t), intent(in) :: process
    end function process_get_real_component
    module function process_extract_active_component_mci &
         (process) result (i_active)
      integer :: i_active
      class(process_t), intent(in) :: process
    end function process_extract_active_component_mci
    module function process_uses_real_partition (process) result (val)
      logical :: val
      class(process_t), intent(in) :: process
    end function process_uses_real_partition
    module function process_get_md5sum_prc &
         (process, i_component) result (md5sum)
      character(32) :: md5sum
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
    end function process_get_md5sum_prc
    module function process_get_md5sum_mci (process, i_mci) result (md5sum)
      character(32) :: md5sum
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
    end function process_get_md5sum_mci
    module function process_get_md5sum_cfg (process) result (md5sum)
      character(32) :: md5sum
      class(process_t), intent(in) :: process
    end function process_get_md5sum_cfg
    module function process_get_n_cores (process) result (n)
      integer :: n
      class(process_t), intent(in) :: process
    end function process_get_n_cores
    module function process_get_base_i_term &
         (process, i_component) result (i_term)
      integer :: i_term
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
    end function process_get_base_i_term
    module function process_get_core_term (process, i_term) result (core)
      class(prc_core_t), pointer :: core
      class(process_t), intent(in), target :: process
      integer, intent(in) :: i_term
    end function process_get_core_term
    module function process_get_core_ptr (process, i_core) result (core)
      class(prc_core_t), pointer :: core
      class(process_t), intent(in), target :: process
      integer, intent(in) :: i_core
    end function process_get_core_ptr
    module function process_get_term_ptr (process, i) result (term)
      type(process_term_t), pointer :: term
      class(process_t), intent(in), target :: process
      integer, intent(in) :: i
    end function process_get_term_ptr
    module function process_get_i_term (process, i_core) result (i_term)
      integer :: i_term
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_core
    end function process_get_i_term
    module function process_get_i_core (process, i_term) result (i_core)
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_term
      integer :: i_core
    end function process_get_i_core
    module subroutine process_set_i_mci_work (process, i_mci)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
    end subroutine process_set_i_mci_work
    pure module function process_get_i_mci_work &
         (process, i_mci) result (i_mci_work)
      integer :: i_mci_work
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
    end function process_get_i_mci_work
    elemental module function process_get_i_sub (process, i_term) result (i_sub)
      integer :: i_sub
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_term
    end function process_get_i_sub
    elemental module function process_get_i_term_virtual &
         (process) result (i_term)
      integer :: i_term
      class(process_t), intent(in) :: process
    end function process_get_i_term_virtual
    elemental module function process_component_is_active_single &
         (process, i_comp) result (val)
      logical :: val
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_comp
    end function process_component_is_active_single
    pure module function process_component_is_active_all (process) result (val)
      logical, dimension(:), allocatable :: val
      class(process_t), intent(in) :: process
    end function process_component_is_active_all
    module function process_get_n_pass_default (process) result (n_pass)
      class(process_t), intent(in) :: process
      integer :: n_pass
    end function process_get_n_pass_default
    module function process_adapt_grids_default (process, pass) result (flag)
      class(process_t), intent(in) :: process
      integer, intent(in) :: pass
      logical :: flag
    end function process_adapt_grids_default
    module function process_adapt_weights_default (process, pass) result (flag)
      class(process_t), intent(in) :: process
      integer, intent(in) :: pass
      logical :: flag
    end function process_adapt_weights_default
    module function process_get_n_it_default (process, pass) result (n_it)
      class(process_t), intent(in) :: process
      integer, intent(in) :: pass
      integer :: n_it
    end function process_get_n_it_default
    module function process_get_n_calls_default (process, pass) result (n_calls)
      class(process_t), intent(in) :: process
      integer, intent(in) :: pass
      integer :: n_calls
    end function process_get_n_calls_default
    module subroutine process_set_run_id (process, run_id)
      class(process_t), intent(inout) :: process
      type(string_t), intent(in) :: run_id
    end subroutine process_set_run_id
    module function process_get_id (process) result (id)
      class(process_t), intent(in) :: process
      type(string_t) :: id
    end function process_get_id
    module function process_get_num_id (process) result (id)
      class(process_t), intent(in) :: process
      integer :: id
    end function process_get_num_id
    module function process_get_run_id (process) result (id)
      class(process_t), intent(in) :: process
      type(string_t) :: id
    end function process_get_run_id
    module function process_get_library_name (process) result (id)
      class(process_t), intent(in) :: process
      type(string_t) :: id
    end function process_get_library_name
    module function process_get_n_in (process) result (n)
      class(process_t), intent(in) :: process
      integer :: n
    end function process_get_n_in
    module function process_get_n_mci (process) result (n)
      class(process_t), intent(in) :: process
      integer :: n
    end function process_get_n_mci
    module function process_get_n_components (process) result (n)
      class(process_t), intent(in) :: process
      integer :: n
    end function process_get_n_components
    module function process_get_n_terms (process) result (n)
      class(process_t), intent(in) :: process
      integer :: n
    end function process_get_n_terms
    module subroutine process_get_i_component (process, i_mci, i_component)
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
      integer, dimension(:), intent(out), allocatable :: i_component
    end subroutine process_get_i_component
    module function process_get_component_id (process, i_component) result (id)
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
      type(string_t) :: id
    end function process_get_component_id
    module function process_get_component_def_ptr &
         (process, i_component) result (ptr)
      type(process_component_def_t), pointer :: ptr
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
    end function process_get_component_def_ptr
    module subroutine process_extract_core (process, i_term, core)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_term
      class(prc_core_t), intent(inout), allocatable :: core
    end subroutine process_extract_core
    module subroutine process_restore_core (process, i_term, core)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_term
      class(prc_core_t), intent(inout), allocatable :: core
    end subroutine process_restore_core
    module function process_get_constants (process, i_core) result (data)
      type(process_constants_t) :: data
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_core
    end function process_get_constants
    module function process_get_config (process) result (config)
      type(process_config_data_t) :: config
      class(process_t), intent(in) :: process
    end function process_get_config
    module function process_get_md5sum_constants (process, i_component, &
       type_string, nlo_type) result (this_md5sum)
      character(32) :: this_md5sum
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
      type(string_t), intent(in) :: type_string
      integer, intent(in) :: nlo_type
    end function process_get_md5sum_constants
    module subroutine process_get_term_flv_out (process, i_term, flv)
      class(process_t), intent(in), target :: process
      integer, intent(in) :: i_term
      type(flavor_t), dimension(:,:), allocatable, intent(out) :: flv
    end subroutine process_get_term_flv_out
    module function process_contains_unstable (process, model) result (flag)
      class(process_t), intent(in) :: process
      class(model_data_t), intent(in), target :: model
      logical :: flag
    end function process_contains_unstable
    module function process_get_sqrts (process) result (sqrts)
      class(process_t), intent(in) :: process
      real(default) :: sqrts
    end function process_get_sqrts
    module function process_get_energy (process) result (e)
      class(process_t), intent(in) :: process
      real(default), dimension(:), allocatable :: e
    end function process_get_energy
    module function process_get_polarization (process) result (pol)
      class(process_t), intent(in) :: process
      real(default), dimension(process%beam_config%data%n) :: pol
    end function process_get_polarization
    module function process_get_meta (process) result (meta)
      type(process_metadata_t) :: meta
      class(process_t), intent(in) :: process
    end function process_get_meta
    module function process_has_matrix_element &
         (process, i, is_term_index) result (active)
      logical :: active
      class(process_t), intent(in) :: process
      integer, intent(in), optional :: i
      logical, intent(in), optional :: is_term_index
    end function process_has_matrix_element
    module function process_get_beam_data_ptr (process) result (beam_data)
      class(process_t), intent(in), target :: process
      type(beam_data_t), pointer :: beam_data
    end function process_get_beam_data_ptr
    module function process_get_beam_config (process) result (beam_config)
      type(process_beam_config_t) :: beam_config
      class(process_t), intent(in) :: process
    end function process_get_beam_config
    module function process_get_beam_config_ptr (process) result (beam_config)
      type(process_beam_config_t), pointer :: beam_config
      class(process_t), intent(in), target :: process
    end function process_get_beam_config_ptr
    module function process_get_pdf_set (process) result (pdf_set)
      class(process_t), intent(in) :: process
      integer :: pdf_set
    end function process_get_pdf_set
    module function process_pcm_contains_pdfs (process) result (has_pdfs)
      logical :: has_pdfs
      class(process_t), intent(in) :: process
    end function process_pcm_contains_pdfs
    module function process_get_beam_file (process) result (file)
      class(process_t), intent(in) :: process
      type(string_t) :: file
    end function process_get_beam_file
    module function process_get_var_list_ptr (process) result (ptr)
      class(process_t), intent(in), target :: process
      type(var_list_t), pointer :: ptr
    end function process_get_var_list_ptr
    module function process_get_model_ptr (process) result (ptr)
      class(process_t), intent(in) :: process
      class(model_data_t), pointer :: ptr
    end function process_get_model_ptr
    module subroutine process_make_rng (process, rng)
      class(process_t), intent(inout) :: process
      class(rng_t), intent(out), allocatable :: rng
    end subroutine process_make_rng
    module function process_compute_amplitude (process, i_core, i, j, p, &
         f, h, c, fac_scale, ren_scale, alpha_qcd_forced) result (amp)
      class(process_t), intent(in), target :: process
      integer, intent(in) :: i_core
      integer, intent(in) :: i, j
      type(vector4_t), dimension(:), intent(in) :: p
      integer, intent(in) :: f, h, c
      real(default), intent(in), optional :: fac_scale, ren_scale
      real(default), intent(in), allocatable, optional :: alpha_qcd_forced
      complex(default) :: amp
    end function process_compute_amplitude
    module subroutine process_check_library_sanity (process)
      class(process_t), intent(in) :: process
    end subroutine process_check_library_sanity
    module subroutine process_reset_library_ptr (process)
      class(process_t), intent(inout) :: process
    end subroutine process_reset_library_ptr
    module subroutine process_set_counter_mci_entry (process, i_mci, counter)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: i_mci
      type(process_counter_t), intent(in) :: counter
    end subroutine process_set_counter_mci_entry
    module subroutine process_pacify (process, efficiency_reset, error_reset)
      class(process_t), intent(inout) :: process
      logical, intent(in), optional :: efficiency_reset, error_reset
    end subroutine process_pacify
    module subroutine test_allocate_sf_channels (process, n)
      class(process_t), intent(inout) :: process
      integer, intent(in) :: n
    end subroutine test_allocate_sf_channels
    module subroutine test_set_component_sf_channel (process, c)
      class(process_t), intent(inout) :: process
      integer, dimension(:), intent(in) :: c
    end subroutine test_set_component_sf_channel
    module subroutine test_get_mci_ptr (process, mci)
      class(process_t), intent(in), target :: process
      class(mci_t), intent(out), pointer :: mci
    end subroutine test_get_mci_ptr
    module subroutine process_init_mci_work (process, mci_work, i)
      class(process_t), intent(in), target :: process
      type(mci_work_t), intent(out) :: mci_work
      integer, intent(in) :: i
    end subroutine process_init_mci_work
    module function process_get_connected_states (process, i_component, &
         connected_terms) result (connected)
      type(connected_state_t), dimension(:), allocatable :: connected
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
      type(connected_state_t), dimension(:), intent(in) :: connected_terms
    end function process_get_connected_states
    module subroutine process_init_nlo_settings (process, var_list)
      class(process_t), intent(inout) :: process
      type(var_list_t), intent(in), target :: var_list
    end subroutine process_init_nlo_settings
    elemental module function process_get_nlo_type_component_single &
         (process, i_component) result (val)
      integer :: val
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
    end function process_get_nlo_type_component_single
    pure module function process_get_nlo_type_component_all &
         (process) result (val)
      integer, dimension(:), allocatable :: val
      class(process_t), intent(in) :: process
    end function process_get_nlo_type_component_all
    module function process_is_nlo_calculation (process) result (nlo)
      logical :: nlo
      class(process_t), intent(in) :: process
    end function process_is_nlo_calculation
    module function process_get_negative_sf (process) result (neg_sf)
      logical :: neg_sf
      class(process_t), intent(in) :: process
    end function process_get_negative_sf
    module function process_is_combined_nlo_integration &
         (process) result (combined)
      logical :: combined
      class(process_t), intent(in) :: process
    end function process_is_combined_nlo_integration
    pure module function process_component_is_real_finite &
         (process, i_component) result (val)
      logical :: val
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
    end function process_component_is_real_finite
    elemental module function process_get_component_nlo_type &
         (process, i_component) result (nlo_type)
      integer :: nlo_type
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
    end function process_get_component_nlo_type
    module function process_get_component_core_ptr &
         (process, i_component) result (core)
      class(process_t), intent(in), target :: process
      integer, intent(in) :: i_component
      class(prc_core_t), pointer :: core
    end function process_get_component_core_ptr
    module function process_get_component_associated_born &
         (process, i_component) result (i_born)
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_component
      integer :: i_born
    end function process_get_component_associated_born
    module function process_get_first_real_component (process) result (i_real)
      integer :: i_real
      class(process_t), intent(in) :: process
    end function process_get_first_real_component
    module function process_get_first_real_term (process) result (i_real)
      integer :: i_real
      class(process_t), intent(in) :: process
      integer :: i_component, i_term
    end function process_get_first_real_term
    elemental module function process_get_associated_real_fin &
         (process, i_component) result (i_real)
       integer :: i_real
       class(process_t), intent(in) :: process
       integer, intent(in) :: i_component
    end function process_get_associated_real_fin
    pure module function process_select_i_term (process, i_mci) result (i_term)
      integer :: i_term
      class(process_t), intent(in) :: process
      integer, intent(in) :: i_mci
    end function process_select_i_term
    module subroutine process_prepare_any_external_code (process)
      class(process_t), intent(inout), target :: process
    end subroutine process_prepare_any_external_code
  end interface

contains

  subroutine process_init &
       (process, proc_id, lib, os_data, model, var_list, beam_structure)
    class(process_t), intent(out) :: process
    type(string_t), intent(in) :: proc_id
    type(process_library_t), intent(in), target :: lib
    type(os_data_t), intent(in) :: os_data
    class(model_t), intent(in), target :: model
    type(var_list_t), intent(inout), target, optional :: var_list
    type(beam_structure_t), intent(in), optional :: beam_structure
    integer :: next_rng_seed
    if (debug_on) call msg_debug (D_PROCESS_INTEGRATION, "process_init")
    associate &
         (meta => process%meta, env => process%env, config => process%config)
      call env%init &
           (model, lib, os_data, var_list, beam_structure)
      call meta%init &
           (proc_id, lib, env%get_var_list_ptr ())
      call config%init &
           (meta, env)
      call dispatch_rng_factory &
           (process%rng_factory, env%get_var_list_ptr (), next_rng_seed)
      call update_rng_seed_in_var_list (var_list, next_rng_seed)
      call dispatch_pcm &
           (process%pcm, config%process_def%is_nlo ())
      associate (pcm => process%pcm)
        call pcm%init (env, meta)
        call pcm%allocate_components (process%component, meta)
        call pcm%categorize_components (config)
      end associate
    end associate
  end subroutine process_init

  subroutine dispatch_pcm (pcm, is_nlo)
    class(pcm_t), allocatable, intent(out) :: pcm
    logical, intent(in) :: is_nlo
    if (.not. is_nlo) then
       allocate (pcm_default_t :: pcm)
    else
       allocate (pcm_nlo_t :: pcm)
    end if
  end subroutine dispatch_pcm

  subroutine dispatch_test_me_core (core, core_def, model, &
       helicity_selection, qcd, use_color_factors, has_beam_pol)
    use prc_test_core, only: test_t
    class(prc_core_t), allocatable, intent(inout) :: core
    class(prc_core_def_t), intent(in) :: core_def
    class(model_data_t), intent(in), target, optional :: model
    type(helicity_selection_t), intent(in), optional :: helicity_selection
    type(qcd_t), intent(in), optional :: qcd
    logical, intent(in), optional :: use_color_factors
    logical, intent(in), optional :: has_beam_pol
    allocate (test_t :: core)
  end subroutine dispatch_test_me_core

  subroutine dispatch_template_core (core, core_def, model, &
       helicity_selection, qcd, use_color_factors, has_beam_pol)
    use prc_template_me, only: prc_template_me_t
    class(prc_core_t), allocatable, intent(inout) :: core
    class(prc_core_def_t), intent(in) :: core_def
    class(model_data_t), intent(in), target, optional :: model
    type(helicity_selection_t), intent(in), optional :: helicity_selection
    type(qcd_t), intent(in), optional :: qcd
    logical, intent(in), optional :: use_color_factors
    logical, intent(in), optional :: has_beam_pol
    allocate (prc_template_me_t :: core)
    select type (core)
    type is (prc_template_me_t)
       call core%set_parameters (model)
    end select
  end subroutine dispatch_template_core

  subroutine process_setup_test_cores (process, type_string)
    class(process_t), intent(inout) :: process
    class(prc_core_t), allocatable :: core
    type(string_t), intent(in), optional :: type_string
    if (present (type_string)) then
       select case (char (type_string))
       case ("template")
          call process%setup_cores (dispatch_template_core)
       case ("test_me")
          call process%setup_cores (dispatch_test_me_core)
       case default
          call msg_bug ("process setup test cores: unsupported type string")
       end select
    else
       call process%setup_cores (dispatch_test_me_core)
    end if
  end subroutine process_setup_test_cores


end module process
