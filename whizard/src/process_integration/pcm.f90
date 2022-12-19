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

module pcm

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use model_data, only: model_data_t
  use models, only: model_t
  use quantum_numbers, only: quantum_numbers_t, quantum_numbers_mask_t
  use variables, only: var_list_t
  use nlo_data, only: nlo_settings_t
  use nlo_data, only: fks_template_t
  use nlo_data, only: FKS_DEFAULT, FKS_RESONANCES
  use mci_base, only: mci_t
  use phs_base, only: phs_config_t
  use mappings, only: mapping_defaults_t
  use phs_forests, only: phs_parameters_t
  use phs_fks, only: isr_kinematics_t, real_kinematics_t
  use phs_fks, only: phs_identifier_t
  use fks_regions, only: region_data_t
  use phs_fks, only: phs_fks_generator_t
  use phs_fks, only: dalitz_plot_t
  use phs_fks, only: phs_fks_config_t, get_filtered_resonance_histories
  use dispatch_phase_space, only: dispatch_phs
  use real_subtraction, only: real_subtraction_t, soft_mismatch_t
  use real_subtraction, only: INTEGRATION, FIXED_ORDER_EVENTS
  use real_subtraction, only: real_partition_t, powheg_damping_simple_t
  use real_subtraction, only: real_partition_fixed_order_t
  use virtual, only: virtual_t
  use dglap_remnant, only: dglap_remnant_t
  use blha_config, only: blha_master_t

  use pcm_base
  use process_config
  use process_mci, only: process_mci_entry_t
  use process_mci, only: REAL_SINGULAR, REAL_FINITE

  implicit none
  private

  public :: pcm_default_t
  public :: pcm_nlo_t
  public :: pcm_nlo_workspace_t

  type, extends (pcm_t) :: pcm_default_t
   contains
     procedure :: allocate_workspace => pcm_default_allocate_workspace
     procedure :: final => pcm_default_final
     procedure :: is_nlo => pcm_default_is_nlo
     procedure :: init => pcm_default_init
     procedure :: categorize_components => pcm_default_categorize_components
     procedure :: init_phs_config => pcm_default_init_phs_config
     procedure :: allocate_cores => pcm_default_allocate_cores
     procedure :: prepare_any_external_code => &
          pcm_default_prepare_any_external_code
     procedure :: setup_blha => pcm_default_setup_blha
     procedure :: prepare_blha_core => pcm_default_prepare_blha_core
     procedure :: set_blha_methods => pcm_default_set_blha_methods
     procedure :: get_blha_flv_states => pcm_default_get_blha_flv_states
     procedure :: setup_mci => pcm_default_setup_mci
     procedure :: call_dispatch_mci => pcm_default_call_dispatch_mci
     procedure :: complete_setup => pcm_default_complete_setup
     procedure :: init_component => pcm_default_init_component
  end type pcm_default_t

  type, extends (pcm_workspace_t) :: pcm_default_workspace_t
  contains
    procedure :: final => pcm_default_workspace_final
    procedure :: is_nlo => pcm_default_workspace_is_nlo
  end type pcm_default_workspace_t

  type, extends (pcm_t) :: pcm_nlo_t
     type(string_t) :: id
     logical :: combined_integration = .false.
     logical :: vis_fks_regions = .false.
     integer, dimension(:), allocatable :: nlo_type
     integer, dimension(:), allocatable :: nlo_type_core
     integer, dimension(:), allocatable :: component_type
     integer :: i_born = 0
     integer :: i_real = 0
     integer :: i_sub = 0
     type(nlo_settings_t) :: settings
     type(region_data_t) :: region_data
     logical :: use_real_partition = .false.
     logical :: use_real_singular = .false.
     real(default) :: real_partition_scale = 0
     class(real_partition_t), allocatable :: real_partition
     type(dalitz_plot_t) :: dalitz_plot
     type(quantum_numbers_t), dimension(:,:), allocatable :: qn_real, qn_born
  contains
    procedure :: init => pcm_nlo_init
    procedure :: init_nlo_settings => pcm_nlo_init_nlo_settings
    procedure :: categorize_components => pcm_nlo_categorize_components
    procedure :: init_phs_config => pcm_nlo_init_phs_config
    procedure :: allocate_cores => pcm_nlo_allocate_cores
    procedure :: prepare_any_external_code => &
         pcm_nlo_prepare_any_external_code
    procedure :: setup_blha => pcm_nlo_setup_blha
    procedure :: complete_setup => pcm_nlo_complete_setup
    procedure :: prepare_blha_core => pcm_nlo_prepare_blha_core
    procedure :: set_blha_methods => pcm_nlo_set_blha_methods
    procedure :: get_blha_flv_states => pcm_nlo_get_blha_flv_states
    procedure :: setup_mci => pcm_nlo_setup_mci
    procedure :: call_dispatch_mci => pcm_nlo_call_dispatch_mci
    procedure :: handle_threshold_core => pcm_nlo_handle_threshold_core
    procedure :: setup_region_data => pcm_nlo_setup_region_data
    procedure :: setup_real_partition => pcm_nlo_setup_real_partition
    procedure :: init_component => pcm_nlo_init_component
    procedure :: record_inactive_components => pcm_nlo_record_inactive_components
    procedure :: core_is_radiation => pcm_nlo_core_is_radiation
    procedure :: get_n_flv_born => pcm_nlo_get_n_flv_born
    procedure :: get_n_flv_real => pcm_nlo_get_n_flv_real
    procedure :: get_n_alr => pcm_nlo_get_n_alr
    procedure :: get_flv_states => pcm_nlo_get_flv_states
    procedure :: get_qn => pcm_nlo_get_qn
    procedure :: has_massive_emitter => pcm_nlo_has_massive_emitter
    procedure :: get_mass_info => pcm_nlo_get_mass_info
    procedure :: allocate_workspace => pcm_nlo_allocate_workspace
    procedure :: init_qn => pcm_nlo_init_qn
    procedure :: allocate_ps_matching => pcm_nlo_allocate_ps_matching
    procedure :: activate_dalitz_plot => pcm_nlo_activate_dalitz_plot
    procedure :: register_dalitz_plot => pcm_nlo_register_dalitz_plot
    procedure :: setup_phs_generator => pcm_nlo_setup_phs_generator
    procedure :: final => pcm_nlo_final
    procedure :: is_nlo => pcm_nlo_is_nlo
  end type pcm_nlo_t

  type, extends (pcm_workspace_t) :: pcm_nlo_workspace_t
     type(real_kinematics_t), pointer :: real_kinematics => null ()
     type(isr_kinematics_t), pointer :: isr_kinematics => null ()
     type(real_subtraction_t) :: real_sub
     type(virtual_t) :: virtual
     type(soft_mismatch_t) :: soft_mismatch
     type(dglap_remnant_t) :: dglap_remnant
     integer, dimension(:), allocatable :: i_mci_to_real_component
  contains
    procedure :: set_radiation_event => pcm_nlo_workspace_set_radiation_event
    procedure :: set_subtraction_event => pcm_nlo_workspace_set_subtraction_event
    procedure :: disable_subtraction => pcm_nlo_workspace_disable_subtraction
    procedure :: init_config => pcm_nlo_workspace_init_config
    procedure :: setup_real_component => pcm_nlo_workspace_setup_real_component
    procedure :: init_real_and_isr_kinematics => &
         pcm_nlo_workspace_init_real_and_isr_kinematics
    procedure :: set_real_and_isr_kinematics => &
        pcm_nlo_workspace_set_real_and_isr_kinematics
    procedure :: init_real_subtraction => pcm_nlo_workspace_init_real_subtraction
    procedure :: set_momenta_and_scales_virtual => &
       pcm_nlo_workspace_set_momenta_and_scales_virtual
    procedure :: set_fac_scale => pcm_nlo_workspace_set_fac_scale
    procedure :: set_momenta => pcm_nlo_workspace_set_momenta
    procedure :: get_momenta => pcm_nlo_workspace_get_momenta
    procedure :: get_xi_max => pcm_nlo_workspace_get_xi_max
    procedure :: set_x_rad => pcm_nlo_workspace_set_x_rad
    procedure :: init_virtual => pcm_nlo_workspace_init_virtual
    procedure :: disable_virtual_subtraction => &
         pcm_nlo_workspace_disable_virtual_subtraction
    procedure :: compute_sqme_virt => pcm_nlo_workspace_compute_sqme_virt
    procedure :: compute_sqme_mismatch => pcm_nlo_workspace_compute_sqme_mismatch
    procedure :: compute_sqme_dglap_remnant => &
         pcm_nlo_workspace_compute_sqme_dglap_remnant
    procedure :: set_fixed_order_event_mode => &
         pcm_nlo_workspace_set_fixed_order_event_mode
    procedure :: init_soft_mismatch => pcm_nlo_workspace_init_soft_mismatch
    procedure :: init_dglap_remnant => pcm_nlo_workspace_init_dglap_remnant
    procedure :: is_fixed_order_nlo_events &
         => pcm_nlo_workspace_is_fixed_order_nlo_events
    procedure :: final => pcm_nlo_workspace_final
    procedure :: is_nlo => pcm_nlo_workspace_is_nlo
  procedure :: powheg_kinematic_factors_real => &
       pcm_nlo_workspace_powheg_kinematic_factors_real
  end type pcm_nlo_workspace_t


  interface
    module subroutine pcm_default_final (pcm)
      class(pcm_default_t), intent(inout) :: pcm
    end subroutine pcm_default_final
    module function pcm_default_is_nlo (pcm) result (is_nlo)
      logical :: is_nlo
      class(pcm_default_t), intent(in) :: pcm
    end function pcm_default_is_nlo
    module subroutine pcm_default_init (pcm, env, meta)
      class(pcm_default_t), intent(out) :: pcm
      type(process_environment_t), intent(in) :: env
      type(process_metadata_t), intent(in) :: meta
    end subroutine pcm_default_init
    module subroutine pcm_default_workspace_final (pcm_work)
      class(pcm_default_workspace_t), intent(inout) :: pcm_work
    end subroutine pcm_default_workspace_final
    module function pcm_default_workspace_is_nlo (pcm_work) result (is_nlo)
      logical :: is_nlo
      class(pcm_default_workspace_t), intent(inout) :: pcm_work
    end function pcm_default_workspace_is_nlo
    module subroutine pcm_default_categorize_components (pcm, config)
      class(pcm_default_t), intent(inout) :: pcm
      type(process_config_data_t), intent(in) :: config
    end subroutine pcm_default_categorize_components
    module subroutine pcm_default_init_phs_config &
         (pcm, phs_entry, meta, env, phs_par, mapping_defs)
      class(pcm_default_t), intent(inout) :: pcm
      type(process_phs_config_t), &
           dimension(:), allocatable, intent(out) :: phs_entry
      type(process_metadata_t), intent(in) :: meta
      type(process_environment_t), intent(in) :: env
      type(mapping_defaults_t), intent(in) :: mapping_defs
      type(phs_parameters_t), intent(in) :: phs_par
    end subroutine pcm_default_init_phs_config
    module subroutine pcm_default_allocate_cores (pcm, config, core_entry)
      class(pcm_default_t), intent(inout) :: pcm
      type(process_config_data_t), intent(in) :: config
      type(core_entry_t), dimension(:), allocatable, intent(out) :: core_entry
    end subroutine pcm_default_allocate_cores
    module subroutine pcm_default_prepare_any_external_code &
         (pcm, core_entry, i_core, libname, model, var_list)
      class(pcm_default_t), intent(in) :: pcm
      type(core_entry_t), intent(inout) :: core_entry
      integer, intent(in) :: i_core
      type(string_t), intent(in) :: libname
      type(model_data_t), intent(in), target :: model
      type(var_list_t), intent(in) :: var_list
    end subroutine pcm_default_prepare_any_external_code
    module subroutine pcm_default_setup_blha (pcm, core_entry)
      class(pcm_default_t), intent(in) :: pcm
      type(core_entry_t), intent(inout) :: core_entry
    end subroutine pcm_default_setup_blha
    module subroutine pcm_default_prepare_blha_core (pcm, core_entry, model)
      class(pcm_default_t), intent(in) :: pcm
      type(core_entry_t), intent(inout) :: core_entry
      class(model_data_t), intent(in), target :: model
    end subroutine pcm_default_prepare_blha_core
    module subroutine pcm_default_set_blha_methods (pcm, blha_master, var_list)
      class(pcm_default_t), intent(inout) :: pcm
      type(blha_master_t), intent(inout) :: blha_master
      type(var_list_t), intent(in) :: var_list
    end subroutine pcm_default_set_blha_methods
    module subroutine pcm_default_get_blha_flv_states &
         (pcm, core_entry, flv_born, flv_real)
      class(pcm_default_t), intent(in) :: pcm
      type(core_entry_t), dimension(:), intent(in) :: core_entry
      integer, dimension(:,:), allocatable, intent(out) :: flv_born
      integer, dimension(:,:), allocatable, intent(out) :: flv_real
    end subroutine pcm_default_get_blha_flv_states
    module subroutine pcm_default_setup_mci (pcm, mci_entry)
      class(pcm_default_t), intent(inout) :: pcm
      type(process_mci_entry_t), &
           dimension(:), allocatable, intent(out) :: mci_entry
    end subroutine pcm_default_setup_mci
    module subroutine pcm_default_call_dispatch_mci (pcm, &
         dispatch_mci, var_list, process_id, mci_template)
      class(pcm_default_t), intent(inout) :: pcm
      procedure(dispatch_mci_proc) :: dispatch_mci
      type(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: process_id
      class(mci_t), allocatable, intent(out) :: mci_template
    end subroutine pcm_default_call_dispatch_mci
    module subroutine pcm_default_complete_setup &
         (pcm, core_entry, component, model)
      class(pcm_default_t), intent(inout) :: pcm
      type(core_entry_t), dimension(:), intent(in) :: core_entry
      type(process_component_t), dimension(:), intent(inout) :: component
      type(model_t), intent(in), target :: model
    end subroutine pcm_default_complete_setup
    module subroutine pcm_default_init_component (pcm, component, i, active, &
         phs_config, env, meta, config)
      class(pcm_default_t), intent(in) :: pcm
      type(process_component_t), intent(out) :: component
      integer, intent(in) :: i
      logical, intent(in) :: active
      class(phs_config_t), allocatable, intent(in) :: phs_config
      type(process_environment_t), intent(in) :: env
      type(process_metadata_t), intent(in) :: meta
      type(process_config_data_t), intent(in) :: config
    end subroutine pcm_default_init_component
    module subroutine pcm_nlo_init (pcm, env, meta)
      class(pcm_nlo_t), intent(out) :: pcm
      type(process_metadata_t), intent(in) :: meta
      type(process_environment_t), intent(in) :: env
    end subroutine pcm_nlo_init
    module subroutine pcm_nlo_init_nlo_settings (pcm, var_list)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(var_list_t), intent(in), target :: var_list
    end subroutine pcm_nlo_init_nlo_settings
    module subroutine pcm_nlo_categorize_components (pcm, config)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(process_config_data_t), intent(in) :: config
    end subroutine pcm_nlo_categorize_components
    module subroutine pcm_nlo_init_phs_config &
         (pcm, phs_entry, meta, env, phs_par, mapping_defs)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(process_phs_config_t), &
           dimension(:), allocatable, intent(out) :: phs_entry
      type(process_metadata_t), intent(in) :: meta
      type(process_environment_t), intent(in) :: env
      type(mapping_defaults_t), intent(in) :: mapping_defs
      type(phs_parameters_t), intent(in) :: phs_par
    end subroutine pcm_nlo_init_phs_config
    module subroutine pcm_nlo_allocate_cores (pcm, config, core_entry)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(process_config_data_t), intent(in) :: config
      type(core_entry_t), dimension(:), allocatable, intent(out) :: core_entry
    end subroutine pcm_nlo_allocate_cores
    module subroutine pcm_nlo_prepare_any_external_code &
         (pcm, core_entry, i_core, libname, model, var_list)
      class(pcm_nlo_t), intent(in) :: pcm
      type(core_entry_t), intent(inout) :: core_entry
      integer, intent(in) :: i_core
      type(string_t), intent(in) :: libname
      type(model_data_t), intent(in), target :: model
      type(var_list_t), intent(in) :: var_list
    end subroutine pcm_nlo_prepare_any_external_code
    module subroutine pcm_nlo_setup_blha (pcm, core_entry)
      class(pcm_nlo_t), intent(in) :: pcm
      type(core_entry_t), intent(inout) :: core_entry
    end subroutine pcm_nlo_setup_blha
    module subroutine pcm_nlo_complete_setup (pcm, core_entry, component, model)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(core_entry_t), dimension(:), intent(in) :: core_entry
      type(process_component_t), dimension(:), intent(inout) :: component
      type(model_t), intent(in), target :: model
    end subroutine pcm_nlo_complete_setup
    module subroutine pcm_nlo_prepare_blha_core (pcm, core_entry, model)
      class(pcm_nlo_t), intent(in) :: pcm
      type(core_entry_t), intent(inout) :: core_entry
      class(model_data_t), intent(in), target :: model
    end subroutine pcm_nlo_prepare_blha_core
    module subroutine pcm_nlo_set_blha_methods (pcm, blha_master, var_list)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(blha_master_t), intent(inout) :: blha_master
      type(var_list_t), intent(in) :: var_list
    end subroutine pcm_nlo_set_blha_methods
    module subroutine pcm_nlo_get_blha_flv_states &
         (pcm, core_entry, flv_born, flv_real)
      class(pcm_nlo_t), intent(in) :: pcm
      type(core_entry_t), dimension(:), intent(in) :: core_entry
      integer, dimension(:,:), allocatable, intent(out) :: flv_born
      integer, dimension(:,:), allocatable, intent(out) :: flv_real
    end subroutine pcm_nlo_get_blha_flv_states
    module subroutine pcm_nlo_setup_mci (pcm, mci_entry)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(process_mci_entry_t), &
           dimension(:), allocatable, intent(out) :: mci_entry
    end subroutine pcm_nlo_setup_mci
    module subroutine pcm_nlo_call_dispatch_mci (pcm, &
         dispatch_mci, var_list, process_id, mci_template)
      class(pcm_nlo_t), intent(inout) :: pcm
      procedure(dispatch_mci_proc) :: dispatch_mci
      type(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: process_id
      class(mci_t), allocatable, intent(out) :: mci_template
    end subroutine pcm_nlo_call_dispatch_mci
    module subroutine pcm_nlo_handle_threshold_core (pcm, core_entry)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(core_entry_t), dimension(:), intent(in) :: core_entry
    end subroutine pcm_nlo_handle_threshold_core
    module subroutine pcm_nlo_setup_region_data &
         (pcm, core_entry, model, alpha_power, alphas_power, phs_config)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(core_entry_t), dimension(:), intent(in) :: core_entry
      type(model_t), intent(in), target :: model
      integer, intent(in) :: alpha_power, alphas_power
      class(phs_config_t), intent(inout), optional :: phs_config
    end subroutine pcm_nlo_setup_region_data
    module subroutine pcm_nlo_init_component (pcm, component, i, active, &
         phs_config, env, meta, config)
      class(pcm_nlo_t), intent(in) :: pcm
      type(process_component_t), intent(out) :: component
      integer, intent(in) :: i
      logical, intent(in) :: active
      class(phs_config_t), allocatable, intent(in) :: phs_config
      type(process_environment_t), intent(in) :: env
      type(process_metadata_t), intent(in) :: meta
      type(process_config_data_t), intent(in) :: config
    end subroutine pcm_nlo_init_component
    module subroutine pcm_nlo_record_inactive_components (pcm, component, meta)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(process_component_t), dimension(:), intent(in) :: component
      type(process_metadata_t), intent(inout) :: meta
    end subroutine pcm_nlo_record_inactive_components
    module function pcm_nlo_core_is_radiation (pcm, i_core) result (is_rad)
      logical :: is_rad
      class(pcm_nlo_t), intent(in) :: pcm
      integer, intent(in) :: i_core
    end function pcm_nlo_core_is_radiation
    module function pcm_nlo_get_n_flv_born (pcm_nlo) result (n_flv)
      integer :: n_flv
      class(pcm_nlo_t), intent(in) :: pcm_nlo
    end function pcm_nlo_get_n_flv_born
    module function pcm_nlo_get_n_flv_real (pcm_nlo) result (n_flv)
      integer :: n_flv
      class(pcm_nlo_t), intent(in) :: pcm_nlo
    end function pcm_nlo_get_n_flv_real
    module function pcm_nlo_get_n_alr (pcm) result (n_alr)
      integer :: n_alr
      class(pcm_nlo_t), intent(in) :: pcm
    end function pcm_nlo_get_n_alr
    module function pcm_nlo_get_flv_states (pcm, born) result (flv)
      integer, dimension(:,:), allocatable :: flv
      class(pcm_nlo_t), intent(in) :: pcm
      logical, intent(in) :: born
    end function pcm_nlo_get_flv_states
    module function pcm_nlo_get_qn (pcm, born) result (qn)
      type(quantum_numbers_t), dimension(:,:), allocatable :: qn
      class(pcm_nlo_t), intent(in) :: pcm
      logical, intent(in) :: born
    end function pcm_nlo_get_qn
    module function pcm_nlo_has_massive_emitter (pcm) result (val)
      logical :: val
      class(pcm_nlo_t), intent(in) :: pcm
    end function pcm_nlo_has_massive_emitter
    module function pcm_nlo_get_mass_info (pcm, i_flv) result (massive)
      class(pcm_nlo_t), intent(in) :: pcm
      integer, intent(in) :: i_flv
      logical, dimension(:), allocatable :: massive
    end function pcm_nlo_get_mass_info
    module subroutine pcm_nlo_init_qn (pcm, model)
      class(pcm_nlo_t), intent(inout) :: pcm
      class(model_data_t), intent(in) :: model
    end subroutine pcm_nlo_init_qn
    module subroutine pcm_nlo_activate_dalitz_plot (pcm, filename)
      class(pcm_nlo_t), intent(inout) :: pcm
      type(string_t), intent(in) :: filename
    end subroutine pcm_nlo_activate_dalitz_plot
    module subroutine pcm_nlo_register_dalitz_plot (pcm, emitter, p)
      class(pcm_nlo_t), intent(inout) :: pcm
      integer, intent(in) :: emitter
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine pcm_nlo_register_dalitz_plot
    module subroutine pcm_nlo_setup_phs_generator (pcm, pcm_work, generator, &
       sqrts, mode, singular_jacobian)
      class(pcm_nlo_t), intent(in) :: pcm
      type(phs_fks_generator_t), intent(inout) :: generator
      type(pcm_nlo_workspace_t), intent(in), target :: pcm_work
      real(default), intent(in) :: sqrts
      integer, intent(in), optional:: mode
      logical, intent(in), optional :: singular_jacobian
    end subroutine pcm_nlo_setup_phs_generator
    module subroutine pcm_nlo_final (pcm)
      class(pcm_nlo_t), intent(inout) :: pcm
    end subroutine pcm_nlo_final
    module function pcm_nlo_is_nlo (pcm) result (is_nlo)
      logical :: is_nlo
      class(pcm_nlo_t), intent(in) :: pcm
    end function pcm_nlo_is_nlo
    module subroutine pcm_nlo_workspace_set_radiation_event (pcm_work)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
    end subroutine pcm_nlo_workspace_set_radiation_event
    module subroutine pcm_nlo_workspace_set_subtraction_event (pcm_work)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
    end subroutine pcm_nlo_workspace_set_subtraction_event
    module subroutine pcm_nlo_workspace_disable_subtraction (pcm_work)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
    end subroutine pcm_nlo_workspace_disable_subtraction
    module subroutine pcm_nlo_workspace_init_config (pcm_work, pcm, &
         active_components, nlo_types, energy, i_real_fin, model)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      class(pcm_t), intent(in) :: pcm
      logical, intent(in), dimension(:) :: active_components
      integer, intent(in), dimension(:) :: nlo_types
      real(default), intent(in), dimension(:) :: energy
      integer, intent(in) :: i_real_fin
      class(model_data_t), intent(in) :: model
    end subroutine pcm_nlo_workspace_init_config
    module subroutine pcm_nlo_workspace_setup_real_component (pcm_work, pcm, &
       subtraction_disabled)
      class(pcm_nlo_workspace_t), intent(inout), target :: pcm_work
      class(pcm_t), intent(in) :: pcm
      logical, intent(in) :: subtraction_disabled
    end subroutine pcm_nlo_workspace_setup_real_component
    module subroutine pcm_nlo_workspace_init_real_and_isr_kinematics &
         (pcm_work, pcm, energy)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      class(pcm_t), intent(in) :: pcm
      real(default), dimension(:), intent(in) :: energy
    end subroutine pcm_nlo_workspace_init_real_and_isr_kinematics
    module subroutine pcm_nlo_workspace_set_real_and_isr_kinematics &
         (pcm_work, phs_identifiers, sqrts)
      class(pcm_nlo_workspace_t), intent(inout), target :: pcm_work
      type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
      real(default), intent(in) :: sqrts
    end subroutine pcm_nlo_workspace_set_real_and_isr_kinematics
    module subroutine pcm_nlo_workspace_init_real_subtraction (pcm_work, pcm)
      class(pcm_nlo_workspace_t), intent(inout), target :: pcm_work
      class(pcm_t), intent(in) :: pcm
    end subroutine pcm_nlo_workspace_init_real_subtraction
    module subroutine pcm_nlo_workspace_set_momenta_and_scales_virtual &
         (pcm_work, p, ren_scale, fac_scale, es_scale)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), allocatable, intent(in) :: ren_scale
      real(default), intent(in) :: fac_scale
      real(default), allocatable, intent(in) :: es_scale
    end subroutine pcm_nlo_workspace_set_momenta_and_scales_virtual
    module subroutine pcm_nlo_workspace_set_fac_scale (pcm_work, fac_scale)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      real(default), intent(in) :: fac_scale
    end subroutine pcm_nlo_workspace_set_fac_scale
    module subroutine pcm_nlo_workspace_set_momenta (pcm_work, &
         p_born, p_real, i_phs, cms)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      type(vector4_t), dimension(:), intent(in) :: p_born, p_real
      integer, intent(in) :: i_phs
      logical, intent(in), optional :: cms
    end subroutine pcm_nlo_workspace_set_momenta
    module function pcm_nlo_workspace_get_momenta (pcm_work, pcm, &
         i_phs, born_phsp, cms) result (p)
      type(vector4_t), dimension(:), allocatable :: p
      class(pcm_nlo_workspace_t), intent(in) :: pcm_work
      class(pcm_t), intent(in) :: pcm
      integer, intent(in) :: i_phs
      logical, intent(in) :: born_phsp
      logical, intent(in), optional :: cms
    end function pcm_nlo_workspace_get_momenta
    module function pcm_nlo_workspace_get_xi_max (pcm_work, alr) result (xi_max)
      real(default) :: xi_max
      class(pcm_nlo_workspace_t), intent(in) :: pcm_work
      integer, intent(in) :: alr
    end function pcm_nlo_workspace_get_xi_max
    module subroutine pcm_nlo_workspace_set_x_rad (pcm_work, x_tot)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      real(default), intent(in), dimension(:) :: x_tot
    end subroutine pcm_nlo_workspace_set_x_rad
    module subroutine pcm_nlo_workspace_init_virtual (pcm_work, pcm, model)
      class(pcm_nlo_workspace_t), intent(inout), target :: pcm_work
      class(pcm_t), intent(in) :: pcm
      class(model_data_t), intent(in) :: model
    end subroutine pcm_nlo_workspace_init_virtual
    module subroutine pcm_nlo_workspace_disable_virtual_subtraction (pcm_work)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
    end subroutine pcm_nlo_workspace_disable_virtual_subtraction
    module subroutine pcm_nlo_workspace_compute_sqme_virt (pcm_work, pcm, p, &
         alpha_coupling, separate_uborns, sqme_virt)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      class(pcm_t), intent(in) :: pcm
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), dimension(2), intent(in) :: alpha_coupling
      logical, intent(in) :: separate_uborns
      real(default), dimension(:), allocatable, intent(inout) :: sqme_virt
    end subroutine pcm_nlo_workspace_compute_sqme_virt
    module subroutine pcm_nlo_workspace_compute_sqme_mismatch (pcm_work, pcm, &
         alpha_s, separate_uborns, sqme_mism)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      class(pcm_t), intent(in) :: pcm
      real(default), intent(in) :: alpha_s
      logical, intent(in) :: separate_uborns
      real(default), dimension(:), allocatable, intent(inout) :: sqme_mism
    end subroutine pcm_nlo_workspace_compute_sqme_mismatch
    module subroutine pcm_nlo_workspace_compute_sqme_dglap_remnant (pcm_work, &
         pcm, alpha_coupling, separate_uborns, sqme_dglap)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      class(pcm_t), intent(in) :: pcm
      real(default), dimension(2), intent(in) :: alpha_coupling
      logical, intent(in) :: separate_uborns
      real(default), dimension(:), allocatable, intent(inout) :: sqme_dglap
    end subroutine pcm_nlo_workspace_compute_sqme_dglap_remnant
    module subroutine pcm_nlo_workspace_set_fixed_order_event_mode (pcm_work)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
    end subroutine pcm_nlo_workspace_set_fixed_order_event_mode
    module subroutine pcm_nlo_workspace_init_soft_mismatch (pcm_work, pcm)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      class(pcm_t), intent(in) :: pcm
    end subroutine pcm_nlo_workspace_init_soft_mismatch
    module subroutine pcm_nlo_workspace_init_dglap_remnant (pcm_work, pcm)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      class(pcm_t), intent(in) :: pcm
    end subroutine pcm_nlo_workspace_init_dglap_remnant
    module function pcm_nlo_workspace_is_fixed_order_nlo_events &
         (pcm_work) result (is_fnlo)
      logical :: is_fnlo
      class(pcm_nlo_workspace_t), intent(in) :: pcm_work
    end function pcm_nlo_workspace_is_fixed_order_nlo_events
    module subroutine pcm_nlo_workspace_final (pcm_work)
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
    end subroutine pcm_nlo_workspace_final
    module function pcm_nlo_workspace_is_nlo (pcm_work) result (is_nlo)
      logical :: is_nlo
      class(pcm_nlo_workspace_t), intent(inout) :: pcm_work
    end function pcm_nlo_workspace_is_nlo
    module function pcm_nlo_workspace_powheg_kinematic_factors_real &
         (pcm_work, sqme_real, alr) result (sqme_real_corr)
      real(default) :: sqme_real_corr
      class(pcm_nlo_workspace_t), intent(in) :: pcm_work
      real(default), intent(in) :: sqme_real
      integer, intent(in) :: alr
    end function pcm_nlo_workspace_powheg_kinematic_factors_real
  end interface

contains

  subroutine pcm_default_allocate_workspace (pcm, work)
    class(pcm_default_t), intent(in) :: pcm
    class(pcm_workspace_t), intent(inout), allocatable :: work
    allocate (pcm_default_workspace_t :: work)
  end subroutine pcm_default_allocate_workspace

  subroutine pcm_nlo_setup_real_partition (pcm)
    class(pcm_nlo_t), intent(inout) :: pcm
    if (pcm%use_real_partition) then
       if (.not. allocated (pcm%real_partition)) then
          allocate (real_partition_fixed_order_t :: pcm%real_partition)
          select type (partition => pcm%real_partition)
          type is (real_partition_fixed_order_t)
             call pcm%region_data%get_all_ftuples (partition%fks_pairs)
             partition%scale = pcm%real_partition_scale
          end select
       end if
    end if
  end subroutine pcm_nlo_setup_real_partition

  subroutine pcm_nlo_allocate_workspace (pcm, work)
    class(pcm_nlo_t), intent(in) :: pcm
    class(pcm_workspace_t), intent(inout), allocatable :: work
    allocate (pcm_nlo_workspace_t :: work)
  end subroutine pcm_nlo_allocate_workspace

  subroutine pcm_nlo_allocate_ps_matching (pcm)
    class(pcm_nlo_t), intent(inout) :: pcm
    if (.not. allocated (pcm%real_partition)) then
       allocate (powheg_damping_simple_t :: pcm%real_partition)
    end if
  end subroutine pcm_nlo_allocate_ps_matching


end module pcm
