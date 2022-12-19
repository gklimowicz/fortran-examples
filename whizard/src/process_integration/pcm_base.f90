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

module pcm_base

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface, only: os_data_t

  use process_libraries, only: process_library_t

  use prc_core_def
  use prc_core

  use variables, only: var_list_t
  use mappings, only: mapping_defaults_t
  use phs_base, only: phs_config_t
  use phs_forests, only: phs_parameters_t
  use mci_base, only: mci_t
  use model_data, only: model_data_t
  use models, only: model_t

  use blha_config, only: blha_master_t
  use blha_olp_interfaces, only: blha_template_t
  use process_config
  use process_mci, only: process_mci_entry_t

  implicit none
  private

  public :: core_entry_t
  public :: pcm_t
  public :: dispatch_mci_proc
  public :: pcm_workspace_t



  type :: core_entry_t
     integer :: i_component = 0
     logical :: active = .false.
     class(prc_core_def_t), pointer :: core_def => null ()
     type(blha_template_t), allocatable :: blha_config
     class(prc_core_t), allocatable :: core
   contains
     procedure :: get_core_ptr => core_entry_get_core_ptr
     procedure :: configure => core_entry_configure
  end type core_entry_t

  type, abstract :: pcm_t
     logical :: initialized = .false.
     logical :: has_pdfs = .false.
     integer :: n_components = 0
     integer :: n_cores = 0
     integer :: n_mci = 0
     logical, dimension(:), allocatable :: component_selected
     logical, dimension(:), allocatable :: component_active
     integer, dimension(:), allocatable :: i_phs_config
     integer, dimension(:), allocatable :: i_core
     integer, dimension(:), allocatable :: i_mci
     type(blha_template_t) :: blha_defaults
     logical :: uses_blha = .false.
     type(os_data_t) :: os_data
  contains
    procedure(pcm_allocate_workspace), deferred :: allocate_workspace
    procedure(pcm_is_nlo), deferred :: is_nlo
    procedure(pcm_final), deferred :: final
    procedure(pcm_init), deferred :: init
    procedure :: set_blha_defaults => pcm_set_blha_defaults
    procedure(pcm_set_blha_methods), deferred :: set_blha_methods
    procedure(pcm_get_blha_flv_states), deferred :: get_blha_flv_states
    procedure :: allocate_components => pcm_allocate_components
    procedure(pcm_categorize_components), deferred :: categorize_components
    procedure(pcm_allocate_cores), deferred :: allocate_cores
    procedure(pcm_prepare_any_external_code), deferred :: &
         prepare_any_external_code
    procedure(pcm_setup_blha), deferred :: setup_blha
    procedure(pcm_prepare_blha_core), deferred :: prepare_blha_core
    procedure(pcm_setup_mci), deferred :: setup_mci
    procedure(pcm_call_dispatch_mci), deferred :: call_dispatch_mci
    procedure(pcm_complete_setup), deferred :: complete_setup
    procedure :: get_i_core => pcm_get_i_core
    procedure(pcm_init_phs_config), deferred :: init_phs_config
    procedure(pcm_init_component), deferred :: init_component
    procedure :: record_inactive_components => pcm_record_inactive_components
  end type pcm_t

  type, abstract :: pcm_workspace_t
!     class(pcm_t), pointer :: config => null ()
    logical :: bad_point = .false.
  contains
    procedure(pcm_work_final), deferred :: final
    procedure(pcm_work_is_nlo), deferred :: is_nlo
    procedure :: is_valid => pcm_work_is_valid
    procedure :: set_bad_point => pcm_work_set_bad_point
  end type pcm_workspace_t


  abstract interface
     subroutine pcm_allocate_workspace (pcm, work)
       import
       class(pcm_t), intent(in) :: pcm
       class(pcm_workspace_t), intent(inout), allocatable :: work
     end subroutine pcm_allocate_workspace
  end interface

  abstract interface
     function pcm_is_nlo (pcm) result (is_nlo)
        import
        logical :: is_nlo
        class(pcm_t), intent(in) :: pcm
     end function pcm_is_nlo
  end interface

  abstract interface
     subroutine pcm_final (pcm)
        import
        class(pcm_t), intent(inout) :: pcm
     end subroutine pcm_final
  end interface

  abstract interface
     subroutine pcm_init (pcm, env, meta)
       import
       class(pcm_t), intent(out) :: pcm
       type(process_environment_t), intent(in) :: env
       type(process_metadata_t), intent(in) :: meta
     end subroutine pcm_init
  end interface

  abstract interface
     subroutine pcm_set_blha_methods (pcm, blha_master, var_list)
       import
       class(pcm_t), intent(inout) :: pcm
       type(blha_master_t), intent(inout) :: blha_master
       type(var_list_t), intent(in) :: var_list
     end subroutine pcm_set_blha_methods
  end interface

  abstract interface
     subroutine pcm_get_blha_flv_states (pcm, core_entry, flv_born, flv_real)
       import
       class(pcm_t), intent(in) :: pcm
       type(core_entry_t), dimension(:), intent(in) :: core_entry
       integer, dimension(:,:), allocatable, intent(out) :: flv_born
       integer, dimension(:,:), allocatable, intent(out) :: flv_real
     end subroutine pcm_get_blha_flv_states
  end interface

  abstract interface
     subroutine pcm_categorize_components (pcm, config)
       import
       class(pcm_t), intent(inout) :: pcm
       type(process_config_data_t), intent(in) :: config
     end subroutine pcm_categorize_components
  end interface

  abstract interface
     subroutine pcm_allocate_cores (pcm, config, core_entry)
       import
       class(pcm_t), intent(inout) :: pcm
       type(process_config_data_t), intent(in) :: config
       type(core_entry_t), dimension(:), allocatable, intent(out) :: core_entry
     end subroutine pcm_allocate_cores
  end interface

  abstract interface
     subroutine pcm_prepare_any_external_code &
          (pcm, core_entry, i_core, libname, model, var_list)
       import
       class(pcm_t), intent(in) :: pcm
       type(core_entry_t), intent(inout) :: core_entry
       integer, intent(in) :: i_core
       type(string_t), intent(in) :: libname
       type(model_data_t), intent(in), target :: model
       type(var_list_t), intent(in) :: var_list
     end subroutine pcm_prepare_any_external_code
  end interface

  abstract interface
     subroutine pcm_setup_blha (pcm, core_entry)
       import
       class(pcm_t), intent(in) :: pcm
       type(core_entry_t), intent(inout) :: core_entry
     end subroutine pcm_setup_blha
  end interface

  abstract interface
     subroutine pcm_prepare_blha_core (pcm, core_entry, model)
       import
       class(pcm_t), intent(in) :: pcm
       type(core_entry_t), intent(inout) :: core_entry
       class(model_data_t), intent(in), target :: model
     end subroutine pcm_prepare_blha_core
  end interface

  abstract interface
     subroutine dispatch_mci_proc (mci, var_list, process_id, is_nlo)
       import
       class(mci_t), allocatable, intent(out) :: mci
       type(var_list_t), intent(in) :: var_list
       type(string_t), intent(in) :: process_id
       logical, intent(in), optional :: is_nlo
     end subroutine dispatch_mci_proc
  end interface

  abstract interface
     subroutine pcm_setup_mci (pcm, mci_entry)
       import
       class(pcm_t), intent(inout) :: pcm
       type(process_mci_entry_t), &
            dimension(:), allocatable, intent(out) :: mci_entry
     end subroutine pcm_setup_mci
  end interface

  abstract interface
     subroutine pcm_call_dispatch_mci (pcm, &
          dispatch_mci, var_list, process_id, mci_template)
       import
       class(pcm_t), intent(inout) :: pcm
       procedure(dispatch_mci_proc) :: dispatch_mci
       type(var_list_t), intent(in) :: var_list
       type(string_t), intent(in) :: process_id
       class(mci_t), intent(out), allocatable :: mci_template
     end subroutine pcm_call_dispatch_mci
  end interface

  abstract interface
     subroutine pcm_complete_setup (pcm, core_entry, component, model)
       import
       class(pcm_t), intent(inout) :: pcm
       type(core_entry_t), dimension(:), intent(in) :: core_entry
       type(process_component_t), dimension(:), intent(inout) :: component
       type(model_t), intent(in), target :: model
     end subroutine pcm_complete_setup
  end interface

  abstract interface
     subroutine pcm_init_phs_config &
          (pcm, phs_entry, meta, env, phs_par, mapping_defs)
       import
       class(pcm_t), intent(inout) :: pcm
       type(process_phs_config_t), &
            dimension(:), allocatable, intent(out) :: phs_entry
       type(process_metadata_t), intent(in) :: meta
       type(process_environment_t), intent(in) :: env
       type(mapping_defaults_t), intent(in) :: mapping_defs
       type(phs_parameters_t), intent(in) :: phs_par
     end subroutine pcm_init_phs_config
  end interface

  abstract interface
     subroutine pcm_init_component &
          (pcm, component, i, active, phs_config, env, meta, config)
       import
       class(pcm_t), intent(in) :: pcm
       type(process_component_t), intent(out) :: component
       integer, intent(in) :: i
       logical, intent(in) :: active
       class(phs_config_t), allocatable, intent(in) :: phs_config
       type(process_environment_t), intent(in) :: env
       type(process_metadata_t), intent(in) :: meta
       type(process_config_data_t), intent(in) :: config
     end subroutine pcm_init_component
  end interface

  abstract interface
     subroutine pcm_work_final (pcm_work)
        import
        class(pcm_workspace_t), intent(inout) :: pcm_work
     end subroutine pcm_work_final
  end interface

  abstract interface
     function pcm_work_is_nlo (pcm_work) result (is_nlo)
        import
        logical :: is_nlo
        class(pcm_workspace_t), intent(inout) :: pcm_work
      end function pcm_work_is_nlo
  end interface


  interface
    module function core_entry_get_core_ptr (core_entry) result (core)
      class(core_entry_t), intent(in), target :: core_entry
      class(prc_core_t), pointer :: core
    end function core_entry_get_core_ptr
    module subroutine core_entry_configure (core_entry, lib, id)
      class(core_entry_t), intent(inout) :: core_entry
      type(process_library_t), intent(in), target :: lib
      type(string_t), intent(in) :: id
    end subroutine core_entry_configure
    module subroutine pcm_set_blha_defaults (pcm, polarized_beams, var_list)
      class(pcm_t), intent(inout) :: pcm
      type(var_list_t), intent(in) :: var_list
      logical, intent(in) :: polarized_beams
    end subroutine pcm_set_blha_defaults
    module subroutine pcm_allocate_components (pcm, comp, meta)
      class(pcm_t), intent(inout) :: pcm
      type(process_component_t), dimension(:), allocatable, intent(out) :: comp
      type(process_metadata_t), intent(in) :: meta
    end subroutine pcm_allocate_components
    module function pcm_get_i_core (pcm, i_component) result (i_core)
      class(pcm_t), intent(in) :: pcm
      integer, intent(in) :: i_component
      integer :: i_core
    end function pcm_get_i_core
    module subroutine pcm_record_inactive_components (pcm, component, meta)
      class(pcm_t), intent(inout) :: pcm
      type(process_component_t), dimension(:), intent(in) :: component
      type(process_metadata_t), intent(inout) :: meta
    end subroutine pcm_record_inactive_components
    module function pcm_work_is_valid (pcm_work) result (valid)
      logical :: valid
      class(pcm_workspace_t), intent(in) :: pcm_work
    end function pcm_work_is_valid
    pure module subroutine pcm_work_set_bad_point (pcm_work, bad_point)
      class(pcm_workspace_t), intent(inout) :: pcm_work
      logical, intent(in) :: bad_point
    end subroutine pcm_work_set_bad_point
  end interface

end module pcm_base
