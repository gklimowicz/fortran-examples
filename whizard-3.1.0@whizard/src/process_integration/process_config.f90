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

module process_config

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use sf_base
  use sf_mappings
  use mappings, only: mapping_defaults_t
  use phs_forests, only: phs_parameters_t
  use sm_qcd
  use integration_results
  use flavors
  use interactions
  use model_data
  use models
  use process_libraries
  use process_constants
  use prc_core
  use beams
  use mci_base
  use beam_structures
  use dispatch_beams, only: dispatch_qcd
  use phs_base
  use expr_base
  use variables

  implicit none
  private

  public :: flagged
  public :: set_flag
  public :: process_config_data_t
  public :: process_environment_t
  public :: process_metadata_t
  public :: process_phs_config_t
  public :: process_beam_config_t
  public :: process_component_t
  public :: process_term_t

  integer, parameter, public :: COMP_DEFAULT = 0
  integer, parameter, public :: COMP_REAL_FIN = 1
  integer, parameter, public :: COMP_MASTER = 2
  integer, parameter, public :: COMP_VIRT = 3
  integer, parameter, public :: COMP_REAL = 4
  integer, parameter, public :: COMP_REAL_SING = 5
  integer, parameter, public :: COMP_MISMATCH = 6
  integer, parameter, public :: COMP_PDF = 7
  integer, parameter, public :: COMP_SUB = 8
  integer, parameter, public :: COMP_RESUM = 9

  integer, parameter, public :: F_PACIFY = 1
  integer, parameter, public :: F_SHOW_VAR_LIST = 11
  integer, parameter, public :: F_SHOW_EXPRESSIONS = 12
  integer, parameter, public :: F_SHOW_LIB = 13
  integer, parameter, public :: F_SHOW_MODEL = 14
  integer, parameter, public :: F_SHOW_QCD = 15
  integer, parameter, public :: F_SHOW_OS_DATA = 16
  integer, parameter, public :: F_SHOW_RNG = 17
  integer, parameter, public :: F_SHOW_BEAMS = 18

  type :: process_config_data_t
     class(process_def_t), pointer :: process_def => null ()
     integer :: n_in = 0
     integer :: n_components = 0
     integer :: n_terms = 0
     integer :: n_mci = 0
     type(string_t) :: model_name
     class(model_data_t), pointer :: model => null ()
     type(qcd_t) :: qcd
     class(expr_factory_t), allocatable :: ef_cuts
     class(expr_factory_t), allocatable :: ef_scale
     class(expr_factory_t), allocatable :: ef_fac_scale
     class(expr_factory_t), allocatable :: ef_ren_scale
     class(expr_factory_t), allocatable :: ef_weight
     character(32) :: md5sum = ""
   contains
     procedure :: write => process_config_data_write
     procedure :: init => process_config_data_init
     procedure :: get_qcd => process_config_data_get_qcd
     procedure :: compute_md5sum => process_config_data_compute_md5sum
     procedure :: get_md5sum => process_config_data_get_md5sum
  end type process_config_data_t

  type :: process_environment_t
     private
     type(model_t), pointer :: model => null ()
     type(var_list_t), pointer :: var_list => null ()
     logical :: var_list_is_set = .false.
     type(process_library_t), pointer :: lib => null ()
     type(beam_structure_t) :: beam_structure
     type(os_data_t) :: os_data
   contains
     procedure :: final => process_environment_final
     procedure :: write => process_environment_write
     procedure :: write_formatted => process_environment_write_formatted
     ! generic :: write (formatted) => write_formatted
     procedure :: init => process_environment_init
     procedure :: got_var_list => process_environment_got_var_list
     procedure :: get_var_list_ptr => process_environment_get_var_list_ptr
     procedure :: get_model_ptr => process_environment_get_model_ptr
     procedure :: get_lib_ptr => process_environment_get_lib_ptr
     procedure :: reset_lib_ptr => process_environment_reset_lib_ptr
     procedure :: check_lib_sanity => process_environment_check_lib_sanity
     procedure :: fill_process_constants => &
          process_environment_fill_process_constants
     procedure :: get_beam_structure => process_environment_get_beam_structure
     procedure :: has_pdfs => process_environment_has_pdfs
     procedure :: has_polarized_beams => process_environment_has_polarized_beams
     procedure :: get_os_data => process_environment_get_os_data
  end type process_environment_t

  type :: process_metadata_t
     integer :: type = PRC_UNKNOWN
     type(string_t) :: id
     integer :: num_id = 0
     type(string_t) :: run_id
     type(string_t), allocatable :: lib_name
     integer :: lib_update_counter = 0
     integer :: lib_index = 0
     integer :: n_components = 0
     type(string_t), dimension(:), allocatable :: component_id
     type(string_t), dimension(:), allocatable :: component_description
     logical, dimension(:), allocatable :: active
   contains
     procedure :: write => process_metadata_write
     procedure :: show => process_metadata_show
     procedure :: init => process_metadata_init
     procedure :: deactivate_component => process_metadata_deactivate_component
  end type process_metadata_t

  type :: process_phs_config_t
     type(phs_parameters_t) :: phs_par
     type(mapping_defaults_t) :: mapping_defs
     class(phs_config_t), allocatable :: phs_config
   contains
     procedure :: write => process_phs_config_write
     procedure :: write_formatted => process_phs_config_write_formatted
     ! generic :: write (formatted) => write_formatted
  end type process_phs_config_t

  type :: process_beam_config_t
     type(beam_data_t) :: data
     integer :: n_strfun = 0
     integer :: n_channel = 1
     integer :: n_sfpar = 0
     type(sf_config_t), dimension(:), allocatable :: sf
     type(sf_channel_t), dimension(:), allocatable :: sf_channel
     logical :: azimuthal_dependence = .false.
     logical :: lab_is_cm = .true.
     character(32) :: md5sum = ""
     logical :: sf_trace = .false.
     type(string_t) :: sf_trace_file
   contains
     procedure :: write => process_beam_config_write
     procedure :: final => process_beam_config_final
     procedure :: init_beam_structure => process_beam_config_init_beam_structure
     procedure :: init_scattering => process_beam_config_init_scattering
     procedure :: init_decay => process_beam_config_init_decay
     procedure :: startup_message => process_beam_config_startup_message
     procedure :: init_sf_chain => process_beam_config_init_sf_chain
     procedure :: allocate_sf_channels => process_beam_config_allocate_sf_channels
     procedure :: set_sf_channel => process_beam_config_set_sf_channel
     procedure :: sf_startup_message => process_beam_config_sf_startup_message
     procedure :: get_pdf_set => process_beam_config_get_pdf_set
     procedure :: get_beam_file => process_beam_config_get_beam_file
     procedure :: compute_md5sum => process_beam_config_compute_md5sum
     procedure :: get_md5sum => process_beam_config_get_md5sum
     procedure :: has_structure_function => &
          process_beam_config_has_structure_function
  end type process_beam_config_t

  type :: process_component_t
     type(process_component_def_t), pointer :: config => null ()
     integer :: index = 0
     logical :: active = .false.
     integer, dimension(:), allocatable :: i_term
     integer :: i_mci = 0
     class(phs_config_t), allocatable :: phs_config
     character(32) :: md5sum_phs = ""
     integer :: component_type = COMP_DEFAULT
   contains
     procedure :: final => process_component_final
     procedure :: write => process_component_write
     procedure :: init => process_component_init
     procedure :: is_active => process_component_is_active
     procedure :: configure_phs => process_component_configure_phs
     procedure :: compute_md5sum => process_component_compute_md5sum
     procedure :: collect_channels => process_component_collect_channels
     procedure :: get_config => process_component_get_config
     procedure :: get_md5sum => process_component_get_md5sum
     procedure :: get_n_phs_par => process_component_get_n_phs_par
     procedure :: get_phs_config => process_component_get_phs_config
     procedure :: get_nlo_type => process_component_get_nlo_type
     procedure :: needs_mci_entry => process_component_needs_mci_entry
     procedure :: can_be_integrated => process_component_can_be_integrated
  end type process_component_t

  type :: process_term_t
     integer :: i_term_global = 0
     integer :: i_component = 0
     integer :: i_term = 0
     integer :: i_sub = 0
     integer :: i_core = 0
     integer :: n_allowed = 0
     type(process_constants_t) :: data
     real(default) :: alpha_s = 0
     integer, dimension(:), allocatable :: flv, hel, col
     integer :: n_sub, n_sub_color, n_sub_spin
     type(interaction_t) :: int
     type(interaction_t), pointer :: int_eff => null ()
   contains
     procedure :: write => process_term_write
     procedure :: write_state_summary => process_term_write_state_summary
     procedure :: final => process_term_final
     procedure :: init => process_term_init
     procedure :: setup_interaction => process_term_setup_interaction
     procedure :: get_process_constants => process_term_get_process_constants
  end type process_term_t


  interface
    module function flagged (v_list, id, def) result (flag)
      logical :: flag
      integer, dimension(:), intent(in) :: v_list
      integer, intent(in) :: id
      logical, intent(in), optional :: def
    end function flagged
    module subroutine set_flag (v_list, value, flag)
      integer, dimension(:), intent(inout), allocatable :: v_list
      integer, intent(in) :: value
      logical, intent(in), optional :: flag
    end subroutine set_flag
    module subroutine process_config_data_write &
         (config, u, counters, model, expressions)
      class(process_config_data_t), intent(in) :: config
      integer, intent(in) :: u
      logical, intent(in) :: counters
      logical, intent(in) :: model
      logical, intent(in) :: expressions
    end subroutine process_config_data_write
    module function process_config_data_get_qcd (config) result (qcd)
      class(process_config_data_t), intent(in) :: config
      type(qcd_t) :: qcd
    end function process_config_data_get_qcd
    module subroutine process_config_data_compute_md5sum (config)
      class(process_config_data_t), intent(inout) :: config
    end subroutine process_config_data_compute_md5sum
    pure module function process_config_data_get_md5sum (config) result (md5)
      character(32) :: md5
      class(process_config_data_t), intent(in) :: config
    end function process_config_data_get_md5sum
    module subroutine process_environment_final (env)
      class(process_environment_t), intent(inout) :: env
    end subroutine process_environment_final
    module subroutine process_environment_write (env, unit, &
         show_var_list, show_model, show_lib, show_beams, show_os_data)
      class(process_environment_t), intent(in) :: env
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_var_list
      logical, intent(in), optional :: show_model
      logical, intent(in), optional :: show_lib
      logical, intent(in), optional :: show_beams
      logical, intent(in), optional :: show_os_data
    end subroutine process_environment_write
    module subroutine process_environment_write_formatted &
         (dtv, unit, iotype, v_list, iostat, iomsg)
      class(process_environment_t), intent(in) :: dtv
      integer, intent(in) :: unit
      character(*), intent(in) :: iotype
      integer, dimension(:), intent(in) :: v_list
      integer, intent(out) :: iostat
      character(*), intent(inout) :: iomsg
    end subroutine process_environment_write_formatted
    module subroutine process_environment_init &
         (env, model, lib, os_data, var_list, beam_structure)
      class(process_environment_t), intent(out) :: env
      type(model_t), intent(in), target :: model
      type(process_library_t), intent(in), target :: lib
      type(os_data_t), intent(in) :: os_data
      type(var_list_t), intent(in), target, optional :: var_list
      type(beam_structure_t), intent(in), optional :: beam_structure
    end subroutine process_environment_init
    module function process_environment_got_var_list (env) result (flag)
      class(process_environment_t), intent(in) :: env
      logical :: flag
    end function process_environment_got_var_list
    module function process_environment_get_var_list_ptr (env) result (var_list)
      class(process_environment_t), intent(in) :: env
      type(var_list_t), pointer :: var_list
    end function process_environment_get_var_list_ptr
    module function process_environment_get_model_ptr (env) result (model)
      class(process_environment_t), intent(in) :: env
      type(model_t), pointer :: model
    end function process_environment_get_model_ptr
    module function process_environment_get_lib_ptr (env) result (lib)
      class(process_environment_t), intent(inout) :: env
      type(process_library_t), pointer :: lib
    end function process_environment_get_lib_ptr
    module subroutine process_environment_reset_lib_ptr (env)
      class(process_environment_t), intent(inout) :: env
    end subroutine process_environment_reset_lib_ptr
    module subroutine process_environment_check_lib_sanity (env, meta)
      class(process_environment_t), intent(in) :: env
      type(process_metadata_t), intent(in) :: meta
    end subroutine process_environment_check_lib_sanity
    module subroutine process_environment_fill_process_constants &
         (env, id, i_component, data)
      class(process_environment_t), intent(in) :: env
      type(string_t), intent(in) :: id
      integer, intent(in) :: i_component
      type(process_constants_t), intent(out) :: data
    end subroutine process_environment_fill_process_constants
    module function process_environment_get_beam_structure &
         (env) result (beam_structure)
      class(process_environment_t), intent(in) :: env
      type(beam_structure_t) :: beam_structure
    end function process_environment_get_beam_structure
    module function process_environment_has_pdfs (env) result (flag)
      class(process_environment_t), intent(in) :: env
      logical :: flag
    end function process_environment_has_pdfs
    module function process_environment_has_polarized_beams (env) result (flag)
      class(process_environment_t), intent(in) :: env
      logical :: flag
    end function process_environment_has_polarized_beams
    module function process_environment_get_os_data (env) result (os_data)
      class(process_environment_t), intent(in) :: env
      type(os_data_t) :: os_data
    end function process_environment_get_os_data
    module subroutine process_metadata_write (meta, u, screen)
      class(process_metadata_t), intent(in) :: meta
      integer, intent(in) :: u
      logical, intent(in) :: screen
    end subroutine process_metadata_write
    module subroutine process_metadata_show (meta, u, model_name)
      class(process_metadata_t), intent(in) :: meta
      integer, intent(in) :: u
      type(string_t), intent(in) :: model_name
    end subroutine process_metadata_show
    module subroutine process_metadata_init (meta, id, lib, var_list)
      class(process_metadata_t), intent(out) :: meta
      type(string_t), intent(in) :: id
      type(process_library_t), intent(in), target :: lib
      type(var_list_t), intent(in) :: var_list
    end subroutine process_metadata_init
    module subroutine process_metadata_deactivate_component (meta, i)
      class(process_metadata_t), intent(inout) :: meta
      integer, intent(in) :: i
    end subroutine process_metadata_deactivate_component
    module subroutine process_phs_config_write (phs_config, unit)
      class(process_phs_config_t), intent(in) :: phs_config
      integer, intent(in), optional :: unit
    end subroutine process_phs_config_write
    module subroutine process_phs_config_write_formatted &
         (dtv, unit, iotype, v_list, iostat, iomsg)
      class(process_phs_config_t), intent(in) :: dtv
      integer, intent(in) :: unit
      character(*), intent(in) :: iotype
      integer, dimension(:), intent(in) :: v_list
      integer, intent(out) :: iostat
      character(*), intent(inout) :: iomsg
    end subroutine process_phs_config_write_formatted
    module subroutine process_beam_config_write (object, unit, verbose)
      class(process_beam_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine process_beam_config_write
    module subroutine process_beam_config_final (object)
      class(process_beam_config_t), intent(inout) :: object
    end subroutine process_beam_config_final
    module subroutine process_beam_config_init_beam_structure &
         (beam_config, beam_structure, sqrts, model, decay_rest_frame)
      class(process_beam_config_t), intent(out) :: beam_config
      type(beam_structure_t), intent(in) :: beam_structure
      logical, intent(in), optional :: decay_rest_frame
      real(default), intent(in) :: sqrts
      class(model_data_t), intent(in), target :: model
    end subroutine process_beam_config_init_beam_structure
    module subroutine process_beam_config_init_scattering &
         (beam_config, flv_in, sqrts, beam_structure)
      class(process_beam_config_t), intent(out) :: beam_config
      type(flavor_t), dimension(2), intent(in) :: flv_in
      real(default), intent(in) :: sqrts
      type(beam_structure_t), intent(in), optional :: beam_structure
    end subroutine process_beam_config_init_scattering
    module subroutine process_beam_config_init_decay &
         (beam_config, flv_in, rest_frame, beam_structure)
      class(process_beam_config_t), intent(out) :: beam_config
      type(flavor_t), dimension(1), intent(in) :: flv_in
      logical, intent(in), optional :: rest_frame
      type(beam_structure_t), intent(in), optional :: beam_structure
    end subroutine process_beam_config_init_decay
    module subroutine process_beam_config_startup_message &
         (beam_config, unit, beam_structure)
      class(process_beam_config_t), intent(in) :: beam_config
      integer, intent(in), optional :: unit
      type(beam_structure_t), intent(in), optional :: beam_structure
    end subroutine process_beam_config_startup_message
    module subroutine process_beam_config_init_sf_chain &
         (beam_config, sf_config, sf_trace_file)
      class(process_beam_config_t), intent(inout) :: beam_config
      type(sf_config_t), dimension(:), intent(in) :: sf_config
      type(string_t), intent(in), optional :: sf_trace_file
    end subroutine process_beam_config_init_sf_chain
    module subroutine process_beam_config_allocate_sf_channels &
         (beam_config, n_channel)
      class(process_beam_config_t), intent(inout) :: beam_config
      integer, intent(in) :: n_channel
    end subroutine process_beam_config_allocate_sf_channels
    module subroutine process_beam_config_set_sf_channel &
         (beam_config, c, sf_channel)
      class(process_beam_config_t), intent(inout) :: beam_config
      integer, intent(in) :: c
      type(sf_channel_t), intent(in) :: sf_channel
    end subroutine process_beam_config_set_sf_channel
    module subroutine process_beam_config_sf_startup_message &
         (beam_config, sf_string, unit)
      class(process_beam_config_t), intent(in) :: beam_config
      type(string_t), intent(in) :: sf_string
      integer, intent(in), optional :: unit
    end subroutine process_beam_config_sf_startup_message
    module function process_beam_config_get_pdf_set &
         (beam_config) result (pdf_set)
      class(process_beam_config_t), intent(in) :: beam_config
      integer :: pdf_set
    end function process_beam_config_get_pdf_set
    module function process_beam_config_get_beam_file &
         (beam_config) result (file)
      class(process_beam_config_t), intent(in) :: beam_config
      type(string_t) :: file
    end function process_beam_config_get_beam_file
    module subroutine process_beam_config_compute_md5sum (beam_config)
      class(process_beam_config_t), intent(inout) :: beam_config
    end subroutine process_beam_config_compute_md5sum
    pure module function process_beam_config_get_md5sum &
         (beam_config) result (md5)
      character(32) :: md5
      class(process_beam_config_t), intent(in) :: beam_config
    end function process_beam_config_get_md5sum
    pure module function process_beam_config_has_structure_function &
         (beam_config) result (has_sf)
      logical :: has_sf
      class(process_beam_config_t), intent(in) :: beam_config
    end function process_beam_config_has_structure_function
    module subroutine process_component_final (object)
      class(process_component_t), intent(inout) :: object
    end subroutine process_component_final
    module subroutine process_component_write (object, unit)
      class(process_component_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine process_component_write
    module subroutine process_component_init (component, &
         i_component, env, meta, config, &
         active, &
         phs_config_template)
      class(process_component_t), intent(out) :: component
      integer, intent(in) :: i_component
      type(process_environment_t), intent(in) :: env
      type(process_metadata_t), intent(in) :: meta
      type(process_config_data_t), intent(in) :: config
      logical, intent(in) :: active
      class(phs_config_t), intent(in), allocatable :: phs_config_template
    end subroutine process_component_init
    elemental module function process_component_is_active &
         (component) result (active)
      logical :: active
      class(process_component_t), intent(in) :: component
    end function process_component_is_active
    module subroutine process_component_configure_phs &
         (component, sqrts, beam_config, rebuild, &
          ignore_mismatch, subdir)
      class(process_component_t), intent(inout) :: component
      real(default), intent(in) :: sqrts
      type(process_beam_config_t), intent(in) :: beam_config
      logical, intent(in), optional :: rebuild
      logical, intent(in), optional :: ignore_mismatch
      type(string_t), intent(in), optional :: subdir
    end subroutine process_component_configure_phs
    module subroutine process_component_compute_md5sum (component)
      class(process_component_t), intent(inout) :: component
    end subroutine process_component_compute_md5sum
    module subroutine process_component_collect_channels (component, coll)
      class(process_component_t), intent(inout) :: component
      type(phs_channel_collection_t), intent(inout) :: coll
    end subroutine process_component_collect_channels
    module function process_component_get_config (component) &
           result (config)
      type(process_component_def_t) :: config
      class(process_component_t), intent(in) :: component
    end function process_component_get_config
    pure module function process_component_get_md5sum (component) result (md5)
      type(string_t) :: md5
      class(process_component_t), intent(in) :: component
    end function process_component_get_md5sum
    module function process_component_get_n_phs_par (component) result (n_par)
      class(process_component_t), intent(in) :: component
      integer :: n_par
    end function process_component_get_n_phs_par
    module subroutine process_component_get_phs_config (component, phs_config)
      class(process_component_t), intent(in), target :: component
      class(phs_config_t), intent(out), pointer :: phs_config
    end subroutine process_component_get_phs_config
    elemental module function process_component_get_nlo_type &
         (component) result (nlo_type)
      integer :: nlo_type
      class(process_component_t), intent(in) :: component
    end function process_component_get_nlo_type
    module function process_component_needs_mci_entry &
         (component, combined_integration) result (value)
      logical :: value
      class(process_component_t), intent(in) :: component
      logical, intent(in), optional :: combined_integration
    end function process_component_needs_mci_entry
    elemental module function process_component_can_be_integrated &
         (component) result (active)
      logical :: active
      class(process_component_t), intent(in) :: component
    end function process_component_can_be_integrated
    module subroutine process_term_write (term, unit)
      class(process_term_t), intent(in) :: term
      integer, intent(in), optional :: unit
    end subroutine process_term_write
    module subroutine process_term_write_state_summary (term, core, unit)
      class(process_term_t), intent(in) :: term
      class(prc_core_t), intent(in) :: core
      integer, intent(in), optional :: unit
    end subroutine process_term_write_state_summary
    module subroutine process_term_final (term)
      class(process_term_t), intent(inout) :: term
    end subroutine process_term_final
    module subroutine process_term_init &
         (term, i_term_global, i_component, i_term, core, model, &
          nlo_type, use_beam_pol, subtraction_method, &
          has_pdfs, n_emitters)
      class(process_term_t), intent(inout), target :: term
      integer, intent(in) :: i_term_global
      integer, intent(in) :: i_component
      integer, intent(in) :: i_term
      class(prc_core_t), intent(inout) :: core
      class(model_data_t), intent(in), target :: model
      integer, intent(in), optional :: nlo_type
      logical, intent(in), optional :: use_beam_pol
      type(string_t), intent(in), optional :: subtraction_method
      logical, intent(in), optional :: has_pdfs
      integer, intent(in), optional :: n_emitters
    end subroutine process_term_init
    module subroutine process_term_setup_interaction (term, core, model, &
       nlo_type, pol_beams, has_pdfs, use_internal_color, n_emitters)
      class(process_term_t), intent(inout) :: term
      class(prc_core_t), intent(inout) :: core
      class(model_data_t), intent(in), target :: model
      logical, intent(in), optional :: pol_beams
      logical, intent(in), optional :: has_pdfs
      integer, intent(in), optional :: nlo_type
      logical, intent(in), optional :: use_internal_color
      integer, intent(in), optional :: n_emitters
    end subroutine process_term_setup_interaction
    module subroutine process_term_get_process_constants &
         (term, prc_constants)
      class(process_term_t), intent(inout) :: term
      type(process_constants_t), intent(out) :: prc_constants
    end subroutine process_term_get_process_constants
  end interface

contains

  subroutine process_config_data_init (config, meta, env)
    class(process_config_data_t), intent(out) :: config
    type(process_metadata_t), intent(in) :: meta
    type(process_environment_t), intent(in) :: env
    config%process_def => env%lib%get_process_def_ptr (meta%id)
    config%n_in = config%process_def%get_n_in ()
    config%n_components = size (meta%component_id)
    config%model => env%get_model_ptr ()
    config%model_name = config%model%get_name ()
    if (env%got_var_list ()) then
       call dispatch_qcd &
            (config%qcd, env%get_var_list_ptr (), env%get_os_data ())
    end if
  end subroutine process_config_data_init


end module process_config
