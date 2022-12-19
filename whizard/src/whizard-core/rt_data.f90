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

module rt_data

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use lexers
  use parser
  use models
  use subevents
  use pdg_arrays
  use variables, only: var_list_t
  use process_libraries
  use prclib_stacks
  use prc_core, only: helicity_selection_t
  use beam_structures
  use event_base, only: event_callback_t
  use user_files
  use process_stacks
  use iterations

  implicit none
  private

  public :: rt_data_t
  public :: fix_system_dependencies
  public :: show_description_of_string
  public :: show_tex_descriptions

  type :: rt_parse_nodes_t
     type(parse_node_t), pointer :: cuts_lexpr => null ()
     type(parse_node_t), pointer :: scale_expr => null ()
     type(parse_node_t), pointer :: fac_scale_expr => null ()
     type(parse_node_t), pointer :: ren_scale_expr => null ()
     type(parse_node_t), pointer :: weight_expr => null ()
     type(parse_node_t), pointer :: selection_lexpr => null ()
     type(parse_node_t), pointer :: reweight_expr => null ()
     type(parse_node_t), pointer :: analysis_lexpr => null ()
     type(parse_node_p), dimension(:), allocatable :: alt_setup
   contains
     procedure :: clear => rt_parse_nodes_clear
     procedure :: write => rt_parse_nodes_write
     procedure :: show => rt_parse_nodes_show
  end type rt_parse_nodes_t

  type :: rt_data_t
     type(lexer_t), pointer :: lexer => null ()
     type(rt_data_t), pointer :: context => null ()
     type(string_t), dimension(:), allocatable :: export
     type(var_list_t) :: var_list
     type(iterations_list_t) :: it_list
     type(os_data_t) :: os_data
     type(model_list_t) :: model_list
     type(model_t), pointer :: model => null ()
     logical :: model_is_copy = .false.
     type(model_t), pointer :: preload_model => null ()
     type(model_t), pointer :: fallback_model => null ()
     type(prclib_stack_t) :: prclib_stack
     type(process_library_t), pointer :: prclib => null ()
     type(beam_structure_t) :: beam_structure
     type(rt_parse_nodes_t) :: pn
     type(process_stack_t) :: process_stack
     type(string_t), dimension(:), allocatable :: sample_fmt
     class(event_callback_t), allocatable :: event_callback
     type(file_list_t), pointer :: out_files => null ()
     logical :: quit = .false.
     integer :: quit_code = 0
     type(string_t) :: logfile
     logical :: nlo_fixed_order = .false.
     logical, dimension(0:5) :: selected_nlo_parts = .false.
     integer, dimension(:), allocatable :: nlo_component
   contains
     procedure :: write => rt_data_write
     procedure :: write_vars => rt_data_write_vars
     procedure :: write_model_list => rt_data_write_model_list
     procedure :: write_libraries => rt_data_write_libraries
     procedure :: write_beams => rt_data_write_beams
     procedure :: write_expr => rt_data_write_expr
     procedure :: write_process_stack => rt_data_write_process_stack
     procedure :: write_var_descriptions => rt_data_write_var_descriptions
     procedure :: show_description_of_string => rt_data_show_description_of_string
     procedure :: clear_beams => rt_data_clear_beams
     procedure :: global_init => rt_data_global_init
     procedure :: local_init => rt_data_local_init
     procedure :: init_pointer_variables => rt_data_init_pointer_variables
     procedure :: activate => rt_data_activate
     procedure :: deactivate => rt_data_deactivate
     procedure :: copy_globals => rt_data_copy_globals
     procedure :: restore_globals => rt_data_restore_globals
     procedure :: write_exports => rt_data_write_exports
     procedure :: get_n_export => rt_data_get_n_export
     procedure :: append_exports => rt_data_append_exports
     procedure :: handle_exports => rt_data_handle_exports
     procedure :: transfer_process_stack => rt_data_transfer_process_stack
     procedure :: final => rt_data_global_final
     procedure :: local_final => rt_data_local_final
     procedure :: read_model => rt_data_read_model
     procedure :: read_ufo_model => rt_data_read_ufo_model
     procedure :: init_fallback_model => rt_data_init_fallback_model
     procedure :: select_model => rt_data_select_model
     procedure :: unselect_model => rt_data_unselect_model
     procedure :: ensure_model_copy => rt_data_ensure_model_copy
     procedure :: model_set_real => rt_data_model_set_real
     procedure :: modify_particle => rt_data_modify_particle
     procedure :: get_var_list_ptr => rt_data_get_var_list_ptr
     procedure :: append_log => rt_data_append_log
     procedure :: append_int => rt_data_append_int
     procedure :: append_real => rt_data_append_real
     procedure :: append_cmplx => rt_data_append_cmplx
     procedure :: append_subevt => rt_data_append_subevt
     procedure :: append_pdg_array => rt_data_append_pdg_array
     procedure :: append_string => rt_data_append_string
     procedure :: import_values => rt_data_import_values
     procedure :: unset_values => rt_data_unset_values
     procedure :: set_log => rt_data_set_log
     procedure :: set_int => rt_data_set_int
     procedure :: set_real => rt_data_set_real
     procedure :: set_cmplx => rt_data_set_cmplx
     procedure :: set_subevt => rt_data_set_subevt
     procedure :: set_pdg_array => rt_data_set_pdg_array
     procedure :: set_string => rt_data_set_string
     procedure :: get_lval => rt_data_get_lval
     procedure :: get_ival => rt_data_get_ival
     procedure :: get_rval => rt_data_get_rval
     procedure :: get_cval => rt_data_get_cval
     procedure :: get_pval => rt_data_get_pval
     procedure :: get_aval => rt_data_get_aval
     procedure :: get_sval => rt_data_get_sval
     procedure :: contains => rt_data_contains
     procedure :: is_known => rt_data_is_known
     procedure :: add_prclib => rt_data_add_prclib
     procedure :: update_prclib => rt_data_update_prclib
     procedure :: get_helicity_selection => rt_data_get_helicity_selection
     procedure :: show_beams => rt_data_show_beams
     procedure :: get_sqrts => rt_data_get_sqrts
     procedure :: pacify => rt_data_pacify
     procedure :: set_event_callback => rt_data_set_event_callback
     procedure :: has_event_callback => rt_data_has_event_callback
     procedure :: get_event_callback => rt_data_get_event_callback
  end type rt_data_t


  interface
    module subroutine rt_parse_nodes_clear (rt_pn, name)
      class(rt_parse_nodes_t), intent(inout) :: rt_pn
      type(string_t), intent(in) :: name
    end subroutine rt_parse_nodes_clear
    module subroutine rt_parse_nodes_write (object, unit)
      class(rt_parse_nodes_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine rt_parse_nodes_write
    module subroutine rt_parse_nodes_show (rt_pn, name, unit)
      class(rt_parse_nodes_t), intent(in) :: rt_pn
      type(string_t), intent(in) :: name
      integer, intent(in), optional :: unit
    end subroutine rt_parse_nodes_show
    module subroutine rt_data_write (object, unit, vars, pacify)
      class(rt_data_t), intent(in) :: object
      integer, intent(in), optional :: unit
      type(string_t), dimension(:), intent(in), optional :: vars
      logical, intent(in), optional :: pacify
    end subroutine rt_data_write
    module subroutine rt_data_write_vars (object, unit, vars)
      class(rt_data_t), intent(in), target :: object
      integer, intent(in), optional :: unit
      type(string_t), dimension(:), intent(in) :: vars
    end subroutine rt_data_write_vars
    module subroutine rt_data_write_model_list (object, unit)
      class(rt_data_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine rt_data_write_model_list
    module subroutine rt_data_write_libraries (object, unit, libpath)
      class(rt_data_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: libpath
    end subroutine rt_data_write_libraries
    module subroutine rt_data_write_beams (object, unit)
      class(rt_data_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine rt_data_write_beams
    module subroutine rt_data_write_expr (object, unit)
      class(rt_data_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine rt_data_write_expr
    module subroutine rt_data_write_process_stack (object, unit)
      class(rt_data_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine rt_data_write_process_stack
    module subroutine rt_data_write_var_descriptions &
         (rt_data, unit, ascii_output)
      class(rt_data_t), intent(in) :: rt_data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: ascii_output
    end subroutine rt_data_write_var_descriptions
    module subroutine rt_data_show_description_of_string (rt_data, string, &
         unit, ascii_output)
      class(rt_data_t), intent(in) :: rt_data
      type(string_t), intent(in) :: string
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: ascii_output
    end subroutine rt_data_show_description_of_string
    module subroutine rt_data_clear_beams (global)
      class(rt_data_t), intent(inout) :: global
    end subroutine rt_data_clear_beams
    module subroutine rt_data_global_init (global, paths, logfile)
      class(rt_data_t), intent(out), target :: global
      type(paths_t), intent(in), optional :: paths
      type(string_t), intent(in), optional :: logfile
    end subroutine rt_data_global_init
    module subroutine rt_data_local_init (local, global, env)
      class(rt_data_t), intent(inout), target :: local
      type(rt_data_t), intent(in), target :: global
      integer, intent(in), optional :: env
    end subroutine rt_data_local_init
    module subroutine rt_data_init_pointer_variables (local)
      class(rt_data_t), intent(inout), target :: local
    end subroutine rt_data_init_pointer_variables
    module subroutine rt_data_deactivate (local, global, keep_local)
      class(rt_data_t), intent(inout), target :: local
      class(rt_data_t), intent(inout), optional, target :: global
      logical, intent(in), optional :: keep_local
    end subroutine rt_data_deactivate
    module subroutine rt_data_copy_globals (global, local)
      class(rt_data_t), intent(in) :: global
      class(rt_data_t), intent(inout) :: local
    end subroutine rt_data_copy_globals
    module subroutine rt_data_restore_globals (global, local)
      class(rt_data_t), intent(inout) :: global
      class(rt_data_t), intent(inout) :: local
    end subroutine rt_data_restore_globals
    module subroutine rt_data_write_exports (rt_data, unit)
      class(rt_data_t), intent(in) :: rt_data
      integer, intent(in), optional :: unit
    end subroutine rt_data_write_exports
    module function rt_data_get_n_export (rt_data) result (n)
      class(rt_data_t), intent(in) :: rt_data
      integer :: n
    end function rt_data_get_n_export
    module subroutine rt_data_append_exports (rt_data, export)
      class(rt_data_t), intent(inout) :: rt_data
      type(string_t), dimension(:), intent(in) :: export
    end subroutine rt_data_append_exports
    module subroutine rt_data_handle_exports (local, global)
      class(rt_data_t), intent(inout), target :: local
      class(rt_data_t), intent(inout), target :: global
    end subroutine rt_data_handle_exports
    module subroutine rt_data_transfer_process_stack (local, global)
      class(rt_data_t), intent(inout), target :: local
      class(rt_data_t), intent(inout), target :: global
    end subroutine rt_data_transfer_process_stack
    module subroutine rt_data_global_final (global)
      class(rt_data_t), intent(inout) :: global
    end subroutine rt_data_global_final
    module subroutine rt_data_local_final (local)
      class(rt_data_t), intent(inout) :: local
    end subroutine rt_data_local_final
    module subroutine rt_data_read_model (global, name, model, scheme)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: scheme
      type(model_t), pointer, intent(out) :: model
    end subroutine rt_data_read_model
    module subroutine rt_data_read_ufo_model (global, name, model, ufo_path)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      type(model_t), pointer, intent(out) :: model
      type(string_t), intent(in), optional :: ufo_path
    end subroutine rt_data_read_ufo_model
    module subroutine rt_data_init_fallback_model (global, name, filename)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name, filename
    end subroutine rt_data_init_fallback_model
    module subroutine rt_data_select_model (global, name, scheme, ufo, ufo_path)
      class(rt_data_t), intent(inout), target :: global
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: scheme
      logical, intent(in), optional :: ufo
      type(string_t), intent(in), optional :: ufo_path
    end subroutine rt_data_select_model
    module subroutine rt_data_unselect_model (global)
      class(rt_data_t), intent(inout), target :: global
    end subroutine rt_data_unselect_model
    module subroutine rt_data_ensure_model_copy (global)
      class(rt_data_t), intent(inout), target :: global
    end subroutine rt_data_ensure_model_copy
    module subroutine rt_data_model_set_real &
         (global, name, rval, verbose, pacified)
      class(rt_data_t), intent(inout), target :: global
      type(string_t), intent(in) :: name
      real(default), intent(in) :: rval
      logical, intent(in), optional :: verbose, pacified
    end subroutine rt_data_model_set_real
    module subroutine rt_data_modify_particle &
         (global, pdg, polarized, stable, decay, &
         isotropic_decay, diagonal_decay, decay_helicity)
      class(rt_data_t), intent(inout), target :: global
      integer, intent(in) :: pdg
      logical, intent(in), optional :: polarized, stable
      logical, intent(in), optional :: isotropic_decay, diagonal_decay
      integer, intent(in), optional :: decay_helicity
      type(string_t), dimension(:), intent(in), optional :: decay
    end subroutine rt_data_modify_particle
    module function rt_data_get_var_list_ptr (global) result (var_list)
      class(rt_data_t), intent(in), target :: global
      type(var_list_t), pointer :: var_list
    end function rt_data_get_var_list_ptr
    module subroutine rt_data_append_log (local, name, lval, intrinsic, user)
      class(rt_data_t), intent(inout) :: local
      type(string_t), intent(in) :: name
      logical, intent(in), optional :: lval
      logical, intent(in), optional :: intrinsic, user
    end subroutine rt_data_append_log
    module subroutine rt_data_append_int (local, name, ival, intrinsic, user)
      class(rt_data_t), intent(inout) :: local
      type(string_t), intent(in) :: name
      integer, intent(in), optional :: ival
      logical, intent(in), optional :: intrinsic, user
    end subroutine rt_data_append_int
    module subroutine rt_data_append_real (local, name, rval, intrinsic, user)
      class(rt_data_t), intent(inout) :: local
      type(string_t), intent(in) :: name
      real(default), intent(in), optional :: rval
      logical, intent(in), optional :: intrinsic, user
    end subroutine rt_data_append_real
    module subroutine rt_data_append_cmplx (local, name, cval, intrinsic, user)
      class(rt_data_t), intent(inout) :: local
      type(string_t), intent(in) :: name
      complex(default), intent(in), optional :: cval
      logical, intent(in), optional :: intrinsic, user
    end subroutine rt_data_append_cmplx
    module subroutine rt_data_append_subevt (local, name, pval, intrinsic, user)
      class(rt_data_t), intent(inout) :: local
      type(string_t), intent(in) :: name
      type(subevt_t), intent(in), optional :: pval
      logical, intent(in) :: intrinsic, user
    end subroutine rt_data_append_subevt
    module subroutine rt_data_append_pdg_array &
         (local, name, aval, intrinsic, user)
      class(rt_data_t), intent(inout) :: local
      type(string_t), intent(in) :: name
      type(pdg_array_t), intent(in), optional :: aval
      logical, intent(in), optional :: intrinsic, user
    end subroutine rt_data_append_pdg_array
    module subroutine rt_data_append_string (local, name, sval, intrinsic, user)
      class(rt_data_t), intent(inout) :: local
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: sval
      logical, intent(in), optional :: intrinsic, user
    end subroutine rt_data_append_string
    module subroutine rt_data_import_values (local)
      class(rt_data_t), intent(inout) :: local
    end subroutine rt_data_import_values
    module subroutine rt_data_unset_values (global)
      class(rt_data_t), intent(inout) :: global
    end subroutine rt_data_unset_values
    module subroutine rt_data_set_log &
         (global, name, lval, is_known, force, verbose)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      logical, intent(in) :: lval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: force, verbose
    end subroutine rt_data_set_log
    module subroutine rt_data_set_int &
         (global, name, ival, is_known, force, verbose)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      integer, intent(in) :: ival
      logical, intent(in) :: is_known
      logical, intent(in), optional :: force, verbose
    end subroutine rt_data_set_int
    module subroutine rt_data_set_real &
         (global, name, rval, is_known, force, verbose, pacified)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      real(default), intent(in) :: rval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: force, verbose, pacified
    end subroutine rt_data_set_real
    module subroutine rt_data_set_cmplx &
         (global, name, cval, is_known, force, verbose, pacified)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      complex(default), intent(in) :: cval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: force, verbose, pacified
    end subroutine rt_data_set_cmplx
    module subroutine rt_data_set_subevt &
         (global, name, pval, is_known, force, verbose)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      type(subevt_t), intent(in) :: pval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: force, verbose
    end subroutine rt_data_set_subevt
    module subroutine rt_data_set_pdg_array &
         (global, name, aval, is_known, force, verbose)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      type(pdg_array_t), intent(in) :: aval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: force, verbose
    end subroutine rt_data_set_pdg_array
    module subroutine rt_data_set_string &
         (global, name, sval, is_known, force, verbose)
      class(rt_data_t), intent(inout) :: global
      type(string_t), intent(in) :: name
      type(string_t), intent(in) :: sval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: force, verbose
    end subroutine rt_data_set_string
    module function rt_data_get_lval (global, name) result (lval)
      logical :: lval
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_get_lval
    module function rt_data_get_ival (global, name) result (ival)
      integer :: ival
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_get_ival
    module function rt_data_get_rval (global, name) result (rval)
      real(default) :: rval
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_get_rval
    module function rt_data_get_cval (global, name) result (cval)
      complex(default) :: cval
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_get_cval
    module function rt_data_get_aval (global, name) result (aval)
      type(pdg_array_t) :: aval
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_get_aval
    module function rt_data_get_pval (global, name) result (pval)
      type(subevt_t) :: pval
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_get_pval
    module function rt_data_get_sval (global, name) result (sval)
      type(string_t) :: sval
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_get_sval
    module function rt_data_contains (global, name) result (lval)
      logical :: lval
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_contains
    module  function rt_data_is_known (global, name) result (lval)
      logical :: lval
      class(rt_data_t), intent(in), target :: global
      type(string_t), intent(in) :: name
    end function rt_data_is_known
    module subroutine rt_data_add_prclib (global, prclib_entry)
      class(rt_data_t), intent(inout) :: global
      type(prclib_entry_t), intent(inout), pointer :: prclib_entry
    end subroutine rt_data_add_prclib
    module subroutine rt_data_update_prclib (global, lib)
      class(rt_data_t), intent(inout) :: global
      type(process_library_t), intent(in), target :: lib
    end subroutine rt_data_update_prclib
    module function rt_data_get_helicity_selection &
         (rt_data) result (helicity_selection)
      class(rt_data_t), intent(in) :: rt_data
      type(helicity_selection_t) :: helicity_selection
    end function rt_data_get_helicity_selection
    module subroutine rt_data_show_beams (rt_data, unit)
      class(rt_data_t), intent(in) :: rt_data
      integer, intent(in), optional :: unit
    end subroutine rt_data_show_beams
    module function rt_data_get_sqrts (rt_data) result (sqrts)
      class(rt_data_t), intent(in) :: rt_data
      real(default) :: sqrts
    end function rt_data_get_sqrts
    module subroutine rt_data_pacify (rt_data, efficiency_reset, error_reset)
      class(rt_data_t), intent(inout) :: rt_data
      logical, intent(in), optional :: efficiency_reset, error_reset
    end subroutine rt_data_pacify
    module subroutine rt_data_set_event_callback (global, callback)
      class(rt_data_t), intent(inout) :: global
      class(event_callback_t), intent(in) :: callback
    end subroutine rt_data_set_event_callback
    module function rt_data_has_event_callback (global) result (flag)
      class(rt_data_t), intent(in) :: global
      logical :: flag
    end function rt_data_has_event_callback
    module function rt_data_get_event_callback (global) result (callback)
      class(rt_data_t), intent(in) :: global
      class(event_callback_t), allocatable :: callback
    end function rt_data_get_event_callback
    module subroutine fix_system_dependencies (global)
      class(rt_data_t), intent(inout), target :: global
    end subroutine fix_system_dependencies
    module subroutine show_description_of_string (string)
      type(string_t), intent(in) :: string
    end subroutine show_description_of_string
    module subroutine show_tex_descriptions ()
    end subroutine show_tex_descriptions
  end interface

contains

  subroutine rt_data_activate (local)
    class(rt_data_t), intent(inout), target :: local
    class(rt_data_t), pointer :: global
    global => local%context
    if (associated (global)) then
       local%lexer => global%lexer
       call global%copy_globals (local)
       local%os_data = global%os_data
       local%logfile = global%logfile
       if (associated (global%prclib)) then
          local%prclib => &
               local%prclib_stack%get_library_ptr (global%prclib%get_name ())
       end if
       call local%import_values ()
       call local%process_stack%link (global%process_stack)
       local%it_list = global%it_list
       local%beam_structure = global%beam_structure
       local%pn = global%pn
       if (allocated (local%sample_fmt))  deallocate (local%sample_fmt)
       if (allocated (global%sample_fmt)) then
          allocate (local%sample_fmt (size (global%sample_fmt)), &
               source = global%sample_fmt)
       end if
       local%out_files => global%out_files
       local%model => global%model
       local%model_is_copy = .false.
    else if (.not. associated (local%model)) then
       local%model => local%preload_model
       local%model_is_copy = .false.
    end if
    if (associated (local%model)) then
       call local%model%link_var_list (local%var_list)
       call local%var_list%set_string (var_str ("$model_name"), &
            local%model%get_name (), is_known = .true.)
    else
       call local%var_list%set_string (var_str ("$model_name"), &
            var_str (""), is_known = .false.)
    end if
  end subroutine rt_data_activate


end module rt_data
