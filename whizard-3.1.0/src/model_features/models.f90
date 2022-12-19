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

module models

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds, only: default
  use kinds, only: c_default_float
  use iso_varying_string, string_t => varying_string
  use os_interface
  use model_data

  use syntax_rules
  use parser
  use variables
  use expr_base
  use eval_trees

  implicit none
  private

  public :: model_t
  public :: syntax_model_file_init
  public :: syntax_model_file_final
  public :: syntax_model_file_write
  public :: create_test_model
  public :: model_list_t

  integer, parameter :: PAR_NONE = 0, PAR_UNUSED = -1
  integer, parameter :: PAR_INDEPENDENT = 1, PAR_DERIVED = 2
  integer, parameter :: PAR_EXTERNAL = 3


  type :: parameter_t
     private
     integer :: type  = PAR_NONE
     class(modelpar_data_t), pointer :: data => null ()
     type(string_t) :: block_name
     integer, dimension(:), allocatable :: block_index
     type(parse_node_t), pointer :: pn => null ()
     class(expr_t), allocatable :: expr
   contains
     procedure :: init_independent_value => parameter_init_independent_value
     procedure :: init_independent => parameter_init_independent
     procedure :: init_derived => parameter_init_derived
     procedure :: init_external => parameter_init_external
     procedure :: init_unused => parameter_init_unused
     procedure :: final => parameter_final
     procedure :: reset_derived => parameter_reset_derived
     procedure :: write => parameter_write
     procedure :: show => parameter_show
  end type parameter_t

  type :: slha_entry_t
     integer, dimension(:), allocatable :: block_index
     integer :: i_par = 0
  end type slha_entry_t
  
  type :: slha_block_t
     type(string_t) :: name
     integer :: n_entry = 0
     type(slha_entry_t), dimension(:), allocatable :: entry
  end type slha_block_t
  
  type, extends (model_data_t) :: model_t
     private
     character(32) :: md5sum = ""
     logical :: ufo_model = .false.
     type(string_t) :: ufo_path
     type(string_t), dimension(:), allocatable :: schemes
     type(string_t), allocatable :: selected_scheme
     type(parameter_t), dimension(:), allocatable :: par
     integer :: n_slha_block = 0
     type(slha_block_t), dimension(:), allocatable :: slha_block
     integer :: max_par_name_length = 0
     integer :: max_field_name_length = 0
     type(var_list_t) :: var_list
     type(string_t) :: dlname
     procedure(model_init_external_parameters), nopass, pointer :: &
          init_external_parameters => null ()
     type(dlaccess_t) :: dlaccess
     type(parse_tree_t) :: parse_tree
   contains
     generic :: init => model_init
     procedure, private :: model_init
     procedure, private :: basic_init => model_basic_init
     procedure :: final => model_final
     procedure :: write => model_write
     procedure :: show => model_show
     procedure :: show_fields => model_show_fields
     procedure :: show_stable => model_show_stable
     procedure :: show_unstable => model_show_unstable
     procedure :: show_polarized => model_show_polarized
     procedure :: show_unpolarized => model_show_unpolarized
     procedure :: get_md5sum => model_get_md5sum
     procedure :: &
          set_parameter_constant => model_set_parameter_constant
     procedure, private :: &
          set_parameter_parse_node => model_set_parameter_parse_node
     procedure :: &
          set_parameter_external => model_set_parameter_external
     procedure :: &
          set_parameter_unused => model_set_parameter_unused
     procedure, private :: copy_parameter => model_copy_parameter
     procedure :: update_parameters => model_parameters_update
     procedure, private :: init_field => model_init_field
     procedure, private :: copy_field => model_copy_field
     procedure :: write_var_list => model_write_var_list
     procedure :: link_var_list => model_link_var_list
     procedure :: var_exists => model_var_exists
     procedure :: var_is_locked => model_var_is_locked
     procedure :: set_real => model_var_set_real
     procedure :: get_rval => model_var_get_rval
     procedure :: get_var_list_ptr => model_get_var_list_ptr
     procedure :: is_ufo_model => model_is_ufo_model
     procedure :: get_ufo_path => model_get_ufo_path
     procedure :: has_schemes => model_has_schemes
     procedure :: enable_schemes => model_enable_schemes
     procedure :: set_scheme => model_set_scheme
     procedure :: get_scheme => model_get_scheme
     procedure :: matches => model_matches
     procedure :: supports_custom_slha => model_supports_custom_slha
     procedure :: get_custom_slha_blocks => model_get_custom_slha_blocks
     procedure :: slha_lookup => model_slha_lookup
     procedure :: slha_set_par => model_slha_set_par
     procedure :: read => model_read
     procedure, private :: read_parameter => model_read_parameter
     procedure, private :: read_derived => model_read_derived
     procedure, private :: read_external => model_read_external
     procedure, private :: read_unused => model_read_unused
     procedure, private :: read_field => model_read_field
     procedure, private :: read_vertex => model_read_vertex
     procedure, private :: append_field_vars => model_append_field_vars
     procedure :: init_instance => model_copy
  end type model_t

  type, extends (model_t) :: model_entry_t
     type(model_entry_t), pointer :: next => null ()
  end type model_entry_t

  type :: model_list_t
     type(model_entry_t), pointer :: first => null ()
     type(model_entry_t), pointer :: last => null ()
     type(model_list_t), pointer :: context => null ()
   contains
     procedure :: write => model_list_write
     procedure :: link => model_list_link
     procedure, private :: import => model_list_import
     procedure :: add => model_list_add
     procedure :: read_model => model_list_read_model
     procedure :: append_copy => model_list_append_copy
     procedure :: model_exists => model_list_model_exists
     procedure :: get_model_ptr => model_list_get_model_ptr
     procedure :: final => model_list_final
  end type model_list_t


  abstract interface
     subroutine model_init_external_parameters (par) bind (C)
       import
       real(c_default_float), dimension(*), intent(inout) :: par
     end subroutine model_init_external_parameters
  end interface


  type(syntax_t), target, save :: syntax_model_file


  interface
    module subroutine parameter_init_independent_value &
         (par, par_data, name, value)
      class(parameter_t), intent(out) :: par
      class(modelpar_data_t), intent(in), target :: par_data
      type(string_t), intent(in) :: name
      real(default), intent(in) :: value
    end subroutine parameter_init_independent_value
    module subroutine parameter_init_external (par, par_data, name)
      class(parameter_t), intent(out) :: par
      class(modelpar_data_t), intent(in), target :: par_data
      type(string_t), intent(in) :: name
    end subroutine parameter_init_external
    module subroutine parameter_init_unused (par, par_data, name)
      class(parameter_t), intent(out) :: par
      class(modelpar_data_t), intent(in), target :: par_data
      type(string_t), intent(in) :: name
    end subroutine parameter_init_unused
    module subroutine parameter_final (par)
      class(parameter_t), intent(inout) :: par
    end subroutine parameter_final
    module subroutine parameter_reset_derived (par)
      class(parameter_t), intent(inout) :: par
    end subroutine parameter_reset_derived
    module subroutine parameter_write (par, unit, write_defs)
      class(parameter_t), intent(in) :: par
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: write_defs
    end subroutine parameter_write
    module subroutine parameter_show (par, l, u, partype)
      class(parameter_t), intent(in) :: par
      integer, intent(in) :: l, u
      integer, intent(in) :: partype
    end subroutine parameter_show
    module subroutine model_init &
         (model, name, libname, os_data, n_par, n_prt, n_vtx, ufo)
      class(model_t), intent(inout) :: model
      type(string_t), intent(in) :: name, libname
      type(os_data_t), intent(in) :: os_data
      integer, intent(in) :: n_par, n_prt, n_vtx
      logical, intent(in), optional :: ufo
    end subroutine model_init
    module subroutine model_basic_init (model, name, n_par, n_prt, n_vtx)
      class(model_t), intent(inout) :: model
      type(string_t), intent(in) :: name
      integer, intent(in) :: n_par, n_prt, n_vtx
    end subroutine model_basic_init
    module subroutine model_final (model)
      class(model_t), intent(inout) :: model
    end subroutine model_final
    module subroutine model_write (model, unit, verbose, &
         show_md5sum, show_variables, show_parameters, &
         show_particles, show_vertices, show_scheme)
      class(model_t), intent(in) :: model
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
      logical, intent(in), optional :: show_md5sum
      logical, intent(in), optional :: show_variables
      logical, intent(in), optional :: show_parameters
      logical, intent(in), optional :: show_particles
      logical, intent(in), optional :: show_vertices
      logical, intent(in), optional :: show_scheme
    end subroutine model_write
    module subroutine model_show (model, unit)
      class(model_t), intent(in) :: model
      integer, intent(in), optional :: unit
    end subroutine model_show
    module subroutine model_show_fields (model, l, unit)
      class(model_t), intent(in), target :: model
      integer, intent(in) :: l
      integer, intent(in), optional :: unit
    end subroutine model_show_fields
    module subroutine model_show_stable (model, unit)
      class(model_t), intent(in), target :: model
      integer, intent(in), optional :: unit
    end subroutine model_show_stable
    module subroutine model_show_unstable (model, unit)
      class(model_t), intent(in), target :: model
      integer, intent(in), optional :: unit
    end subroutine model_show_unstable
    module subroutine model_show_polarized (model, unit)
      class(model_t), intent(in), target :: model
      integer, intent(in), optional :: unit
    end subroutine model_show_polarized
    module subroutine model_show_unpolarized (model, unit)
      class(model_t), intent(in), target :: model
      integer, intent(in), optional :: unit
    end subroutine model_show_unpolarized
    module function model_get_md5sum (model) result (md5sum)
      character(32) :: md5sum
      class(model_t), intent(in) :: model
    end function model_get_md5sum
    module subroutine model_set_parameter_constant (model, i, name, value)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(string_t), intent(in) :: name
      real(default), intent(in) :: value
    end subroutine model_set_parameter_constant
    module subroutine model_set_parameter_parse_node &
         (model, i, name, pn, constant)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(string_t), intent(in) :: name
      type(parse_node_t), intent(in), target :: pn
      logical, intent(in) :: constant
    end subroutine model_set_parameter_parse_node
    module subroutine model_set_parameter_external (model, i, name)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(string_t), intent(in) :: name
    end subroutine model_set_parameter_external
    module subroutine model_set_parameter_unused (model, i, name)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(string_t), intent(in) :: name
    end subroutine model_set_parameter_unused
    module subroutine model_copy_parameter (model, i, par)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(parameter_t), intent(in) :: par
    end subroutine model_copy_parameter
    module subroutine model_parameters_update (model)
      class(model_t), intent(inout) :: model
    end subroutine model_parameters_update
    module subroutine model_init_field (model, i, longname, pdg)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(string_t), intent(in) :: longname
      integer, intent(in) :: pdg
    end subroutine model_init_field
    module subroutine model_copy_field (model, i, name_src)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(string_t), intent(in) :: name_src
    end subroutine model_copy_field
    module subroutine model_write_var_list (model, unit, follow_link)
      class(model_t), intent(in) :: model
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: follow_link
    end subroutine model_write_var_list
    module subroutine model_link_var_list (model, var_list)
      class(model_t), intent(inout) :: model
      type(var_list_t), intent(in), target :: var_list
    end subroutine model_link_var_list
    module function model_var_exists (model, name) result (flag)
      class(model_t), intent(in) :: model
      type(string_t), intent(in) :: name
      logical :: flag
    end function model_var_exists
    module function model_var_is_locked (model, name) result (flag)
      class(model_t), intent(in) :: model
      type(string_t), intent(in) :: name
      logical :: flag
    end function model_var_is_locked  
    module subroutine model_var_set_real (model, name, rval, verbose, pacified)
      class(model_t), intent(inout) :: model
      type(string_t), intent(in) :: name
      real(default), intent(in) :: rval
      logical, intent(in), optional :: verbose, pacified
    end subroutine model_var_set_real
    module function model_var_get_rval (model, name) result (rval)
      class(model_t), intent(in) :: model
      type(string_t), intent(in) :: name
      real(default) :: rval
    end function model_var_get_rval
    module function model_get_var_list_ptr (model) result (var_list)
      type(var_list_t), pointer :: var_list
      class(model_t), intent(in), target :: model
    end function model_get_var_list_ptr
    module function model_is_ufo_model (model) result (flag)
      class(model_t), intent(in) :: model
      logical :: flag
    end function model_is_ufo_model
    module function model_get_ufo_path (model) result (path)
      class(model_t), intent(in) :: model
      type(string_t) :: path
    end function model_get_ufo_path
    module function model_has_schemes (model) result (flag)
      logical :: flag
      class(model_t), intent(in) :: model
    end function model_has_schemes
    module subroutine model_enable_schemes (model, scheme)
      class(model_t), intent(inout) :: model
      type(string_t), dimension(:), intent(in) :: scheme
    end subroutine model_enable_schemes
    module subroutine model_set_scheme (model, scheme)
      class(model_t), intent(inout) :: model
      type(string_t), intent(in), optional :: scheme
    end subroutine model_set_scheme
    module function model_get_scheme (model) result (scheme)
      class(model_t), intent(in) :: model
      type(string_t) :: scheme
    end function model_get_scheme
    module function model_matches &
         (model, name, scheme, ufo, ufo_path) result (flag)
      logical :: flag
      class(model_t), intent(in) :: model
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: scheme
      logical, intent(in), optional :: ufo
      type(string_t), intent(in), optional :: ufo_path
    end function model_matches
    module function model_supports_custom_slha (model) result (flag)
      class(model_t), intent(in) :: model
      logical :: flag    
    end function model_supports_custom_slha
    module subroutine model_get_custom_slha_blocks (model, block_name)
      class(model_t), intent(in) :: model
      type(string_t), dimension(:), allocatable :: block_name
    end subroutine model_get_custom_slha_blocks
    module subroutine model_slha_lookup &
         (model, block_name, block_index, par_data)
      class(model_t), intent(in) :: model
      type(string_t), intent(in) :: block_name
      integer, dimension(:), intent(in) :: block_index
      class(modelpar_data_t), pointer, intent(out) :: par_data
    end subroutine model_slha_lookup
    module subroutine model_slha_set_par (model, block_name, block_index, value)
      class(model_t), intent(inout) :: model
      type(string_t), intent(in) :: block_name
      integer, dimension(:), intent(in) :: block_index
      real(default), intent(in) :: value
    end subroutine model_slha_set_par
    module subroutine syntax_model_file_init ()
    end subroutine syntax_model_file_init
    module subroutine syntax_model_file_final ()
    end subroutine syntax_model_file_final
    module subroutine syntax_model_file_write (unit)
      integer, intent(in), optional :: unit
    end subroutine syntax_model_file_write
    module subroutine model_read (model, filename, os_data, exist, &
         scheme, ufo, ufo_path_requested, rebuild_mdl)
      class(model_t), intent(out), target :: model
      type(string_t), intent(in) :: filename
      type(os_data_t), intent(in) :: os_data
      logical, intent(out), optional :: exist
      type(string_t), intent(in), optional :: scheme
      logical, intent(in), optional :: ufo
      type(string_t), intent(in), optional :: ufo_path_requested
      logical, intent(in), optional :: rebuild_mdl
    end subroutine model_read
    module subroutine model_read_parameter (model, i, node)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(parse_node_t), intent(in), target :: node
    end subroutine model_read_parameter
    module subroutine model_read_derived (model, i, node)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(parse_node_t), intent(in), target :: node
    end subroutine model_read_derived
    module subroutine model_read_external (model, i, node)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(parse_node_t), intent(in), target :: node
    end subroutine model_read_external
    module subroutine model_read_unused (model, i, node)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(parse_node_t), intent(in), target :: node
    end subroutine model_read_unused
    module subroutine model_read_field (model, i, node)
      class(model_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(parse_node_t), intent(in) :: node
    end subroutine model_read_field
    module subroutine model_read_vertex (model, i, node)
      class(model_t), intent(inout) :: model
      integer, intent(in) :: i
      type(parse_node_t), intent(in) :: node
    end subroutine model_read_vertex
    module subroutine model_append_field_vars (model)
      class(model_t), intent(inout) :: model
    end subroutine model_append_field_vars
    module subroutine create_test_model (model_name, test_model)
      type(string_t), intent(in) :: model_name
      type(model_t), intent(out), pointer :: test_model
    end subroutine create_test_model
    recursive module subroutine model_list_write &
         (object, unit, verbose, follow_link)
      class(model_list_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
      logical, intent(in), optional :: follow_link
    end subroutine model_list_write
    module subroutine model_list_link (model_list, context)
      class(model_list_t), intent(inout) :: model_list
      type(model_list_t), intent(in), target :: context
    end subroutine model_list_link
    module subroutine model_list_import (model_list, current, model)
      class(model_list_t), intent(inout) :: model_list
      type(model_entry_t), pointer, intent(inout) :: current
      type(model_t), optional, pointer, intent(out) :: model
    end subroutine model_list_import
    module subroutine model_list_add (model_list, &
         name, os_data, n_par, n_prt, n_vtx, model)
      class(model_list_t), intent(inout) :: model_list
      type(string_t), intent(in) :: name
      type(os_data_t), intent(in) :: os_data
      integer, intent(in) :: n_par, n_prt, n_vtx
      type(model_t), pointer :: model
    end subroutine model_list_add
    module subroutine model_list_read_model &
         (model_list, name, filename, os_data, model, &
         scheme, ufo, ufo_path, rebuild_mdl)
      class(model_list_t), intent(inout), target :: model_list
      type(string_t), intent(in) :: name, filename
      type(os_data_t), intent(in) :: os_data
      type(model_t), pointer, intent(inout) :: model
      type(string_t), intent(in), optional :: scheme
      logical, intent(in), optional :: ufo
      type(string_t), intent(in), optional :: ufo_path
      logical, intent(in), optional :: rebuild_mdl
    end subroutine model_list_read_model
    module subroutine model_list_append_copy (model_list, orig, model)
      class(model_list_t), intent(inout) :: model_list
      type(model_t), intent(in), target :: orig
      type(model_t), intent(out), pointer, optional :: model
    end subroutine model_list_append_copy
    recursive module function model_list_model_exists &
         (model_list, name, scheme, ufo, ufo_path, follow_link) result (exists)
      class(model_list_t), intent(in) :: model_list
      logical :: exists
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: scheme
      logical, intent(in), optional :: ufo
      type(string_t), intent(in), optional :: ufo_path
      logical, intent(in), optional :: follow_link
    end function model_list_model_exists
    recursive module function model_list_get_model_ptr &
         (model_list, name, scheme, ufo, ufo_path, follow_link) result (model)
      class(model_list_t), intent(in) :: model_list
      type(model_t), pointer :: model
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: scheme
      logical, intent(in), optional :: ufo
      type(string_t), intent(in), optional :: ufo_path
      logical, intent(in), optional :: follow_link
    end function model_list_get_model_ptr
    module subroutine model_list_final (model_list)
      class(model_list_t), intent(inout) :: model_list
    end subroutine model_list_final
    module subroutine model_copy (model, orig)
      class(model_t), intent(out), target :: model
      type(model_t), intent(in) :: orig
    end subroutine model_copy
  end interface

contains

  subroutine parameter_init_independent (par, par_data, name, pn)
    class(parameter_t), intent(out) :: par
    class(modelpar_data_t), intent(in), target :: par_data
    type(string_t), intent(in) :: name
    type(parse_node_t), intent(in), target :: pn
    par%type = PAR_INDEPENDENT
    par%pn => pn
    allocate (eval_tree_t :: par%expr)
    select type (expr => par%expr)
    type is (eval_tree_t)
       call expr%init_numeric_value (pn)
    end select
    par%data => par_data
    call par%data%init (name, par%expr%get_real ())
  end subroutine parameter_init_independent

  subroutine parameter_init_derived (par, par_data, name, pn, var_list)
    class(parameter_t), intent(out) :: par
    class(modelpar_data_t), intent(in), target :: par_data
    type(string_t), intent(in) :: name
    type(parse_node_t), intent(in), target :: pn
    type(var_list_t), intent(in), target :: var_list
    par%type = PAR_DERIVED
    par%pn => pn
    allocate (eval_tree_t :: par%expr)
    select type (expr => par%expr)
    type is (eval_tree_t)
       call expr%init_expr (pn, var_list=var_list)
    end select
    par%data => par_data
!    call par%expr%evaluate ()
    call par%data%init (name, 0._default)
  end subroutine parameter_init_derived


end module models
