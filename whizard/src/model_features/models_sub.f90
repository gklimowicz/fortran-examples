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

submodule (models) models_s

  use io_units
  use diagnostics
  use md5
  use physics_defs, only: UNDEFINED

  use ifiles
  use lexers
  use pdg_arrays
  use ttv_formfactors, only: init_parameters

  implicit none

contains

  module subroutine parameter_init_independent_value &
       (par, par_data, name, value)
    class(parameter_t), intent(out) :: par
    class(modelpar_data_t), intent(in), target :: par_data
    type(string_t), intent(in) :: name
    real(default), intent(in) :: value
    par%type = PAR_INDEPENDENT
    par%data => par_data
    call par%data%init (name, value)
  end subroutine parameter_init_independent_value

  module subroutine parameter_init_external (par, par_data, name)
    class(parameter_t), intent(out) :: par
    class(modelpar_data_t), intent(in), target :: par_data
    type(string_t), intent(in) :: name
    par%type = PAR_EXTERNAL
    par%data => par_data
    call par%data%init (name, 0._default)
  end subroutine parameter_init_external

  module subroutine parameter_init_unused (par, par_data, name)
    class(parameter_t), intent(out) :: par
    class(modelpar_data_t), intent(in), target :: par_data
    type(string_t), intent(in) :: name
    par%type = PAR_UNUSED
    par%data => par_data
    call par%data%init (name, 0._default)
  end subroutine parameter_init_unused

  module subroutine parameter_final (par)
    class(parameter_t), intent(inout) :: par
    if (allocated (par%expr)) then
       call par%expr%final ()
    end if
  end subroutine parameter_final

  module subroutine parameter_reset_derived (par)
    class(parameter_t), intent(inout) :: par
    select case (par%type)
    case (PAR_DERIVED)
       call par%expr%evaluate ()
       par%data = par%expr%get_real ()
    end select
  end subroutine parameter_reset_derived

  module subroutine parameter_write (par, unit, write_defs)
    class(parameter_t), intent(in) :: par
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: write_defs
    logical :: defs
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    defs = .false.;  if (present (write_defs))  defs = write_defs
    select case (par%type)
    case (PAR_INDEPENDENT)
       write (u, "(3x,A)", advance="no")  "parameter"
       call par%data%write (u)
    case (PAR_DERIVED)
       write (u, "(3x,A)", advance="no")  "derived"
       call par%data%write (u)
    case (PAR_EXTERNAL)
       write (u, "(3x,A)", advance="no")  "external"
       call par%data%write (u)
    case (PAR_UNUSED)
       write (u, "(3x,A)", advance="no")  "unused"
       write (u, "(1x,A)", advance="no")  char (par%data%get_name ())
    end select
    select case (par%type)
    case (PAR_INDEPENDENT)
       if (allocated (par%block_index)) then
          write (u, "(1x,A,1x,A,*(1x,I0))") &
               "slha_entry", char (par%block_name), par%block_index
       else
          write (u, "(A)")
       end if
    case (PAR_DERIVED)
       if (defs) then
          call par%expr%write (unit)
       else
          write (u, "(A)")
       end if
    case default
       write (u, "(A)")
    end select
  end subroutine parameter_write

  module subroutine parameter_show (par, l, u, partype)
    class(parameter_t), intent(in) :: par
    integer, intent(in) :: l, u
    integer, intent(in) :: partype
    if (par%type == partype) then
       call par%data%show (l, u)
    end if
  end subroutine parameter_show

  module subroutine model_init &
       (model, name, libname, os_data, n_par, n_prt, n_vtx, ufo)
    class(model_t), intent(inout) :: model
    type(string_t), intent(in) :: name, libname
    type(os_data_t), intent(in) :: os_data
    integer, intent(in) :: n_par, n_prt, n_vtx
    logical, intent(in), optional :: ufo
    type(c_funptr) :: c_fptr
    type(string_t) :: libpath
    integer :: scheme_num
    scheme_num = model%get_scheme_num ()
    call model%basic_init (name, n_par, n_prt, n_vtx)
    if (present (ufo))  model%ufo_model = ufo
    call model%set_scheme_num (scheme_num)
    if (libname /= "") then
       if (.not. os_data%use_testfiles) then
          libpath = os_data%whizard_models_libpath_local
          model%dlname = os_get_dlname ( &
            libpath // "/" // libname, os_data, ignore=.true.)
       end if
       if (model%dlname == "") then
          libpath = os_data%whizard_models_libpath
          model%dlname = os_get_dlname (libpath // "/" // libname, os_data)
       end if
    else
       model%dlname = ""
    end if
    if (model%dlname /= "") then
       if (.not. dlaccess_is_open (model%dlaccess)) then
          if (logging) &
               call msg_message ("Loading model auxiliary library '" &
               // char (libpath) // "/" // char (model%dlname) // "'")
          call dlaccess_init (model%dlaccess, os_data%whizard_models_libpath, &
               model%dlname, os_data)
          if (dlaccess_has_error (model%dlaccess)) then
             call msg_message (char (dlaccess_get_error (model%dlaccess)))
             call msg_fatal ("Loading model auxiliary library '" &
                  // char (model%dlname) // "' failed")
             return
          end if
          c_fptr = dlaccess_get_c_funptr (model%dlaccess, &
               var_str ("init_external_parameters"))
          if (dlaccess_has_error (model%dlaccess)) then
             call msg_message (char (dlaccess_get_error (model%dlaccess)))
             call msg_fatal ("Loading function from auxiliary library '" &
                  // char (model%dlname) // "' failed")
             return
          end if
          call c_f_procpointer (c_fptr, model% init_external_parameters)
       end if
    end if
  end subroutine model_init

  module subroutine model_basic_init (model, name, n_par, n_prt, n_vtx)
    class(model_t), intent(inout) :: model
    type(string_t), intent(in) :: name
    integer, intent(in) :: n_par, n_prt, n_vtx
    allocate (model%par (n_par))
    call model%model_data_t%init (name, n_par, 0, n_prt, n_vtx)
  end subroutine model_basic_init

  module subroutine model_final (model)
    class(model_t), intent(inout) :: model
    integer :: i
    if (allocated (model%par)) then
       do i = 1, size (model%par)
          call model%par(i)%final ()
       end do
    end if
    call model%var_list%final (follow_link=.false.)
    if (model%dlname /= "")  call dlaccess_final (model%dlaccess)
    call parse_tree_final (model%parse_tree)
    call model%model_data_t%final ()
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
    logical :: verb, show_md5, show_par, show_var
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    verb = .false.;  if (present (verbose))  verb = verbose
    show_md5 = .true.;  if (present (show_md5sum)) &
         show_md5 = show_md5sum
    show_par = .true.;  if (present (show_parameters)) &
         show_par = show_parameters
    show_var = verb;  if (present (show_variables)) &
         show_var = show_variables
    write (u, "(A,A,A)") 'model "', char (model%get_name ()), '"'
    if (show_md5 .and. model%md5sum /= "") &
         write (u, "(1x,A,A,A)") "! md5sum = '", model%md5sum, "'"
    if (model%is_ufo_model ()) then
       write (u, "(1x,A)")  "! model derived from UFO source"
    else if (model%has_schemes ()) then
       write (u, "(1x,A)", advance="no")  "! schemes ="
       do i = 1, size (model%schemes)
          if (i > 1)  write (u, "(',')", advance="no")
          write (u, "(1x,A,A,A)", advance="no") &
               "'", char (model%schemes(i)), "'"
       end do
       write (u, *)
       if (allocated (model%selected_scheme)) then
          write (u, "(1x,A,A,A,I0,A)")  &
               "! selected scheme = '", char (model%get_scheme ()), &
               "' (", model%get_scheme_num (), ")"
       end if
    end if
    if (show_par) then
       write (u, "(A)")
       do i = 1, size (model%par)
          call model%par(i)%write (u, write_defs=verbose)
       end do
    end if
    call model%model_data_t%write (unit, verbose, &
         show_md5sum, show_variables, &
         show_parameters=.false., &
         show_particles=show_particles, &
         show_vertices=show_vertices, &
         show_scheme=show_scheme)
    if (show_var) then
       write (u, "(A)")
       call model%var_list%write (unit, follow_link=.false.)
    end if
  end subroutine model_write

  module subroutine model_show (model, unit)
    class(model_t), intent(in) :: model
    integer, intent(in), optional :: unit
    integer :: i, u, l
    u = given_output_unit (unit)
    write (u, "(A,1x,A)")  "Model:", char (model%get_name ())
    if (model%has_schemes ()) then
       write (u, "(2x,A,A,A,I0,A)")  "Scheme: '", &
            char (model%get_scheme ()), "' (", model%get_scheme_num (), ")"
    end if
    l = model%max_field_name_length
    call model%show_fields (l, u)
    l = model%max_par_name_length
    if (any (model%par%type == PAR_INDEPENDENT)) then
       write (u, "(2x,A)")  "Independent parameters:"
       do i = 1, size (model%par)
          call model%par(i)%show (l, u, PAR_INDEPENDENT)
       end do
    end if
    if (any (model%par%type == PAR_DERIVED)) then
       write (u, "(2x,A)")  "Derived parameters:"
       do i = 1, size (model%par)
          call model%par(i)%show (l, u, PAR_DERIVED)
       end do
    end if
    if (any (model%par%type == PAR_EXTERNAL)) then
       write (u, "(2x,A)")  "External parameters:"
       do i = 1, size (model%par)
          call model%par(i)%show (l, u, PAR_EXTERNAL)
       end do
    end if
    if (any (model%par%type == PAR_UNUSED)) then
       write (u, "(2x,A)")  "Unused parameters:"
       do i = 1, size (model%par)
          call model%par(i)%show (l, u, PAR_UNUSED)
       end do
    end if
  end subroutine model_show

  module subroutine model_show_fields (model, l, unit)
    class(model_t), intent(in), target :: model
    integer, intent(in) :: l
    integer, intent(in), optional :: unit
    type(field_data_t), pointer :: field
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(2x,A)")  "Particles:"
    do i = 1, model%get_n_field ()
       field => model%get_field_ptr_by_index (i)
       call field%show (l, u)
    end do
  end subroutine model_show_fields

  module subroutine model_show_stable (model, unit)
    class(model_t), intent(in), target :: model
    integer, intent(in), optional :: unit
    type(field_data_t), pointer :: field
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(A,1x)", advance="no")  "Stable particles:"
    do i = 1, model%get_n_field ()
       field => model%get_field_ptr_by_index (i)
       if (field%is_stable (.false.)) then
          write (u, "(1x,A)", advance="no")  char (field%get_name (.false.))
       end if
       if (field%has_antiparticle ()) then
          if (field%is_stable (.true.)) then
             write (u, "(1x,A)", advance="no")  char (field%get_name (.true.))
          end if
       end if
    end do
    write (u, *)
  end subroutine model_show_stable

  module subroutine model_show_unstable (model, unit)
    class(model_t), intent(in), target :: model
    integer, intent(in), optional :: unit
    type(field_data_t), pointer :: field
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(A,1x)", advance="no")  "Unstable particles:"
    do i = 1, model%get_n_field ()
       field => model%get_field_ptr_by_index (i)
       if (.not. field%is_stable (.false.)) then
          write (u, "(1x,A)", advance="no")  char (field%get_name (.false.))
       end if
       if (field%has_antiparticle ()) then
          if (.not. field%is_stable (.true.)) then
             write (u, "(1x,A)", advance="no")  char (field%get_name (.true.))
          end if
       end if
    end do
    write (u, *)
  end subroutine model_show_unstable

  module subroutine model_show_polarized (model, unit)
    class(model_t), intent(in), target :: model
    integer, intent(in), optional :: unit
    type(field_data_t), pointer :: field
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(A,1x)", advance="no")  "Polarized particles:"
    do i = 1, model%get_n_field ()
       field => model%get_field_ptr_by_index (i)
       if (field%is_polarized (.false.)) then
          write (u, "(1x,A)", advance="no") char (field%get_name (.false.))
       end if
       if (field%has_antiparticle ()) then
          if (field%is_polarized (.true.)) then
             write (u, "(1x,A)", advance="no") char (field%get_name (.true.))
          end if
       end if
    end do
    write (u, *)
  end subroutine model_show_polarized

  module subroutine model_show_unpolarized (model, unit)
    class(model_t), intent(in), target :: model
    integer, intent(in), optional :: unit
    type(field_data_t), pointer :: field
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(A,1x)", advance="no")  "Unpolarized particles:"
    do i = 1, model%get_n_field ()
       field => model%get_field_ptr_by_index (i)
       if (.not. field%is_polarized (.false.)) then
          write (u, "(1x,A)", advance="no") &
               char (field%get_name (.false.))
       end if
       if (field%has_antiparticle ()) then
          if (.not. field%is_polarized (.true.)) then
             write (u, "(1x,A)", advance="no") char (field%get_name (.true.))
          end if
       end if
    end do
    write (u, *)
  end subroutine model_show_unpolarized

  module function model_get_md5sum (model) result (md5sum)
    character(32) :: md5sum
    class(model_t), intent(in) :: model
    md5sum = model%md5sum
  end function model_get_md5sum

  module subroutine model_set_parameter_constant (model, i, name, value)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(string_t), intent(in) :: name
    real(default), intent(in) :: value
    logical, save, target :: known = .true.
    class(modelpar_data_t), pointer :: par_data
    real(default), pointer :: value_ptr
    par_data => model%get_par_real_ptr (i)
    call model%par(i)%init_independent_value (par_data, name, value)
    value_ptr => par_data%get_real_ptr ()
    call model%var_list%append_real_ptr (name, value_ptr, &
         is_known=known, intrinsic=.true.)
    model%max_par_name_length = max (model%max_par_name_length, len (name))
  end subroutine model_set_parameter_constant

  module subroutine model_set_parameter_parse_node &
       (model, i, name, pn, constant)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(string_t), intent(in) :: name
    type(parse_node_t), intent(in), target :: pn
    logical, intent(in) :: constant
    logical, save, target :: known = .true.
    class(modelpar_data_t), pointer :: par_data
    real(default), pointer :: value_ptr
    par_data => model%get_par_real_ptr (i)
    if (constant) then
       call model%par(i)%init_independent (par_data, name, pn)
    else
       call model%par(i)%init_derived (par_data, name, pn, model%var_list)
    end if
    value_ptr => par_data%get_real_ptr ()
    call model%var_list%append_real_ptr (name, value_ptr, &
         is_known=known, locked=.not.constant, intrinsic=.true.)
    model%max_par_name_length = max (model%max_par_name_length, len (name))
  end subroutine model_set_parameter_parse_node

  module subroutine model_set_parameter_external (model, i, name)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(string_t), intent(in) :: name
    logical, save, target :: known = .true.
    class(modelpar_data_t), pointer :: par_data
    real(default), pointer :: value_ptr
    par_data => model%get_par_real_ptr (i)
    call model%par(i)%init_external (par_data, name)
    value_ptr => par_data%get_real_ptr ()
    call model%var_list%append_real_ptr (name, value_ptr, &
         is_known=known, locked=.true., intrinsic=.true.)
    model%max_par_name_length = max (model%max_par_name_length, len (name))
  end subroutine model_set_parameter_external

  module subroutine model_set_parameter_unused (model, i, name)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(string_t), intent(in) :: name
    class(modelpar_data_t), pointer :: par_data
    par_data => model%get_par_real_ptr (i)
    call model%par(i)%init_unused (par_data, name)
    call model%var_list%append_real (name, locked=.true., intrinsic=.true.)
    model%max_par_name_length = max (model%max_par_name_length, len (name))
  end subroutine model_set_parameter_unused

  module subroutine model_copy_parameter (model, i, par)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(parameter_t), intent(in) :: par
    type(string_t) :: name
    real(default) :: value
    name = par%data%get_name ()
    select case (par%type)
    case (PAR_INDEPENDENT)
       if (associated (par%pn)) then
          call model%set_parameter_parse_node (i, name, par%pn, &
               constant = .true.)
       else
          value = par%data%get_real ()
          call model%set_parameter_constant (i, name, value)
       end if
       if (allocated (par%block_index)) then
          model%par(i)%block_name = par%block_name
          model%par(i)%block_index = par%block_index
       end if
    case (PAR_DERIVED)
       call model%set_parameter_parse_node (i, name, par%pn, &
            constant = .false.)
    case (PAR_EXTERNAL)
       call model%set_parameter_external (i, name)
    case (PAR_UNUSED)
       call model%set_parameter_unused (i, name)
    end select
  end subroutine model_copy_parameter

  module subroutine model_parameters_update (model)
    class(model_t), intent(inout) :: model
    integer :: i
    real(default), dimension(:), allocatable :: par
    do i = 1, size (model%par)
       call model%par(i)%reset_derived ()
    end do
    if (associated (model%init_external_parameters)) then
       allocate (par (model%get_n_real ()))
       call model%real_parameters_to_c_array (par)
       call model%init_external_parameters (par)
       call model%real_parameters_from_c_array (par)
       if (model%get_name() == var_str ("SM_tt_threshold")) &
            call set_threshold_parameters ()
    end if
  contains
    subroutine set_threshold_parameters ()
      real(default) :: mpole, wtop
      !!! !!! !!! Workaround for OS-X and BSD which do not load
      !!! !!! !!! the global values created previously. Therefore
      !!! !!! !!! a second initialization for the threshold model,
      !!! !!! !!! where M1S is required to compute the top mass.
      call init_parameters (mpole, wtop, &
           par(20), par(21), par(22), &
           par(19), par(39), par(4), par(1), &
           par(2), par(10), par(24), par(25), &
           par(23), par(26), par(27), par(29), &
           par(30), par(31), par(32), par(33), &
           par(36) > 0._default, par(28))
    end subroutine set_threshold_parameters
  end subroutine model_parameters_update

  module subroutine model_init_field (model, i, longname, pdg)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(string_t), intent(in) :: longname
    integer, intent(in) :: pdg
    type(field_data_t), pointer :: field
    field => model%get_field_ptr_by_index (i)
    call field%init (longname, pdg)
  end subroutine model_init_field

  module subroutine model_copy_field (model, i, name_src)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(string_t), intent(in) :: name_src
    type(field_data_t), pointer :: field_src, field
    field_src => model%get_field_ptr (name_src)
    field => model%get_field_ptr_by_index (i)
    call field%copy_from (field_src)
  end subroutine model_copy_field

  module subroutine model_write_var_list (model, unit, follow_link)
    class(model_t), intent(in) :: model
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: follow_link
    call model%var_list%write (unit, follow_link)
  end subroutine model_write_var_list

  module subroutine model_link_var_list (model, var_list)
    class(model_t), intent(inout) :: model
    type(var_list_t), intent(in), target :: var_list
    call model%var_list%link (var_list)
  end subroutine model_link_var_list

  module function model_var_exists (model, name) result (flag)
    class(model_t), intent(in) :: model
    type(string_t), intent(in) :: name
    logical :: flag
    flag = model%var_list%contains (name, follow_link=.false.)
  end function model_var_exists

  module function model_var_is_locked (model, name) result (flag)
    class(model_t), intent(in) :: model
    type(string_t), intent(in) :: name
    logical :: flag
    flag = model%var_list%is_locked (name, follow_link=.false.)
  end function model_var_is_locked

  module subroutine model_var_set_real (model, name, rval, verbose, pacified)
    class(model_t), intent(inout) :: model
    type(string_t), intent(in) :: name
    real(default), intent(in) :: rval
    logical, intent(in), optional :: verbose, pacified
    call model%var_list%set_real (name, rval, &
         is_known=.true., ignore=.false., &
         verbose=verbose, model_name=model%get_name (), pacified=pacified)
    call model%update_parameters ()
  end subroutine model_var_set_real

  module function model_var_get_rval (model, name) result (rval)
    class(model_t), intent(in) :: model
    type(string_t), intent(in) :: name
    real(default) :: rval
    rval = model%var_list%get_rval (name, follow_link=.false.)
  end function model_var_get_rval

  module function model_get_var_list_ptr (model) result (var_list)
    type(var_list_t), pointer :: var_list
    class(model_t), intent(in), target :: model
    var_list => model%var_list
  end function model_get_var_list_ptr

  module function model_is_ufo_model (model) result (flag)
    class(model_t), intent(in) :: model
    logical :: flag
    flag = model%ufo_model
  end function model_is_ufo_model

  module function model_get_ufo_path (model) result (path)
    class(model_t), intent(in) :: model
    type(string_t) :: path
    if (model%ufo_model) then
       path = model%ufo_path
    else
       path = ""
    end if
  end function model_get_ufo_path

  subroutine model_generate_ufo (filename, os_data, ufo_path, &
       ufo_path_requested)
    type(string_t), intent(in) :: filename
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(out) :: ufo_path
    type(string_t), intent(in), optional :: ufo_path_requested
    type(string_t) :: model_name, omega_path, ufo_dir, ufo_init
    logical :: exist
    call get_model_name (filename, model_name)
    call msg_message ("Model: Generating model '" // char (model_name) &
         // "' from UFO sources")
    if (present (ufo_path_requested)) then
       call msg_message ("Model: Searching for UFO sources in '" &
            // char (ufo_path_requested) // "'")
       ufo_path = ufo_path_requested
       ufo_dir = ufo_path_requested // "/" // model_name
       ufo_init = ufo_dir // "/" // "__init__.py"
       inquire (file = char (ufo_init), exist = exist)
    else
       call msg_message ("Model: Searching for UFO sources in &
            &working directory")
       ufo_path = "."
       ufo_dir = ufo_path // "/" // model_name
       ufo_init = ufo_dir // "/" // "__init__.py"
       inquire (file = char (ufo_init), exist = exist)
       if (.not. exist) then
          ufo_path = char (os_data%whizard_modelpath_ufo)
          ufo_dir = ufo_path // "/" // model_name
          ufo_init = ufo_dir // "/" // "__init__.py"
          call msg_message ("Model: Searching for UFO sources in '" &
               // char (os_data%whizard_modelpath_ufo) // "'")
          inquire (file = char (ufo_init), exist = exist)
       end if
    end if
    if (exist) then
       call msg_message ("Model: Found UFO sources for model '" &
            // char (model_name) // "'")
    else
       call msg_fatal ("Model: UFO sources for model '" &
            // char (model_name) // "' not found")
    end if
    omega_path = os_data%whizard_omega_binpath // "/omega_UFO.opt"
    call os_system_call (omega_path &
         // " -model:UFO_dir " // ufo_dir &
         // " -model:exec" &
         // " -model:write_WHIZARD" &
         // " > " // filename)
    inquire (file = char (filename), exist = exist)
    if (exist) then
        call msg_message ("Model: Model file '" // char (filename) //&
             "' generated")
    else
        call msg_fatal ("Model: Model file '" // char (filename) &
             // "' could not be generated")
    end if
  contains
    subroutine get_model_name (filename, model_name)
       type(string_t), intent(in) :: filename
       type(string_t), intent(out) :: model_name
       type(string_t) :: string
       string = filename
       call split (string, model_name, ".")
    end subroutine get_model_name
  end subroutine model_generate_ufo

  module function model_has_schemes (model) result (flag)
    logical :: flag
    class(model_t), intent(in) :: model
    flag = allocated (model%schemes)
  end function model_has_schemes

  module subroutine model_enable_schemes (model, scheme)
    class(model_t), intent(inout) :: model
    type(string_t), dimension(:), intent(in) :: scheme
    allocate (model%schemes (size (scheme)), source = scheme)
  end subroutine model_enable_schemes

  module subroutine model_set_scheme (model, scheme)
    class(model_t), intent(inout) :: model
    type(string_t), intent(in), optional :: scheme
    logical :: ok
    integer :: i
    if (model%has_schemes ()) then
       if (present (scheme)) then
          ok = .false.
          CHECK_SCHEME: do i = 1, size (model%schemes)
             if (scheme == model%schemes(i)) then
                allocate (model%selected_scheme, source = scheme)
                call model%set_scheme_num (i)
                ok = .true.
                exit CHECK_SCHEME
             end if
          end do CHECK_SCHEME
          if (.not. ok) then
             call msg_fatal &
                  ("Model '" // char (model%get_name ()) &
                  // "': scheme '" // char (scheme) // "' not supported")
          end if
       else
          allocate (model%selected_scheme, source = model%schemes(1))
          call model%set_scheme_num (1)
       end if
    else
       if (present (scheme)) then
          call msg_error &
                  ("Model '" // char (model%get_name ()) &
                  // "' does not support schemes")
       end if
    end if
  end subroutine model_set_scheme

  module function model_get_scheme (model) result (scheme)
    class(model_t), intent(in) :: model
    type(string_t) :: scheme
    if (allocated (model%selected_scheme)) then
       scheme = model%selected_scheme
    else
       scheme = ""
    end if
  end function model_get_scheme

  module function model_matches &
       (model, name, scheme, ufo, ufo_path) result (flag)
    logical :: flag
    class(model_t), intent(in) :: model
    type(string_t), intent(in) :: name
    type(string_t), intent(in), optional :: scheme
    logical, intent(in), optional :: ufo
    type(string_t), intent(in), optional :: ufo_path
    logical :: ufo_model
    ufo_model = .false.;  if (present (ufo))  ufo_model = ufo
    if (name /= model%get_name ()) then
       flag = .false.
    else if (ufo_model .neqv. model%is_ufo_model ()) then
       flag = .false.
    else if (ufo_model) then
       if (present (ufo_path)) then
          flag = model%get_ufo_path () == ufo_path
       else
          flag = .true.
       end if
    else if (model%has_schemes ()) then
       if (present (scheme)) then
          flag = model%get_scheme () == scheme
       else
          flag = model%get_scheme_num () == 1
       end if
    else if (present (scheme)) then
       flag = .false.
    else
       flag = .true.
    end if
  end function model_matches

  module function model_supports_custom_slha (model) result (flag)
    class(model_t), intent(in) :: model
    logical :: flag
    
    flag = model%n_slha_block > 0
    
  end function model_supports_custom_slha

  module subroutine model_get_custom_slha_blocks (model, block_name)
    class(model_t), intent(in) :: model
    type(string_t), dimension(:), allocatable :: block_name
    
    integer :: i
    
    allocate (block_name (model%n_slha_block))
    do i = 1, size (block_name)
       block_name(i) = model%slha_block(i)%name
    end do

  end subroutine model_get_custom_slha_blocks
  
  subroutine model_record_slha_block_entry (model, i_par, node)
    class(model_t), intent(inout) :: model
    integer, intent(in) :: i_par
    type(parse_node_t), intent(in), target :: node

    type(parse_node_t), pointer :: node_block_name, node_index
    type(string_t) :: block_name
    integer :: n_index, i, i_block
    integer, dimension(:), allocatable :: block_index

    node_block_name => node%get_sub_ptr (2)
    select case (char (node_block_name%get_rule_key ()))
    case ("block_name")
       block_name = node_block_name%get_string ()
    case ("QNUMBERS")
       block_name = "QNUMBERS"
    case default
       block_name = node_block_name%get_string ()
    end select
    n_index = node%get_n_sub () - 2
    allocate (block_index (n_index))
    node_index => node_block_name%get_next_ptr ()
    do i = 1, n_index
       block_index(i) = node_index%get_integer ()
       node_index => node_index%get_next_ptr ()
    end do
    i_block = 0
    FIND_BLOCK: do i = 1, model%n_slha_block
       if (model%slha_block(i)%name == block_name) then
          i_block = i
          exit FIND_BLOCK
       end if
    end do FIND_BLOCK
    if (i_block == 0) then
       call model_add_slha_block (model, block_name)
       i_block = model%n_slha_block
    end if
    associate (b => model%slha_block(i_block))
      call add_slha_block_entry (b, block_index, i_par)
    end associate
    model%par(i_par)%block_name = block_name
    model%par(i_par)%block_index = block_index
  end subroutine model_record_slha_block_entry
  
  subroutine model_add_slha_block (model, block_name)
    class(model_t), intent(inout) :: model
    type(string_t), intent(in) :: block_name
    
    if (.not. allocated (model%slha_block))  allocate (model%slha_block (1))
    if (model%n_slha_block == size (model%slha_block))  call grow
    model%n_slha_block = model%n_slha_block + 1
    associate (b => model%slha_block(model%n_slha_block))
      b%name = block_name
      allocate (b%entry (1))
    end associate
    
  contains
    
    subroutine grow
      type(slha_block_t), dimension(:), allocatable :: b_tmp
      call move_alloc (model%slha_block, b_tmp)
      allocate (model%slha_block (2 * size (b_tmp)))
      model%slha_block(:size (b_tmp)) = b_tmp(:)
    end subroutine grow
     
  end subroutine model_add_slha_block
  
  subroutine add_slha_block_entry (b, block_index, i_par)
    type(slha_block_t), intent(inout) :: b
    integer, dimension(:), intent(in) :: block_index
    integer, intent(in) :: i_par
    
    if (b%n_entry == size (b%entry))  call grow
    b%n_entry = b%n_entry + 1
    associate (entry => b%entry(b%n_entry))
      entry%block_index = block_index
      entry%i_par = i_par
    end associate
    
  contains
    
    subroutine grow
      type(slha_entry_t), dimension(:), allocatable :: entry_tmp
      call move_alloc (b%entry, entry_tmp)
      allocate (b%entry (2 * size (entry_tmp)))
      b%entry(:size (entry_tmp)) = entry_tmp(:)
    end subroutine grow
    
  end subroutine add_slha_block_entry

  module subroutine model_slha_lookup &
       (model, block_name, block_index, par_data)
    class(model_t), intent(in) :: model
    type(string_t), intent(in) :: block_name
    integer, dimension(:), intent(in) :: block_index
    class(modelpar_data_t), pointer, intent(out) :: par_data
    
    integer :: i, j
    
    par_data => null ()
    if (allocated (model%slha_block)) then
       do i = 1, model%n_slha_block
          associate (block => model%slha_block(i))
            if (block%name == block_name) then
               do j = 1, block%n_entry
                  associate (entry => block%entry(j))
                    if (size (entry%block_index) == size (block_index)) then
                       if (all (entry%block_index == block_index)) then
                          par_data => model%par(entry%i_par)%data
                          return
                       end if
                    end if
                  end associate
               end do
            end if
          end associate
       end do
    end if
    
  end subroutine model_slha_lookup
    
  module subroutine model_slha_set_par (model, block_name, block_index, value)
    class(model_t), intent(inout) :: model
    type(string_t), intent(in) :: block_name
    integer, dimension(:), intent(in) :: block_index
    real(default), intent(in) :: value
    
    class(modelpar_data_t), pointer :: par_data
    
    call model%slha_lookup (block_name, block_index, par_data)
    if (associated (par_data)) then
       par_data = value
    end if

  end subroutine model_slha_set_par
  
  subroutine define_model_file_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "SEQ model_def = model_name_def " // &
         "scheme_header parameters external_pars particles vertices")
    call ifile_append (ifile, "SEQ model_name_def = model model_name")
    call ifile_append (ifile, "KEY model")
    call ifile_append (ifile, "QUO model_name = '""'...'""'")
    call ifile_append (ifile, "SEQ scheme_header = scheme_decl?")
    call ifile_append (ifile, "SEQ scheme_decl = schemes '=' scheme_list")
    call ifile_append (ifile, "KEY schemes")
    call ifile_append (ifile, "LIS scheme_list = scheme_name+")
    call ifile_append (ifile, "QUO scheme_name = '""'...'""'")
    call ifile_append (ifile, "SEQ parameters = generic_par_def*")
    call ifile_append (ifile, "ALT generic_par_def = &
         &parameter_def | derived_def | unused_def | scheme_block")
    call ifile_append (ifile, "SEQ parameter_def = parameter par_name " // &
         "'=' any_real_value slha_annotation?")
    call ifile_append (ifile, "ALT any_real_value = " &
         // "neg_real_value | pos_real_value | real_value")
    call ifile_append (ifile, "SEQ neg_real_value = '-' real_value")
    call ifile_append (ifile, "SEQ pos_real_value = '+' real_value")
    call ifile_append (ifile, "KEY parameter")
    call ifile_append (ifile, "IDE par_name")
    ! call ifile_append (ifile, "KEY '='")          !!! Key already exists
    call ifile_append (ifile, "SEQ slha_annotation = " // &
         "slha_entry slha_block_name slha_entry_index*")
    call ifile_append (ifile, "KEY slha_entry")
    call ifile_append (ifile, "IDE slha_block_name")
    call ifile_append (ifile, "INT slha_entry_index")
    call ifile_append (ifile, "SEQ derived_def = derived par_name " // &
         "'=' expr")
    call ifile_append (ifile, "KEY derived")
    call ifile_append (ifile, "SEQ unused_def = unused par_name")
    call ifile_append (ifile, "KEY unused")
    call ifile_append (ifile, "SEQ external_pars = external_def*")
    call ifile_append (ifile, "SEQ external_def = external par_name")
    call ifile_append (ifile, "KEY external")
    call ifile_append (ifile, "SEQ scheme_block = &
         &scheme_block_beg scheme_block_body scheme_block_end")
    call ifile_append (ifile, "SEQ scheme_block_beg = select scheme")
    call ifile_append (ifile, "SEQ scheme_block_body = scheme_block_case*")
    call ifile_append (ifile, "SEQ scheme_block_case = &
         &scheme scheme_id parameters")
    call ifile_append (ifile, "ALT scheme_id = scheme_list | other")
    call ifile_append (ifile, "SEQ scheme_block_end = end select")
    call ifile_append (ifile, "KEY select")
    call ifile_append (ifile, "KEY scheme")
    call ifile_append (ifile, "KEY other")
    call ifile_append (ifile, "KEY end")
    call ifile_append (ifile, "SEQ particles = particle_def*")
    call ifile_append (ifile, "SEQ particle_def = particle name_def " // &
         "prt_pdg prt_details")
    call ifile_append (ifile, "KEY particle")
    call ifile_append (ifile, "SEQ prt_pdg = signed_int")
    call ifile_append (ifile, "ALT prt_details = prt_src | prt_properties")
    call ifile_append (ifile, "SEQ prt_src = like name_def prt_properties")
    call ifile_append (ifile, "KEY like")
    call ifile_append (ifile, "SEQ prt_properties = prt_property*")
    call ifile_append (ifile, "ALT prt_property = " // &
         "parton | invisible | gauge | left | right | " // &
         "prt_name | prt_anti | prt_tex_name | prt_tex_anti | " // &
         "prt_spin | prt_isospin | prt_charge | " // &
         "prt_color | prt_mass | prt_width")
    call ifile_append (ifile, "KEY parton")
    call ifile_append (ifile, "KEY invisible")
    call ifile_append (ifile, "KEY gauge")
    call ifile_append (ifile, "KEY left")
    call ifile_append (ifile, "KEY right")
    call ifile_append (ifile, "SEQ prt_name = name name_def+")
    call ifile_append (ifile, "SEQ prt_anti = anti name_def+")
    call ifile_append (ifile, "SEQ prt_tex_name = tex_name name_def")
    call ifile_append (ifile, "SEQ prt_tex_anti = tex_anti name_def")
    call ifile_append (ifile, "KEY name")
    call ifile_append (ifile, "KEY anti")
    call ifile_append (ifile, "KEY tex_name")
    call ifile_append (ifile, "KEY tex_anti")
    call ifile_append (ifile, "ALT name_def = name_string | name_id")
    call ifile_append (ifile, "QUO name_string = '""'...'""'")
    call ifile_append (ifile, "IDE name_id")
    call ifile_append (ifile, "SEQ prt_spin = spin frac")
    call ifile_append (ifile, "KEY spin")
    call ifile_append (ifile, "SEQ prt_isospin = isospin frac")
    call ifile_append (ifile, "KEY isospin")
    call ifile_append (ifile, "SEQ prt_charge = charge frac")
    call ifile_append (ifile, "KEY charge")
    call ifile_append (ifile, "SEQ prt_color = color integer_literal")
    call ifile_append (ifile, "KEY color")
    call ifile_append (ifile, "SEQ prt_mass = mass par_name")
    call ifile_append (ifile, "KEY mass")
    call ifile_append (ifile, "SEQ prt_width = width par_name")
    call ifile_append (ifile, "KEY width")
    call ifile_append (ifile, "SEQ vertices = vertex_def*")
    call ifile_append (ifile, "SEQ vertex_def = vertex name_def+")
    call ifile_append (ifile, "KEY vertex")
    call define_expr_syntax (ifile, particles=.false., analysis=.false.)
  end subroutine define_model_file_syntax

  module subroutine syntax_model_file_init ()
    type(ifile_t) :: ifile
    call define_model_file_syntax (ifile)
    call syntax_init (syntax_model_file, ifile)
    call ifile_final (ifile)
  end subroutine syntax_model_file_init

  subroutine lexer_init_model_file (lexer)
    type(lexer_t), intent(out) :: lexer
    call lexer_init (lexer, &
         comment_chars = "#!", &
         quote_chars = '"{', &
         quote_match = '"}', &
         single_chars = ":(),", &
         special_class = [ "+-*/^", "<>=  " ] , &
         keyword_list = syntax_get_keyword_list_ptr (syntax_model_file))
  end subroutine lexer_init_model_file

  module subroutine syntax_model_file_final ()
    call syntax_final (syntax_model_file)
  end subroutine syntax_model_file_final

  module subroutine syntax_model_file_write (unit)
    integer, intent(in), optional :: unit
    call syntax_write (syntax_model_file, unit)
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
    type(string_t) :: file
    type(stream_t), target :: stream
    type(lexer_t) :: lexer
    integer :: unit
    character(32) :: model_md5sum
    type(parse_node_t), pointer :: nd_model_def, nd_model_name_def
    type(parse_node_t), pointer :: nd_schemes, nd_scheme_decl
    type(parse_node_t), pointer :: nd_parameters
    type(parse_node_t), pointer :: nd_external_pars
    type(parse_node_t), pointer :: nd_particles, nd_vertices
    type(string_t) :: model_name, lib_name
    integer :: n_parblock, n_par, i_par, n_ext, n_prt, n_vtx
    type(parse_node_t), pointer :: nd_par_def
    type(parse_node_t), pointer :: nd_ext_def
    type(parse_node_t), pointer :: nd_prt
    type(parse_node_t), pointer :: nd_vtx
    logical :: ufo_model, model_exist, rebuild
    ufo_model = .false.;  if (present (ufo))  ufo_model = ufo
    rebuild = .true.; if (present (rebuild_mdl))  rebuild = rebuild_mdl
    file = filename
    inquire (file=char(file), exist=model_exist)
    if ((.not. model_exist) .and. (.not. os_data%use_testfiles)) then
       file = os_data%whizard_modelpath_local // "/" // filename
       inquire (file = char (file), exist = model_exist)
    end if
    if (.not. model_exist) then
       file = os_data%whizard_modelpath // "/" // filename
       inquire (file = char (file), exist = model_exist)
    end if
    if (ufo_model .and. rebuild) then
       file = filename
       call model_generate_ufo (filename, os_data, model%ufo_path, &
            ufo_path_requested=ufo_path_requested)
       inquire (file = char (file), exist = model_exist)
    end if
    if (.not. model_exist) then
       call msg_fatal ("Model file '" // char (filename) // "' not found")
       if (present (exist))  exist = .false.
       return
    end if
    if (present (exist))  exist = .true.

    if (logging) call msg_message ("Reading model file '" // char (file) // "'")

    unit = free_unit ()
    open (file=char(file), unit=unit, action="read", status="old")
    model_md5sum = md5sum (unit)
    close (unit)

    call lexer_init_model_file (lexer)
    call stream_init (stream, char (file))
    call lexer_assign_stream (lexer, stream)
    call parse_tree_init (model%parse_tree, syntax_model_file, lexer)
    call stream_final (stream)
    call lexer_final (lexer)

    nd_model_def => model%parse_tree%get_root_ptr ()
    nd_model_name_def => parse_node_get_sub_ptr (nd_model_def)
    model_name = parse_node_get_string &
         (parse_node_get_sub_ptr (nd_model_name_def, 2))
    nd_schemes => nd_model_name_def%get_next_ptr ()

    call find_block &
         ("scheme_header", nd_schemes, nd_scheme_decl, nd_next=nd_parameters)
    call find_block &
         ("parameters", nd_parameters, nd_par_def, n_parblock, nd_external_pars)
    call find_block &
         ("external_pars", nd_external_pars, nd_ext_def, n_ext, nd_particles)
    call find_block &
         ("particles", nd_particles, nd_prt, n_prt, nd_vertices)
    call find_block &
         ("vertices", nd_vertices, nd_vtx, n_vtx)
    
    if (associated (nd_external_pars)) then
       lib_name = "external." // model_name
    else
       lib_name = ""
    end if

    if (associated (nd_scheme_decl)) then
       call handle_schemes (nd_scheme_decl, scheme)
    end if

    n_par = 0
    call count_parameters (nd_par_def, n_parblock, n_par)
    
    call model%init &
         (model_name, lib_name, os_data, n_par + n_ext, n_prt, n_vtx, ufo)
    model%md5sum = model_md5sum

    if (associated (nd_par_def)) then
       i_par = 0
       call handle_parameters (nd_par_def, n_parblock, i_par)
    end if
    if (associated (nd_ext_def)) then
       call handle_external (nd_ext_def, n_par, n_ext)
    end if
    call model%update_parameters ()
    if (associated (nd_prt)) then
       call handle_fields (nd_prt, n_prt)
    end if
    if (associated (nd_vtx)) then
       call handle_vertices (nd_vtx, n_vtx)
    end if

    call model%freeze_vertices ()
    call model%append_field_vars ()

  contains

    subroutine find_block (key, nd, nd_item, n_item, nd_next)
      character(*), intent(in) :: key
      type(parse_node_t), pointer, intent(inout) :: nd
      type(parse_node_t), pointer, intent(out) :: nd_item
      integer, intent(out), optional :: n_item
      type(parse_node_t), pointer, intent(out), optional :: nd_next
      if (associated (nd)) then
         if (nd%get_rule_key () == key) then
            nd_item => nd%get_sub_ptr ()
            if (present (n_item))  n_item = nd%get_n_sub ()
            if (present (nd_next))  nd_next => nd%get_next_ptr ()
         else
            nd_item => null ()
            if (present (n_item))  n_item = 0
            if (present (nd_next))  nd_next => nd
            nd => null ()
         end if
      else
         nd_item => null ()
         if (present (n_item))  n_item = 0
         if (present (nd_next))  nd_next => null ()
      end if
    end subroutine find_block

    subroutine handle_schemes (nd_scheme_decl, scheme)
      type(parse_node_t), pointer, intent(in) :: nd_scheme_decl
      type(string_t), intent(in), optional :: scheme
      type(parse_node_t), pointer :: nd_list, nd_entry
      type(string_t), dimension(:), allocatable :: schemes
      integer :: i, n_schemes
      nd_list => nd_scheme_decl%get_sub_ptr (3)
      nd_entry => nd_list%get_sub_ptr ()
      n_schemes = nd_list%get_n_sub ()
      allocate (schemes (n_schemes))
      do i = 1, n_schemes
         schemes(i) = nd_entry%get_string ()
         nd_entry => nd_entry%get_next_ptr ()
      end do
      if (present (scheme)) then
         do i = 1, n_schemes
            if (schemes(i) == scheme)  goto 10   ! block exit
         end do
         call msg_fatal ("Scheme '" // char (scheme) &
              // "' is not supported by model '" // char (model_name) // "'")
      end if
10    continue
      call model%enable_schemes (schemes)
      call model%set_scheme (scheme)
    end subroutine handle_schemes

    subroutine select_scheme (nd_scheme_block, n_parblock_sub, nd_par_def)
      type(parse_node_t), pointer, intent(in) :: nd_scheme_block
      integer, intent(out) :: n_parblock_sub
      type(parse_node_t), pointer, intent(out) :: nd_par_def
      type(parse_node_t), pointer :: nd_scheme_body
      type(parse_node_t), pointer :: nd_scheme_case, nd_scheme_id, nd_scheme
      type(string_t) :: scheme
      integer :: n_cases, i
      scheme = model%get_scheme ()
      nd_scheme_body => nd_scheme_block%get_sub_ptr (2)
      nd_parameters => null ()
      select case (char (nd_scheme_body%get_rule_key ()))
      case ("scheme_block_body")
         n_cases = nd_scheme_body%get_n_sub ()
         FIND_SCHEME: do i = 1, n_cases
            nd_scheme_case => nd_scheme_body%get_sub_ptr (i)
            nd_scheme_id => nd_scheme_case%get_sub_ptr (2)
            select case (char (nd_scheme_id%get_rule_key ()))
            case ("scheme_list")
               nd_scheme => nd_scheme_id%get_sub_ptr ()
               do while (associated (nd_scheme))
                  if (scheme == nd_scheme%get_string ()) then
                     nd_parameters => nd_scheme_id%get_next_ptr ()
                     exit FIND_SCHEME
                  end if
                  nd_scheme => nd_scheme%get_next_ptr ()
               end do
            case ("other")
               nd_parameters => nd_scheme_id%get_next_ptr ()
               exit FIND_SCHEME
            case default
               print *, "'", char (nd_scheme_id%get_rule_key ()), "'"
               call msg_bug ("Model read: impossible scheme rule")
            end select
         end do FIND_SCHEME
      end select
      if (associated (nd_parameters)) then
         select case (char (nd_parameters%get_rule_key ()))
         case ("parameters")
            n_parblock_sub = nd_parameters%get_n_sub ()
            if (n_parblock_sub > 0) then
               nd_par_def => nd_parameters%get_sub_ptr ()
            else
               nd_par_def => null ()
            end if
         case default
            n_parblock_sub = 0
            nd_par_def => null ()
         end select
      else
         n_parblock_sub = 0
         nd_par_def => null ()
      end if
    end subroutine select_scheme

    recursive subroutine count_parameters (nd_par_def_in, n_parblock, n_par)
      type(parse_node_t), pointer, intent(in) :: nd_par_def_in
      integer, intent(in) :: n_parblock
      integer, intent(inout) :: n_par
      type(parse_node_t), pointer :: nd_par_def, nd_par_key
      type(parse_node_t), pointer :: nd_par_def_sub
      integer :: n_parblock_sub
      integer :: i
      nd_par_def => nd_par_def_in
      do i = 1, n_parblock
         nd_par_key => nd_par_def%get_sub_ptr ()
         select case (char (nd_par_key%get_rule_key ()))
         case ("parameter", "derived", "unused")
            n_par = n_par + 1
         case ("scheme_block_beg")
            call select_scheme (nd_par_def, n_parblock_sub, nd_par_def_sub)
            if (n_parblock_sub > 0) then
               call count_parameters (nd_par_def_sub, n_parblock_sub, n_par)
            end if
         case default
            print *, "'", char (nd_par_key%get_rule_key ()), "'"
            call msg_bug ("Model read: impossible parameter rule")
         end select
         nd_par_def => parse_node_get_next_ptr (nd_par_def)
      end do
    end subroutine count_parameters

    recursive subroutine handle_parameters (nd_par_def_in, n_parblock, i_par)
      type(parse_node_t), pointer, intent(in) :: nd_par_def_in
      integer, intent(in) :: n_parblock
      integer, intent(inout) :: i_par
      type(parse_node_t), pointer :: nd_par_def, nd_par_key
      type(parse_node_t), pointer :: nd_par_def_sub
      integer :: n_parblock_sub
      integer :: i
      nd_par_def => nd_par_def_in
      do i = 1, n_parblock
         nd_par_key => nd_par_def%get_sub_ptr ()
         select case (char (nd_par_key%get_rule_key ()))
         case ("parameter")
            i_par = i_par + 1
            call model%read_parameter (i_par, nd_par_def)
         case ("derived")
            i_par = i_par + 1
            call model%read_derived (i_par, nd_par_def)
         case ("unused")
            i_par = i_par + 1
            call model%read_unused (i_par, nd_par_def)
         case ("scheme_block_beg")
            call select_scheme (nd_par_def, n_parblock_sub, nd_par_def_sub)
            if (n_parblock_sub > 0) then
               call handle_parameters (nd_par_def_sub, n_parblock_sub, i_par)
            end if
         end select
         nd_par_def => parse_node_get_next_ptr (nd_par_def)
      end do
    end subroutine handle_parameters

    subroutine handle_external (nd_ext_def, n_par, n_ext)
      type(parse_node_t), pointer, intent(inout) :: nd_ext_def
      integer, intent(in) :: n_par, n_ext
      integer :: i
      do i = n_par + 1, n_par + n_ext
         call model%read_external (i, nd_ext_def)
         nd_ext_def => parse_node_get_next_ptr (nd_ext_def)
      end do
!     real(c_default_float), dimension(:), allocatable :: par
!       if (associated (model%init_external_parameters)) then
!          allocate (par (model%get_n_real ()))
!          call model%real_parameters_to_c_array (par)
!          call model%init_external_parameters (par)
!          call model%real_parameters_from_c_array (par)
!       end if
    end subroutine handle_external

    subroutine handle_fields (nd_prt, n_prt)
      type(parse_node_t), pointer, intent(inout) :: nd_prt
      integer, intent(in) :: n_prt
      integer :: i
      do i = 1, n_prt
         call model%read_field (i, nd_prt)
         nd_prt => parse_node_get_next_ptr (nd_prt)
      end do
    end subroutine handle_fields

    subroutine handle_vertices (nd_vtx, n_vtx)
      type(parse_node_t), pointer, intent(inout) :: nd_vtx
      integer, intent(in) :: n_vtx
      integer :: i
      do i = 1, n_vtx
         call model%read_vertex (i, nd_vtx)
         nd_vtx => parse_node_get_next_ptr (nd_vtx)
      end do
    end subroutine handle_vertices

  end subroutine model_read

  module subroutine model_read_parameter (model, i, node)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(parse_node_t), intent(in), target :: node
    type(parse_node_t), pointer :: node_name, node_val, node_slha_entry
    type(string_t) :: name
    node_name => parse_node_get_sub_ptr (node, 2)
    name = parse_node_get_string (node_name)
    node_val => parse_node_get_next_ptr (node_name, 2)
    call model%set_parameter_parse_node (i, name, node_val, constant=.true.)
    node_slha_entry => parse_node_get_next_ptr (node_val)
    if (associated (node_slha_entry)) then
       call model_record_slha_block_entry (model, i, node_slha_entry)
    end if
  end subroutine model_read_parameter

  module subroutine model_read_derived (model, i, node)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(parse_node_t), intent(in), target :: node
    type(string_t) :: name
    type(parse_node_t), pointer :: pn_expr
    name = parse_node_get_string (parse_node_get_sub_ptr (node, 2))
    pn_expr => parse_node_get_sub_ptr (node, 4)
    call model%set_parameter_parse_node (i, name, pn_expr, constant=.false.)
  end subroutine model_read_derived

  module subroutine model_read_external (model, i, node)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(parse_node_t), intent(in), target :: node
    type(string_t) :: name
    name = parse_node_get_string (parse_node_get_sub_ptr (node, 2))
    call model%set_parameter_external (i, name)
  end subroutine model_read_external

  module subroutine model_read_unused (model, i, node)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(parse_node_t), intent(in), target :: node
    type(string_t) :: name
    name = parse_node_get_string (parse_node_get_sub_ptr (node, 2))
    call model%set_parameter_unused (i, name)
  end subroutine model_read_unused

  module subroutine model_read_field (model, i, node)
    class(model_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(parse_node_t), intent(in) :: node
    type(parse_node_t), pointer :: nd_src, nd_props, nd_prop
    type(string_t) :: longname
    integer :: pdg
    type(string_t) :: name_src
    type(string_t), dimension(:), allocatable :: name
    type(field_data_t), pointer :: field, field_src
    longname = parse_node_get_string (parse_node_get_sub_ptr (node, 2))
    pdg = read_frac (parse_node_get_sub_ptr (node, 3))
    field => model%get_field_ptr_by_index (i)
    call field%init (longname, pdg)
    nd_src => parse_node_get_sub_ptr (node, 4)
    if (associated (nd_src)) then
       if (parse_node_get_rule_key (nd_src) == "prt_src") then
          name_src = parse_node_get_string (parse_node_get_sub_ptr (nd_src, 2))
          field_src => model%get_field_ptr (name_src, check=.true.)
          call field%copy_from (field_src)
          nd_props => parse_node_get_sub_ptr (nd_src, 3)
       else
          nd_props => nd_src
       end if
       nd_prop => parse_node_get_sub_ptr (nd_props)
       do while (associated (nd_prop))
          select case (char (parse_node_get_rule_key (nd_prop)))
          case ("invisible")
             call field%set (is_visible=.false.)
          case ("parton")
             call field%set (is_parton=.true.)
          case ("gauge")
             call field%set (is_gauge=.true.)
          case ("left")
             call field%set (is_left_handed=.true.)
          case ("right")
             call field%set (is_right_handed=.true.)
          case ("prt_name")
             call read_names (nd_prop, name)
             call field%set (name=name)
          case ("prt_anti")
             call read_names (nd_prop, name)
             call field%set (anti=name)
          case ("prt_tex_name")
             call field%set ( &
                  tex_name = parse_node_get_string &
                  (parse_node_get_sub_ptr (nd_prop, 2)))
          case ("prt_tex_anti")
             call field%set ( &
                  tex_anti = parse_node_get_string &
                  (parse_node_get_sub_ptr (nd_prop, 2)))
          case ("prt_spin")
             call field%set ( &
                  spin_type = read_frac &
                  (parse_node_get_sub_ptr (nd_prop, 2), 2))
          case ("prt_isospin")
             call field%set ( &
                  isospin_type = read_frac &
                  (parse_node_get_sub_ptr (nd_prop, 2), 2))
          case ("prt_charge")
             call field%set ( &
                  charge_type = read_frac &
                  (parse_node_get_sub_ptr (nd_prop, 2), 3))
          case ("prt_color")
             call field%set ( &
                  color_type = parse_node_get_integer &
                  (parse_node_get_sub_ptr (nd_prop, 2)))
          case ("prt_mass")
             call field%set ( &
                  mass_data = model%get_par_data_ptr &
                  (parse_node_get_string &
                  (parse_node_get_sub_ptr (nd_prop, 2))))
          case ("prt_width")
             call field%set ( &
                  width_data = model%get_par_data_ptr &
                  (parse_node_get_string &
                  (parse_node_get_sub_ptr (nd_prop, 2))))
          case default
             call msg_bug (" Unknown particle property '" &
                  // char (parse_node_get_rule_key (nd_prop)) // "'")
          end select
          if (allocated (name))  deallocate (name)
          nd_prop => parse_node_get_next_ptr (nd_prop)
       end do
    end if
    call field%freeze ()
  end subroutine model_read_field

  module subroutine model_read_vertex (model, i, node)
    class(model_t), intent(inout) :: model
    integer, intent(in) :: i
    type(parse_node_t), intent(in) :: node
    type(string_t), dimension(:), allocatable :: name
    call read_names (node, name)
    call model%set_vertex (i, name)
  end subroutine model_read_vertex

  subroutine read_names (node, name)
    type(parse_node_t), intent(in) :: node
    type(string_t), dimension(:), allocatable, intent(inout) :: name
    type(parse_node_t), pointer :: nd_name
    integer :: n_names, i
    n_names = parse_node_get_n_sub (node) - 1
    allocate (name (n_names))
    nd_name => parse_node_get_sub_ptr (node, 2)
    do i = 1, n_names
       name(i) = parse_node_get_string (nd_name)
       nd_name => parse_node_get_next_ptr (nd_name)
    end do
  end subroutine read_names

  function read_frac (nd_frac, base) result (qn_type)
    integer :: qn_type
    type(parse_node_t), intent(in) :: nd_frac
    integer, intent(in), optional :: base
    type(parse_node_t), pointer :: nd_num, nd_den
    integer :: num, den
    nd_num => parse_node_get_sub_ptr (nd_frac)
    nd_den => parse_node_get_next_ptr (nd_num)
    select case (char (parse_node_get_rule_key (nd_num)))
    case ("integer_literal")
       num = parse_node_get_integer (nd_num)
    case ("neg_int")
       num = - parse_node_get_integer (parse_node_get_sub_ptr (nd_num, 2))
    case ("pos_int")
       num = parse_node_get_integer (parse_node_get_sub_ptr (nd_num, 2))
    case default
       call parse_tree_bug (nd_num, "int|neg_int|pos_int")
    end select
    if (associated (nd_den)) then
       den = parse_node_get_integer (parse_node_get_sub_ptr (nd_den, 2))
    else
       den = 1
    end if
    if (present (base)) then
       if (den == 1) then
          qn_type = sign (1 + abs (num) * base, num)
       else if (den == base) then
          qn_type = sign (abs (num) + 1, num)
       else
          call parse_node_write_rec (nd_frac)
          call msg_fatal (" Fractional quantum number: wrong denominator")
       end if
    else
       if (den == 1) then
          qn_type = num
       else
          call parse_node_write_rec (nd_frac)
          call msg_fatal (" Wrong type: no fraction expected")
       end if
    end if
  end function read_frac

  module subroutine model_append_field_vars (model)
    class(model_t), intent(inout) :: model
    type(pdg_array_t) :: aval
    type(field_data_t), dimension(:), pointer :: field_array
    type(field_data_t), pointer :: field
    type(string_t) :: name
    type(string_t), dimension(:), allocatable :: name_array
    integer, dimension(:), allocatable :: pdg
    logical, dimension(:), allocatable :: mask
    integer :: i, j
    field_array => model%get_field_array_ptr ()
    aval = UNDEFINED
    call model%var_list%append_pdg_array (var_str ("particle"), &
         aval, locked = .true., intrinsic=.true.)
    do i = 1, size (field_array)
       aval = field_array(i)%get_pdg ()
       name = field_array(i)%get_longname ()
       call model%var_list%append_pdg_array &
            (name, aval, locked=.true., intrinsic=.true.)
       call field_array(i)%get_name_array (.false., name_array)
       do j = 1, size (name_array)
          call model%var_list%append_pdg_array (name_array(j), &
               aval, locked=.true., intrinsic=.true.)
       end do
       model%max_field_name_length = &
            max (model%max_field_name_length, len (name_array(1)))
       aval = - field_array(i)%get_pdg ()
       call field_array(i)%get_name_array (.true., name_array)
       do j = 1, size (name_array)
          call model%var_list%append_pdg_array (name_array(j), &
               aval, locked=.true., intrinsic=.true.)
       end do
       if (size (name_array) > 0) then
          model%max_field_name_length = &
               max (model%max_field_name_length, len (name_array(1)))
       end if
    end do
    call model%get_all_pdg (pdg)
    allocate (mask (size (pdg)))
    do i = 1, size (pdg)
       field => model%get_field_ptr (pdg(i))
       mask(i) = field%get_charge_type () /= 1
    end do
    aval = pack (pdg, mask)
    call model%var_list%append_pdg_array (var_str ("charged"), &
         aval, locked = .true., intrinsic=.true.)
    do i = 1, size (pdg)
       field => model%get_field_ptr (pdg(i))
       mask(i) = field%get_charge_type () == 1
    end do
    aval = pack (pdg, mask)
    call model%var_list%append_pdg_array (var_str ("neutral"), &
         aval, locked = .true., intrinsic=.true.)
    do i = 1, size (pdg)
       field => model%get_field_ptr (pdg(i))
       mask(i) = field%get_color_type () /= 1
    end do
    aval = pack (pdg, mask)
    call model%var_list%append_pdg_array (var_str ("colored"), &
         aval, locked = .true., intrinsic=.true.)
  end subroutine model_append_field_vars

  module subroutine create_test_model (model_name, test_model)
    type(string_t), intent(in) :: model_name
    type(model_t), intent(out), pointer :: test_model
    type(os_data_t) :: os_data
    type(model_list_t) :: model_list
    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
       (model_name, model_name // var_str (".mdl"), os_data, test_model)
  end subroutine create_test_model

  recursive module subroutine model_list_write &
       (object, unit, verbose, follow_link)
    class(model_list_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    logical, intent(in), optional :: follow_link
    type(model_entry_t), pointer :: current
    logical :: rec
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    rec = .true.;  if (present (follow_link))  rec = follow_link
    if (rec .and. associated (object%context)) then
       call object%context%write (unit, verbose, follow_link)
    end if
    current => object%first
    if (associated (current)) then
       do while (associated (current))
          call current%write (unit, verbose)
          current => current%next
          if (associated (current))  write (u, *)
       end do
    end if
  end subroutine model_list_write

  module subroutine model_list_link (model_list, context)
    class(model_list_t), intent(inout) :: model_list
    type(model_list_t), intent(in), target :: context
    model_list%context => context
  end subroutine model_list_link

  module subroutine model_list_import (model_list, current, model)
    class(model_list_t), intent(inout) :: model_list
    type(model_entry_t), pointer, intent(inout) :: current
    type(model_t), optional, pointer, intent(out) :: model
    if (associated (current)) then
       if (associated (model_list%first)) then
          model_list%last%next => current
       else
          model_list%first => current
       end if
       model_list%last => current
       if (present (model))  model => current%model_t
       current => null ()
    end if
  end subroutine model_list_import

  module subroutine model_list_add (model_list, &
       name, os_data, n_par, n_prt, n_vtx, model)
    class(model_list_t), intent(inout) :: model_list
    type(string_t), intent(in) :: name
    type(os_data_t), intent(in) :: os_data
    integer, intent(in) :: n_par, n_prt, n_vtx
    type(model_t), pointer :: model
    type(model_entry_t), pointer :: current
    if (model_list%model_exists (name, follow_link=.false.)) then
       model => null ()
    else
       allocate (current)
       call current%init (name, var_str (""), os_data, &
            n_par, n_prt, n_vtx)
       call model_list%import (current, model)
    end if
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
    class(model_list_t), pointer :: global_model_list
    type(model_entry_t), pointer :: current
    logical :: exist
    if (.not. model_list%model_exists (name, &
         scheme, ufo, ufo_path, follow_link=.true.)) then
       allocate (current)
       call current%read (filename, os_data, exist, &
            scheme=scheme, ufo=ufo, ufo_path_requested=ufo_path, &
            rebuild_mdl=rebuild_mdl)
       if (.not. exist)  return
       if (current%get_name () /= name) then
          call msg_fatal ("Model file '" // char (filename) // &
               "' contains model '" // char (current%get_name ()) // &
               "' instead of '" // char (name) // "'")
          call current%final ();  deallocate (current)
          return
       end if
       global_model_list => model_list
       do while (associated (global_model_list%context))
          global_model_list => global_model_list%context
       end do
       call global_model_list%import (current, model)
    else
       model => model_list%get_model_ptr (name, scheme, ufo, ufo_path)
    end if
  end subroutine model_list_read_model

  module subroutine model_list_append_copy (model_list, orig, model)
    class(model_list_t), intent(inout) :: model_list
    type(model_t), intent(in), target :: orig
    type(model_t), intent(out), pointer, optional :: model
    type(model_entry_t), pointer :: copy
    allocate (copy)
    call copy%init_instance (orig)
    call model_list%import (copy, model)
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
    type(model_entry_t), pointer :: current
    logical :: rec
    rec = .true.;  if (present (follow_link))  rec = follow_link
    current => model_list%first
    do while (associated (current))
       if (current%matches (name, scheme, ufo, ufo_path)) then
          exists = .true.
          return
       end if
       current => current%next
    end do
    if (rec .and. associated (model_list%context)) then
       exists = model_list%context%model_exists (name, &
            scheme, ufo, ufo_path, follow_link)
    else
       exists = .false.
    end if
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
    type(model_entry_t), pointer :: current
    logical :: rec
    rec = .true.;  if (present (follow_link))  rec = follow_link
    current => model_list%first
    do while (associated (current))
       if (current%matches (name, scheme, ufo, ufo_path)) then
          model => current%model_t
          return
       end if
       current => current%next
    end do
    if (rec .and. associated (model_list%context)) then
       model => model_list%context%get_model_ptr (name, &
            scheme, ufo, ufo_path, follow_link)
    else
       model => null ()
    end if
  end function model_list_get_model_ptr

  module subroutine model_list_final (model_list)
    class(model_list_t), intent(inout) :: model_list
    type(model_entry_t), pointer :: current
    model_list%last => null ()
    do while (associated (model_list%first))
       current => model_list%first
       model_list%first => model_list%first%next
       call current%final ()
       deallocate (current)
    end do
  end subroutine model_list_final

  module subroutine model_copy (model, orig)
    class(model_t), intent(out), target :: model
    type(model_t), intent(in) :: orig
    integer :: n_par, n_prt, n_vtx
    integer :: i
    n_par = orig%get_n_real ()
    n_prt = orig%get_n_field ()
    n_vtx = orig%get_n_vtx ()
    call model%basic_init (orig%get_name (), n_par, n_prt, n_vtx)
    if (allocated (orig%schemes)) then
       model%schemes = orig%schemes
       if (allocated (orig%selected_scheme)) then
          model%selected_scheme = orig%selected_scheme
          call model%set_scheme_num (orig%get_scheme_num ())
       end if
    end if
    if (allocated (orig%slha_block)) then
       model%slha_block = orig%slha_block
    end if
    model%md5sum = orig%md5sum
    model%ufo_model = orig%ufo_model
    model%ufo_path = orig%ufo_path
    if (allocated (orig%par)) then
       do i = 1, n_par
          call model%copy_parameter (i, orig%par(i))
       end do
    end if
    model%init_external_parameters => orig%init_external_parameters
    call model%model_data_t%copy_from (orig)
    model%max_par_name_length = orig%max_par_name_length
    call model%append_field_vars ()
  end subroutine model_copy


end submodule models_s

