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

submodule (rt_data) rt_data_s

  use io_units
  use format_utils, only: write_separator
  use format_defs, only: FMT_19, FMT_12
  use system_dependencies
  use diagnostics

  implicit none

contains

  module subroutine rt_parse_nodes_clear (rt_pn, name)
    class(rt_parse_nodes_t), intent(inout) :: rt_pn
    type(string_t), intent(in) :: name
    select case (char (name))
    case ("cuts")
       rt_pn%cuts_lexpr => null ()
    case ("scale")
       rt_pn%scale_expr => null ()
    case ("factorization_scale")
       rt_pn%fac_scale_expr => null ()
    case ("renormalization_scale")
       rt_pn%ren_scale_expr => null ()
    case ("weight")
       rt_pn%weight_expr => null ()
    case ("selection")
       rt_pn%selection_lexpr => null ()
    case ("reweight")
       rt_pn%reweight_expr => null ()
    case ("analysis")
       rt_pn%analysis_lexpr => null ()
    end select
  end subroutine rt_parse_nodes_clear

  module subroutine rt_parse_nodes_write (object, unit)
    class(rt_parse_nodes_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    call wrt ("Cuts", object%cuts_lexpr)
    call write_separator (u)
    call wrt ("Scale", object%scale_expr)
    call write_separator (u)
    call wrt ("Factorization scale", object%fac_scale_expr)
    call write_separator (u)
    call wrt ("Renormalization scale", object%ren_scale_expr)
    call write_separator (u)
    call wrt ("Weight", object%weight_expr)
    call write_separator (u, 2)
    call wrt ("Event selection", object%selection_lexpr)
    call write_separator (u)
    call wrt ("Event reweighting factor", object%reweight_expr)
    call write_separator (u)
    call wrt ("Event analysis", object%analysis_lexpr)
    if (allocated (object%alt_setup)) then
       call write_separator (u, 2)
       write (u, "(1x,A,':')")  "Alternative setups"
       do i = 1, size (object%alt_setup)
          call write_separator (u)
          call wrt ("Commands", object%alt_setup(i)%ptr)
       end do
    end if
  contains
    subroutine wrt (title, pn)
      character(*), intent(in) :: title
      type(parse_node_t), intent(in), pointer :: pn
      if (associated (pn)) then
         write (u, "(1x,A,':')")  title
         call write_separator (u)
         call parse_node_write_rec (pn, u)
      else
         write (u, "(1x,A,':',1x,A)")  title, "[undefined]"
      end if
    end subroutine wrt
  end subroutine rt_parse_nodes_write

  module subroutine rt_parse_nodes_show (rt_pn, name, unit)
    class(rt_parse_nodes_t), intent(in) :: rt_pn
    type(string_t), intent(in) :: name
    integer, intent(in), optional :: unit
    type(parse_node_t), pointer :: pn
    integer :: u
    u = given_output_unit (unit)
    select case (char (name))
    case ("cuts")
       pn => rt_pn%cuts_lexpr
    case ("scale")
       pn => rt_pn%scale_expr
    case ("factorization_scale")
       pn => rt_pn%fac_scale_expr
    case ("renormalization_scale")
       pn => rt_pn%ren_scale_expr
    case ("weight")
       pn => rt_pn%weight_expr
    case ("selection")
       pn => rt_pn%selection_lexpr
    case ("reweight")
       pn => rt_pn%reweight_expr
    case ("analysis")
       pn => rt_pn%analysis_lexpr
    end select
    if (associated (pn)) then
       write (u, "(A,1x,A,1x,A)")  "Expression:", char (name), "(parse tree):"
       call parse_node_write_rec (pn, u)
    else
       write (u, "(A,1x,A,A)")  "Expression:", char (name), ": [undefined]"
    end if
  end subroutine rt_parse_nodes_show

  module subroutine rt_data_write (object, unit, vars, pacify)
    class(rt_data_t), intent(in) :: object
    integer, intent(in), optional :: unit
    type(string_t), dimension(:), intent(in), optional :: vars
    logical, intent(in), optional :: pacify
    integer :: u, i
    u = given_output_unit (unit)
    call write_separator (u, 2)
    write (u, "(1x,A)")  "Runtime data:"
    if (object%get_n_export () > 0) then
       call write_separator (u, 2)
       write (u, "(1x,A)")  "Exported objects and variables:"
       call write_separator (u)
       call object%write_exports (u)
    end if
    if (present (vars)) then
       if (size (vars) /= 0) then
          call write_separator (u, 2)
          write (u, "(1x,A)")  "Selected variables:"
          call write_separator (u)
          call object%write_vars (u, vars)
       end if
    else
       call write_separator (u, 2)
       if (associated (object%model)) then
          call object%model%write_var_list (u, follow_link=.true.)
       else
          call object%var_list%write (u, follow_link=.true.)
       end if
    end if
    if (object%it_list%get_n_pass () > 0) then
       call write_separator (u, 2)
       write (u, "(1x)", advance="no")
       call object%it_list%write (u)
    end if
    if (associated (object%model)) then
       call write_separator (u, 2)
       call object%model%write (u)
    end if
    call object%prclib_stack%write (u)
    call object%beam_structure%write (u)
    call write_separator (u, 2)
    call object%pn%write (u)
    if (allocated (object%sample_fmt)) then
       call write_separator (u)
       write (u, "(1x,A)", advance="no")  "Event sample formats = "
       do i = 1, size (object%sample_fmt)
          if (i > 1)  write (u, "(A,1x)", advance="no")  ","
          write (u, "(A)", advance="no")  char (object%sample_fmt(i))
       end do
       write (u, "(A)")
    end if
    call write_separator (u)
    write (u, "(1x,A)", advance="no")  "Event callback:"
    if (allocated (object%event_callback)) then
       call object%event_callback%write (u)
    else
       write (u, "(1x,A)")  "[undefined]"
    end if
    call object%process_stack%write (u, pacify)
    write (u, "(1x,A,1x,L1)")  "quit     :", object%quit
    write (u, "(1x,A,1x,I0)")  "quit_code:", object%quit_code
    call write_separator (u, 2)
    write (u, "(1x,A,1x,A)")   "Logfile  :", "'" // trim (char (object%logfile)) // "'"
    call write_separator (u, 2)
  end subroutine rt_data_write

  module subroutine rt_data_write_vars (object, unit, vars)
    class(rt_data_t), intent(in), target :: object
    integer, intent(in), optional :: unit
    type(string_t), dimension(:), intent(in) :: vars
    type(var_list_t), pointer :: var_list
    integer :: u, i
    u = given_output_unit (unit)
    var_list => object%get_var_list_ptr ()
    do i = 1, size (vars)
       associate (var => vars(i))
         if (var_list%contains (var, follow_link=.true.)) then
            call var_list%write_var (var, unit = u, &
                 follow_link = .true., defined=.true.)
         end if
       end associate
    end do
  end subroutine rt_data_write_vars

  module subroutine rt_data_write_model_list (object, unit)
    class(rt_data_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    call object%model_list%write (u)
  end subroutine rt_data_write_model_list

  module subroutine rt_data_write_libraries (object, unit, libpath)
    class(rt_data_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: libpath
    integer :: u
    u = given_output_unit (unit)
    call object%prclib_stack%write (u, libpath)
  end subroutine rt_data_write_libraries

  module subroutine rt_data_write_beams (object, unit)
    class(rt_data_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    call write_separator (u, 2)
    call object%beam_structure%write (u)
    call write_separator (u, 2)
  end subroutine rt_data_write_beams

  module subroutine rt_data_write_expr (object, unit)
    class(rt_data_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    call write_separator (u, 2)
    call object%pn%write (u)
    call write_separator (u, 2)
  end subroutine rt_data_write_expr

  module subroutine rt_data_write_process_stack (object, unit)
    class(rt_data_t), intent(in) :: object
    integer, intent(in), optional :: unit
    call object%process_stack%write (unit)
  end subroutine rt_data_write_process_stack

  module subroutine rt_data_write_var_descriptions (rt_data, unit, ascii_output)
    class(rt_data_t), intent(in) :: rt_data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: ascii_output
    integer :: u
    logical :: ao
    u = given_output_unit (unit)
    ao = .false.;  if (present (ascii_output))  ao = ascii_output
    call rt_data%var_list%write (u, follow_link=.true., &
         descriptions=.true., ascii_output=ao)
  end subroutine rt_data_write_var_descriptions

  module subroutine rt_data_show_description_of_string (rt_data, string, &
       unit, ascii_output)
    class(rt_data_t), intent(in) :: rt_data
    type(string_t), intent(in) :: string
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: ascii_output
    integer :: u
    logical :: ao
    u = given_output_unit (unit)
    ao = .false.;  if (present (ascii_output))  ao = ascii_output
    call rt_data%var_list%write_var (string, unit=u, follow_link=.true., &
         defined=.false., descriptions=.true., ascii_output=ao)
  end subroutine rt_data_show_description_of_string

  module subroutine rt_data_clear_beams (global)
    class(rt_data_t), intent(inout) :: global
    call global%beam_structure%final_sf ()
    call global%beam_structure%final_pol ()
    call global%beam_structure%final_mom ()
  end subroutine rt_data_clear_beams

  module subroutine rt_data_global_init (global, paths, logfile)
    class(rt_data_t), intent(out), target :: global
    type(paths_t), intent(in), optional :: paths
    type(string_t), intent(in), optional :: logfile
    integer :: seed
    call global%os_data%init (paths)
    if (present (logfile)) then
       global%logfile = logfile
    else
       global%logfile = ""
    end if
    allocate (global%out_files)
    call system_clock (seed)
    call global%var_list%init_defaults (seed, paths)
    call global%init_pointer_variables ()
    call global%process_stack%init_var_list (global%var_list)
  end subroutine rt_data_global_init

  module subroutine rt_data_local_init (local, global, env)
    class(rt_data_t), intent(inout), target :: local
    type(rt_data_t), intent(in), target :: global
    integer, intent(in), optional :: env
    local%context => global
    call local%process_stack%link (global%process_stack)
    call local%process_stack%init_var_list (local%var_list)
    call local%process_stack%link_var_list (global%var_list)
    call local%var_list%append_string (var_str ("$model_name"), &
         var_str (""), intrinsic=.true.)
    call local%init_pointer_variables ()
    local%fallback_model => global%fallback_model
    local%os_data = global%os_data
    local%logfile = global%logfile
    call local%model_list%link (global%model_list)
    local%model => global%model
    if (associated (local%model)) then
       call local%model%link_var_list (local%var_list)
    end if
    if (allocated (global%event_callback)) then
       allocate (local%event_callback, source = global%event_callback)
    end if
  end subroutine rt_data_local_init

  module subroutine rt_data_init_pointer_variables (local)
    class(rt_data_t), intent(inout), target :: local
    logical, target, save :: known = .true.
    call local%var_list%append_string_ptr (var_str ("$fc"), &
         local%os_data%fc, known, intrinsic=.true., &
         description=var_str('This string variable gives the ' // &
         '\ttt{Fortran} compiler used within \whizard. It can ' // &
         'only be accessed, not set by the user. (cf. also ' // &
         '\ttt{\$fcflags}, \ttt{\$fclibs})'))
    call local%var_list%append_string_ptr (var_str ("$fcflags"), &
         local%os_data%fcflags, known, intrinsic=.true., &
         description=var_str('This string variable gives the ' // &
         'compiler flags for the \ttt{Fortran} compiler used ' // &
         'within \whizard. It can only be accessed, not set by ' // &
         'the user. (cf. also \ttt{\$fc}, \ttt{\$fclibs})'))
    call local%var_list%append_string_ptr (var_str ("$fclibs"), &
         local%os_data%fclibs, known, intrinsic=.true., &
         description=var_str('This string variable gives the ' // &
         'linked libraries for the \ttt{Fortran} compiler used ' // &
         'within \whizard. It can only be accessed, not set by ' // &
         'the user. (cf. also \ttt{\$fc}, \ttt{\$fcflags})'))
  end subroutine rt_data_init_pointer_variables

  module subroutine rt_data_deactivate (local, global, keep_local)
    class(rt_data_t), intent(inout), target :: local
    class(rt_data_t), intent(inout), optional, target :: global
    logical, intent(in), optional :: keep_local
    type(string_t) :: local_model, local_scheme
    logical :: same_model, delete
    delete = .true.;  if (present (keep_local))  delete = .not. keep_local
    if (present (global)) then
       if (associated (global%model) .and. associated (local%model)) then
          local_model = local%model%get_name ()
          if (global%model%has_schemes ()) then
             local_scheme = local%model%get_scheme ()
             same_model = &
                  global%model%matches (local_model, local_scheme)
          else
             same_model = global%model%matches (local_model)
          end if
       else
          same_model = .false.
       end if
       if (delete) then
          call local%process_stack%clear ()
          call local%unselect_model ()
          call local%unset_values ()
       else if (associated (local%model)) then
          call local%ensure_model_copy ()
       end if
       if (.not. same_model .and. associated (global%model)) then
          if (global%model%has_schemes ()) then
             call msg_message ("Restoring model '" // &
                  char (global%model%get_name ()) // "', scheme '" // &
                  char (global%model%get_scheme ()) // "'")
          else
             call msg_message ("Restoring model '" // &
                  char (global%model%get_name ()) // "'")
          end if
       end if
       if (associated (global%model)) then
          call global%model%link_var_list (global%var_list)
       end if
       call global%restore_globals (local)
    else
       call local%unselect_model ()
    end if
  end subroutine rt_data_deactivate

  module subroutine rt_data_copy_globals (global, local)
    class(rt_data_t), intent(in) :: global
    class(rt_data_t), intent(inout) :: local
    local%prclib_stack = global%prclib_stack
  end subroutine rt_data_copy_globals

  module subroutine rt_data_restore_globals (global, local)
    class(rt_data_t), intent(inout) :: global
    class(rt_data_t), intent(inout) :: local
    global%prclib_stack = local%prclib_stack
    call local%handle_exports (global)
  end subroutine rt_data_restore_globals

  module subroutine rt_data_write_exports (rt_data, unit)
    class(rt_data_t), intent(in) :: rt_data
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    do i = 1, rt_data%get_n_export ()
       write (u, "(A)")  char (rt_data%export(i))
    end do
  end subroutine rt_data_write_exports

  module function rt_data_get_n_export (rt_data) result (n)
    class(rt_data_t), intent(in) :: rt_data
    integer :: n
    if (allocated (rt_data%export)) then
       n = size (rt_data%export)
    else
       n = 0
    end if
  end function rt_data_get_n_export

  module subroutine rt_data_append_exports (rt_data, export)
    class(rt_data_t), intent(inout) :: rt_data
    type(string_t), dimension(:), intent(in) :: export
    logical, dimension(:), allocatable :: mask
    type(string_t), dimension(:), allocatable :: tmp
    integer :: i, j, n
    if (.not. allocated (rt_data%export))  allocate (rt_data%export (0))
    n = size (rt_data%export)
    allocate (mask (size (export)), source=.false.)
    do i = 1, size (export)
       mask(i) = all (export(i) /= rt_data%export) &
            .and. all (export(i) /= export(:i-1))
    end do
    if (count (mask) > 0) then
       allocate (tmp (n + count (mask)))
       tmp(1:n) = rt_data%export(:)
       j = n
       do i = 1, size (export)
          if (mask(i)) then
             j = j + 1
             tmp(j) = export(i)
          end if
       end do
       call move_alloc (from=tmp, to=rt_data%export)
    end if
  end subroutine rt_data_append_exports

  module subroutine rt_data_handle_exports (local, global)
    class(rt_data_t), intent(inout), target :: local
    class(rt_data_t), intent(inout), target :: global
    type(string_t) :: export
    integer :: i
    if (local%get_n_export () > 0) then
       do i = 1, local%get_n_export ()
          export = local%export(i)
          select case (char (export))
          case ("results")
             call msg_message ("Exporting integration results &
                  &to outer environment")
             call local%transfer_process_stack (global)
          case default
             call msg_bug ("handle exports: '" &
                  // char (export) // "' unsupported")
          end select
       end do
    end if
  end subroutine rt_data_handle_exports

  module subroutine rt_data_transfer_process_stack (local, global)
    class(rt_data_t), intent(inout), target :: local
    class(rt_data_t), intent(inout), target :: global
    type(process_entry_t), pointer :: process
    type(string_t) :: process_id
    do
       call local%process_stack%pop_last (process)
       if (.not. associated (process))  exit
       process_id = process%get_id ()
       call global%process_stack%push (process)
       call global%process_stack%fill_result_vars (process_id)
       call global%process_stack%update_result_vars &
            (process_id, global%var_list)
    end do
  end subroutine rt_data_transfer_process_stack

  module subroutine rt_data_global_final (global)
    class(rt_data_t), intent(inout) :: global
    call global%process_stack%final ()
    call global%prclib_stack%final ()
    call global%model_list%final ()
    call global%var_list%final (follow_link=.false.)
    if (associated (global%out_files)) then
       call file_list_final (global%out_files)
       deallocate (global%out_files)
    end if
  end subroutine rt_data_global_final

  module subroutine rt_data_local_final (local)
    class(rt_data_t), intent(inout) :: local
    call local%process_stack%clear ()
    call local%model_list%final ()
    call local%var_list%final (follow_link=.false.)
  end subroutine rt_data_local_final

  module subroutine rt_data_read_model (global, name, model, scheme)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    type(string_t), intent(in), optional :: scheme
    type(model_t), pointer, intent(out) :: model
    type(string_t) :: filename
    filename = name // ".mdl"
    call global%model_list%read_model &
         (name, filename, global%os_data, model, scheme)
  end subroutine rt_data_read_model

  module subroutine rt_data_read_ufo_model (global, name, model, ufo_path)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    type(model_t), pointer, intent(out) :: model
    type(string_t), intent(in), optional :: ufo_path
    type(string_t) :: filename
    filename = name // ".ufo.mdl"
    call global%model_list%read_model &
         (name, filename, global%os_data, model, ufo=.true., ufo_path=ufo_path)
  end subroutine rt_data_read_ufo_model

  module subroutine rt_data_init_fallback_model (global, name, filename)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name, filename
    call global%model_list%read_model &
         (name, filename, global%os_data, global%fallback_model)
  end subroutine rt_data_init_fallback_model

  module subroutine rt_data_select_model (global, name, scheme, ufo, ufo_path)
    class(rt_data_t), intent(inout), target :: global
    type(string_t), intent(in) :: name
    type(string_t), intent(in), optional :: scheme
    logical, intent(in), optional :: ufo
    type(string_t), intent(in), optional :: ufo_path
    logical :: same_model, ufo_model
    ufo_model = .false.;  if (present (ufo))  ufo_model = ufo
    if (associated (global%model)) then
       same_model = global%model%matches (name, scheme, ufo)
    else
       same_model = .false.
    end if
    if (.not. same_model) then
       global%model => global%model_list%get_model_ptr (name, scheme, ufo)
       if (.not. associated (global%model)) then
          if (ufo_model) then
             call global%read_ufo_model (name, global%model, ufo_path)
          else
             call global%read_model (name, global%model)
          end if
          global%model_is_copy = .false.
       else if (associated (global%context)) then
          global%model_is_copy = &
               global%model_list%model_exists (name, scheme, ufo, &
               follow_link=.false.)
       else
          global%model_is_copy = .false.
       end if
    end if
    if (associated (global%model)) then
       call global%model%link_var_list (global%var_list)
       call global%var_list%set_string (var_str ("$model_name"), &
            name, is_known = .true.)
       if (global%model%is_ufo_model ()) then
          call msg_message ("Switching to model '" // char (name) // "' " &
               // "(generated from UFO source)")
       else if (global%model%has_schemes ()) then
          call msg_message ("Switching to model '" // char (name) // "', " &
               // "scheme '" // char (global%model%get_scheme ()) // "'")
       else
          call msg_message ("Switching to model '" // char (name) // "'")
       end if
    else
       call global%var_list%set_string (var_str ("$model_name"), &
            var_str (""), is_known = .false.)
    end if
  end subroutine rt_data_select_model

  module subroutine rt_data_unselect_model (global)
    class(rt_data_t), intent(inout), target :: global
    if (associated (global%model)) then
       global%model => null ()
       global%model_is_copy = .false.
    end if
  end subroutine rt_data_unselect_model

  module subroutine rt_data_ensure_model_copy (global)
    class(rt_data_t), intent(inout), target :: global
    if (associated (global%context)) then
       if (.not. global%model_is_copy) then
          call global%model_list%append_copy (global%model, global%model)
          global%model_is_copy = .true.
          call global%model%link_var_list (global%var_list)
       end if
    end if
  end subroutine rt_data_ensure_model_copy

  module subroutine rt_data_model_set_real &
       (global, name, rval, verbose, pacified)
    class(rt_data_t), intent(inout), target :: global
    type(string_t), intent(in) :: name
    real(default), intent(in) :: rval
    logical, intent(in), optional :: verbose, pacified
    call global%ensure_model_copy ()
    call global%model%set_real (name, rval, verbose, pacified)
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
    call global%ensure_model_copy ()
    if (present (polarized)) then
       if (polarized) then
          call global%model%set_polarized (pdg)
       else
          call global%model%set_unpolarized (pdg)
       end if
    end if
    if (present (stable)) then
       if (stable) then
          call global%model%set_stable (pdg)
       else if (present (decay)) then
          call global%model%set_unstable &
               (pdg, decay, isotropic_decay, diagonal_decay, decay_helicity)
       else
          call msg_bug ("Setting particle unstable: missing decay processes")
       end if
    end if
  end subroutine rt_data_modify_particle

  module function rt_data_get_var_list_ptr (global) result (var_list)
    class(rt_data_t), intent(in), target :: global
    type(var_list_t), pointer :: var_list
    if (associated (global%model)) then
       var_list => global%model%get_var_list_ptr ()
    else
       var_list => global%var_list
    end if
  end function rt_data_get_var_list_ptr

  module subroutine rt_data_append_log (local, name, lval, intrinsic, user)
    class(rt_data_t), intent(inout) :: local
    type(string_t), intent(in) :: name
    logical, intent(in), optional :: lval
    logical, intent(in), optional :: intrinsic, user
    call local%var_list%append_log (name, lval, &
         intrinsic = intrinsic, user = user)
  end subroutine rt_data_append_log

  module subroutine rt_data_append_int (local, name, ival, intrinsic, user)
    class(rt_data_t), intent(inout) :: local
    type(string_t), intent(in) :: name
    integer, intent(in), optional :: ival
    logical, intent(in), optional :: intrinsic, user
    call local%var_list%append_int (name, ival, &
         intrinsic = intrinsic, user = user)
  end subroutine rt_data_append_int

  module subroutine rt_data_append_real (local, name, rval, intrinsic, user)
    class(rt_data_t), intent(inout) :: local
    type(string_t), intent(in) :: name
    real(default), intent(in), optional :: rval
    logical, intent(in), optional :: intrinsic, user
    call local%var_list%append_real (name, rval, &
         intrinsic = intrinsic, user = user)
  end subroutine rt_data_append_real

  module subroutine rt_data_append_cmplx (local, name, cval, intrinsic, user)
    class(rt_data_t), intent(inout) :: local
    type(string_t), intent(in) :: name
    complex(default), intent(in), optional :: cval
    logical, intent(in), optional :: intrinsic, user
    call local%var_list%append_cmplx (name, cval, &
         intrinsic = intrinsic, user = user)
  end subroutine rt_data_append_cmplx

  module subroutine rt_data_append_subevt (local, name, pval, intrinsic, user)
    class(rt_data_t), intent(inout) :: local
    type(string_t), intent(in) :: name
    type(subevt_t), intent(in), optional :: pval
    logical, intent(in) :: intrinsic, user
    call local%var_list%append_subevt (name, &
         intrinsic = intrinsic, user = user)
  end subroutine rt_data_append_subevt

  module subroutine rt_data_append_pdg_array &
       (local, name, aval, intrinsic, user)
    class(rt_data_t), intent(inout) :: local
    type(string_t), intent(in) :: name
    type(pdg_array_t), intent(in), optional :: aval
    logical, intent(in), optional :: intrinsic, user
    call local%var_list%append_pdg_array (name, aval, &
         intrinsic = intrinsic, user = user)
  end subroutine rt_data_append_pdg_array

  module subroutine rt_data_append_string (local, name, sval, intrinsic, user)
    class(rt_data_t), intent(inout) :: local
    type(string_t), intent(in) :: name
    type(string_t), intent(in), optional :: sval
    logical, intent(in), optional :: intrinsic, user
    call local%var_list%append_string (name, sval, &
         intrinsic = intrinsic, user = user)
  end subroutine rt_data_append_string

  module subroutine rt_data_import_values (local)
    class(rt_data_t), intent(inout) :: local
    type(rt_data_t), pointer :: global
    global => local%context
    if (associated (global)) then
       call local%var_list%import (global%var_list)
    end if
  end subroutine rt_data_import_values

  module subroutine rt_data_unset_values (global)
    class(rt_data_t), intent(inout) :: global
    call global%var_list%undefine (follow_link=.false.)
  end subroutine rt_data_unset_values

  module subroutine rt_data_set_log &
       (global, name, lval, is_known, force, verbose)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    logical, intent(in) :: lval
    logical, intent(in) :: is_known
    logical, intent(in), optional :: force, verbose
    call global%var_list%set_log (name, lval, is_known, &
         force=force, verbose=verbose)
  end subroutine rt_data_set_log

  module subroutine rt_data_set_int &
       (global, name, ival, is_known, force, verbose)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    integer, intent(in) :: ival
    logical, intent(in) :: is_known
    logical, intent(in), optional :: force, verbose
    call global%var_list%set_int (name, ival, is_known, &
         force=force, verbose=verbose)
  end subroutine rt_data_set_int

  module subroutine rt_data_set_real &
       (global, name, rval, is_known, force, verbose, pacified)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    real(default), intent(in) :: rval
    logical, intent(in) :: is_known
    logical, intent(in), optional :: force, verbose, pacified
    call global%var_list%set_real (name, rval, is_known, &
         force=force, verbose=verbose, pacified=pacified)
  end subroutine rt_data_set_real

  module subroutine rt_data_set_cmplx &
       (global, name, cval, is_known, force, verbose, pacified)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    complex(default), intent(in) :: cval
    logical, intent(in) :: is_known
    logical, intent(in), optional :: force, verbose, pacified
    call global%var_list%set_cmplx (name, cval, is_known, &
         force=force, verbose=verbose, pacified=pacified)
  end subroutine rt_data_set_cmplx

  module subroutine rt_data_set_subevt &
       (global, name, pval, is_known, force, verbose)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    type(subevt_t), intent(in) :: pval
    logical, intent(in) :: is_known
    logical, intent(in), optional :: force, verbose
    call global%var_list%set_subevt (name, pval, is_known, &
         force=force, verbose=verbose)
  end subroutine rt_data_set_subevt

  module subroutine rt_data_set_pdg_array &
       (global, name, aval, is_known, force, verbose)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    type(pdg_array_t), intent(in) :: aval
    logical, intent(in) :: is_known
    logical, intent(in), optional :: force, verbose
    call global%var_list%set_pdg_array (name, aval, is_known, &
         force=force, verbose=verbose)
  end subroutine rt_data_set_pdg_array

  module subroutine rt_data_set_string &
       (global, name, sval, is_known, force, verbose)
    class(rt_data_t), intent(inout) :: global
    type(string_t), intent(in) :: name
    type(string_t), intent(in) :: sval
    logical, intent(in) :: is_known
    logical, intent(in), optional :: force, verbose
    call global%var_list%set_string (name, sval, is_known, &
         force=force, verbose=verbose)
  end subroutine rt_data_set_string

  module function rt_data_get_lval (global, name) result (lval)
    logical :: lval
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    lval = var_list%get_lval (name)
  end function rt_data_get_lval

  module function rt_data_get_ival (global, name) result (ival)
    integer :: ival
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    ival = var_list%get_ival (name)
  end function rt_data_get_ival

  module function rt_data_get_rval (global, name) result (rval)
    real(default) :: rval
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    rval = var_list%get_rval (name)
  end function rt_data_get_rval

  module function rt_data_get_cval (global, name) result (cval)
    complex(default) :: cval
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    cval = var_list%get_cval (name)
  end function rt_data_get_cval

  module function rt_data_get_aval (global, name) result (aval)
    type(pdg_array_t) :: aval
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    aval = var_list%get_aval (name)
  end function rt_data_get_aval

  module function rt_data_get_pval (global, name) result (pval)
    type(subevt_t) :: pval
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    pval = var_list%get_pval (name)
  end function rt_data_get_pval

  module function rt_data_get_sval (global, name) result (sval)
    type(string_t) :: sval
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    sval = var_list%get_sval (name)
  end function rt_data_get_sval

  module function rt_data_contains (global, name) result (lval)
    logical :: lval
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    lval = var_list%contains (name)
  end function rt_data_contains

  module  function rt_data_is_known (global, name) result (lval)
    logical :: lval
    class(rt_data_t), intent(in), target :: global
    type(string_t), intent(in) :: name
    type(var_list_t), pointer :: var_list
    var_list => global%get_var_list_ptr ()
    lval = var_list%is_known (name)
  end function rt_data_is_known

  module subroutine rt_data_add_prclib (global, prclib_entry)
    class(rt_data_t), intent(inout) :: global
    type(prclib_entry_t), intent(inout), pointer :: prclib_entry
    call global%prclib_stack%push (prclib_entry)
    call global%update_prclib (global%prclib_stack%get_first_ptr ())
  end subroutine rt_data_add_prclib

  module subroutine rt_data_update_prclib (global, lib)
    class(rt_data_t), intent(inout) :: global
    type(process_library_t), intent(in), target :: lib
    global%prclib => lib
    if (global%var_list%contains (&
         var_str ("$library_name"), follow_link = .false.)) then
       call global%var_list%set_string (var_str ("$library_name"), &
            global%prclib%get_name (), is_known=.true.)
    else
       call global%var_list%append_string ( &
            var_str ("$library_name"), global%prclib%get_name (), &
            intrinsic = .true.)
    end if
  end subroutine rt_data_update_prclib

  module function rt_data_get_helicity_selection &
       (rt_data) result (helicity_selection)
    class(rt_data_t), intent(in) :: rt_data
    type(helicity_selection_t) :: helicity_selection
    associate (var_list => rt_data%var_list)
      helicity_selection%active = var_list%get_lval (&
           var_str ("?helicity_selection_active"))
      if (helicity_selection%active) then
         helicity_selection%threshold = var_list%get_rval (&
              var_str ("helicity_selection_threshold"))
         helicity_selection%cutoff = var_list%get_ival (&
              var_str ("helicity_selection_cutoff"))
      end if
    end associate
  end function rt_data_get_helicity_selection

  module subroutine rt_data_show_beams (rt_data, unit)
    class(rt_data_t), intent(in) :: rt_data
    integer, intent(in), optional :: unit
    type(string_t) :: s
    integer :: u
    u = given_output_unit (unit)
    associate (beams => rt_data%beam_structure, var_list => rt_data%var_list)
      call beams%write (u)
      if (.not. beams%asymmetric () .and. beams%get_n_beam () == 2) then
         write (u, "(2x,A," // FMT_19 // ",1x,'GeV')") "sqrts =", &
              var_list%get_rval (var_str ("sqrts"))
      end if
      if (beams%contains ("pdf_builtin")) then
         s = var_list%get_sval (var_str ("$pdf_builtin_set"))
         if (s /= "") then
            write (u, "(2x,A,1x,3A)")  "PDF set =", '"', char (s), '"'
         else
            write (u, "(2x,A,1x,A)")  "PDF set =", "[undefined]"
         end if
      end if
      if (beams%contains ("lhapdf")) then
         s = var_list%get_sval (var_str ("$lhapdf_dir"))
         if (s /= "") then
            write (u, "(2x,A,1x,3A)")  "LHAPDF dir    =", '"', char (s), '"'
         end if
         s = var_list%get_sval (var_str ("$lhapdf_file"))
         if (s /= "") then
            write (u, "(2x,A,1x,3A)")  "LHAPDF file   =", '"', char (s), '"'
            write (u, "(2x,A,1x,I0)") "LHAPDF member =", &
                 var_list%get_ival (var_str ("lhapdf_member"))
         else
            write (u, "(2x,A,1x,A)")  "LHAPDF file   =", "[undefined]"
         end if
      end if
      if (beams%contains ("lhapdf_photon")) then
         s = var_list%get_sval (var_str ("$lhapdf_dir"))
         if (s /= "") then
            write (u, "(2x,A,1x,3A)")  "LHAPDF dir    =", '"', char (s), '"'
         end if
         s = var_list%get_sval (var_str ("$lhapdf_photon_file"))
         if (s /= "") then
            write (u, "(2x,A,1x,3A)")  "LHAPDF file   =", '"', char (s), '"'
            write (u, "(2x,A,1x,I0)") "LHAPDF member =", &
                 var_list%get_ival (var_str ("lhapdf_member"))
            write (u, "(2x,A,1x,I0)") "LHAPDF scheme =", &
                 var_list%get_ival (&
                 var_str ("lhapdf_photon_scheme"))
         else
            write (u, "(2x,A,1x,A)")  "LHAPDF file   =", "[undefined]"
         end if
      end if
      if (beams%contains ("isr")) then
         write (u, "(2x,A," // FMT_19 // ")") "ISR alpha        =", &
              var_list%get_rval (var_str ("isr_alpha"))
         write (u, "(2x,A," // FMT_19 // ")") "ISR Q max        =", &
              var_list%get_rval (var_str ("isr_q_max"))
         write (u, "(2x,A," // FMT_19 // ")") "ISR mass         =", &
              var_list%get_rval (var_str ("isr_mass"))
         write (u, "(2x,A,1x,I0)") "ISR order        =", &
              var_list%get_ival (var_str ("isr_order"))
         write (u, "(2x,A,1x,L1)") "ISR recoil       =", &
              var_list%get_lval (var_str ("?isr_recoil"))
         write (u, "(2x,A,1x,L1)") "ISR energy cons. =", &
              var_list%get_lval (var_str ("?isr_keep_energy"))
      end if
      if (beams%contains ("epa")) then
         write (u, "(2x,A," // FMT_19 // ")") "EPA alpha         =", &
              var_list%get_rval (var_str ("epa_alpha"))
         write (u, "(2x,A," // FMT_19 // ")") "EPA x min         =", &
              var_list%get_rval (var_str ("epa_x_min"))
         write (u, "(2x,A," // FMT_19 // ")") "EPA Q min         =", &
              var_list%get_rval (var_str ("epa_q_min"))
         write (u, "(2x,A," // FMT_19 // ")") "EPA Q max         =", &
              var_list%get_rval (var_str ("epa_q_max"))
         write (u, "(2x,A," // FMT_19 // ")") "EPA mass          =", &
              var_list%get_rval (var_str ("epa_mass"))
         write (u, "(2x,A,1x,L1)") "EPA recoil        =", &
              var_list%get_lval (var_str ("?epa_recoil"))
         write (u, "(2x,A,1x,L1)") "EPA  energy cons. =", &
              var_list%get_lval (var_str ("?epa_keep_energy"))
      end if
      if (beams%contains ("ewa")) then
         write (u, "(2x,A," // FMT_19 // ")") "EWA x min       =", &
              var_list%get_rval (var_str ("ewa_x_min"))
         write (u, "(2x,A," // FMT_19 // ")") "EWA Pt max      =", &
              var_list%get_rval (var_str ("ewa_pt_max"))
         write (u, "(2x,A," // FMT_19 // ")") "EWA mass        =", &
              var_list%get_rval (var_str ("ewa_mass"))
         write (u, "(2x,A,1x,L1)") "EWA recoil       =", &
              var_list%get_lval (var_str ("?ewa_recoil"))
         write (u, "(2x,A,1x,L1)") "EWA energy cons. =", &
              var_list%get_lval (var_str ("ewa_keep_energy"))
      end if
      if (beams%contains ("circe1")) then
         write (u, "(2x,A,1x,I0)") "CIRCE1 version    =", &
              var_list%get_ival (var_str ("circe1_ver"))
         write (u, "(2x,A,1x,I0)") "CIRCE1 revision   =", &
              var_list%get_ival (var_str ("circe1_rev"))
         s = var_list%get_sval (var_str ("$circe1_acc"))
         write (u, "(2x,A,1x,A)") "CIRCE1 acceler.   =", char (s)
         write (u, "(2x,A,1x,I0)") "CIRCE1 chattin.   =", &
              var_list%get_ival (var_str ("circe1_chat"))
         write (u, "(2x,A," // FMT_19 // ")") "CIRCE1 sqrts      =", &
              var_list%get_rval (var_str ("circe1_sqrts"))
         write (u, "(2x,A," // FMT_19 // ")") "CIRCE1 epsil.     =", &
              var_list%get_rval (var_str ("circe1_eps"))
         write (u, "(2x,A,1x,L1)") "CIRCE1 phot. 1  =", &
              var_list%get_lval (var_str ("?circe1_photon1"))
         write (u, "(2x,A,1x,L1)") "CIRCE1 phot. 2  =", &
              var_list%get_lval (var_str ("?circe1_photon2"))
         write (u, "(2x,A,1x,L1)") "CIRCE1 generat. =", &
              var_list%get_lval (var_str ("?circe1_generate"))
         write (u, "(2x,A,1x,L1)") "CIRCE1 mapping  =", &
              var_list%get_lval (var_str ("?circe1_map"))
         write (u, "(2x,A," // FMT_19 // ")") "CIRCE1 map. slope =", &
              var_list%get_rval (var_str ("circe1_mapping_slope"))
         write (u, "(2x,A,1x,L1)") "CIRCE recoil photon =", &
              var_list%get_lval (var_str ("?circe1_with_radiation"))
      end if
      if (beams%contains ("circe2")) then
         s = var_list%get_sval (var_str ("$circe2_design"))
         write (u, "(2x,A,1x,A)") "CIRCE2 design   =", char (s)
         s = var_list%get_sval (var_str ("$circe2_file"))
         write (u, "(2x,A,1x,A)") "CIRCE2 file     =", char (s)
         write (u, "(2x,A,1x,L1)") "CIRCE2 polarized =", &
              var_list%get_lval (var_str ("?circe2_polarized"))
      end if
      if (beams%contains ("gaussian")) then
         write (u, "(2x,A,1x," // FMT_12 // ")") "Gaussian spread 1    =", &
              var_list%get_rval (var_str ("gaussian_spread1"))
         write (u, "(2x,A,1x," // FMT_12 // ")") "Gaussian spread 2    =", &
              var_list%get_rval (var_str ("gaussian_spread2"))
      end if
      if (beams%contains ("beam_events")) then
         s = var_list%get_sval (var_str ("$beam_events_file"))
         write (u, "(2x,A,1x,A)") "Beam events file     =", char (s)
         write (u, "(2x,A,1x,L1)") "Beam events EOF warn =", &
              var_list%get_lval (var_str ("?beam_events_warn_eof"))
      end if
    end associate
  end subroutine rt_data_show_beams

  module function rt_data_get_sqrts (rt_data) result (sqrts)
    class(rt_data_t), intent(in) :: rt_data
    real(default) :: sqrts
    sqrts = rt_data%var_list%get_rval (var_str ("sqrts"))
  end function rt_data_get_sqrts

  module subroutine rt_data_pacify (rt_data, efficiency_reset, error_reset)
    class(rt_data_t), intent(inout) :: rt_data
    logical, intent(in), optional :: efficiency_reset, error_reset
    type(process_entry_t), pointer :: process
    process => rt_data%process_stack%first
    do while (associated (process))
       call process%pacify (efficiency_reset, error_reset)
       process => process%next
    end do
  end subroutine rt_data_pacify

  module subroutine rt_data_set_event_callback (global, callback)
    class(rt_data_t), intent(inout) :: global
    class(event_callback_t), intent(in) :: callback
    if (allocated (global%event_callback))  deallocate (global%event_callback)
    allocate (global%event_callback, source = callback)
  end subroutine rt_data_set_event_callback

  module function rt_data_has_event_callback (global) result (flag)
    class(rt_data_t), intent(in) :: global
    logical :: flag
    flag = allocated (global%event_callback)
  end function rt_data_has_event_callback

  module function rt_data_get_event_callback (global) result (callback)
    class(rt_data_t), intent(in) :: global
    class(event_callback_t), allocatable :: callback
    if (allocated (global%event_callback)) then
       allocate (callback, source = global%event_callback)
    end if
  end function rt_data_get_event_callback

  module subroutine fix_system_dependencies (global)
    class(rt_data_t), intent(inout), target :: global
    type(var_list_t), pointer :: var_list

    var_list => global%get_var_list_ptr ()
    call var_list%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true., force=.true.)
    call var_list%set_log (var_str ("?openmp_is_active"), &
         .false., is_known = .true., force=.true.)
    call var_list%set_int (var_str ("openmp_num_threads_default"), &
         1, is_known = .true., force=.true.)
    call var_list%set_int (var_str ("openmp_num_threads"), &
         1, is_known = .true., force=.true.)
    call var_list%set_int (var_str ("real_range"), &
         307, is_known = .true., force=.true.)
    call var_list%set_int (var_str ("real_precision"), &
         15, is_known = .true., force=.true.)
    call var_list%set_real (var_str ("real_epsilon"), &
         1.e-16_default, is_known = .true., force=.true.)
    call var_list%set_real (var_str ("real_tiny"), &
         1.e-300_default, is_known = .true., force=.true.)

    global%os_data%fc = "Fortran-compiler"
    global%os_data%fcflags = "Fortran-flags"
    global%os_data%fclibs = "Fortran-libs"

  end subroutine fix_system_dependencies

  module subroutine show_description_of_string (string)
    type(string_t), intent(in) :: string
    type(rt_data_t), target :: global
    call global%global_init ()
    call global%show_description_of_string (string, ascii_output=.true.)
  end subroutine show_description_of_string

  module subroutine show_tex_descriptions ()
    type(rt_data_t), target :: global
    call global%global_init ()
    call fix_system_dependencies (global)
    call global%set_int (var_str ("seed"), 0, is_known=.true.)
    call global%var_list%sort ()
    call global%write_var_descriptions ()
  end subroutine show_tex_descriptions


end submodule rt_data_s

