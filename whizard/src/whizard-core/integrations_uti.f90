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

module integrations_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use ifiles
  use lexers
  use parser
  use flavors
  use interactions, only: reset_interaction_counter
  use phs_forests
  use eval_trees
  use models
  use rt_data
  use process_configurations_ut, only: prepare_test_library
  use compilations, only: compile_library

  use integrations

  use phs_wood_ut, only: write_test_phs_file

  implicit none
  private

  public :: integrations_1
  public :: integrations_2
  public :: integrations_3
  public :: integrations_4
  public :: integrations_5
  public :: integrations_6
  public :: integrations_7
  public :: integrations_8
  public :: integrations_9
  public :: integrations_history_1

contains

  subroutine integrations_1 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global

    write (u, "(A)")  "* Test output: integrations_1"
    write (u, "(A)")  "*   Purpose: integrate test process"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()

    libname = "integration_1"
    procname = "prc_config_a"

    call prepare_test_library (global, libname, 1)
    call compile_library (libname, global)

    call global%set_string (var_str ("$run_id"), &
         var_str ("integrations1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true.)

    call global%write (u, vars = [ &
         var_str ("$method"), &
         var_str ("sqrts"), &
         var_str ("$integration_method"), &
         var_str ("$phs_method"), &
         var_str ("$run_id")])

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_1"

  end subroutine integrations_1

  subroutine integrations_2 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global

    type(string_t) :: cut_expr_text
    type(ifile_t) :: ifile
    type(stream_t) :: stream
    type(parse_tree_t) :: parse_tree

    type(string_t), dimension(0) :: empty_string_array

    write (u, "(A)")  "* Test output: integrations_2"
    write (u, "(A)")  "*   Purpose: integrate test process with cut"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()

    write (u, "(A)")  "* Prepare a cut expression"
    write (u, "(A)")

    call syntax_pexpr_init ()
    cut_expr_text = "all Pt > 100 [s]"
    call ifile_append (ifile, cut_expr_text)
    call stream_init (stream, ifile)
    call parse_tree_init_lexpr (parse_tree, stream, .true.)
    global%pn%cuts_lexpr => parse_tree%get_root_ptr ()

    write (u, "(A)")  "* Build and initialize a test process"
    write (u, "(A)")

    libname = "integration_3"
    procname = "prc_config_a"

    call prepare_test_library (global, libname, 1)
    call compile_library (libname, global)

    call global%set_string (var_str ("$run_id"), &
         var_str ("integrations1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true.)

    call global%write (u, vars = empty_string_array)

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_2"

  end subroutine integrations_2

  subroutine integrations_3 (u)
    use kinds, only: default
    use iso_varying_string, string_t => varying_string
    use interactions, only: reset_interaction_counter
    use models
    use rt_data
    use process_configurations_ut, only: prepare_test_library
    use compilations, only: compile_library
    use integrations

    implicit none

    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global
    integer :: u_phs

    write (u, "(A)")  "* Test output: integrations_3"
    write (u, "(A)")  "*   Purpose: integrate test process"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and parameters"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()

    libname = "integration_3"
    procname = "prc_config_a"

    call prepare_test_library (global, libname, 1)
    call compile_library (libname, global)

    call global%set_string (var_str ("$run_id"), &
         var_str ("integrations1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("default"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?phs_s_mapping"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    write (u, "(A)")  "* Create a scratch phase-space file"
    write (u, "(A)")

    u_phs = free_unit ()
    open (u_phs, file = "integrations_3.phs", &
         status = "replace", action = "write")
    call write_test_phs_file (u_phs, var_str ("prc_config_a_i1"))
    close (u_phs)

    call global%set_string (var_str ("$phs_file"),&
         var_str ("integrations_3.phs"), is_known = .true.)

    call global%it_list%init ([1], [1000])

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true.)

    call global%write (u, vars = [ &
         var_str ("$phs_method"), &
         var_str ("$phs_file")])

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_phs_forest_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_3"

  end subroutine integrations_3

  subroutine integrations_4 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global

    write (u, "(A)")  "* Test output: integrations_4"
    write (u, "(A)")  "*   Purpose: integrate test process using VAMP"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and parameters"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()

    libname = "integrations_4_lib"
    procname = "integrations_4"

    call prepare_test_library (global, libname, 1, [procname])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true.)

    call global%pacify (efficiency_reset = .true., error_reset = .true.)
    call global%write (u, vars = [var_str ("$integration_method")], &
            pacify = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_4"

  end subroutine integrations_4

  subroutine integrations_5 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global

    write (u, "(A)")  "* Test output: integrations_5"
    write (u, "(A)")  "*   Purpose: integrate test process using VAMP"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and parameters"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()

    libname = "integrations_5_lib"
    procname = "integrations_5"

    call prepare_test_library (global, libname, 1, [procname])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([3], [1000])

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true.)

    call global%pacify (efficiency_reset = .true., error_reset = .true.)
    call global%write (u, vars = [var_str ("$integration_method")], &
            pacify = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_5"

  end subroutine integrations_5

  subroutine integrations_6 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global
    type(string_t), dimension(0) :: no_vars

    write (u, "(A)")  "* Test output: integrations_6"
    write (u, "(A)")  "*   Purpose: integrate test process using VAMP"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and parameters"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()

    libname = "integrations_6_lib"
    procname = "integrations_6"

    call prepare_test_library (global, libname, 1, [procname])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([3, 3, 3], [1000, 1000, 1000], &
         adapt = [.true., .true., .false.], &
         adapt_code = [var_str ("wg"), var_str ("g"), var_str ("")])

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true.)

    call global%pacify (efficiency_reset = .true., error_reset = .true.)
    call global%write (u, vars = no_vars, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_6"

  end subroutine integrations_6

  subroutine integrations_7 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global
    type(string_t), dimension(0) :: no_vars
    integer :: iostat, u_phs
    character(95) :: buffer
    type(string_t) :: phs_file
    logical :: exist

    write (u, "(A)")  "* Test output: integrations_7"
    write (u, "(A)")  "*   Purpose: integrate test process using VAMP"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and parameters"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()

    libname = "integrations_7_lib"
    procname = "integrations_7"

    call prepare_test_library (global, libname, 1, [procname])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .true., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?phs_s_mapping"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([3, 3, 3], [1000, 1000, 1000], &
         adapt = [.true., .true., .false.], &
         adapt_code = [var_str ("wg"), var_str ("g"), var_str ("")])

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true.)

    call global%pacify (efficiency_reset = .true., error_reset = .true.)
    call global%write (u, vars = no_vars, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_phs_forest_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Generated phase-space file"
    write (u, "(A)")

    phs_file = procname // ".r1.i1.phs"
    inquire (file = char (phs_file), exist = exist)
    if (exist) then
       u_phs = free_unit ()
       open (u_phs, file = char (phs_file), action = "read", status = "old")
       iostat = 0
       do while (iostat == 0)
          read (u_phs, "(A)", iostat = iostat)  buffer
          if (iostat == 0)  write (u, "(A)")  trim (buffer)
       end do
       close (u_phs)
    else
       write (u, "(A)")  "[file is missing]"
    end if

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_7"

  end subroutine integrations_7

  subroutine integrations_8 (u)
    use kinds, only: default
    use iso_varying_string, string_t => varying_string
    use interactions, only: reset_interaction_counter
    use phs_forests
    use models
    use rt_data
    use process_configurations_ut, only: prepare_test_library
    use compilations, only: compile_library
    use integrations

    implicit none

    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global
    type(flavor_t) :: flv
    type(string_t) :: name

    write (u, "(A)")  "* Test output: integrations_8"
    write (u, "(A)")  "*   Purpose: integrate test process using VAMP &
         &with structure function"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and parameters"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()

    libname = "integrations_8_lib"
    procname = "integrations_8"

    call prepare_test_library (global, libname, 1, [procname])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .true., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?phs_s_mapping"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%model_set_real (var_str ("ms"), 0._default)

    call reset_interaction_counter ()

    call flv%init (25, global%model)

    name = flv%get_name ()
    call global%beam_structure%init_sf ([name, name], [1])
    call global%beam_structure%set_sf (1, 1, var_str ("sf_test_1"))

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call global%it_list%init ([1], [1000])
    call integrate_process (procname, global, local_stack=.true.)

    call global%write (u, vars = [var_str ("ms")])

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_phs_forest_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_8"

  end subroutine integrations_8

  subroutine integrations_9 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global

    type(string_t) :: wgt_expr_text
    type(ifile_t) :: ifile
    type(stream_t) :: stream
    type(parse_tree_t) :: parse_tree

    write (u, "(A)")  "* Test output: integrations_9"
    write (u, "(A)")  "*   Purpose: integrate test process"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()

    write (u, "(A)")  "* Prepare a weight expression"
    write (u, "(A)")

    call syntax_pexpr_init ()
    wgt_expr_text = "eval 2 * sgn (Pz) - 1 [s]"
    call ifile_append (ifile, wgt_expr_text)
    call stream_init (stream, ifile)
    call parse_tree_init_expr (parse_tree, stream, .true.)
    global%pn%weight_expr => parse_tree%get_root_ptr ()

    write (u, "(A)")  "* Build and evaluate a test process"
    write (u, "(A)")

    libname = "integration_9"
    procname = "prc_config_a"

    call prepare_test_library (global, libname, 1)
    call compile_library (libname, global)

    call global%set_string (var_str ("$run_id"), &
         var_str ("integrations1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true.)

    call global%write (u, vars = [ &
         var_str ("$method"), &
         var_str ("sqrts"), &
         var_str ("$integration_method"), &
         var_str ("$phs_method"), &
         var_str ("$run_id")])

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_9"

  end subroutine integrations_9

  subroutine integrations_history_1 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global
    type(string_t), dimension(0) :: no_vars
    integer :: iostat, u_his
    character(91) :: buffer
    type(string_t) :: his_file, ps_file, pdf_file
    logical :: exist, exist_ps, exist_pdf

    write (u, "(A)")  "* Test output: integrations_history_1"
    write (u, "(A)")  "*   Purpose: test integration history files"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and parameters"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()

    libname = "integrations_history_1_lib"
    procname = "integrations_history_1"

    call global%set_log (var_str ("?vis_history"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?phs_s_mapping"),&
         .false., is_known = .true.)

    call prepare_test_library (global, libname, 1, [procname])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .true., is_known = .true.)
    call global%set_real (var_str ("error_threshold"),&
         5E-6_default, is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known=.true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([2, 2, 2], [1000, 1000, 1000], &
         adapt = [.true., .true., .false.], &
         adapt_code = [var_str ("wg"), var_str ("g"), var_str ("")])

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call reset_interaction_counter ()
    call integrate_process (procname, global, local_stack=.true., &
         eff_reset = .true.)

    call global%pacify (efficiency_reset = .true., error_reset = .true.)
    call global%write (u, vars = no_vars, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Generated history files"
    write (u, "(A)")

    his_file = procname // ".r1.history.tex"
    ps_file  = procname // ".r1.history.ps"
    pdf_file = procname // ".r1.history.pdf"
    inquire (file = char (his_file), exist = exist)
    if (exist) then
       u_his = free_unit ()
       open (u_his, file = char (his_file), action = "read", status = "old")
       iostat = 0
       do while (iostat == 0)
          read (u_his, "(A)", iostat = iostat)  buffer
          if (iostat == 0)  write (u, "(A)")  trim (buffer)
       end do
       close (u_his)
    else
       write (u, "(A)")  "[History LaTeX file is missing]"
    end if
    inquire (file = char (ps_file), exist = exist_ps)
    if (exist_ps) then
       write (u, "(A)")  "[History Postscript file exists and is nonempty]"
    else
       write (u, "(A)")  "[History Postscript file is missing/non-regular]"
    end if
    inquire (file = char (pdf_file), exist = exist_pdf)
    if (exist_pdf) then
       write (u, "(A)")  "[History PDF file exists and is nonempty]"
    else
       write (u, "(A)")  "[History PDF file is missing/non-regular]"
    end if

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_phs_forest_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integrations_history_1"

  end subroutine integrations_history_1


end module integrations_uti

