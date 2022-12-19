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

module commands_uti

    use kinds, only: default
    use kinds, only: i64
    use iso_varying_string, string_t => varying_string
    use io_units
    use ifiles
    use parser
    use interactions, only: reset_interaction_counter
    use prclib_stacks
    use analysis
    use variables, only: var_list_t
    use models
    use slha_interface
    use rt_data
    use event_base, only: generic_event_t, event_callback_t
    use commands

  implicit none
  private

  public :: commands_1
  public :: commands_2
  public :: commands_3
  public :: commands_4
  public :: commands_5
  public :: commands_6
  public :: commands_7
  public :: commands_8
  public :: commands_9
  public :: commands_10
  public :: commands_11
  public :: commands_12
  public :: commands_13
  public :: commands_14
  public :: commands_15
  public :: commands_16
  public :: commands_17
  public :: commands_18
  public :: commands_19
  public :: commands_20
  public :: commands_21
  public :: commands_22
  public :: commands_23
  public :: commands_24
  public :: commands_25
  public :: commands_26
  public :: commands_27
  public :: commands_28
  public :: commands_29
  public :: commands_30
  public :: commands_31
  public :: commands_32
  public :: commands_33
  public :: commands_34

  type, extends (event_callback_t) :: event_callback_34_t
     private
     integer :: u = 0
   contains
     procedure :: write => event_callback_34_write
     procedure :: proc => event_callback_34
  end type event_callback_34_t


contains

  subroutine commands_1 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_1"
    write (u, "(A)")  "*   Purpose: compile and execute empty command list"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Parse empty file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"

    if (associated (pn_root)) then
       call command_list%compile (pn_root, global)
    end if

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"

    call global%activate ()
    call command_list%execute (global)
    call global%deactivate ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call syntax_cmd_list_final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_1"

  end subroutine commands_1

  subroutine commands_2 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_2"
    write (u, "(A)")  "*   Purpose: set model"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')

    call ifile_write (ifile, u)

    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_2"

  end subroutine commands_2

  subroutine commands_3 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: commands_3"
    write (u, "(A)")  "*   Purpose: define process"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()
    call global%var_list%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd3"))
    call global%add_prclib (lib)

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process t3 = s, s => s, s')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%prclib_stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_3"

  end subroutine commands_3

  subroutine commands_4 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: commands_4"
    write (u, "(A)")  "*   Purpose: define process and compile library"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()
    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd4"))
    call global%add_prclib (lib)

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process t4 = s, s => s, s')
    call ifile_append (ifile, 'compile ("lib_cmd4")')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%prclib_stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_4"

  end subroutine commands_4

  subroutine commands_5 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: commands_5"
    write (u, "(A)")  "*   Purpose: define process, iterations, and integrate"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()
    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)
    call global%var_list%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known=.true.)
    call global%var_list%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known=.true.)
    call global%var_list%set_log (var_str ("?vis_history"),&
         .false., is_known=.true.)
    call global%var_list%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%var_list%set_real (var_str ("sqrts"), &
         1000._default, is_known=.true.)
    call global%var_list%set_int (var_str ("seed"), 0, is_known=.true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd5"))
    call global%add_prclib (lib)

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process t5 = s, s => s, s')
    call ifile_append (ifile, 'compile')
    call ifile_append (ifile, 'iterations = 1:1000')
    call ifile_append (ifile, 'integrate (t5)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call reset_interaction_counter ()
    call command_list%execute (global)

    call global%it_list%write (u)
    write (u, "(A)")
    call global%process_stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_5"

  end subroutine commands_5

  subroutine commands_6 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_6"
    write (u, "(A)")  "*   Purpose: define and set variables"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()
    call global%write_vars (u, [ &
         var_str ("$run_id"), &
         var_str ("?unweighted"), &
         var_str ("sqrts")])

    write (u, "(A)")
    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, '$run_id = "run1"')
    call ifile_append (ifile, '?unweighted = false')
    call ifile_append (ifile, 'sqrts = 1000')
    call ifile_append (ifile, 'int j = 10')
    call ifile_append (ifile, 'real x = 1000.')
    call ifile_append (ifile, 'complex z = 5')
    call ifile_append (ifile, 'string $text = "abcd"')
    call ifile_append (ifile, 'logical ?flag = true')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write_vars (u, [ &
         var_str ("$run_id"), &
         var_str ("?unweighted"), &
         var_str ("sqrts"), &
         var_str ("j"), &
         var_str ("x"), &
         var_str ("z"), &
         var_str ("$text"), &
         var_str ("?flag")])


    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call syntax_cmd_list_final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_6"

  end subroutine commands_6

  subroutine commands_7 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_7"
    write (u, "(A)")  "*   Purpose: declare process libraries"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()
    call global%var_list%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    global%os_data%fc = "Fortran-compiler"
    global%os_data%fcflags = "Fortran-flags"
    global%os_data%fclibs = "Fortran-libs"

    write (u, "(A)")
    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'library = "lib_cmd7_1"')
    call ifile_append (ifile, 'library = "lib_cmd7_2"')
    call ifile_append (ifile, 'library = "lib_cmd7_1"')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write_libraries (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call syntax_cmd_list_final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_7"

  end subroutine commands_7

  subroutine commands_8 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: commands_8"
    write (u, "(A)")  "*   Purpose: define process, integrate, generate events"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)
    call global%var_list%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known=.true.)
    call global%var_list%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known=.true.)
    call global%var_list%set_log (var_str ("?vis_history"),&
         .false., is_known=.true.)
    call global%var_list%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%var_list%set_real (var_str ("sqrts"), &
         1000._default, is_known=.true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd8"))
    call global%add_prclib (lib)

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process commands_8_p = s, s => s, s')
    call ifile_append (ifile, 'compile')
    call ifile_append (ifile, 'iterations = 1:1000')
    call ifile_append (ifile, 'integrate (commands_8_p)')
    call ifile_append (ifile, '?unweighted = false')
    call ifile_append (ifile, 'n_events = 3')
    call ifile_append (ifile, '?read_raw = false')
    call ifile_append (ifile, 'simulate (commands_8_p)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"

    call command_list%execute (global)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_8"

  end subroutine commands_8

  subroutine commands_9 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(string_t), dimension(0) :: no_vars

    write (u, "(A)")  "* Test output: commands_9"
    write (u, "(A)")  "*   Purpose: define cuts"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'cuts = all Pt > 0 [particle]')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write (u, vars = no_vars)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_9"

  end subroutine commands_9

  subroutine commands_10 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_10"
    write (u, "(A)")  "*   Purpose: define beams"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = QCD')
    call ifile_append (ifile, 'sqrts = 1000')
    call ifile_append (ifile, 'beams = p, p')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write_beams (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_10"

  end subroutine commands_10

  subroutine commands_11 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_11"
    write (u, "(A)")  "*   Purpose: define beams with structure functions"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = QCD')
    call ifile_append (ifile, 'sqrts = 1100')
    call ifile_append (ifile, 'beams = p, p => lhapdf => pdf_builtin, isr')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write_beams (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_11"

  end subroutine commands_11

  subroutine commands_12 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: commands_12"
    write (u, "(A)")  "*   Purpose: generate events and rescan"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()

    call global%global_init ()
    call global%var_list%append_log (&
         var_str ("?rebuild_phase_space"), .false., &
         intrinsic=.true.)
    call global%var_list%append_log (&
         var_str ("?rebuild_grids"), .false., &
         intrinsic=.true.)
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)
    call global%var_list%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known=.true.)
    call global%var_list%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known=.true.)
    call global%var_list%set_log (var_str ("?vis_history"),&
         .false., is_known=.true.)
    call global%var_list%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%var_list%set_real (var_str ("sqrts"), &
         1000._default, is_known=.true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd12"))
    call global%add_prclib (lib)

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process commands_12_p = s, s => s, s')
    call ifile_append (ifile, 'compile')
    call ifile_append (ifile, 'iterations = 1:1000')
    call ifile_append (ifile, 'integrate (commands_12_p)')
    call ifile_append (ifile, '?unweighted = false')
    call ifile_append (ifile, 'n_events = 3')
    call ifile_append (ifile, '?read_raw = false')
    call ifile_append (ifile, 'simulate (commands_12_p)')
    call ifile_append (ifile, '?write_raw = false')
    call ifile_append (ifile, 'rescan "commands_12_p" (commands_12_p)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"

    call command_list%execute (global)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_12"

  end subroutine commands_12

  subroutine commands_13 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib
    logical :: exist

    write (u, "(A)")  "* Test output: commands_13"
    write (u, "(A)")  "*   Purpose: generate events and rescan"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)
    call global%var_list%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known=.true.)
    call global%var_list%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known=.true.)
    call global%var_list%set_real (var_str ("sqrts"), &
         1000._default, is_known=.true.)
    call global%var_list%set_log (var_str ("?vis_history"),&
         .false., is_known=.true.)
    call global%var_list%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd13"))
    call global%add_prclib (lib)

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process commands_13_p = s, s => s, s')
    call ifile_append (ifile, 'compile')
    call ifile_append (ifile, 'iterations = 1:1000')
    call ifile_append (ifile, 'integrate (commands_13_p)')
    call ifile_append (ifile, '?unweighted = false')
    call ifile_append (ifile, 'n_events = 1')
    call ifile_append (ifile, '?read_raw = false')
    call ifile_append (ifile, 'sample_format = weight_stream')
    call ifile_append (ifile, 'simulate (commands_13_p)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"

    call command_list%execute (global)

    write (u, "(A)")
    write (u, "(A)")  "* Verify output files"
    write (u, "(A)")

    inquire (file = "commands_13_p.evx", exist = exist)
    if (exist)  write (u, "(1x,A)")  "raw"

    inquire (file = "commands_13_p.weights.dat", exist = exist)
    if (exist)  write (u, "(1x,A)")  "weight_stream"

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_13"

  end subroutine commands_13

  subroutine commands_14 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_14"
    write (u, "(A)")  "*   Purpose: define and compile empty libraries"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_cmd_list_init ()

    call global%global_init ()

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'library = "lib1"')
    call ifile_append (ifile, 'library = "lib2"')
    call ifile_append (ifile, 'compile ()')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%prclib_stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()

    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_14"

  end subroutine commands_14

  subroutine commands_15 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: commands_15"
    write (u, "(A)")  "*   Purpose: define process and compile library"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()
    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)
    call global%var_list%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known=.true.)
    call global%var_list%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known=.true.)
    call global%var_list%set_real (var_str ("sqrts"), &
         1000._default, is_known=.true.)
    call global%var_list%set_log (var_str ("?vis_history"),&
         .false., is_known=.true.)
    call global%var_list%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd15"))
    call global%add_prclib (lib)

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process t15 = s, s => s, s')
    call ifile_append (ifile, 'iterations = 1:1000')
    call ifile_append (ifile, 'integrate (t15)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%prclib_stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_15"

  end subroutine commands_15

  subroutine commands_16 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_16"
    write (u, "(A)")  "*   Purpose: declare an observable"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, '$obs_label = "foo"')
    call ifile_append (ifile, '$obs_unit = "cm"')
    call ifile_append (ifile, '$title = "Observable foo"')
    call ifile_append (ifile, '$description = "This is observable foo"')
    call ifile_append (ifile, 'observable foo')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Record two data items"
    write (u, "(A)")

    call analysis_record_data (var_str ("foo"), 1._default)
    call analysis_record_data (var_str ("foo"), 3._default)

    write (u, "(A)")  "* Display analysis store"
    write (u, "(A)")

    call analysis_write (u, verbose=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_16"

  end subroutine commands_16

  subroutine commands_17 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(string_t), dimension(3) :: name
    integer :: i

    write (u, "(A)")  "* Test output: commands_17"
    write (u, "(A)")  "*   Purpose: declare histograms"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, '$obs_label = "foo"')
    call ifile_append (ifile, '$obs_unit = "cm"')
    call ifile_append (ifile, '$title = "Histogram foo"')
    call ifile_append (ifile, '$description = "This is histogram foo"')
    call ifile_append (ifile, 'histogram foo (0,5,1)')
    call ifile_append (ifile, '$title = "Histogram bar"')
    call ifile_append (ifile, '$description = "This is histogram bar"')
    call ifile_append (ifile, 'n_bins = 2')
    call ifile_append (ifile, 'histogram bar (0,5)')
    call ifile_append (ifile, '$title = "Histogram gee"')
    call ifile_append (ifile, '$description = "This is histogram gee"')
    call ifile_append (ifile, '?normalize_bins = true')
    call ifile_append (ifile, 'histogram gee (0,5)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Record two data items"
    write (u, "(A)")

    name(1) = "foo"
    name(2) = "bar"
    name(3) = "gee"

    do i = 1, 3
       call analysis_record_data (name(i), 0.1_default, &
            weight = 0.25_default)
       call analysis_record_data (name(i), 3.1_default)
       call analysis_record_data (name(i), 4.1_default, &
            excess = 0.5_default)
       call analysis_record_data (name(i), 7.1_default)
    end do

    write (u, "(A)")  "* Display analysis store"
    write (u, "(A)")

    call analysis_write (u, verbose=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_17"

  end subroutine commands_17

  subroutine commands_18 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_18"
    write (u, "(A)")  "*   Purpose: declare a plot"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, '$obs_label = "foo"')
    call ifile_append (ifile, '$obs_unit = "cm"')
    call ifile_append (ifile, '$title = "Plot foo"')
    call ifile_append (ifile, '$description = "This is plot foo"')
    call ifile_append (ifile, '$x_label = "x axis"')
    call ifile_append (ifile, '$y_label = "y axis"')
    call ifile_append (ifile, '?x_log = false')
    call ifile_append (ifile, '?y_log = true')
    call ifile_append (ifile, 'x_min = -1')
    call ifile_append (ifile, 'x_max = 1')
    call ifile_append (ifile, 'y_min = 0.1')
    call ifile_append (ifile, 'y_max = 1000')
    call ifile_append (ifile, 'plot foo')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Record two data items"
    write (u, "(A)")

    call analysis_record_data (var_str ("foo"), 0._default, 20._default, &
         xerr = 0.25_default)
    call analysis_record_data (var_str ("foo"), 0.5_default, 0.2_default, &
         yerr = 0.07_default)
    call analysis_record_data (var_str ("foo"), 3._default, 2._default)

    write (u, "(A)")  "* Display analysis store"
    write (u, "(A)")

    call analysis_write (u, verbose=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_18"

  end subroutine commands_18

  subroutine commands_19 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_19"
    write (u, "(A)")  "*   Purpose: combine two plots to a graph"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'plot a')
    call ifile_append (ifile, 'plot b')
    call ifile_append (ifile, '$title = "Graph foo"')
    call ifile_append (ifile, '$description = "This is graph foo"')
    call ifile_append (ifile, 'graph foo = a & b')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Display analysis object"
    write (u, "(A)")

    call analysis_write (var_str ("foo"), u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_19"

  end subroutine commands_19

  subroutine commands_20 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_20"
    write (u, "(A)")  "*   Purpose: record data"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization: create observable, histogram, plot"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    call analysis_init_observable (var_str ("o"))
    call analysis_init_histogram (var_str ("h"), 0._default, 1._default, 3, &
         normalize_bins = .false.)
    call analysis_init_plot (var_str ("p"))

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'record o (1.234)')
    call ifile_append (ifile, 'record h (0.5)')
    call ifile_append (ifile, 'record p (1, 2)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Display analysis object"
    write (u, "(A)")

    call analysis_write (u, verbose = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_20"

  end subroutine commands_20

  subroutine commands_21 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: commands_21"
    write (u, "(A)")  "*   Purpose: create and use analysis expression"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization: create observable"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)
    call global%var_list%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known=.true.)
    call global%var_list%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known=.true.)
    call global%var_list%set_log (var_str ("?vis_history"),&
         .false., is_known=.true.)
    call global%var_list%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%var_list%set_real (var_str ("sqrts"), &
         1000._default, is_known=.true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd8"))
    call global%add_prclib (lib)

    call analysis_init_observable (var_str ("m"))

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process commands_21_p = s, s => s, s')
    call ifile_append (ifile, 'compile')
    call ifile_append (ifile, 'iterations = 1:100')
    call ifile_append (ifile, 'integrate (commands_21_p)')
    call ifile_append (ifile, '?unweighted = true')
    call ifile_append (ifile, 'n_events = 3')
    call ifile_append (ifile, '?read_raw = false')
    call ifile_append (ifile, 'observable m')
    call ifile_append (ifile, 'analysis = record m (eval M [s])')
    call ifile_append (ifile, 'simulate (commands_21_p)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Display analysis object"
    write (u, "(A)")

    call analysis_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_21"

  end subroutine commands_21

  subroutine commands_22 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    integer :: u_file, iostat
    logical :: exist
    character(80) :: buffer

    write (u, "(A)")  "* Test output: commands_22"
    write (u, "(A)")  "*   Purpose: write analysis data"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization: create observable"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    call analysis_init_observable (var_str ("m"))
    call analysis_record_data (var_str ("m"), 125._default)

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, '$out_file = "commands_22.dat"')
    call ifile_append (ifile, 'write_analysis')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Display analysis data"
    write (u, "(A)")

    inquire (file = "commands_22.dat", exist = exist)
    if (.not. exist) then
       write (u, "(A)")  "ERROR: File commands_22.dat not found"
       return
    end if

    u_file = free_unit ()
    open (u_file, file = "commands_22.dat", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_22"

  end subroutine commands_22

  subroutine commands_23 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    integer :: u_file, iostat
    character(256) :: buffer
    logical :: exist
    type(graph_options_t) :: graph_options

    write (u, "(A)")  "* Test output: commands_23"
    write (u, "(A)")  "*   Purpose: write and compile analysis data"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization: create and fill histogram"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    call graph_options%init ()
    call graph_options%set (title = var_str ("Histogram for test: commands 23"), &
         description = var_str ("This is a test."), &
         width_mm = 125, height_mm = 85)
    call analysis_init_histogram (var_str ("h"), &
         0._default, 10._default, 2._default, .false., &
         graph_options = graph_options)
    call analysis_record_data (var_str ("h"), 1._default)
    call analysis_record_data (var_str ("h"), 1._default)
    call analysis_record_data (var_str ("h"), 1._default)
    call analysis_record_data (var_str ("h"), 1._default)
    call analysis_record_data (var_str ("h"), 3._default)
    call analysis_record_data (var_str ("h"), 3._default)
    call analysis_record_data (var_str ("h"), 3._default)
    call analysis_record_data (var_str ("h"), 5._default)
    call analysis_record_data (var_str ("h"), 7._default)
    call analysis_record_data (var_str ("h"), 7._default)
    call analysis_record_data (var_str ("h"), 7._default)
    call analysis_record_data (var_str ("h"), 7._default)
    call analysis_record_data (var_str ("h"), 9._default)
    call analysis_record_data (var_str ("h"), 9._default)
    call analysis_record_data (var_str ("h"), 9._default)
    call analysis_record_data (var_str ("h"), 9._default)
    call analysis_record_data (var_str ("h"), 9._default)
    call analysis_record_data (var_str ("h"), 9._default)
    call analysis_record_data (var_str ("h"), 9._default)

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, '$out_file = "commands_23.dat"')
    call ifile_append (ifile, 'compile_analysis')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Delete Postscript output"
    write (u, "(A)")

    inquire (file = "commands_23.ps", exist = exist)
    if (exist) then
       u_file = free_unit ()
       open (u_file, file = "commands_23.ps", action = "write", status = "old")
       close (u_file, status = "delete")
    end if
    inquire (file = "commands_23.ps", exist = exist)
    write (u, "(1x,A,L1)")  "Postcript output exists = ", exist

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* TeX file"
    write (u, "(A)")

    inquire (file = "commands_23.tex", exist = exist)
    if (.not. exist) then
       write (u, "(A)")  "ERROR: File commands_23.tex not found"
       return
    end if

    u_file = free_unit ()
    open (u_file, file = "commands_23.tex", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)
    write (u, *)

    inquire (file = "commands_23.ps", exist = exist)
    write (u, "(1x,A,L1)")  "Postcript output exists = ", exist

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_23"

  end subroutine commands_23

  subroutine commands_24 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_24"
    write (u, "(A)")  "*   Purpose: check graph and drawing options"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, '$title = "Title"')
    call ifile_append (ifile, '$description = "Description"')
    call ifile_append (ifile, '$x_label = "X Label"')
    call ifile_append (ifile, '$y_label = "Y Label"')
    call ifile_append (ifile, 'graph_width_mm = 111')
    call ifile_append (ifile, 'graph_height_mm = 222')
    call ifile_append (ifile, 'x_min = -11')
    call ifile_append (ifile, 'x_max = 22')
    call ifile_append (ifile, 'y_min = -33')
    call ifile_append (ifile, 'y_max = 44')
    call ifile_append (ifile, '$gmlcode_bg = "GML Code BG"')
    call ifile_append (ifile, '$gmlcode_fg = "GML Code FG"')
    call ifile_append (ifile, '$fill_options = "Fill Options"')
    call ifile_append (ifile, '$draw_options = "Draw Options"')
    call ifile_append (ifile, '$err_options = "Error Options"')
    call ifile_append (ifile, '$symbol = "Symbol"')
    call ifile_append (ifile, 'histogram foo (0,1)')
    call ifile_append (ifile, 'plot bar')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Display analysis store"
    write (u, "(A)")

    call analysis_write (u, verbose=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_24"

  end subroutine commands_24

  subroutine commands_25 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_25"
    write (u, "(A)")  "*   Purpose: declare local environment for process"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_cmd_list_init ()
    call global%global_init ()
    call global%var_list%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'library = "commands_25_lib"')
    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process commands_25_p1 = g, g => g, g &
         &{ model = "QCD" }')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)
    call global%write_libraries (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_25"

  end subroutine commands_25

  subroutine commands_26 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_26"
    write (u, "(A)")  "*   Purpose: declare alternative setups for simulation"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'int i = 0')
    call ifile_append (ifile, 'alt_setup = ({ i = 1 }, { i = 2 })')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write_expr (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_26"

  end subroutine commands_26

  subroutine commands_27 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: commands_27"
    write (u, "(A)")  "*   Purpose: modify particle properties"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call global%global_init ()
    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)
    call global%var_list%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known=.true.)
    call global%var_list%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known=.true.)
    call global%var_list%set_log (var_str ("?vis_history"),&
         .false., is_known=.true.)
    call global%var_list%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)

    allocate (lib)
    call lib%init (var_str ("commands_27_lib"))
    call global%add_prclib (lib)

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'ff = 0.4')
    call ifile_append (ifile, 'process d1 = s => f, fbar')
    call ifile_append (ifile, 'unstable s (d1)')
    call ifile_append (ifile, 'polarized f, fbar')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Show model"
    write (u, "(A)")

    call global%model%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extra Input"
    write (u, "(A)")

    call ifile_final (ifile)
    call ifile_append (ifile, '?diagonal_decay = true')
    call ifile_append (ifile, 'unstable s (d1)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%final ()
    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Show model"
    write (u, "(A)")

    call global%model%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extra Input"
    write (u, "(A)")

    call ifile_final (ifile)
    call ifile_append (ifile, '?isotropic_decay = true')
    call ifile_append (ifile, 'unstable s (d1)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%final ()
    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Show model"
    write (u, "(A)")

    call global%model%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extra Input"
    write (u, "(A)")

    call ifile_final (ifile)
    call ifile_append (ifile, 'stable s')
    call ifile_append (ifile, 'unpolarized f')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%final ()
    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Show model"
    write (u, "(A)")

    call global%model%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_model_file_init ()
    call syntax_cmd_list_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_27"

  end subroutine commands_27

  subroutine commands_28 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root1, pn_root2
    type(string_t), dimension(0) :: no_vars

    write (u, "(A)")  "* Test output: commands_28"
    write (u, "(A)")  "*   Purpose: quit the program"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file: quit without code"
    write (u, "(A)")

    call ifile_append (ifile, 'quit')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root1, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root1, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write (u, vars = no_vars)

    write (u, "(A)")
    write (u, "(A)")  "*  Input file: quit with code"
    write (u, "(A)")

    call ifile_final (ifile)
    call command_list%final ()
    call ifile_append (ifile, 'quit ( 3 + 4 )')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root2, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root2, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write (u, vars = no_vars)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_28"

  end subroutine commands_28

  subroutine commands_29 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(var_list_t), pointer :: model_vars
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_29"
    write (u, "(A)")  "*   Purpose: test SLHA interface"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call syntax_model_file_init ()
    call syntax_slha_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Model MSSM, read SLHA file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "MSSM"')
    call ifile_append (ifile, '?slha_read_decays = true')
    call ifile_append (ifile, 'read_slha ("sps1ap_decays.slha")')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Model MSSM, default values:"
    write (u, "(A)")

    call global%model%write (u, verbose = .false., &
         show_vertices = .false., show_particles = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Selected global variables"
    write (u, "(A)")

    model_vars => global%model%get_var_list_ptr ()

    call model_vars%write_var (var_str ("mch1"), u)
    call model_vars%write_var (var_str ("wch1"), u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")  "* Model MSSM, values from SLHA file"
    write (u, "(A)")

    call global%model%write (u, verbose = .false., &
         show_vertices = .false., show_particles = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Selected global variables"
    write (u, "(A)")

    model_vars => global%model%get_var_list_ptr ()

    call model_vars%write_var (var_str ("mch1"), u)
    call model_vars%write_var (var_str ("wch1"), u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_slha_final ()
    call syntax_model_file_final ()
    call syntax_cmd_list_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_29"

  end subroutine commands_29

  subroutine commands_30 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_30"
    write (u, "(A)")  "*   Purpose: define scales"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'scale = 200 GeV')
    call ifile_append (ifile, &
         'factorization_scale = eval Pt [particle]')
    call ifile_append (ifile, &
         'renormalization_scale = eval E [particle]')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write_expr (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_30"

  end subroutine commands_30

  subroutine commands_31 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_31"
    write (u, "(A)")  "*   Purpose: define weight/reweight"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'weight = eval Pz [particle]')
    call ifile_append (ifile, 'reweight = eval M2 [particle]')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write_expr (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_31"

  end subroutine commands_31

  subroutine commands_32 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root

    write (u, "(A)")  "* Test output: commands_32"
    write (u, "(A)")  "*   Purpose: define selection"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'selection = any PDG == 13 [particle]')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    call global%write_expr (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_32"

  end subroutine commands_32

  subroutine commands_33 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    integer :: u_file, iostat
    character(3) :: buffer

    write (u, "(A)")  "* Test output: commands_33"
    write (u, "(A)")  "*   Purpose: execute shell command"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    write (u, "(A)")  "*  Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'exec ("echo foo >> bar")')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "*  Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root, u)

    write (u, "(A)")
    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)
    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)
    u_file = free_unit ()
    open (u_file, file = "bar", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (iostat /= 0) exit
    end do
    write (u, "(A,A)")  "should be 'foo': ", trim (buffer)
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_33"

  end subroutine commands_33

  subroutine commands_34 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(command_list_t), target :: command_list
    type(rt_data_t), target :: global
    type(parse_node_t), pointer :: pn_root
    type(prclib_entry_t), pointer :: lib
    type(event_callback_34_t) :: event_callback

    write (u, "(A)")  "* Test output: commands_34"
    write (u, "(A)")  "*   Purpose: write analysis data"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization: create observable"
    write (u, "(A)")

    call syntax_cmd_list_init ()
    call global%global_init ()

    call syntax_model_file_init ()
    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%var_list%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known=.true.)
    call global%var_list%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known=.true.)
    call global%var_list%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known=.true.)
    call global%var_list%set_real (var_str ("sqrts"), &
         1000._default, is_known=.true.)
    call global%var_list%set_log (var_str ("?vis_history"),&
         .false., is_known=.true.)
    call global%var_list%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%var_list%set_int (var_str ("seed"), 0, is_known=.true.)

    allocate (lib)
    call lib%init (var_str ("lib_cmd34"))
    call global%add_prclib (lib)

    write (u, "(A)")  "* Prepare callback for writing analysis to I/O unit"
    write (u, "(A)")

    event_callback%u = u
    call global%set_event_callback (event_callback)

    write (u, "(A)")  "* Input file"
    write (u, "(A)")

    call ifile_append (ifile, 'model = "Test"')
    call ifile_append (ifile, 'process commands_34_p = s, s => s, s')
    call ifile_append (ifile, 'compile')
    call ifile_append (ifile, 'iterations = 1:1000')
    call ifile_append (ifile, 'integrate (commands_34_p)')
    call ifile_append (ifile, 'observable sq')
    call ifile_append (ifile, 'analysis = record sq (sqrts)')
    call ifile_append (ifile, 'n_events = 4')
    call ifile_append (ifile, 'event_callback_interval = 3')
    call ifile_append (ifile, 'simulate (commands_34_p)')

    call ifile_write (ifile, u)

    write (u, "(A)")
    write (u, "(A)")  "* Parse file"
    write (u, "(A)")

    call parse_ifile (ifile, pn_root)

    write (u, "(A)")  "* Compile command list"
    write (u, "(A)")

    call command_list%compile (pn_root, global)

    call command_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Execute command list"
    write (u, "(A)")

    call command_list%execute (global)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call ifile_final (ifile)

    call analysis_final ()
    call command_list%final ()
    call global%final ()
    call syntax_cmd_list_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: commands_34"

  end subroutine commands_34


  subroutine parse_ifile (ifile, pn_root, u)
    use ifiles
    use lexers
    use parser
    use commands
    type(ifile_t), intent(in) :: ifile
    type(parse_node_t), pointer, intent(out) :: pn_root
    integer, intent(in), optional :: u
    type(stream_t), target :: stream
    type(lexer_t), target :: lexer
    type(parse_tree_t) :: parse_tree

    call lexer_init_cmd_list (lexer)
    call stream_init (stream, ifile)
    call lexer_assign_stream (lexer, stream)

    call parse_tree_init (parse_tree, syntax_cmd_list, lexer)
    if (present (u))  call parse_tree_write (parse_tree, u)
    pn_root => parse_tree%get_root_ptr ()

    call stream_final (stream)
    call lexer_final (lexer)
  end subroutine parse_ifile

  subroutine event_callback_34_write (event_callback, unit)
    class(event_callback_34_t), intent(in) :: event_callback
    integer, intent(in), optional :: unit
  end subroutine event_callback_34_write

  subroutine event_callback_34 (event_callback, i, event)
    class(event_callback_34_t), intent(in) :: event_callback
    integer(i64), intent(in) :: i
    class(generic_event_t), intent(in) :: event
    call analysis_write (event_callback%u)
  end subroutine event_callback_34


end module commands_uti

