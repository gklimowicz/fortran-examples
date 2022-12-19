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

module process_stacks_uti

  use iso_varying_string, string_t => varying_string
  use os_interface
  use sm_qcd
  use models
  use model_data
  use variables, only: var_list_t
  use process_libraries
  use rng_base
  use prc_test, only: prc_test_create_library
  use process, only: process_t
  use instances, only: process_instance_t
  use processes_ut, only: prepare_test_process

  use process_stacks

  use rng_base_ut, only: rng_test_factory_t

  implicit none
  private

  public :: process_stacks_1
  public :: process_stacks_2
  public :: process_stacks_3
  public :: process_stacks_4

contains

  subroutine process_stacks_1 (u)
    integer, intent(in) :: u
    type(process_stack_t) :: stack

    write (u, "(A)")  "* Test output: process_stacks_1"
    write (u, "(A)")  "*   Purpose: display an empty process stack"
    write (u, "(A)")

    call stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_stacks_1"

  end subroutine process_stacks_1

  subroutine process_stacks_2 (u)
    integer, intent(in) :: u
    type(process_stack_t) :: stack
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(var_list_t) :: var_list
    type(process_entry_t), pointer :: process => null ()

    write (u, "(A)")  "* Test output: process_stacks_2"
    write (u, "(A)")  "*   Purpose: fill a process stack"
    write (u, "(A)")

    write (u, "(A)")  "* Build, initialize and store two test processes"
    write (u, "(A)")

    libname = "process_stacks2"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    call model%init_test ()
    call var_list%append_string (var_str ("$run_id"))
    call var_list%append_log (var_str ("?alphas_is_fixed"), .true.)
    call var_list%append_int (var_str ("seed"), 0)

    allocate (process)

    call var_list%set_string &
         (var_str ("$run_id"), var_str ("run1"), is_known=.true.)
    call process%init (procname, lib, os_data, model, var_list)
    call stack%push (process)

    allocate (process)

    call var_list%set_string &
         (var_str ("$run_id"), var_str ("run2"), is_known=.true.)
    call process%init (procname, lib, os_data, model, var_list)
    call stack%push (process)

    call stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call stack%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_stacks_2"

  end subroutine process_stacks_2

  subroutine process_stacks_3 (u)
    integer, intent(in) :: u
    type(process_stack_t) :: stack
    type(model_t), target :: model
    type(string_t) :: procname
    type(process_entry_t), pointer :: process => null ()
    type(process_instance_t), target :: process_instance

    write (u, "(A)")  "* Test output: process_stacks_3"
    write (u, "(A)")  "*   Purpose: setup process variables"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process variables"
    write (u, "(A)")

    procname = "processes_test"
    call model%init_test ()

    write (u, "(A)")  "* Initialize process variables"
    write (u, "(A)")

    call stack%init_var_list ()
    call stack%init_result_vars (procname)
    call stack%write_var_list (u)

    write (u, "(A)")
    write (u, "(A)")  "* Build and integrate a test process"
    write (u, "(A)")

    allocate (process)
    call prepare_test_process (process%process_t, process_instance, model)
    call process_instance%integrate (1, 1, 1000)
    call process_instance%final ()
    call process%final_integration (1)
    call stack%push (process)

    write (u, "(A)")  "* Fill process variables"
    write (u, "(A)")

    call stack%fill_result_vars (procname)
    call stack%write_var_list (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call stack%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_stacks_3"

  end subroutine process_stacks_3

  subroutine process_stacks_4 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(process_stack_t), target :: stack1, stack2
    type(model_t), target :: model
    type(string_t) :: libname
    type(string_t) :: procname1, procname2
    type(os_data_t) :: os_data
    type(process_entry_t), pointer :: process => null ()

    write (u, "(A)")  "* Test output: process_stacks_4"
    write (u, "(A)")  "*   Purpose: link process stacks"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process variables"
    write (u, "(A)")

    libname = "process_stacks_4_lib"
    procname1 = "process_stacks_4a"
    procname2 = "process_stacks_4b"

    call os_data%init ()

    write (u, "(A)")  "* Initialize first process"
    write (u, "(A)")

    call prc_test_create_library (procname1, lib)

    call model%init_test ()

    allocate (process)
    call process%init (procname1, lib, os_data, model)
    call stack1%push (process)

    write (u, "(A)")  "* Initialize second process"
    write (u, "(A)")

    call stack2%link (stack1)

    call prc_test_create_library (procname2, lib)

    allocate (process)

    call process%init (procname2, lib, os_data, model)
    call stack2%push (process)

    write (u, "(A)")  "* Show linked stacks"
    write (u, "(A)")

    call stack2%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call stack2%final ()
    call stack1%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_stacks_4"

  end subroutine process_stacks_4


end module process_stacks_uti

