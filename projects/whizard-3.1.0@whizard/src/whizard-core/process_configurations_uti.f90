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

module process_configurations_uti

  use iso_varying_string, string_t => varying_string
  use particle_specifiers, only: new_prt_spec
  use prclib_stacks
  use models
  use rt_data

  use process_configurations

  implicit none
  private

  public :: process_configurations_1
  public :: process_configurations_2

  public :: prepare_test_library

contains

  subroutine prepare_test_library (global, libname, mode, procname)
    type(rt_data_t), intent(inout), target :: global
    type(string_t), intent(in) :: libname
    integer, intent(in) :: mode
    type(string_t), intent(in), dimension(:), optional :: procname
    type(prclib_entry_t), pointer :: lib
    type(string_t) :: prc_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    integer :: n_components
    type(process_configuration_t) :: prc_config

    if (.not. associated (global%prclib_stack%get_first_ptr ())) then
       allocate (lib)
       call lib%init (libname)
       call global%add_prclib (lib)
    end if

    if (btest (mode, 0)) then

       call global%select_model (var_str ("Test"))

       if (present (procname)) then
          prc_name = procname(1)
       else
          prc_name = "prc_config_a"
       end if
       n_components = 1
       allocate (prt_in (2), prt_out (2))
       prt_in = [var_str ("s"), var_str ("s")]
       prt_out = [var_str ("s"), var_str ("s")]

       call global%set_string (var_str ("$method"),&
            var_str ("unit_test"), is_known = .true.)

       call prc_config%init (prc_name, &
            size (prt_in), n_components, &
            global%model, global%var_list)
       call prc_config%setup_component (1, &
            new_prt_spec (prt_in), new_prt_spec (prt_out), &
            global%model, global%var_list)
       call prc_config%record (global)

       deallocate (prt_in, prt_out)

    end if

    if (btest (mode, 1)) then

       call global%select_model (var_str ("QED"))

       if (present (procname)) then
          prc_name = procname(2)
       else
          prc_name = "prc_config_b"
       end if
       n_components = 1
       allocate (prt_in (2), prt_out (2))
       prt_in = [var_str ("e+"), var_str ("e-")]
       prt_out = [var_str ("m+"), var_str ("m-")]

       call global%set_string (var_str ("$method"),&
            var_str ("omega"), is_known = .true.)

       call prc_config%init (prc_name, &
            size (prt_in), n_components, &
            global%model, global%var_list)
       call prc_config%setup_component (1, &
            new_prt_spec (prt_in), new_prt_spec (prt_out), &
            global%model, global%var_list)
       call prc_config%record (global)

       deallocate (prt_in, prt_out)

    end if

    if (btest (mode, 2)) then

       call global%select_model (var_str ("Test"))

       if (present (procname)) then
          prc_name = procname(1)
       else
          prc_name = "prc_config_a"
       end if
       n_components = 1
       allocate (prt_in (1), prt_out (2))
       prt_in = [var_str ("s")]
       prt_out = [var_str ("f"), var_str ("fbar")]

       call global%set_string (var_str ("$method"),&
            var_str ("unit_test"), is_known = .true.)

       call prc_config%init (prc_name, &
            size (prt_in), n_components, &
            global%model, global%var_list)
       call prc_config%setup_component (1, &
            new_prt_spec (prt_in), new_prt_spec (prt_out), &
            global%model, global%var_list)
       call prc_config%record (global)

       deallocate (prt_in, prt_out)

    end if

  end subroutine prepare_test_library


  subroutine process_configurations_1 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global

    write (u, "(A)")  "* Test output: process_configurations_1"
    write (u, "(A)")  "*   Purpose: configure test processes"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    write (u, "(A)")  "* Configure processes as prc_test, model Test"
    write (u, "(A)")  "*                     and omega, model QED"
    write (u, *)

    call global%set_int (var_str ("process_num_id"), &
         42, is_known = .true.)
    call prepare_test_library (global, var_str ("prc_config_lib_1"), 3)

    global%os_data%fc = "Fortran-compiler"
    global%os_data%fcflags = "Fortran-flags"
    global%os_data%fclibs = "Fortran-libs"

    call global%write_libraries (u)

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_configurations_1"

  end subroutine process_configurations_1

  subroutine process_configurations_2 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global

    type(string_t) :: libname
    type(prclib_entry_t), pointer :: lib
    type(string_t) :: prc_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    integer :: n_components
    type(process_configuration_t) :: prc_config

    write (u, "(A)")  "* Test output: process_configurations_2"
    write (u, "(A)")  "*   Purpose: configure test processes with options"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()

    write (u, "(A)")  "* Configure processes as omega, model QED"
    write (u, *)

    libname = "prc_config_lib_2"

    allocate (lib)
    call lib%init (libname)
    call global%add_prclib (lib)

    call global%select_model (var_str ("QED"))

    prc_name = "prc_config_c"
    n_components = 2
    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("e+"), var_str ("e-")]
    prt_out = [var_str ("m+"), var_str ("m-")]

    call global%set_string (var_str ("$method"),&
         var_str ("omega"), is_known = .true.)
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    call prc_config%init (prc_name, size (prt_in), n_components, &
         global%model, global%var_list)

    call global%set_log (var_str ("?report_progress"), &
         .true., is_known = .true.)
    call prc_config%setup_component (1, &
         new_prt_spec (prt_in), new_prt_spec (prt_out), global%model, global%var_list)

    call global%set_log (var_str ("?report_progress"), &
         .false., is_known = .true.)
    call global%set_log (var_str ("?omega_openmp"), &
         .true., is_known = .true.)
    call global%set_string (var_str ("$restrictions"),&
         var_str ("3+4~A"), is_known = .true.)
    call global%set_string (var_str ("$omega_flags"), &
         var_str ("-fusion:progress_file omega_prc_config.log"), &
         is_known = .true.)
    call prc_config%setup_component (2, &
         new_prt_spec (prt_in), new_prt_spec (prt_out), global%model, global%var_list)

    call prc_config%record (global)

    deallocate (prt_in, prt_out)

    global%os_data%fc = "Fortran-compiler"
    global%os_data%fcflags = "Fortran-flags"
    global%os_data%fclibs = "Fortran-libs"

    call global%write_vars (u, [ &
         var_str ("$model_name"), &
         var_str ("$method"), &
         var_str ("?report_progress"), &
         var_str ("$restrictions"), &
         var_str ("$omega_flags")])
    write (u, "(A)")
    call global%write_libraries (u)

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_configurations_2"

  end subroutine process_configurations_2


end module process_configurations_uti

