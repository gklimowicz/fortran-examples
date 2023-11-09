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

module compilations_uti

  use iso_varying_string, string_t => varying_string
  use io_units
  use models
  use rt_data
  use process_configurations_ut, only: prepare_test_library

  use compilations

  implicit none
  private

  public :: compilations_1
  public :: compilations_2
  public :: compilations_3
  public :: compilations_static_1
  public :: compilations_static_2

contains

  subroutine compilations_1 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global

    write (u, "(A)")  "* Test output: compilations_1"
    write (u, "(A)")  "*   Purpose: configure and compile test process"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()

    libname = "compilation_1"
    procname = "prc_comp_1"
    call prepare_test_library (global, libname, 1, [procname])

    call compile_library (libname, global)

    call global%write_libraries (u)

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: compilations_1"

  end subroutine compilations_1

  subroutine compilations_2 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname
    type(rt_data_t), target :: global

    write (u, "(A)")  "* Test output: compilations_2"
    write (u, "(A)")  "*   Purpose: configure and compile test process"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    libname = "compilation_2"
    procname = "prc_comp_2"
    call prepare_test_library (global, libname, 2, [procname,procname])

    call compile_library (libname, global)

    call global%write_libraries (u, libpath = .false.)

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: compilations_2"

  end subroutine compilations_2

  subroutine compilations_3 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname, exename
    type(rt_data_t), target :: global
    type(compilation_t) :: compilation
    integer :: u_file
    character(80) :: buffer

    write (u, "(A)")  "* Test output: compilations_3"
    write (u, "(A)")  "*   Purpose: make static executable"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize library"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    libname = "compilations_3_lib"
    procname = "prc_comp_3"
    exename = "compilations_3"

    call prepare_test_library (global, libname, 2, [procname,procname])

    call compilation%init (exename, [libname])
    call compilation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write dispatcher"
    write (u, "(A)")

    call compilation%write_dispatcher ()

    u_file = free_unit ()
    open (u_file, file = char (exename) // "_prclib_dispatcher.f90", &
         status = "old", action = "read")
    do
       read (u_file, "(A)", end = 1)  buffer
       write (u, "(A)")  trim (buffer)
    end do
1   close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Write Makefile"
    write (u, "(A)")

    associate (os_data => global%os_data)
      os_data%fc = "fortran-compiler"
      os_data%cxx = "c++-compiler"
      os_data%whizard_includes = "my-includes"
      os_data%fcflags = "my-fcflags"
      os_data%fclibs = "my-fclibs"
      os_data%cxxflags = "my-cxxflags"
      os_data%cxxlibs = "my-cxxlibs"
      os_data%ldflags = "my-ldflags"
      os_data%ldflags_static = "my-ldflags-static"
      os_data%ldflags_hepmc = "my-ldflags-hepmc"
      os_data%ldflags_lcio = "my-ldflags-lcio"
      os_data%ldflags_hoppet = "my-ldflags-hoppet"
      os_data%ldflags_looptools = "my-ldflags-looptools"
      os_data%whizard_ldflags = "my-ldwhizard"
      os_data%whizard_libtool = "my-libtool"
    end associate

    call compilation%write_makefile &
         (global%os_data, verbose = .true., overwrite_os = .true.)

    open (u_file, file = char (exename) // ".makefile", &
         status = "old", action = "read")
    do
       read (u_file, "(A)", end = 2)  buffer
       write (u, "(A)")  trim (buffer)
    end do
2   close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: compilations_3"

  end subroutine compilations_3

  subroutine compilations_static_1 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname, exename
    type(rt_data_t), target :: global
    type(compilation_item_t) :: item
    type(compilation_t) :: compilation
    logical :: exist

    write (u, "(A)")  "* Test output: compilations_static_1"
    write (u, "(A)")  "*   Purpose: make static executable"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize library"

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    libname = "compilations_static_1_lib"
    procname = "prc_comp_stat_1"
    exename = "compilations_static_1"

    call prepare_test_library (global, libname, 2, [procname,procname])

    call compilation%init (exename, [libname])

    write (u, "(A)")
    write (u, "(A)")  "* Write dispatcher"

    call compilation%write_dispatcher ()

    write (u, "(A)")
    write (u, "(A)")  "* Write Makefile"

    call compilation%write_makefile (global%os_data, verbose = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Build libraries"

    call item%init (libname, global%prclib_stack, global%var_list)
    call item%compile &
         (global%model, global%os_data, force=.true., recompile=.false.)
    call item%success ()

    write (u, "(A)")
    write (u, "(A)")  "* Check executable (should be absent)"
    write (u, "(A)")

    call compilation%make_clean_exe (global%os_data)
    inquire (file = char (exename), exist = exist)
    write (u, "(A,A,L1)")  char (exename), " exists = ", exist

    write (u, "(A)")
    write (u, "(A)")  "* Build executable"
    write (u, "(A)")

    call compilation%make_compile (global%os_data)
    call compilation%make_link (global%os_data)

    write (u, "(A)")  "* Check executable (should be present)"
    write (u, "(A)")

    inquire (file = char (exename), exist = exist)
    write (u, "(A,A,L1)")  char (exename), " exists = ", exist

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call compilation%make_clean_exe (global%os_data)

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: compilations_static_1"

  end subroutine compilations_static_1

  subroutine compilations_static_2 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname, exename
    type(rt_data_t), target :: global
    logical :: exist
    integer :: u_file

    write (u, "(A)")  "* Test output: compilations_static_2"
    write (u, "(A)")  "*   Purpose: make static executable"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize library and compile"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    libname = "compilations_static_2_lib"
    procname = "prc_comp_stat_2"
    exename = "compilations_static_2"

    call prepare_test_library (global, libname, 2, [procname,procname])

    call compile_executable (exename, [libname], global)

    write (u, "(A)")  "* Check executable (should be present)"
    write (u, "(A)")

    inquire (file = char (exename), exist = exist)
    write (u, "(A,A,L1)")  char (exename), " exists = ", exist

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    u_file = free_unit ()
    open (u_file, file = char (exename), status = "old", action = "write")
    close (u_file, status = "delete")

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: compilations_static_2"

  end subroutine compilations_static_2


end module compilations_uti

