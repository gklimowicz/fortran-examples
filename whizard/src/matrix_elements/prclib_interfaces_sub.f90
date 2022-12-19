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

submodule (prclib_interfaces) prclib_interfaces_s

  use io_units
  use system_defs, only: TAB
  use string_utils, only: lower_case
  use diagnostics

  implicit none

contains

  module subroutine prc_writer_write_use_line (writer, unit, id, feature)
    class(prc_writer_f_module_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t) :: id, feature
    write (unit, "(2x,9A)")  "use ", char (writer%get_module_name (id)), &
         ", only: ", char (writer%get_c_procname (id, feature)), &
         " => ", char (writer%get_procname (feature))
  end subroutine prc_writer_write_use_line

  module subroutine prc_writer_init_test (writer)
    class(prc_writer_t), intent(out) :: writer
    writer%md5sum = "1234567890abcdef1234567890abcdef"
  end subroutine prc_writer_init_test

  module subroutine prclib_driver_record_write (object, unit)
    class(prclib_driver_record_t), intent(in) :: object
    integer, intent(in) :: unit
    integer :: j
    class(prc_writer_t), pointer :: writer
    write (unit, "(3x,A,2x,'[',A,']')")  &
         char (object%id), char (object%model_name)
    if (allocated (object%feature)) then
       writer => object%writer
       write (unit, "(5x,A,A)", advance="no") &
            char (writer%type_name ()), ":"
       do j = 1, size (object%feature)
          write (unit, "(1x,A)", advance="no") &
               char (object%feature(j))
       end do
       write (unit, *)
    end if
  end subroutine prclib_driver_record_write

  module subroutine prclib_driver_record_write_use_line (record, unit, feature)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    type(string_t), intent(in) :: feature
    select type (writer => record%writer)
    class is (prc_writer_f_module_t)
       call writer%write_use_line (unit, record%id, feature)
    end select
  end subroutine prclib_driver_record_write_use_line

  module subroutine prclib_driver_record_write_interface (record, unit, feature)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    type(string_t), intent(in) :: feature
    select type (writer => record%writer)
    class is (prc_writer_f_module_t)
    class default
       call writer%write_interface (unit, record%id, feature)
    end select
  end subroutine prclib_driver_record_write_interface

  module subroutine prclib_driver_record_write_interfaces (record, unit)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    integer :: i
    do i = 1, size (record%feature)
       call record%writer%write_interface (unit, record%id, record%feature(i))
    end do
  end subroutine prclib_driver_record_write_interfaces

  module subroutine prclib_driver_record_write_wrappers (record, unit)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    integer :: i
    select type (writer => record%writer)
    class is (prc_writer_f_module_t)
       do i = 1, size (record%feature)
          call writer%write_wrapper (unit, record%id, record%feature(i))
       end do
    end select
  end subroutine prclib_driver_record_write_wrappers

  module subroutine prclib_driver_record_write_makefile_code &
       (record, unit, os_data, verbose, testflag)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    call record%writer%write_makefile_code &
         (unit, record%id, os_data, verbose, testflag)
  end subroutine prclib_driver_record_write_makefile_code

  module subroutine prclib_driver_record_write_source_code (record)
    class(prclib_driver_record_t), intent(in) :: record
    call record%writer%write_source_code (record%id)
  end subroutine prclib_driver_record_write_source_code

  module subroutine prclib_driver_record_before_compile (record)
    class(prclib_driver_record_t), intent(in) :: record
    call record%writer%before_compile (record%id)
  end subroutine prclib_driver_record_before_compile

  module subroutine prclib_driver_record_after_compile (record)
    class(prclib_driver_record_t), intent(in) :: record
    call record%writer%after_compile (record%id)
  end subroutine prclib_driver_record_after_compile

  module subroutine prclib_driver_write (object, unit, libpath)
    class(prclib_driver_t), intent(in) :: object
    integer, intent(in) :: unit
    logical, intent(in), optional :: libpath
    logical :: write_lib
    integer :: i
    write_lib = .true.
    if (present (libpath))  write_lib = libpath
    write (unit, "(1x,A,A)")  &
         "External matrix-element code library: ", char (object%basename)
    select type (object)
    type is (prclib_driver_dynamic_t)
       write (unit, "(3x,A,L1)")  "static    = F"
    class default
       write (unit, "(3x,A,L1)")  "static    = T"
    end select
    write (unit, "(3x,A,L1)")  "loaded    = ", object%loaded
    write (unit, "(3x,A,A,A)") "MD5 sum   = '", object%md5sum, "'"
    if (write_lib) then
       write (unit, "(3x,A,A,A)") "Mdl flags = '", &
            char (object%modellibs_ldflags), "'"
    end if
    select type (object)
    type is (prclib_driver_dynamic_t)
       write (unit, *)
       call object%dlaccess%write (unit)
    end select
    write (unit, *)
    if (allocated (object%record)) then
       write (unit, "(1x,A)")  "Matrix-element code entries:"
       do i = 1, object%n_processes
          call object%record(i)%write (unit)
       end do
    else
       write (unit, "(1x,A)")  "Matrix-element code entries: [undefined]"
    end if
  end subroutine prclib_driver_write

  module subroutine prclib_driver_init (driver, n_processes)
    class(prclib_driver_t), intent(inout) :: driver
    integer, intent(in) :: n_processes
    driver%n_processes = n_processes
    allocate (driver%record (n_processes))
  end subroutine prclib_driver_init

  module subroutine prclib_driver_set_md5sum (driver, md5sum)
    class(prclib_driver_t), intent(inout) :: driver
    character(32), intent(in) :: md5sum
    driver%md5sum = md5sum
  end subroutine prclib_driver_set_md5sum

  module subroutine prclib_driver_set_record (driver, i, &
       id, model_name, features, writer)
    class(prclib_driver_t), intent(inout) :: driver
    integer, intent(in) :: i
    type(string_t), intent(in) :: id
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: features
    class(prc_writer_t), intent(in), pointer :: writer
    if (i > 0) then
       associate (record => driver%record(i))
         record%id = id
         record%model_name = model_name
         allocate (record%feature (size (features)))
         record%feature = features
         record%writer => writer
       end associate
    end if
  end subroutine prclib_driver_set_record

  module subroutine prclib_driver_write_interfaces (driver, unit, feature)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: feature
    integer :: i
    do i = 1, driver%n_processes
       call driver%record(i)%write_use_line (unit, feature)
    end do
    write (unit, "(2x,9A)")  "implicit none"
    do i = 1, driver%n_processes
       call driver%record(i)%write_interface (unit, feature)
    end do
  end subroutine prclib_driver_write_interfaces

  module subroutine prclib_driver_generate_makefile (driver, unit, os_data, verbose, testflag)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    integer :: i
    write (unit, "(A)")  "# WHIZARD: Makefile for process library '" &
         // char (driver%basename) // "'"
    write (unit, "(A)")  "# Automatically generated file, do not edit"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Integrity check (don't modify the following line!)"
    write (unit, "(A)")  "MD5SUM = '" // driver%md5sum // "'"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Library name"
    write (unit, "(A)")  "BASE = " // char (driver%basename)
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Compiler"
    write (unit, "(A)")  "FC = " // char (os_data%fc)
    write (unit, "(A)")  "CC = " // char (os_data%cc)
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Included libraries"
    write (unit, "(A)")  "FCINCL = " // char (os_data%whizard_includes)
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Compiler flags"
    write (unit, "(A)")  "FCFLAGS = " // char (os_data%fcflags)
    write (unit, "(A)")  "FCFLAGS_PIC = " // char (os_data%fcflags_pic)
    write (unit, "(A)")  "CFLAGS = " // char (os_data%cflags)
    write (unit, "(A)")  "CFLAGS_PIC = " // char (os_data%cflags_pic)
    write (unit, "(A)")  "LDFLAGS = " // char (os_data%whizard_ldflags) &
         // " " // char (os_data%ldflags) // " " // &
         char (driver%modellibs_ldflags)
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# LaTeX setup"
    write (unit, "(A)")  "LATEX = " // char (os_data%latex)
    write (unit, "(A)")  "MPOST = " // char (os_data%mpost)
    write (unit, "(A)")  "DVIPS = " // char (os_data%dvips)
    write (unit, "(A)")  "PS2PDF = " // char (os_data%ps2pdf)
    write (unit, "(A)")  'TEX_FLAGS = "$$TEXINPUTS:' // &
         char(os_data%whizard_texpath) // '"'
    write (unit, "(A)")  'MP_FLAGS  = "$$MPINPUTS:' // &
         char(os_data%whizard_texpath) // '"'
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Libtool"
    write (unit, "(A)")  "LIBTOOL = " // char (os_data%whizard_libtool)
    if (verbose) then
       write (unit, "(A)")  "FCOMPILE = $(LIBTOOL) --tag=FC --mode=compile"
       write (unit, "(A)")  "CCOMPILE = $(LIBTOOL) --tag=CC --mode=compile"
       write (unit, "(A)")  "LINK = $(LIBTOOL) --tag=FC --mode=link"
    else
       write (unit, "(A)")  "FCOMPILE = @$(LIBTOOL) --silent --tag=FC --mode=compile"
       write (unit, "(A)")  "CCOMPILE = @$(LIBTOOL) --silent --tag=CC --mode=compile"
       write (unit, "(A)")  "LINK = @$(LIBTOOL) --silent --tag=FC --mode=link"
    end if
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Compile commands (default)"
    write (unit, "(A)")  "LTFCOMPILE = $(FCOMPILE) $(FC) -c &
         &$(FCINCL) $(FCFLAGS) $(FCFLAGS_PIC)"
    write (unit, "(A)")  "LTCCOMPILE = $(CCOMPILE) $(CC) -c &
         &$(CFLAGS) $(CFLAGS_PIC)"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Default target"
    write (unit, "(A)")  "all: link diags"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Matrix-element code files"
    do i = 1, size (driver%record)
       call driver%record(i)%write_makefile_code (unit, os_data, verbose, testflag)
    end do
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Library driver"
    write (unit, "(A)")  "$(BASE).lo: $(BASE).f90 $(OBJECTS)"
    write (unit, "(A)")  TAB // "$(LTFCOMPILE) $<"
    if (.not. verbose) then
       write (unit, "(A)")  TAB // '@echo  "  FC       " $@'
    end if
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Library"
    write (unit, "(A)")  "$(BASE).la: $(BASE).lo $(OBJECTS)"
    if (.not. verbose) then
       write (unit, "(A)")  TAB // '@echo  "  FCLD     " $@'
    end if
    write (unit, "(A)")  TAB // "$(LINK) $(FC) -module -rpath /dev/null &
         &$(FCFLAGS) $(LDFLAGS) -o $(BASE).la $^"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Main targets"
    write (unit, "(A)")  "link: compile $(BASE).la"
    write (unit, "(A)")  "compile: source $(OBJECTS) $(TEX_OBJECTS) $(BASE).lo"
    write (unit, "(A)")  "compile_tex: $(TEX_OBJECTS)"
    write (unit, "(A)")  "source: $(SOURCES) $(BASE).f90 $(TEX_SOURCES)"
    write (unit, "(A)")  ".PHONY: link diags compile compile_tex source"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Specific cleanup targets"
    do i = 1, size (driver%record)
       write (unit, "(A)")  "clean-" // char (driver%record(i)%id) // ":"
       write (unit, "(A)")  ".PHONY: clean-" // char (driver%record(i)%id)
    end do
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# Generic cleanup targets"
    write (unit, "(A)")  "clean-library:"
    if (verbose) then
       write (unit, "(A)")  TAB // "rm -f $(BASE).la"
    else
       write (unit, "(A)")  TAB // '@echo  "  RM        $(BASE).la"'
       write (unit, "(A)")  TAB // "@rm -f $(BASE).la"
    end if
    write (unit, "(A)")  "clean-objects:"
    if (verbose) then
       write (unit, "(A)")  TAB // "rm -f $(BASE).lo $(BASE)_driver.mod &
            &$(CLEAN_OBJECTS)"
    else
       write (unit, "(A)")  TAB // '@echo  "  RM        $(BASE).lo &
            &$(BASE)_driver.mod $(CLEAN_OBJECTS)"'
       write (unit, "(A)")  TAB // "@rm -f $(BASE).lo $(BASE)_driver.mod &
            &$(CLEAN_OBJECTS)"
    end if
    write (unit, "(A)")  "clean-source:"
    if (verbose) then
       write (unit, "(A)")  TAB // "rm -f $(CLEAN_SOURCES)"
    else
       write (unit, "(A)")  TAB // '@echo  "  RM        $(CLEAN_SOURCES)"'
       write (unit, "(A)")  TAB // "@rm -f $(CLEAN_SOURCES)"
    end if
    write (unit, "(A)")  "clean-driver:"
    if (verbose) then
       write (unit, "(A)")  TAB // "rm -f $(BASE).f90"
    else
       write (unit, "(A)")  TAB // '@echo  "  RM        $(BASE).f90"'
       write (unit, "(A)")  TAB // "@rm -f $(BASE).f90"
    end if
    write (unit, "(A)")  "clean-makefile:"
    if (verbose) then
       write (unit, "(A)")  TAB // "rm -f $(BASE).makefile"
    else
       write (unit, "(A)")  TAB // '@echo  "  RM        $(BASE).makefile"'
       write (unit, "(A)")  TAB // "@rm -f $(BASE).makefile"
    end if
    write (unit, "(A)")  ".PHONY: clean-library clean-objects &
         &clean-source clean-driver clean-makefile"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "clean: clean-library clean-objects clean-source"
    write (unit, "(A)")  "distclean: clean clean-driver clean-makefile"
    write (unit, "(A)")  ".PHONY: clean distclean"
  end subroutine prclib_driver_generate_makefile

  module subroutine prclib_driver_generate_code (driver, unit)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t) :: prefix
    integer :: i

    prefix = driver%basename // "_"

    write (unit, "(A)")  "! WHIZARD matrix-element code interface"
    write (unit, "(A)")  "!"
    write (unit, "(A)")  "! Automatically generated file, do not edit"
    call driver%write_module (unit, prefix)
    call driver%write_lib_md5sum_fun (unit, prefix)
    call driver%write_get_n_processes_fun (unit, prefix)
    call driver%write_get_process_id_fun (unit, prefix)
    call driver%write_get_model_name_fun (unit, prefix)
    call driver%write_get_md5sum_fun (unit, prefix)
    call driver%write_string_to_array_fun (unit, prefix)
    call driver%write_get_openmp_status_fun (unit, prefix)
    call driver%write_get_int_fun (unit, prefix, var_str ("n_in"))
    call driver%write_get_int_fun (unit, prefix, var_str ("n_out"))
    call driver%write_get_int_fun (unit, prefix, var_str ("n_flv"))
    call driver%write_get_int_fun (unit, prefix, var_str ("n_hel"))
    call driver%write_get_int_fun (unit, prefix, var_str ("n_col"))
    call driver%write_get_int_fun (unit, prefix, var_str ("n_cin"))
    call driver%write_get_int_fun (unit, prefix, var_str ("n_cf"))
    call driver%write_set_int_sub (unit, prefix, var_str ("flv_state"))
    call driver%write_set_int_sub (unit, prefix, var_str ("hel_state"))
    call driver%write_set_col_state_sub (unit, prefix)
    call driver%write_set_color_factors_sub (unit, prefix)
    call driver%write_get_fptr_sub (unit, prefix)
    do i = 1, driver%n_processes
       call driver%record(i)%write_wrappers (unit)
    end do
  end subroutine prclib_driver_generate_code

  module subroutine prclib_driver_write_module (unit, prefix)
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Module: define library driver as an extension &
         &of the abstract driver type."
    write (unit, "(A)")  "! This is used _only_ by the library dispatcher &
         &of a static executable."
    write (unit, "(A)")  "! For a dynamical library, the stand-alone proce&
         &dures are linked via libdl."
    write (unit, "(A)")  ""
    write (unit, "(A)")  "module " &
         // char (prefix) // "driver"
    write (unit, "(A)")  "  use iso_c_binding"
    write (unit, "(A)")  "  use iso_varying_string, string_t => varying_string"
    write (unit, "(A)")  "  use diagnostics"
    write (unit, "(A)")  "  use prclib_interfaces"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "  implicit none"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "  type, extends (prclib_driver_t) :: " &
         // char (prefix) // "driver_t"
    write (unit, "(A)")  "   contains"
    write (unit, "(A)")  "     procedure :: get_c_funptr => " &
         // char (prefix) // "driver_get_c_funptr"
    write (unit, "(A)")  "  end type " &
         // char (prefix) // "driver_t"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "contains"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "  function " &
         // char (prefix) // "driver_get_c_funptr (driver, feature) result &
         &(c_fptr)"
    write (unit, "(A)")  "    class(" &
         // char (prefix) // "driver_t), intent(inout) :: driver"
    write (unit, "(A)")  "    type(string_t), intent(in) :: feature"
    write (unit, "(A)")  "    type(c_funptr) :: c_fptr"
    call write_decl ("get_n_processes", "get_n_processes")
    call write_decl ("get_stringptr", "get_process_id_ptr")
    call write_decl ("get_stringptr", "get_model_name_ptr")
    call write_decl ("get_stringptr", "get_md5sum_ptr")
    call write_decl ("get_log", "get_openmp_status")
    call write_decl ("get_int", "get_n_in")
    call write_decl ("get_int", "get_n_out")
    call write_decl ("get_int", "get_n_flv")
    call write_decl ("get_int", "get_n_hel")
    call write_decl ("get_int", "get_n_col")
    call write_decl ("get_int", "get_n_cin")
    call write_decl ("get_int", "get_n_cf")
    call write_decl ("set_int_tab1", "set_flv_state_ptr")
    call write_decl ("set_int_tab1", "set_hel_state_ptr")
    call write_decl ("set_col_state", "set_col_state_ptr")
    call write_decl ("set_color_factors", "set_color_factors_ptr")
    call write_decl ("get_fptr", "get_fptr")
    write (unit, "(A)")  "    select case (char (feature))"
    call write_case ("get_n_processes")
    call write_case ("get_process_id_ptr")
    call write_case ("get_model_name_ptr")
    call write_case ("get_md5sum_ptr")
    call write_case ("get_openmp_status")
    call write_case ("get_n_in")
    call write_case ("get_n_out")
    call write_case ("get_n_flv")
    call write_case ("get_n_hel")
    call write_case ("get_n_col")
    call write_case ("get_n_cin")
    call write_case ("get_n_cf")
    call write_case ("set_flv_state_ptr")
    call write_case ("set_hel_state_ptr")
    call write_case ("set_col_state_ptr")
    call write_case ("set_color_factors_ptr")
    call write_case ("get_fptr")
    write (unit, "(A)")  "    case default"
    write (unit, "(A)")  "       call msg_bug ('prclib2 driver setup: unknown &
         &function name')"
    write (unit, "(A)")  "    end select"
    write (unit, "(A)")  "  end function " &
         // char (prefix) // "driver_get_c_funptr"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "end module " &
         // char (prefix) // "driver"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Stand-alone external procedures: used for both &
         &static and dynamic linkage"
  contains
    subroutine write_decl (template, feature)
      character(*), intent(in) :: template, feature
      write (unit, "(A)")  "    procedure(prc_" // template // ") &"
      write (unit, "(A)")  "         :: " &
           // char (prefix) // feature
    end subroutine write_decl
    subroutine write_case (feature)
      character(*), intent(in) :: feature
      write (unit, "(A)")  "    case ('" // feature // "')"
      write (unit, "(A)")  "       c_fptr = c_funloc (" &
           // char (prefix) // feature // ")"
    end subroutine write_case
  end subroutine prclib_driver_write_module

  module subroutine prclib_driver_write_lib_md5sum_fun (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! The MD5 sum of the library"
    write (unit, "(A)")  "function " // char (prefix) &
         // "md5sum () result (md5sum)"
    write (unit, "(A)")  "  implicit none"
    write (unit, "(A)")  "  character(32) :: md5sum"
    write (unit, "(A)")  "  md5sum = '" // driver%md5sum // "'"
    write (unit, "(A)")  "end function " // char (prefix) // "md5sum"
  end subroutine prclib_driver_write_lib_md5sum_fun

  module subroutine write_get_n_processes_fun (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Return the number of processes in this library"
    write (unit, "(A)")  "function " // char (prefix) &
         // "get_n_processes () result (n) bind(C)"
    write (unit, "(A)")  "  use iso_c_binding"
    write (unit, "(A)")  "  implicit none"
    write (unit, "(A)")  "  integer(c_int) :: n"
    write (unit, "(A,I0)")  "  n = ", driver%n_processes
    write (unit, "(A)")  "end function " // char (prefix) &
         // "get_n_processes"
  end subroutine write_get_n_processes_fun

  subroutine get_string_via_cptr (string, i, get_stringptr)
    type(string_t), intent(out) :: string
    integer, intent(in) :: i
    procedure(prc_get_stringptr) :: get_stringptr
    type(c_ptr) :: cptr
    integer(c_int) :: pid, len
    character(kind=c_char), dimension(:), pointer :: c_array
    pid = i
    call get_stringptr (pid, cptr, len)
    if (c_associated (cptr)) then
       call c_f_pointer (cptr, c_array, shape = [len])
       call set_string (c_array)
       call get_stringptr (0_c_int, cptr, len)
    else
       string = ""
    end if
  contains
    subroutine set_string (buffer)
      character(len, kind=c_char), dimension(1), intent(in) :: buffer
      string = buffer(1)
    end subroutine set_string
  end subroutine get_string_via_cptr

  module subroutine write_string_to_array_fun (unit, prefix)
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Auxiliary: convert character string &
         &to array pointer"
    write (unit, "(A)")  "subroutine " // char (prefix) &
         // "string_to_array (string, a)"
    write (unit, "(A)")  "  use iso_c_binding"
    write (unit, "(A)")  "  implicit none"
    write (unit, "(A)")  "  character(*), intent(in) :: string"
    write (unit, "(A)")  "  character(kind=c_char), dimension(:), &
         &allocatable, intent(out) :: a"
    write (unit, "(A)")  "  integer :: i"
    write (unit, "(A)")  "  allocate (a (len (string)))"
    write (unit, "(A)")  "  do i = 1, size (a)"
    write (unit, "(A)")  "     a(i) = string(i:i)"
    write (unit, "(A)")  "  end do"
    write (unit, "(A)")  "end subroutine " // char (prefix) &
         // "string_to_array"
  end subroutine write_string_to_array_fun

  subroutine write_string_to_array_interface (unit, prefix)
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    write (unit, "(2x,A)")  "interface"
    write (unit, "(2x,A)")  "   subroutine " // char (prefix) &
         // "string_to_array (string, a)"
    write (unit, "(2x,A)")  "     use iso_c_binding"
    write (unit, "(2x,A)")  "     implicit none"
    write (unit, "(2x,A)")  "     character(*), intent(in) :: string"
    write (unit, "(2x,A)")  "     character(kind=c_char), dimension(:), &
         &allocatable, intent(out) :: a"
    write (unit, "(2x,A)")  "   end subroutine " // char (prefix) &
         // "string_to_array"
    write (unit, "(2x,A)")  "end interface"
  end subroutine write_string_to_array_interface

  module subroutine write_get_process_id_fun (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    integer :: i
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Return the process ID of process #i &
         &(as a C pointer to a character array)"
    write (unit, "(A)")  "subroutine " // char (prefix) &
         // "get_process_id_ptr (i, cptr, len) bind(C)"
    write (unit, "(A)")  "  use iso_c_binding"
    write (unit, "(A)")  "  implicit none"
    write (unit, "(A)")  "  integer(c_int), intent(in) :: i"
    write (unit, "(A)")  "  type(c_ptr), intent(out) :: cptr"
    write (unit, "(A)")  "  integer(c_int), intent(out) :: len"
    write (unit, "(A)")  "  character(kind=c_char), dimension(:), &
         &allocatable, target, save :: a"
    call write_string_to_array_interface (unit, prefix)
    write (unit, "(A)")  "  select case (i)"
    write (unit, "(A)")  "  case (0);  if (allocated (a))  deallocate (a)"
    do i = 1, driver%n_processes
       write (unit, "(A,I0,9A)")  "  case (", i, ");  ", &
            "call ", char (prefix), "string_to_array ('", &
            char (driver%record(i)%id), "', a)"
    end do
    write (unit, "(A)")  "  end select"
    write (unit, "(A)")  "  if (allocated (a)) then"
    write (unit, "(A)")  "     cptr = c_loc (a)"
    write (unit, "(A)")  "     len = size (a)"
    write (unit, "(A)")  "  else"
    write (unit, "(A)")  "     cptr = c_null_ptr"
    write (unit, "(A)")  "     len = 0"
    write (unit, "(A)")  "  end if"
    write (unit, "(A)")  "end subroutine " // char (prefix) &
         // "get_process_id_ptr"
  end subroutine write_get_process_id_fun

  module subroutine write_get_model_name_fun (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    integer :: i
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Return the model name for process #i &
         &(as a C pointer to a character array)"
    write (unit, "(A)")  "subroutine " // char (prefix) &
         // "get_model_name_ptr (i, cptr, len) bind(C)"
    write (unit, "(A)")  "  use iso_c_binding"
    write (unit, "(A)")  "  implicit none"
    write (unit, "(A)")  "  integer(c_int), intent(in) :: i"
    write (unit, "(A)")  "  type(c_ptr), intent(out) :: cptr"
    write (unit, "(A)")  "  integer(c_int), intent(out) :: len"
    write (unit, "(A)")  "  character(kind=c_char), dimension(:), &
         &allocatable, target, save :: a"
    call write_string_to_array_interface (unit, prefix)
    write (unit, "(A)")  "  select case (i)"
    write (unit, "(A)")  "  case (0);  if (allocated (a))  deallocate (a)"
    do i = 1, driver%n_processes
       write (unit, "(A,I0,9A)")  "  case (", i, ");  ", &
            "call ", char (prefix), "string_to_array ('" , &
            char (driver%record(i)%model_name), &
            "', a)"
    end do
    write (unit, "(A)")  "  end select"
    write (unit, "(A)")  "  if (allocated (a)) then"
    write (unit, "(A)")  "     cptr = c_loc (a)"
    write (unit, "(A)")  "     len = size (a)"
    write (unit, "(A)")  "  else"
    write (unit, "(A)")  "     cptr = c_null_ptr"
    write (unit, "(A)")  "     len = 0"
    write (unit, "(A)")  "  end if"
    write (unit, "(A)")  "end subroutine " // char (prefix) &
         // "get_model_name_ptr"
  end subroutine write_get_model_name_fun

  module subroutine write_get_md5sum_fun (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    integer :: i
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Return the MD5 sum for the process configuration &
         &(as a C pointer to a character array)"
    write (unit, "(A)")  "subroutine " // char (prefix) &
         // "get_md5sum_ptr (i, cptr, len) bind(C)"
    write (unit, "(A)")  "  use iso_c_binding"
    call driver%write_interfaces (unit, var_str ("md5sum"))
    write (unit, "(A)")  "  interface"
    write (unit, "(A)")  "     function " // char (prefix) &
         // "md5sum () result (md5sum)"
    write (unit, "(A)")  "       character(32) :: md5sum"
    write (unit, "(A)")  "     end function " // char (prefix) // "md5sum"
    write (unit, "(A)")  "  end interface"
    write (unit, "(A)")  "  integer(c_int), intent(in) :: i"
    write (unit, "(A)")  "  type(c_ptr), intent(out) :: cptr"
    write (unit, "(A)")  "  integer(c_int), intent(out) :: len"
    write (unit, "(A)")  "  character(kind=c_char), dimension(32), &
         &target, save :: md5sum"
    write (unit, "(A)")  "  select case (i)"
    write (unit, "(A)")  "  case (0)"
    write (unit, "(A)")  "     call copy (" // char (prefix) // "md5sum ())"
    write (unit, "(A)")  "     cptr = c_loc (md5sum)"
    do i = 1, driver%n_processes
       write (unit, "(A,I0,A)")  "  case (", i, ")"
       call driver%record(i)%write_md5sum_call (unit)
    end do
    write (unit, "(A)")  "  case default"
    write (unit, "(A)")  "     cptr = c_null_ptr"
    write (unit, "(A)")  "  end select"
    write (unit, "(A)")  "  len = 32"
    write (unit, "(A)")  "contains"
    write (unit, "(A)")  "  subroutine copy (md5sum_tmp)"
    write (unit, "(A)")  "    character, dimension(32), intent(in) :: &
         &md5sum_tmp"
    write (unit, "(A)")  "    md5sum = md5sum_tmp"
    write (unit, "(A)")  "  end subroutine copy"
    write (unit, "(A)")  "end subroutine " // char (prefix) &
           // "get_md5sum_ptr"
  end subroutine write_get_md5sum_fun

  module subroutine prclib_driver_record_write_md5sum_call (record, unit)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    call record%writer%write_md5sum_call (unit, record%id)
  end subroutine prclib_driver_record_write_md5sum_call

  module subroutine prc_writer_f_module_write_md5sum_call (writer, unit, id)
    class(prc_writer_f_module_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(5x,9A)")  "call copy (", &
         char (writer%get_c_procname (id, var_str ("md5sum"))), " ())"
    write (unit, "(5x,9A)")  "cptr = c_loc (md5sum)"
  end subroutine prc_writer_f_module_write_md5sum_call

  module subroutine prc_writer_c_lib_write_md5sum_call (writer, unit, id)
    class(prc_writer_c_lib_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(5x,9A)") &
         "cptr =  ", &
         char (writer%get_c_procname (id, var_str ("get_md5sum"))), " ()"
  end subroutine prc_writer_c_lib_write_md5sum_call

  module function prclib_driver_get_process_id (driver, i) result (string)
    type(string_t) :: string
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    call get_string_via_cptr (string, i, driver%get_process_id_ptr)
  end function prclib_driver_get_process_id

  module function prclib_driver_get_model_name (driver, i) result (string)
    type(string_t) :: string
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    call get_string_via_cptr (string, i, driver%get_model_name_ptr)
  end function prclib_driver_get_model_name

  module function prclib_driver_get_md5sum (driver, i) result (string)
    type(string_t) :: string
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    call get_string_via_cptr (string, i, driver%get_md5sum_ptr)
  end function prclib_driver_get_md5sum

  module subroutine write_get_openmp_status_fun (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    integer :: i
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Return the OpenMP support status"
    write (unit, "(A)")  "function " // char (prefix) &
         // "get_openmp_status (i) result (openmp_status) bind(C)"
    write (unit, "(A)")  "  use iso_c_binding"
    call driver%write_interfaces (unit, var_str ("openmp_supported"))
    write (unit, "(A)")  "  integer(c_int), intent(in) :: i"
    write (unit, "(A)")  "  logical(c_bool) :: openmp_status"
    write (unit, "(A)")  "  select case (i)"
    do i = 1, driver%n_processes
       write (unit, "(A,I0,9A)")  "  case (", i, ");  ", &
            "openmp_status = ", &
            char (driver%record(i)%get_c_procname &
            (var_str ("openmp_supported"))), " ()"
    end do
    write (unit, "(A)")  "  end select"
    write (unit, "(A)")  "end function " // char (prefix) &
         // "get_openmp_status"
  end subroutine write_get_openmp_status_fun

  module subroutine write_get_int_fun (driver, unit, prefix, feature)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    type(string_t), intent(in) :: feature
    integer :: i
    write (unit, "(A)")  ""
    write (unit, "(9A)")  "! Return the value of ", char (feature)
    write (unit, "(9A)")  "function ", char (prefix), &
         "get_", char (feature), " (pid)", &
         " result (", char (feature), ") bind(C)"
    write (unit, "(9A)")  "  use iso_c_binding"
    call driver%write_interfaces (unit, feature)
    write (unit, "(9A)")  "  integer(c_int), intent(in) :: pid"
    write (unit, "(9A)")  "  integer(c_int) :: ", char (feature)
    write (unit, "(9A)")  "  select case (pid)"
    do i = 1, driver%n_processes
       write (unit, "(2x,A,I0,9A)")  "case (", i, ");  ", &
            char (feature), " = ", &
            char (driver%record(i)%get_c_procname (feature)), &
            " ()"
    end do
    write (unit, "(9A)")  "  end select"
    write (unit, "(9A)")  "end function ", char (prefix), &
         "get_", char (feature)
  end subroutine write_get_int_fun

  subroutine write_case_int_fun (record, unit, i, feature)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    integer, intent(in) :: i
    type(string_t), intent(in) :: feature
    write (unit, "(5x,A,I0,9A)")  "case (", i, ");  ", &
         char (feature), " = ", char (record%get_c_procname (feature))
  end subroutine write_case_int_fun

  module subroutine write_set_int_sub (driver, unit, prefix, feature)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    type(string_t), intent(in) :: feature
    integer :: i
    write (unit, "(A)")  ""
    write (unit, "(9A)")  "! Set table: ", char (feature)
    write (unit, "(9A)")  "subroutine ", char (prefix), &
         "set_", char (feature), "_ptr (pid, ", char (feature), &
         ", shape) bind(C)"
    write (unit, "(9A)")  "  use iso_c_binding"
    call driver%write_interfaces (unit, feature)
    write (unit, "(9A)")  "  integer(c_int), intent(in) :: pid"
    write (unit, "(9A)")  "  integer(c_int), dimension(*), intent(out) :: ", &
         char (feature)
    write (unit, "(9A)")  "  integer(c_int), dimension(2), intent(in) :: shape"
    write (unit, "(9A)")  "  integer, dimension(:,:), allocatable :: ", &
         char (feature), "_tmp"
    write (unit, "(9A)")  "  integer :: i, j"
    write (unit, "(9A)")  "  select case (pid)"
    do i = 1, driver%n_processes
       write (unit, "(2x,A,I0,A)")  "case (", i, ")"
       call driver%record(i)%write_int_sub_call (unit, feature)
    end do
    write (unit, "(9A)")  "  end select"
    write (unit, "(9A)")  "end subroutine ", char (prefix), &
         "set_", char (feature), "_ptr"
  end subroutine write_set_int_sub

  module subroutine prclib_driver_record_write_int_sub_call (record, unit, feature)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    type(string_t), intent(in) :: feature
    call record%writer%write_int_sub_call (unit, record%id, feature)
  end subroutine prclib_driver_record_write_int_sub_call

  module subroutine prc_writer_f_module_write_int_sub_call (writer, unit, id, feature)
    class(prc_writer_f_module_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    write (unit, "(5x,9A)")  "allocate (",  char (feature), "_tmp ", &
         "(shape(1), shape(2)))"
    write (unit, "(5x,9A)")  "call ", &
         char (writer%get_c_procname (id, feature)), &
         " (", char (feature), "_tmp)"
    write (unit, "(5x,9A)")  "forall (i=1:shape(1), j=1:shape(2)) "
    write (unit, "(8x,9A)")  char (feature), "(i + shape(1)*(j-1)) = ", &
         char (feature), "_tmp", "(i,j)"
    write (unit, "(5x,9A)")  "end forall"
  end subroutine prc_writer_f_module_write_int_sub_call

  module subroutine prc_writer_c_lib_write_int_sub_call (writer, unit, id, feature)
    class(prc_writer_c_lib_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    write (unit, "(5x,9A)")  "call ", &
         char (writer%get_c_procname (id, feature)), " (", char (feature), ")"
  end subroutine prc_writer_c_lib_write_int_sub_call

  module subroutine write_set_col_state_sub (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    integer :: i
    type(string_t) :: feature
    feature = "col_state"
    write (unit, "(A)")  ""
    write (unit, "(9A)")  "! Set tables: col_state, ghost_flag"
    write (unit, "(9A)")  "subroutine ", char (prefix), &
         "set_col_state_ptr (pid, col_state, ghost_flag, shape) bind(C)"
    write (unit, "(9A)")  "  use iso_c_binding"
    call driver%write_interfaces (unit, feature)
    write (unit, "(9A)")  "  integer(c_int), intent(in) :: pid"
    write (unit, "(9A)") &
         "  integer(c_int), dimension(*), intent(out) :: col_state"
    write (unit, "(9A)") &
         "  logical(c_bool), dimension(*), intent(out) :: ghost_flag"
    write (unit, "(9A)") &
         "  integer(c_int), dimension(3), intent(in) :: shape"
    write (unit, "(9A)") &
         "  integer, dimension(:,:,:), allocatable :: col_state_tmp"
    write (unit, "(9A)") &
         "  logical, dimension(:,:), allocatable :: ghost_flag_tmp"
    write (unit, "(9A)") "  integer :: i, j, k"
    write (unit, "(A)")  "  select case (pid)"
    do i = 1, driver%n_processes
       write (unit, "(A,I0,A)")  "  case (", i, ")"
       call driver%record(i)%write_col_state_call (unit)
    end do
    write (unit, "(A)")  "  end select"
    write (unit, "(9A)")  "end subroutine ", char (prefix), &
         "set_col_state_ptr"
  end subroutine write_set_col_state_sub

  module subroutine prclib_driver_record_write_col_state_call (record, unit)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    call record%writer%write_col_state_call (unit, record%id)
  end subroutine prclib_driver_record_write_col_state_call

  module subroutine prc_writer_f_module_write_col_state_call (writer, unit, id)
    class(prc_writer_f_module_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(9A)")  "  allocate (col_state_tmp ", &
            "(shape(1), shape(2), shape(3)))"
    write (unit, "(5x,9A)")  "allocate (ghost_flag_tmp ", &
               "(shape(2), shape(3)))"
    write (unit, "(5x,9A)")  "call ", &
         char (writer%get_c_procname (id, var_str ("col_state"))), &
         " (col_state_tmp, ghost_flag_tmp)"
    write (unit, "(5x,9A)")  "forall (i = 1:shape(2), j = 1:shape(3))"
    write (unit, "(8x,9A)")  "forall (k = 1:shape(1))"
    write (unit, "(11x,9A)")  &
         "col_state(k + shape(1) * (i + shape(2)*(j-1) - 1)) ", &
         "= col_state_tmp(k,i,j)"
    write (unit, "(8x,9A)")  "end forall"
    write (unit, "(8x,9A)")  &
         "ghost_flag(i + shape(2)*(j-1)) = ghost_flag_tmp(i,j)"
    write (unit, "(5x,9A)")  "end forall"
  end subroutine prc_writer_f_module_write_col_state_call

  module subroutine prc_writer_c_lib_write_col_state_call (writer, unit, id)
    class(prc_writer_c_lib_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(5x,9A)")  "call ", &
         char (writer%get_c_procname (id, var_str ("col_state"))), &
         " (col_state, ghost_flag)"
  end subroutine prc_writer_c_lib_write_col_state_call

  module subroutine write_set_color_factors_sub (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    integer :: i
    type(string_t) :: feature
    feature = "color_factors"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Set tables: color factors"
    write (unit, "(9A)")  "subroutine ", char (prefix), &
         "set_color_factors_ptr (pid, cf_index1, cf_index2, color_factors, ", &
         "shape) bind(C)"
    write (unit, "(A)")  "  use iso_c_binding"
    write (unit, "(A)")  "  use kinds"
    write (unit, "(A)")  "  use omega_color"
    call driver%write_interfaces (unit, feature)
    write (unit, "(A)")  "  integer(c_int), intent(in) :: pid"
    write (unit, "(A)")  "  integer(c_int), dimension(1), intent(in) :: shape"
    write (unit, "(A)")  "  integer(c_int), dimension(*), intent(out) :: &
         &cf_index1, cf_index2"
    write (unit, "(A)")  "  complex(c_default_complex), dimension(*), &
         &intent(out) :: color_factors"
    write (unit, "(A)")  "  type(omega_color_factor), dimension(:), &
         &allocatable :: cf"
    write (unit, "(A)")  "  select case (pid)"
    do i = 1, driver%n_processes
       write (unit, "(2x,A,I0,A)")  "case (", i, ")"
       call driver%record(i)%write_color_factors_call (unit)
    end do
    write (unit, "(A)")  "  end select"
    write (unit, "(A)")  "end subroutine " // char (prefix) &
         // "set_color_factors_ptr"
  end subroutine write_set_color_factors_sub

  module subroutine prclib_driver_record_write_color_factors_call (record, unit)
    class(prclib_driver_record_t), intent(in) :: record
    integer, intent(in) :: unit
    call record%writer%write_color_factors_call (unit, record%id)
  end subroutine prclib_driver_record_write_color_factors_call

  module subroutine prc_writer_f_module_write_color_factors_call (writer, unit, id)
    class(prc_writer_f_module_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(5x,A)")  "allocate (cf (shape(1)))"
    write (unit, "(5x,9A)")  "call ", &
         char (writer%get_c_procname (id, var_str ("color_factors"))), " (cf)"
    write (unit, "(5x,9A)")  "cf_index1(1:shape(1)) = cf%i1"
    write (unit, "(5x,9A)")  "cf_index2(1:shape(1)) = cf%i2"
    write (unit, "(5x,9A)")  "color_factors(1:shape(1)) = cf%factor"
  end subroutine prc_writer_f_module_write_color_factors_call

  module subroutine prc_writer_c_lib_write_color_factors_call (writer, unit, id)
    class(prc_writer_c_lib_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(5x,9A)")  "call ", &
         char (writer%get_c_procname (id, var_str ("color_factors"))), &
         " (cf_index1, cf_index2, color_factors)"
  end subroutine prc_writer_c_lib_write_color_factors_call

  module subroutine prc_writer_c_lib_write_interface (writer, unit, id, feature)
    class(prc_writer_c_lib_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    select case (char (feature))
    case ("md5sum")
       write (unit, "(2x,9A)")  "interface"
       write (unit, "(5x,9A)")  "function ", &
            char (writer%get_c_procname (id, var_str ("get_md5sum"))), &
            " () result (cptr) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "implicit none"
       write (unit, "(7x,9A)")  "type(c_ptr) :: cptr"
       write (unit, "(5x,9A)")  "end function ", &
            char (writer%get_c_procname (id, var_str ("get_md5sum")))
       write (unit, "(2x,9A)")  "end interface"
    case ("openmp_supported")
       write (unit, "(2x,9A)")  "interface"
       write (unit, "(5x,9A)")  "function ", &
            char (writer%get_c_procname (id, feature)), &
            " () result (status) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "implicit none"
       write (unit, "(7x,9A)")  "logical(c_bool) :: status"
       write (unit, "(5x,9A)")  "end function ", &
            char (writer%get_c_procname (id, feature))
       write (unit, "(2x,9A)")  "end interface"
    case ("n_in", "n_out", "n_flv", "n_hel", "n_col", "n_cin", "n_cf")
       write (unit, "(2x,9A)")  "interface"
       write (unit, "(5x,9A)")  "function ", &
            char (writer%get_c_procname (id, feature)), &
            " () result (n) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "implicit none"
       write (unit, "(7x,9A)")  "integer(c_int) :: n"
       write (unit, "(5x,9A)")  "end function ", &
            char (writer%get_c_procname (id, feature))
       write (unit, "(2x,9A)")  "end interface"
    case ("flv_state", "hel_state")
       write (unit, "(2x,9A)")  "interface"
       write (unit, "(5x,9A)")  "subroutine ", &
            char (writer%get_c_procname (id, feature)), &
            " (", char (feature), ") bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "implicit none"
       write (unit, "(7x,9A)")  "integer(c_int), dimension(*), intent(out) ", &
            ":: ", char (feature)
       write (unit, "(5x,9A)")  "end subroutine ", &
            char (writer%get_c_procname (id, feature))
       write (unit, "(2x,9A)")  "end interface"
    case ("col_state")
       write (unit, "(2x,9A)")  "interface"
       write (unit, "(5x,9A)")  "subroutine ", &
            char (writer%get_c_procname (id, feature)), &
            " (col_state, ghost_flag) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "implicit none"
       write (unit, "(7x,9A)")  "integer(c_int), dimension(*), intent(out) ", &
            ":: col_state"
       write (unit, "(7x,9A)")  "logical(c_bool), dimension(*), intent(out) ", &
            ":: ghost_flag"
       write (unit, "(5x,9A)")  "end subroutine ", &
            char (writer%get_c_procname (id, feature))
       write (unit, "(2x,9A)")  "end interface"
    case ("color_factors")
       write (unit, "(2x,9A)")  "interface"
       write (unit, "(5x,9A)")  "subroutine ", &
            char (writer%get_c_procname (id, feature)), &
            " (cf_index1, cf_index2, color_factors) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "implicit none"
       write (unit, "(7x,9A)")  "integer(c_int), dimension(*), &
            &intent(out) :: cf_index1"
       write (unit, "(7x,9A)")  "integer(c_int), dimension(*), &
            &intent(out) :: cf_index2"
       write (unit, "(7x,9A)")  "complex(c_default_complex), dimension(*), &
            &intent(out) :: color_factors"
       write (unit, "(5x,9A)")  "end subroutine ", &
            char (writer%get_c_procname (id, feature))
       write (unit, "(2x,9A)")  "end interface"
    end select
  end subroutine prc_writer_c_lib_write_interface

  module subroutine prclib_driver_set_flv_state (driver, i, flv_state)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    integer, dimension(:,:), allocatable, intent(out) :: flv_state
    integer :: n_tot, n_flv
    integer(c_int) :: pid
    integer(c_int), dimension(:,:), allocatable :: c_flv_state
    pid = i
    n_tot = driver%get_n_in (pid) + driver%get_n_out (pid)
    n_flv = driver%get_n_flv (pid)
    allocate (flv_state (n_tot, n_flv))
    allocate (c_flv_state (n_tot, n_flv))
    call driver%set_flv_state_ptr &
         (pid, c_flv_state, int ([n_tot, n_flv], kind=c_int))
    flv_state = c_flv_state
  end subroutine prclib_driver_set_flv_state

  module subroutine prclib_driver_set_hel_state (driver, i, hel_state)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    integer, dimension(:,:), allocatable, intent(out) :: hel_state
    integer :: n_tot, n_hel
    integer(c_int) :: pid
    integer(c_int), dimension(:,:), allocatable, target :: c_hel_state
    pid = i
    n_tot = driver%get_n_in (pid) + driver%get_n_out (pid)
    n_hel = driver%get_n_hel (pid)
    allocate (hel_state (n_tot, n_hel))
    allocate (c_hel_state (n_tot, n_hel))
    call driver%set_hel_state_ptr &
         (pid, c_hel_state, int ([n_tot, n_hel], kind=c_int))
    hel_state = c_hel_state
  end subroutine prclib_driver_set_hel_state

  module subroutine prclib_driver_set_col_state (driver, i, col_state, ghost_flag)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    integer, dimension(:,:,:), allocatable, intent(out) :: col_state
    logical, dimension(:,:), allocatable, intent(out) :: ghost_flag
    integer :: n_cin, n_tot, n_col
    integer(c_int) :: pid
    integer(c_int), dimension(:,:,:), allocatable :: c_col_state
    logical(c_bool), dimension(:,:), allocatable :: c_ghost_flag
    pid = i
    n_cin = driver%get_n_cin (pid)
    n_tot = driver%get_n_in (pid) + driver%get_n_out (pid)
    n_col = driver%get_n_col (pid)
    allocate (col_state (n_cin, n_tot, n_col))
    allocate (c_col_state (n_cin, n_tot, n_col))
    allocate (ghost_flag (n_tot, n_col))
    allocate (c_ghost_flag (n_tot, n_col))
    call driver%set_col_state_ptr (pid, &
         c_col_state, c_ghost_flag, int ([n_cin, n_tot, n_col], kind=c_int))
    col_state = c_col_state
    ghost_flag = c_ghost_flag
  end subroutine prclib_driver_set_col_state

  module subroutine prclib_driver_set_color_factors (driver, i, color_factors, cf_index)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    complex(default), dimension(:), allocatable, intent(out) :: color_factors
    integer, dimension(:,:), allocatable, intent(out) :: cf_index
    integer :: n_cf
    integer(c_int) :: pid
    complex(c_default_complex), dimension(:), allocatable, target :: c_color_factors
    integer(c_int), dimension(:), allocatable, target :: c_cf_index1
    integer(c_int), dimension(:), allocatable, target :: c_cf_index2
    pid = i
    n_cf = driver%get_n_cf (pid)
    allocate (color_factors (n_cf))
    allocate (c_color_factors (n_cf))
    allocate (c_cf_index1 (n_cf))
    allocate (c_cf_index2 (n_cf))
    call driver%set_color_factors_ptr (pid, &
         c_cf_index1, c_cf_index2, &
         c_color_factors, int ([n_cf], kind=c_int))
    color_factors = c_color_factors
    allocate (cf_index (2, n_cf))
    cf_index(1,:) = c_cf_index1
    cf_index(2,:) = c_cf_index2
  end subroutine prclib_driver_set_color_factors

  module subroutine write_get_fptr_sub (driver, unit, prefix)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: unit
    type(string_t), intent(in) :: prefix
    integer :: i, j
    write (unit, "(A)")  ""
    write (unit, "(A)")  "! Return C pointer to a procedure:"
    write (unit, "(A)")  "! pid = process index;  fid = function index"
    write (unit, "(4A)")  "subroutine ", char (prefix), "get_fptr ", &
         "(pid, fid, fptr) bind(C)"
    write (unit, "(A)")  "  use iso_c_binding"
    write (unit, "(A)")  "  use kinds"
    write (unit, "(A)")  "  implicit none"
    write (unit, "(A)")  "  integer(c_int), intent(in) :: pid"
    write (unit, "(A)")  "  integer(c_int), intent(in) :: fid"
    write (unit, "(A)")  "  type(c_funptr), intent(out) :: fptr"
    do i = 1, driver%n_processes
       call driver%record(i)%write_interfaces (unit)
    end do
    write (unit, "(A)")  "  select case (pid)"
    do i = 1, driver%n_processes
       write (unit, "(2x,A,I0,A)")  "case (", i, ")"
       write (unit, "(5x,A)")  "select case (fid)"
       associate (record => driver%record(i))
         do j = 1, size (record%feature)
            write (unit, "(5x,A,I0,9A)")  "case (", j, ");  ", &
                 "fptr = c_funloc (", &
                 char (record%get_c_procname (record%feature(j))), &
                 ")"
         end do
       end associate
       write (unit, "(5x,A)")  "end select"
    end do
    write (unit, "(A)")  "  end select"
    write (unit, "(3A)")  "end subroutine ", char (prefix), "get_fptr"
  end subroutine write_get_fptr_sub

  module function workspace_prefix (workspace) result (prefix)
    type(string_t), intent(in), optional :: workspace
    type(string_t) :: prefix
    if (present (workspace)) then
       if (workspace /= "") then
          prefix = workspace // "/"
       else
          prefix = ""
       end if
    else
       prefix = ""
    end if
  end function workspace_prefix

  module function workspace_path (workspace) result (path)
    type(string_t), intent(in), optional :: workspace
    type(string_t) :: path
    if (present (workspace)) then
       if (workspace /= "") then
          path = workspace
       else
          path = "."
       end if
    else
       path = "."
    end if
  end function workspace_path

  module function workspace_cmd (workspace) result (cmd)
    type(string_t), intent(in), optional :: workspace
    type(string_t) :: cmd
    if (present (workspace)) then
       if (workspace /= "") then
          cmd = "cd " // workspace // " && "
       else
          cmd = ""
       end if
    else
       cmd = ""
    end if
  end function workspace_cmd

  module subroutine prclib_driver_make_source (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    integer :: i
    do i = 1, driver%n_processes
       call driver%record(i)%write_source_code ()
    end do
    call os_system_call ( &
         workspace_cmd (workspace) &
         // "make source " // os_data%makeflags &
         // " -f " // driver%basename // ".makefile")
  end subroutine prclib_driver_make_source

  module subroutine prclib_driver_make_compile (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    integer :: i
    do i = 1, driver%n_processes
       call driver%record(i)%before_compile ()
    end do
    call os_system_call ( &
         workspace_cmd (workspace) &
         // "make compile " // os_data%makeflags &
         // " -f " // driver%basename // ".makefile")
    do i = 1, driver%n_processes
       call driver%record(i)%after_compile ()
    end do
  end subroutine prclib_driver_make_compile

  module subroutine prclib_driver_make_link (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    integer :: i
    call os_system_call ( &
         workspace_cmd (workspace) &
         // "make link " // os_data%makeflags &
         // " -f " // driver%basename // ".makefile")
  end subroutine prclib_driver_make_link

  module subroutine prclib_driver_clean_library (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    if (driver%makefile_exists ()) then
       call os_system_call ( &
         workspace_cmd (workspace) &
         // "make clean-library " // os_data%makeflags &
            // " -f " // driver%basename // ".makefile")
    end if
  end subroutine prclib_driver_clean_library

  module subroutine prclib_driver_clean_objects (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    if (driver%makefile_exists ()) then
       call os_system_call ( &
         workspace_cmd (workspace) &
         // "make clean-objects " // os_data%makeflags &
            // " -f " // driver%basename // ".makefile")
    end if
  end subroutine prclib_driver_clean_objects

  module subroutine prclib_driver_clean_source (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    if (driver%makefile_exists ()) then
       call os_system_call ( &
         workspace_cmd (workspace) &
         // "make clean-source " // os_data%makeflags &
            // " -f " // driver%basename // ".makefile")
    end if
  end subroutine prclib_driver_clean_source

  module subroutine prclib_driver_clean_driver (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    if (driver%makefile_exists ()) then
       call os_system_call ( &
         workspace_cmd (workspace) &
         // "make clean-driver " // os_data%makeflags &
            // " -f " // driver%basename // ".makefile")
    end if
  end subroutine prclib_driver_clean_driver

  module subroutine prclib_driver_clean_makefile (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    if (driver%makefile_exists ()) then
       call os_system_call ( &
         workspace_cmd (workspace) &
         // "make clean-makefile " // os_data%makeflags &
            // " -f " // driver%basename // ".makefile")
    end if
  end subroutine prclib_driver_clean_makefile

  module subroutine prclib_driver_clean (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    if (driver%makefile_exists ()) then
       call os_system_call ( &
         workspace_cmd (workspace) &
         // "make clean " // os_data%makeflags &
            // " -f " // driver%basename // ".makefile")
    end if
  end subroutine prclib_driver_clean

  module subroutine prclib_driver_distclean (driver, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    if (driver%makefile_exists ()) then
       call os_system_call ( &
         workspace_cmd (workspace) &
         // "make distclean " // os_data%makeflags &
            // " -f " // driver%basename // ".makefile")
    end if
  end subroutine prclib_driver_distclean

  module subroutine prclib_driver_clean_proc (driver, i, os_data, workspace)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    type(string_t) :: id
    if (driver%makefile_exists ()) then
       id = driver%record(i)%id
      call os_system_call ( &
         workspace_cmd (workspace) &
         // "make clean-" // driver%record(i)%id // " " &
           // os_data%makeflags &
           // " -f " // driver%basename // ".makefile")
    end if
  end subroutine prclib_driver_clean_proc

  module function prclib_driver_makefile_exists (driver, workspace) result (flag)
    class(prclib_driver_t), intent(in) :: driver
    type(string_t), intent(in), optional :: workspace
    logical :: flag
    inquire (file = char (workspace_prefix (workspace) &
         &                // driver%basename) // ".makefile", &
         exist = flag)
  end function prclib_driver_makefile_exists

  module subroutine prclib_driver_load (driver, os_data, noerror, workspace)
    class(prclib_driver_t), intent(inout) :: driver
    type(os_data_t), intent(in) :: os_data
    logical, intent(in), optional :: noerror
    type(string_t), intent(in), optional :: workspace
    type(c_funptr) :: c_fptr
    logical :: ignore

    ignore = .false.;  if (present (noerror))  ignore = noerror

    driver%libname = os_get_dlname ( &
         workspace_prefix (workspace) // driver%basename, &
         os_data, noerror, noerror)
    if (driver%libname == "")  return
    select type (driver)
    type is (prclib_driver_dynamic_t)
       if (.not. dlaccess_is_open (driver%dlaccess)) then
          call dlaccess_init &
               (driver%dlaccess, workspace_path (workspace), &
                driver%libname, os_data)
          if (.not. ignore)  call driver%check_dlerror ()
       end if
       driver%loaded = dlaccess_is_open (driver%dlaccess)
    class default
       driver%loaded = .true.
    end select
    if (.not. driver%loaded)  return

    c_fptr = driver%get_c_funptr (var_str ("get_n_processes"))
    call c_f_procpointer (c_fptr, driver%get_n_processes)
    driver%loaded = driver%loaded .and. associated (driver%get_n_processes)

    c_fptr = driver%get_c_funptr (var_str ("get_process_id_ptr"))
    call c_f_procpointer (c_fptr, driver%get_process_id_ptr)
    driver%loaded = driver%loaded .and. associated (driver%get_process_id_ptr)

    c_fptr = driver%get_c_funptr (var_str ("get_model_name_ptr"))
    call c_f_procpointer (c_fptr, driver%get_model_name_ptr)
    driver%loaded = driver%loaded .and. associated (driver%get_model_name_ptr)

    c_fptr = driver%get_c_funptr (var_str ("get_md5sum_ptr"))
    call c_f_procpointer (c_fptr, driver%get_md5sum_ptr)
    driver%loaded = driver%loaded .and. associated (driver%get_md5sum_ptr)

    c_fptr = driver%get_c_funptr (var_str ("get_openmp_status"))
    call c_f_procpointer (c_fptr, driver%get_openmp_status)
    driver%loaded = driver%loaded .and. associated (driver%get_openmp_status)

    c_fptr = driver%get_c_funptr (var_str ("get_n_in"))
    call c_f_procpointer (c_fptr, driver%get_n_in)
    driver%loaded = driver%loaded .and. associated (driver%get_n_in)

    c_fptr = driver%get_c_funptr (var_str ("get_n_out"))
    call c_f_procpointer (c_fptr, driver%get_n_out)
    driver%loaded = driver%loaded .and. associated (driver%get_n_out)

    c_fptr = driver%get_c_funptr (var_str ("get_n_flv"))
    call c_f_procpointer (c_fptr, driver%get_n_flv)
    driver%loaded = driver%loaded .and. associated (driver%get_n_flv)

    c_fptr = driver%get_c_funptr (var_str ("get_n_hel"))
    call c_f_procpointer (c_fptr, driver%get_n_hel)
    driver%loaded = driver%loaded .and. associated (driver%get_n_hel)

    c_fptr = driver%get_c_funptr (var_str ("get_n_col"))
    call c_f_procpointer (c_fptr, driver%get_n_col)
    driver%loaded = driver%loaded .and. associated (driver%get_n_col)

    c_fptr = driver%get_c_funptr (var_str ("get_n_cin"))
    call c_f_procpointer (c_fptr, driver%get_n_cin)
    driver%loaded = driver%loaded .and. associated (driver%get_n_cin)

    c_fptr = driver%get_c_funptr (var_str ("get_n_cf"))
    call c_f_procpointer (c_fptr, driver%get_n_cf)
    driver%loaded = driver%loaded .and. associated (driver%get_n_cf)

    c_fptr = driver%get_c_funptr (var_str ("set_flv_state_ptr"))
    call c_f_procpointer (c_fptr, driver%set_flv_state_ptr)
    driver%loaded = driver%loaded .and. associated (driver%set_flv_state_ptr)

    c_fptr = driver%get_c_funptr (var_str ("set_hel_state_ptr"))
    call c_f_procpointer (c_fptr, driver%set_hel_state_ptr)
    driver%loaded = driver%loaded .and. associated (driver%set_hel_state_ptr)

    c_fptr = driver%get_c_funptr (var_str ("set_col_state_ptr"))
    call c_f_procpointer (c_fptr, driver%set_col_state_ptr)
    driver%loaded = driver%loaded .and. associated (driver%set_col_state_ptr)

    c_fptr = driver%get_c_funptr (var_str ("set_color_factors_ptr"))
    call c_f_procpointer (c_fptr, driver%set_color_factors_ptr)
    driver%loaded = driver%loaded .and. associated (driver%set_color_factors_ptr)

    c_fptr = driver%get_c_funptr (var_str ("get_fptr"))
    call c_f_procpointer (c_fptr, driver%get_fptr)
    driver%loaded = driver%loaded .and. associated (driver%get_fptr)

  end subroutine prclib_driver_load

  module subroutine prclib_driver_unload (driver)
    class(prclib_driver_t), intent(inout) :: driver
    select type (driver)
    type is (prclib_driver_dynamic_t)
       if (dlaccess_is_open (driver%dlaccess)) then
          call dlaccess_final (driver%dlaccess)
          call driver%check_dlerror ()
       end if
    end select
    driver%loaded = .false.
    nullify (driver%get_n_processes)
    nullify (driver%get_process_id_ptr)
    nullify (driver%get_model_name_ptr)
    nullify (driver%get_md5sum_ptr)
    nullify (driver%get_openmp_status)
    nullify (driver%get_n_in)
    nullify (driver%get_n_out)
    nullify (driver%get_n_flv)
    nullify (driver%get_n_hel)
    nullify (driver%get_n_col)
    nullify (driver%get_n_cin)
    nullify (driver%get_n_cf)
    nullify (driver%set_flv_state_ptr)
    nullify (driver%set_hel_state_ptr)
    nullify (driver%set_col_state_ptr)
    nullify (driver%set_color_factors_ptr)
    nullify (driver%get_fptr)
  end subroutine prclib_driver_unload

  module subroutine prclib_driver_check_dlerror (driver)
    class(prclib_driver_dynamic_t), intent(in) :: driver
    if (dlaccess_has_error (driver%dlaccess)) then
       call msg_fatal (char (dlaccess_get_error (driver%dlaccess)))
    end if
  end subroutine prclib_driver_check_dlerror

  module function prclib_driver_dynamic_get_c_funptr &
       (driver, feature) result (c_fptr)
    class(prclib_driver_dynamic_t), intent(inout) :: driver
    type(string_t), intent(in) :: feature
    type(c_funptr) :: c_fptr
    type(string_t) :: prefix, full_name
    prefix = lower_case (driver%basename) // "_"
    full_name = prefix // feature
    c_fptr = dlaccess_get_c_funptr (driver%dlaccess, full_name)
    call driver%check_dlerror ()
  end function prclib_driver_dynamic_get_c_funptr

  module function prclib_driver_get_md5sum_makefile &
       (driver, workspace) result (md5sum)
    class(prclib_driver_t), intent(in) :: driver
    type(string_t), intent(in), optional :: workspace
    character(32) :: md5sum
    type(string_t) :: filename
    character(80) :: buffer
    logical :: exist
    integer :: u, iostat
    md5sum = ""
    filename = workspace_prefix (workspace) // driver%basename // ".makefile"
    inquire (file = char (filename), exist = exist)
    if (exist) then
       u = free_unit ()
       open (u, file = char (filename), action = "read", status = "old")
       iostat = 0
       do
          read (u, "(A)", iostat = iostat)  buffer
          if (iostat /= 0)  exit
          buffer = adjustl (buffer)
          select case (buffer(1:9))
          case ("MD5SUM = ")
             read (buffer(11:), "(A32)")  md5sum
             exit
          end select
       end do
       close (u)
    end if
  end function prclib_driver_get_md5sum_makefile

  module function prclib_driver_get_md5sum_driver (driver, workspace) result (md5sum)
    class(prclib_driver_t), intent(in) :: driver
    type(string_t), intent(in), optional :: workspace
    character(32) :: md5sum
    type(string_t) :: filename
    character(80) :: buffer
    logical :: exist
    integer :: u, iostat
    md5sum = ""
    filename = workspace_prefix (workspace) // driver%basename // ".f90"
    inquire (file = char (filename), exist = exist)
    if (exist) then
       u = free_unit ()
       open (u, file = char (filename), action = "read", status = "old")
       iostat = 0
       do
          read (u, "(A)", iostat = iostat)  buffer
          if (iostat /= 0)  exit
          buffer = adjustl (buffer)
          select case (buffer(1:9))
          case ("md5sum = ")
             read (buffer(11:), "(A32)")  md5sum
             exit
          end select
       end do
       close (u)
    end if
  end function prclib_driver_get_md5sum_driver

  module function prclib_driver_get_md5sum_source &
       (driver, i, workspace) result (md5sum)
    class(prclib_driver_t), intent(in) :: driver
    integer, intent(in) :: i
    type(string_t), intent(in), optional :: workspace
    character(32) :: md5sum
    type(string_t) :: filename
    character(80) :: buffer
    logical :: exist
    integer :: u, iostat
    md5sum = ""

    filename = workspace_prefix (workspace) // driver%record(i)%id // ".f90"
    inquire (file = char (filename), exist = exist)
    if (exist) then
       u = free_unit ()
       open (u, file = char (filename), action = "read", status = "old")
       iostat = 0
       do
          read (u, "(A)", iostat = iostat)  buffer
          if (iostat /= 0)  exit
          buffer = adjustl (buffer)
          select case (buffer(1:9))
          case ("md5sum = ")
             read (buffer(11:), "(A32)")  md5sum
             exit
          end select
       end do
       close (u)
    end if
  end function prclib_driver_get_md5sum_source


end submodule prclib_interfaces_s

