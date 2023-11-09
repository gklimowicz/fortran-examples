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

submodule (compilations) compilations_s

  use io_units
  use system_defs, only: TAB
  use system_dependencies, only: OS_IS_DARWIN
  use diagnostics

  implicit none

contains

  module subroutine compilation_item_init (comp, libname, stack, var_list)
    class(compilation_item_t), intent(out) :: comp
    type(string_t), intent(in) :: libname
    type(prclib_stack_t), intent(inout) :: stack
    type(var_list_t), intent(in) :: var_list
    comp%libname = libname
    comp%lib => stack%get_library_ptr (comp%libname)
    if (.not. associated (comp%lib)) then
       call msg_fatal ("Process library '" // char (comp%libname) &
            // "' has not been declared.")
    end if
    comp%recompile_library = &
         var_list%get_lval (var_str ("?recompile_library"))
    comp%verbose = &
         var_list%get_lval (var_str ("?me_verbose"))
    comp%use_workspace = &
         var_list%is_known (var_str ("$compile_workspace"))
    if (comp%use_workspace) then
       comp%workspace = &
            var_list%get_sval (var_str ("$compile_workspace"))
       if (comp%workspace == "")  comp%use_workspace = .false.
    else
       comp%workspace = ""
    end if
  end subroutine compilation_item_init

  module subroutine compilation_item_compile &
       (comp, model, os_data, force, recompile)
    class(compilation_item_t), intent(inout) :: comp
    class(model_data_t), intent(in), target :: model
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: force, recompile
    if (associated (comp%lib)) then
       if (comp%use_workspace)  call setup_workspace (comp%workspace, os_data)
       call msg_message ("Process library '" &
            // char (comp%libname) // "': compiling ...")
       call comp%lib%configure (os_data)
       if (signal_is_pending ())  return
       call comp%lib%compute_md5sum (model)
       call comp%lib%write_makefile &
            (os_data, force, verbose=comp%verbose, workspace=comp%workspace)
       if (signal_is_pending ())  return
       if (force) then
          call comp%lib%clean &
               (os_data, distclean = .false., workspace=comp%workspace)
          if (signal_is_pending ())  return
       end if
       call comp%lib%write_driver (force, workspace=comp%workspace)
       if (signal_is_pending ())  return
       if (recompile) then
          call comp%lib%load &
               (os_data, keep_old_source = .true., workspace=comp%workspace)
          if (signal_is_pending ())  return
       end if
       call comp%lib%update_status (os_data, workspace=comp%workspace)
    end if
  end subroutine compilation_item_compile

  subroutine setup_workspace (workspace, os_data)
    type(string_t), intent(in) :: workspace
    type(os_data_t), intent(in) :: os_data
    if (verify (workspace, ALLOWED_IN_DIRNAME) == 0) then
       call msg_message ("Compile: preparing workspace directory '" &
            // char (workspace) // "'")
       call os_system_call ("mkdir -p '" // workspace // "'")
    else
       call msg_fatal ("compile: workspace name '" &
            // char (workspace) // "' contains illegal characters")
    end if
  end subroutine setup_workspace

  module subroutine compilation_item_load (comp, os_data)
    class(compilation_item_t), intent(inout) :: comp
    type(os_data_t), intent(in) :: os_data
    if (associated (comp%lib)) then
       call comp%lib%load (os_data, workspace=comp%workspace)
    end if
  end subroutine compilation_item_load

  module subroutine compilation_item_success (comp)
    class(compilation_item_t), intent(in) :: comp
    if (associated (comp%lib)) then
       call msg_message ("Process library '" // char (comp%libname) &
            // "': ... success.")
    else
       call msg_fatal ("Process library '" // char (comp%libname) &
            // "': ... failure.")
    end if
  end subroutine compilation_item_success

  module subroutine compile_library (libname, global)
    type(string_t), intent(in) :: libname
    type(rt_data_t), intent(inout), target :: global
    type(compilation_item_t) :: comp
    logical :: force, recompile
    force = &
         global%var_list%get_lval (var_str ("?rebuild_library"))
    recompile = &
         global%var_list%get_lval (var_str ("?recompile_library"))
    if (associated (global%model)) then
       call comp%init (libname, global%prclib_stack, global%var_list)
       call comp%compile (global%model, global%os_data, force, recompile)
       if (signal_is_pending ())  return
       call comp%load (global%os_data)
       if (signal_is_pending ())  return
    else
       call msg_fatal ("Process library compilation: " &
            // " model is undefined.")
    end if
    call comp%success ()
  end subroutine compile_library

  module subroutine compilation_write (object, unit)
    class(compilation_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Compilation object:"
    write (u, "(3x,3A)")  "executable        = '", &
         char (object%exe_name), "'"
    write (u, "(3x,A)", advance="no")  "process libraries ="
    do i = 1, size (object%lib_name)
       write (u, "(1x,3A)", advance="no")  "'", char (object%lib_name(i)), "'"
    end do
    write (u, *)
  end subroutine compilation_write

  module subroutine compilation_init (compilation, exe_name, lib_name)
    class(compilation_t), intent(out) :: compilation
    type(string_t), intent(in) :: exe_name
    type(string_t), dimension(:), intent(in) :: lib_name
    compilation%exe_name = exe_name
    allocate (compilation%lib_name (size (lib_name)))
    compilation%lib_name = lib_name
  end subroutine compilation_init

  module subroutine compilation_write_dispatcher (compilation)
    class(compilation_t), intent(in) :: compilation
    type(string_t) :: file
    integer :: u, i
    file = compilation%exe_name // "_prclib_dispatcher.f90"
    call msg_message ("Static executable '" // char (compilation%exe_name) &
         // "': writing library dispatcher")
    u = free_unit ()
    open (u, file = char (file), status="replace", action="write")
    write (u, "(3A)")  "! Whizard: process libraries for executable '", &
         char (compilation%exe_name), "'"
    write (u, "(A)")  "! Automatically generated file, do not edit"
    write (u, "(A)")  "subroutine dispatch_prclib_static " // &
         "(driver, basename, modellibs_ldflags)"
    write (u, "(A)")  "  use iso_varying_string, string_t => varying_string"
    write (u, "(A)")  "  use prclib_interfaces"
    do i = 1, size (compilation%lib_name)
       associate (lib_name => compilation%lib_name(i))
         write (u, "(A)")  "  use " // char (lib_name) // "_driver"
       end associate
    end do
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  class(prclib_driver_t), intent(inout), allocatable &
         &:: driver"
    write (u, "(A)")  "  type(string_t), intent(in) :: basename"
    write (u, "(A)")  "  logical, intent(in), optional :: " // &
         "modellibs_ldflags"
    write (u, "(A)")  "  select case (char (basename))"
    do i = 1, size (compilation%lib_name)
       associate (lib_name => compilation%lib_name(i))
         write (u, "(3A)")  "  case ('", char (lib_name), "')"
         write (u, "(3A)")  "     allocate (", char (lib_name), "_driver_t &
              &:: driver)"
       end associate
    end do
    write (u, "(A)")  "  end select"
    write (u, "(A)")  "end subroutine dispatch_prclib_static"
    write (u, *)
    write (u, "(A)")  "subroutine get_prclib_static (libname)"
    write (u, "(A)")  "  use iso_varying_string, string_t => varying_string"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  type(string_t), dimension(:), intent(inout), &
         &allocatable :: libname"
    write (u, "(A,I0,A)")  "  allocate (libname (", &
         size (compilation%lib_name), "))"
    do i = 1, size (compilation%lib_name)
       associate (lib_name => compilation%lib_name(i))
         write (u, "(A,I0,A,A,A)")  "  libname(", i, ") = '", &
              char (lib_name), "'"
       end associate
    end do
    write (u, "(A)")  "end subroutine get_prclib_static"
    close (u)
  end subroutine compilation_write_dispatcher

  module subroutine compilation_write_makefile &
       (compilation, os_data, ext_libtag, verbose, overwrite_os)
    class(compilation_t), intent(in) :: compilation
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: overwrite_os
    type(string_t), intent(in), optional :: ext_libtag
    logical :: overwrite
    type(string_t) :: file, ext_tag
    integer :: u, i
    overwrite = .false.
    if (present (overwrite_os))  overwrite = overwrite_os
    if (present (ext_libtag)) then
       ext_tag = ext_libtag
    else
       ext_tag = ""
    end if
    file = compilation%exe_name // ".makefile"
    call msg_message ("Static executable '" // char (compilation%exe_name) &
         // "': writing makefile")
    u = free_unit ()
    open (u, file = char (file), status="replace", action="write")
    write (u, "(3A)")  "# WHIZARD: Makefile for executable '", &
         char (compilation%exe_name), "'"
    write (u, "(A)")  "# Automatically generated file, do not edit"
    write (u, "(A)") ""
    write (u, "(A)") "# Executable name"
    write (u, "(A)") "EXE = " // char (compilation%exe_name)
    write (u, "(A)") ""
    write (u, "(A)") "# Compiler"
    write (u, "(A)") "FC = " // char (os_data%fc)
    write (u, "(A)") "CXX = " // char (os_data%cxx)
    write (u, "(A)") ""
    write (u, "(A)") "# Included libraries"
    write (u, "(A)") "FCINCL = " // char (os_data%whizard_includes)
    write (u, "(A)") ""
    write (u, "(A)") "# Compiler flags"
    write (u, "(A)") "FCFLAGS = " // char (os_data%fcflags)
    write (u, "(A)") "FCLIBS = " // char (os_data%fclibs)
    write (u, "(A)") "CXXFLAGS = " // char (os_data%cxxflags)
    write (u, "(A)") "CXXLIBSS = " // char (os_data%cxxlibs)
    write (u, "(A)") "LDFLAGS = " // char (os_data%ldflags)
    write (u, "(A)") "LDFLAGS_STATIC = " // char (os_data%ldflags_static)
    write (u, "(A)") "LDFLAGS_HEPMC = " // char (os_data%ldflags_hepmc)
    write (u, "(A)") "LDFLAGS_LCIO = " // char (os_data%ldflags_lcio)
    write (u, "(A)") "LDFLAGS_HOPPET = " // char (os_data%ldflags_hoppet)
    write (u, "(A)") "LDFLAGS_LOOPTOOLS = " // char (os_data%ldflags_looptools)
    write (u, "(A)") "LDWHIZARD = " // char (os_data%whizard_ldflags)
    write (u, "(A)") ""
    write (u, "(A)") "# Libtool"
    write (u, "(A)") "LIBTOOL = " // char (os_data%whizard_libtool)
    if (verbose) then
       write (u, "(A)") "FCOMPILE = $(LIBTOOL) --tag=FC --mode=compile"
       if (OS_IS_DARWIN .and. .not. overwrite) then
          write (u, "(A)") "LINK = $(LIBTOOL) --tag=CXX --mode=link"
       else
          write (u, "(A)") "LINK = $(LIBTOOL) --tag=FC --mode=link"
       end if
    else
       write (u, "(A)") "FCOMPILE = @$(LIBTOOL) --silent --tag=FC --mode=compile"
       if (OS_IS_DARWIN .and. .not. overwrite) then
          write (u, "(A)") "LINK = @$(LIBTOOL) --silent --tag=CXX --mode=link"
       else
          write (u, "(A)") "LINK = @$(LIBTOOL) --silent --tag=FC --mode=link"
       end if
    end if
    write (u, "(A)") ""
    write (u, "(A)") "# Compile commands (default)"
    write (u, "(A)") "LTFCOMPILE = $(FCOMPILE) $(FC) -c $(FCINCL) $(FCFLAGS)"
    write (u, "(A)") ""
    write (u, "(A)") "# Default target"
    write (u, "(A)") "all: link"
    write (u, "(A)") ""
    write (u, "(A)") "# Libraries"
    do i = 1, size (compilation%lib_name)
       associate (lib_name => compilation%lib_name(i))
         write (u, "(A)") "LIBRARIES += " // char (lib_name) // ".la"
         write (u, "(A)") char (lib_name) // ".la:"
         write (u, "(A)") TAB // "$(MAKE) -f " // char (lib_name) // ".makefile"
       end associate
    end do
    write (u, "(A)") ""
    write (u, "(A)") "# Library dispatcher"
    write (u, "(A)") "DISP = $(EXE)_prclib_dispatcher"
    write (u, "(A)") "$(DISP).lo: $(DISP).f90 $(LIBRARIES)"
    if (.not. verbose) then
       write (u, "(A)")  TAB // '@echo  "  FC       " $@'
    end if
    write (u, "(A)") TAB // "$(LTFCOMPILE) $<"
    write (u, "(A)") ""
    write (u, "(A)") "# Executable"
    write (u, "(A)") "$(EXE): $(DISP).lo $(LIBRARIES)"
    if (.not. verbose) then
       if (OS_IS_DARWIN .and. .not. overwrite) then
          write (u, "(A)")  TAB // '@echo  "  CXXLD    " $@'
       else
          write (u, "(A)")  TAB // '@echo  "  FCLD     " $@'
       end if
    end if
    if (OS_IS_DARWIN .and. .not. overwrite) then
       write (u, "(A)") TAB // "$(LINK) $(CXX) -static $(CXXFLAGS) \"
    else
       write (u, "(A)") TAB // "$(LINK) $(FC) -static $(FCFLAGS) \"
    end if
    write (u, "(A)") TAB // "   $(LDWHIZARD) $(LDFLAGS) \"
    write (u, "(A)") TAB // "   -o $(EXE) $^ \"
    write (u, "(A)") TAB // "   $(LDFLAGS_HEPMC) $(LDFLAGS_LCIO) $(LDFLAGS_HOPPET) \"
    if (OS_IS_DARWIN .and. .not. overwrite) then
       write (u, "(A)") TAB // "   $(LDFLAGS_LOOPTOOLS) $(LDFLAGS_STATIC) \"
       write (u, "(A)") TAB // "   $(CXXLIBS) $(FCLIBS)" // char (ext_tag)
    else
       write (u, "(A)") TAB // "   $(LDFLAGS_LOOPTOOLS) $(LDFLAGS_STATIC)" // char (ext_tag)
    end if
    write (u, "(A)") ""
    write (u, "(A)") "# Main targets"
    write (u, "(A)") "link: compile $(EXE)"
    write (u, "(A)") "compile: $(LIBRARIES) $(DISP).lo"
    write (u, "(A)") ".PHONY: link compile"
    write (u, "(A)") ""
    write (u, "(A)") "# Cleanup targets"
    write (u, "(A)") "clean-exe:"
    write (u, "(A)") TAB // "rm -f $(EXE)"
    write (u, "(A)") "clean-objects:"
    write (u, "(A)") TAB // "rm -f $(DISP).lo"
    write (u, "(A)") "clean-source:"
    write (u, "(A)") TAB // "rm -f $(DISP).f90"
    write (u, "(A)") "clean-makefile:"
    write (u, "(A)") TAB // "rm -f $(EXE).makefile"
    write (u, "(A)") ""
    write (u, "(A)") "clean: clean-exe clean-objects clean-source"
    write (u, "(A)") "distclean: clean clean-makefile"
    write (u, "(A)") ".PHONY: clean distclean"
    close (u)
  end subroutine compilation_write_makefile

  module subroutine compilation_make_compile (compilation, os_data)
    class(compilation_t), intent(in) :: compilation
    type(os_data_t), intent(in) :: os_data
    call os_system_call ("make compile " // os_data%makeflags &
         // " -f " // compilation%exe_name // ".makefile")
  end subroutine compilation_make_compile

  module subroutine compilation_make_link (compilation, os_data)
    class(compilation_t), intent(in) :: compilation
    type(os_data_t), intent(in) :: os_data
    call os_system_call ("make link " // os_data%makeflags &
         // " -f " // compilation%exe_name // ".makefile")
  end subroutine compilation_make_link

  module subroutine compilation_make_clean_exe (compilation, os_data)
    class(compilation_t), intent(in) :: compilation
    type(os_data_t), intent(in) :: os_data
    call os_system_call ("make clean-exe " // os_data%makeflags &
         // " -f " // compilation%exe_name // ".makefile")
  end subroutine compilation_make_clean_exe

  module subroutine compile_executable (exename, libname, global)
    type(string_t), intent(in) :: exename
    type(string_t), dimension(:), intent(in) :: libname
    type(rt_data_t), intent(inout), target :: global
    type(compilation_t) :: compilation
    type(compilation_item_t) :: item
    type(string_t) :: ext_libtag
    logical :: force, recompile, verbose
    integer :: i
    ext_libtag = ""
    force = &
         global%var_list%get_lval (var_str ("?rebuild_library"))
    recompile = &
         global%var_list%get_lval (var_str ("?recompile_library"))
    verbose = &
         global%var_list%get_lval (var_str ("?me_verbose"))
    call compilation%init (exename, [libname])
    if (signal_is_pending ())  return
    call compilation%write_dispatcher ()
    if (signal_is_pending ())  return
    do i = 1, size (libname)
       call item%init (libname(i), global%prclib_stack, global%var_list)
       call item%compile (global%model, global%os_data, &
            force=force, recompile=recompile)
       ext_libtag = "" // item%lib%get_static_modelname (global%os_data)
       if (signal_is_pending ())  return
       call item%success ()
    end do
    call compilation%write_makefile &
         (global%os_data, ext_libtag=ext_libtag, verbose=verbose)
    if (signal_is_pending ())  return
    call compilation%make_compile (global%os_data)
    if (signal_is_pending ())  return
    call compilation%make_link (global%os_data)
  end subroutine compile_executable


end submodule compilations_s

