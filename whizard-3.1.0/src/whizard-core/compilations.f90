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

module compilations

  use iso_varying_string, string_t => varying_string
  use os_interface
  use variables, only: var_list_t
  use model_data
  use process_libraries
  use prclib_stacks
  use rt_data

  implicit none
  private

  public :: compilation_item_t
  public :: compile_library
  public :: compilation_t
  public :: compile_executable

  type :: compilation_item_t
     private
     type(string_t) :: libname
     type(string_t) :: static_external_tag
     type(process_library_t), pointer :: lib => null ()
     logical :: recompile_library = .false.
     logical :: verbose = .false.
     logical :: use_workspace = .false.
     type(string_t) :: workspace
   contains
     procedure :: init => compilation_item_init
     procedure :: compile => compilation_item_compile
     procedure :: load => compilation_item_load
     procedure :: success => compilation_item_success
  end type compilation_item_t

  type :: compilation_t
     private
     type(string_t) :: exe_name
     type(string_t), dimension(:), allocatable :: lib_name
   contains
     procedure :: write => compilation_write
     procedure :: init => compilation_init
     procedure :: write_dispatcher => compilation_write_dispatcher
     procedure :: write_makefile => compilation_write_makefile
     procedure :: make_compile => compilation_make_compile
     procedure :: make_link => compilation_make_link
     procedure :: make_clean_exe => compilation_make_clean_exe
  end type compilation_t


  character(*), parameter :: ALLOWED_IN_DIRNAME = &
       "abcdefghijklmnopqrstuvwxyz&
       &ABCDEFGHIJKLMNOPQRSTUVWXYZ&
       &1234567890&
       &.,_-+="

  interface
    module subroutine compilation_item_init (comp, libname, stack, var_list)
      class(compilation_item_t), intent(out) :: comp
      type(string_t), intent(in) :: libname
      type(prclib_stack_t), intent(inout) :: stack
      type(var_list_t), intent(in) :: var_list
    end subroutine compilation_item_init
    module subroutine compilation_item_compile &
         (comp, model, os_data, force, recompile)
      class(compilation_item_t), intent(inout) :: comp
      class(model_data_t), intent(in), target :: model
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: force, recompile
    end subroutine compilation_item_compile
    module subroutine compilation_item_load (comp, os_data)
      class(compilation_item_t), intent(inout) :: comp
      type(os_data_t), intent(in) :: os_data
    end subroutine compilation_item_load
    module subroutine compilation_item_success (comp)
      class(compilation_item_t), intent(in) :: comp
    end subroutine compilation_item_success
    module subroutine compile_library (libname, global)
      type(string_t), intent(in) :: libname
      type(rt_data_t), intent(inout), target :: global
    end subroutine compile_library
    module subroutine compilation_write (object, unit)
      class(compilation_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine compilation_write
    module subroutine compilation_init (compilation, exe_name, lib_name)
      class(compilation_t), intent(out) :: compilation
      type(string_t), intent(in) :: exe_name
      type(string_t), dimension(:), intent(in) :: lib_name
    end subroutine compilation_init
    module subroutine compilation_write_dispatcher (compilation)
      class(compilation_t), intent(in) :: compilation
    end subroutine compilation_write_dispatcher
    module subroutine compilation_write_makefile &
         (compilation, os_data, ext_libtag, verbose, overwrite_os)
      class(compilation_t), intent(in) :: compilation
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      logical, intent(in), optional :: overwrite_os
      type(string_t), intent(in), optional :: ext_libtag
    end subroutine compilation_write_makefile
    module subroutine compilation_make_compile (compilation, os_data)
      class(compilation_t), intent(in) :: compilation
      type(os_data_t), intent(in) :: os_data
    end subroutine compilation_make_compile
    module subroutine compilation_make_link (compilation, os_data)
      class(compilation_t), intent(in) :: compilation
      type(os_data_t), intent(in) :: os_data
    end subroutine compilation_make_link
    module subroutine compilation_make_clean_exe (compilation, os_data)
      class(compilation_t), intent(in) :: compilation
      type(os_data_t), intent(in) :: os_data
    end subroutine compilation_make_clean_exe
    module subroutine compile_executable (exename, libname, global)
      type(string_t), intent(in) :: exename
      type(string_t), dimension(:), intent(in) :: libname
      type(rt_data_t), intent(inout), target :: global
    end subroutine compile_executable
  end interface

end module compilations
