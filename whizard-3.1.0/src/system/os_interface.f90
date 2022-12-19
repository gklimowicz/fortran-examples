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

module os_interface

  use, intrinsic :: iso_c_binding !NODEP!

  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: paths_t
  public :: paths_init
  public :: os_data_t
  public :: dlaccess_t
  public :: dlaccess_init
  public :: dlaccess_final
  public :: dlaccess_has_error
  public :: dlaccess_get_error
  public :: dlaccess_get_c_funptr
  public :: dlaccess_is_open
  public :: os_system_call
  public :: os_dir_exist
  public :: os_file_exist
  public :: os_pack_file
  public :: os_unpack_file
  public :: os_compile_shared
  public :: os_link_shared
  public :: os_link_static
  public :: os_get_dlname
  public :: openmp_set_num_threads_verbose
  public :: mpi_set_logging
  public :: mpi_get_comm_id
  public :: mpi_is_comm_master

  type :: paths_t
     type(string_t) :: prefix
     type(string_t) :: exec_prefix
     type(string_t) :: bindir
     type(string_t) :: libdir
     type(string_t) :: includedir
     type(string_t) :: datarootdir
     type(string_t) :: localprefix
     type(string_t) :: libtool
     type(string_t) :: lhapdfdir
  end type paths_t

  type :: os_data_t
     logical :: use_libtool
     logical :: use_testfiles
     type(string_t) :: fc
     type(string_t) :: fcflags
     type(string_t) :: fcflags_pic
     type(string_t) :: fclibs
     type(string_t) :: fc_src_ext
     type(string_t) :: cc
     type(string_t) :: cflags
     type(string_t) :: cflags_pic
     type(string_t) :: cxx
     type(string_t) :: cxxflags
     type(string_t) :: cxxlibs
     type(string_t) :: obj_ext
     type(string_t) :: ld
     type(string_t) :: ldflags
     type(string_t) :: ldflags_so
     type(string_t) :: ldflags_static
     type(string_t) :: ldflags_hepmc
     type(string_t) :: ldflags_lcio
     type(string_t) :: ldflags_hoppet
     type(string_t) :: ldflags_looptools
     type(string_t) :: shrlib_ext
     type(string_t) :: fc_shrlib_ext
     type(string_t) :: pack_cmd
     type(string_t) :: unpack_cmd
     type(string_t) :: pack_ext
     type(string_t) :: makeflags
     type(string_t) :: prefix
     type(string_t) :: exec_prefix
     type(string_t) :: bindir
     type(string_t) :: libdir
     type(string_t) :: includedir
     type(string_t) :: datarootdir
     type(string_t) :: whizard_omega_binpath
     type(string_t) :: whizard_includes
     type(string_t) :: whizard_ldflags
     type(string_t) :: whizard_libtool
     type(string_t) :: whizard_modelpath
     type(string_t) :: whizard_modelpath_ufo
     type(string_t) :: whizard_models_libpath
     type(string_t) :: whizard_susypath
     type(string_t) :: whizard_gmlpath
     type(string_t) :: whizard_cutspath
     type(string_t) :: whizard_texpath
     type(string_t) :: whizard_sharepath
     type(string_t) :: whizard_testdatapath
     type(string_t) :: whizard_modelpath_local
     type(string_t) :: whizard_models_libpath_local
     type(string_t) :: whizard_omega_binpath_local
     type(string_t) :: whizard_circe2path
     type(string_t) :: whizard_beamsimpath
     type(string_t) :: whizard_mulipath
     type(string_t) :: pdf_builtin_datapath
     logical :: event_analysis = .false.
     logical :: event_analysis_ps  = .false.
     logical :: event_analysis_pdf = .false.
     type(string_t) :: latex
     type(string_t) :: mpost
     type(string_t) :: gml
     type(string_t) :: dvips
     type(string_t) :: ps2pdf
     type(string_t) :: gosampath
     type(string_t) :: golempath
     type(string_t) :: formpath
     type(string_t) :: qgrafpath
     type(string_t) :: ninjapath
     type(string_t) :: samuraipath
   contains
     procedure :: init => os_data_init
     procedure :: write => os_data_write
     procedure :: build_latex_file => os_data_build_latex_file
  end type os_data_t

  type :: dlaccess_t
     private
     type(string_t) :: filename
     type(c_ptr) :: handle = c_null_ptr
     logical :: is_open = .false.
     logical :: has_error = .false.
     type(string_t) :: error
   contains
     procedure :: write => dlaccess_write
     procedure :: init => dlaccess_init
     procedure :: final => dlaccess_final
  end type dlaccess_t


  interface
     function dlopen (filename, flag) result (handle) bind(C)
       import
       character(c_char), dimension(*) :: filename
       integer(c_int), value :: flag
       type(c_ptr) :: handle
     end function dlopen
  end interface

  interface
     function dlclose (handle) result (status) bind(C)
       import
       type(c_ptr), value :: handle
       integer(c_int) :: status
     end function dlclose
  end interface

  interface
     function dlerror () result (str) bind(C)
       import
       type(c_ptr) :: str
     end function dlerror
  end interface

  interface
     function dlsym (handle, symbol) result (fptr) bind(C)
       import
       type(c_ptr), value :: handle
       character(c_char), dimension(*) :: symbol
       type(c_funptr) :: fptr
     end function dlsym
  end interface

  interface
     function system (command) result (status) bind(C)
       import
       integer(c_int) :: status
       character(c_char), dimension(*) :: command
     end function system
  end interface


  interface
    module subroutine paths_init (paths)
      type(paths_t), intent(out) :: paths
    end subroutine paths_init
    module subroutine os_data_init (os_data, paths)
      class(os_data_t), intent(out) :: os_data
      type(paths_t), intent(in), optional :: paths
    end subroutine os_data_init
    module subroutine os_data_write (os_data, unit)
      class(os_data_t), intent(in) :: os_data
      integer, intent(in), optional :: unit
    end subroutine os_data_write
    module subroutine os_data_build_latex_file (os_data, filename, stat_out)
      class(os_data_t), intent(in) :: os_data
      type(string_t), intent(in) :: filename
      integer, intent(out), optional :: stat_out
    end subroutine os_data_build_latex_file
    module subroutine dlaccess_write (object, unit)
      class(dlaccess_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine dlaccess_write
    module subroutine dlaccess_init (dlaccess, prefix, libname, os_data)
      class(dlaccess_t), intent(out) :: dlaccess
      type(string_t), intent(in) :: prefix, libname
      type(os_data_t), intent(in), optional :: os_data
    end subroutine dlaccess_init
    module subroutine dlaccess_final (dlaccess)
      class(dlaccess_t), intent(inout) :: dlaccess
    end subroutine dlaccess_final
   module function dlaccess_has_error (dlaccess) result (flag)
      logical :: flag
      type(dlaccess_t), intent(in) :: dlaccess
    end function dlaccess_has_error
    module function dlaccess_get_error (dlaccess) result (error)
      type(string_t) :: error
      type(dlaccess_t), intent(in) :: dlaccess
    end function dlaccess_get_error
    module function dlaccess_get_c_funptr (dlaccess, fname) result (fptr)
      type(c_funptr) :: fptr
      type(dlaccess_t), intent(inout) :: dlaccess
      type(string_t), intent(in) :: fname
    end function dlaccess_get_c_funptr
    module function dlaccess_is_open (dlaccess) result (flag)
      logical :: flag
      type(dlaccess_t), intent(in) :: dlaccess
    end function dlaccess_is_open
    module subroutine os_system_call (command_string, status, verbose)
      type(string_t), intent(in) :: command_string
      integer, intent(out), optional :: status
      logical, intent(in), optional :: verbose
    end subroutine os_system_call
    module function os_dir_exist (name) result (res)
      type(string_t), intent(in) :: name
      logical :: res
    end function os_dir_exist
    module function os_file_exist (name) result (exist)
      type(string_t), intent(in) :: name
      logical :: exist
    end function os_file_exist
    module subroutine os_pack_file (file, os_data, status)
      type(string_t), intent(in) :: file
      type(os_data_t), intent(in) :: os_data
      integer, intent(out), optional :: status
    end subroutine os_pack_file
    module subroutine os_unpack_file (file, os_data, status)
      type(string_t), intent(in) :: file
      type(os_data_t), intent(in) :: os_data
      integer, intent(out), optional :: status
    end subroutine os_unpack_file
    module subroutine os_compile_shared (src, os_data, status)
      type(string_t), intent(in) :: src
      type(os_data_t), intent(in) :: os_data
      integer, intent(out), optional :: status
    end subroutine os_compile_shared
    module subroutine os_link_shared (objlist, lib, os_data, status)
      type(string_t), intent(in) :: objlist, lib
      type(os_data_t), intent(in) :: os_data
      integer, intent(out), optional :: status
    end subroutine os_link_shared
    module subroutine os_link_static (objlist, exec_name, os_data, status)
      type(string_t), intent(in) :: objlist, exec_name
      type(os_data_t), intent(in) :: os_data
      integer, intent(out), optional :: status
    end subroutine os_link_static
    module function os_get_dlname (lib, os_data, ignore, silent) result (dlname)
      type(string_t) :: dlname
      type(string_t), intent(in) :: lib
      type(os_data_t), intent(in) :: os_data
      logical, intent(in), optional :: ignore, silent
    end function os_get_dlname
    module subroutine openmp_set_num_threads_verbose (num_threads, openmp_logging)
      integer, intent(in) :: num_threads
      logical, intent(in), optional :: openmp_logging
    end subroutine openmp_set_num_threads_verbose
    module subroutine mpi_set_logging (mpi_logging)
      logical, intent(in) :: mpi_logging
    end subroutine mpi_set_logging
    module subroutine mpi_get_comm_id (n_size, rank)
      integer, intent(out) :: n_size
      integer, intent(out) :: rank
    end subroutine mpi_get_comm_id
    module function mpi_is_comm_master () result (flag)
      logical :: flag
    end function mpi_is_comm_master
  end interface

end module os_interface
