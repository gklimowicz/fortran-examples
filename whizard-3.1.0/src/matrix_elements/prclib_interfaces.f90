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

module prclib_interfaces

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use os_interface

  implicit none
  private

  public :: prc_writer_t
  public :: prc_writer_f_module_t
  public :: prc_writer_c_lib_t
  public :: prclib_driver_t
  public :: prclib_driver_dynamic_t
  public :: dispatch_prclib_driver
  public :: prc_get_n_processes
  public :: prc_get_stringptr
  public :: prc_get_log
  public :: prc_get_int
  public :: prc_set_int_tab1
  public :: prc_set_col_state
  public :: prc_set_color_factors
  public :: prc_get_fptr
  public :: write_driver_code
  public :: prclib_unload_hook
  public :: prclib_reload_hook
  public :: workspace_prefix
  public :: workspace_path

  type, abstract :: prc_writer_t
     character(32) :: md5sum = ""
   contains
     procedure(get_const_string), nopass, deferred :: type_name
     procedure, nopass :: get_procname => prc_writer_get_procname
     procedure :: get_c_procname => prc_writer_get_c_procname
     procedure(write_feature_code), deferred :: write_interface
     procedure(write_code_os), deferred :: write_makefile_code
     procedure(write_code_file), deferred :: write_source_code
     procedure(write_code_file), deferred :: before_compile
     procedure(write_code_file), deferred :: after_compile
     procedure :: init_test => prc_writer_init_test
     procedure(write_code), deferred :: write_md5sum_call
     procedure(write_feature_code), deferred :: write_int_sub_call
     procedure(write_code), deferred :: write_col_state_call
     procedure(write_code), deferred :: write_color_factors_call
  end type prc_writer_t

  type, extends (prc_writer_t), abstract :: prc_writer_f_module_t
   contains
     procedure, nopass :: get_module_name => prc_writer_get_module_name
     procedure :: write_use_line => prc_writer_write_use_line
     procedure(prc_write_wrapper), deferred :: write_wrapper
     procedure :: write_md5sum_call => prc_writer_f_module_write_md5sum_call
     procedure :: write_int_sub_call => prc_writer_f_module_write_int_sub_call
     procedure :: write_col_state_call => prc_writer_f_module_write_col_state_call
     procedure :: write_color_factors_call => prc_writer_f_module_write_color_factors_call
  end type prc_writer_f_module_t

  type, extends (prc_writer_t), abstract :: prc_writer_c_lib_t
   contains
     procedure :: write_md5sum_call => prc_writer_c_lib_write_md5sum_call
     procedure :: write_int_sub_call => prc_writer_c_lib_write_int_sub_call
     procedure :: write_col_state_call => prc_writer_c_lib_write_col_state_call
     procedure :: write_color_factors_call => &
          prc_writer_c_lib_write_color_factors_call
     procedure :: write_standard_interface => prc_writer_c_lib_write_interface
  end type prc_writer_c_lib_t

  type :: prclib_driver_record_t
     type(string_t) :: id
     type(string_t) :: model_name
     type(string_t), dimension(:), allocatable :: feature
     class(prc_writer_t), pointer :: writer => null ()
   contains
     procedure :: write => prclib_driver_record_write
     procedure :: get_c_procname => prclib_driver_record_get_c_procname
     procedure :: write_use_line => prclib_driver_record_write_use_line
     procedure :: write_interface => prclib_driver_record_write_interface
     procedure :: write_interfaces => prclib_driver_record_write_interfaces
     procedure :: write_wrappers => prclib_driver_record_write_wrappers
     procedure :: write_makefile_code => prclib_driver_record_write_makefile_code
     procedure :: write_source_code => prclib_driver_record_write_source_code
     procedure :: before_compile => prclib_driver_record_before_compile
     procedure :: after_compile => prclib_driver_record_after_compile
     procedure :: write_md5sum_call => prclib_driver_record_write_md5sum_call
     procedure :: write_int_sub_call => prclib_driver_record_write_int_sub_call
     procedure :: write_col_state_call => prclib_driver_record_write_col_state_call
     procedure :: write_color_factors_call => prclib_driver_record_write_color_factors_call
  end type prclib_driver_record_t

  type, abstract :: prclib_driver_t
     type(string_t) :: basename
     character(32) :: md5sum = ""
     logical :: loaded = .false.
     type(string_t) :: libname
     type(string_t) :: modellibs_ldflags
     integer :: n_processes = 0
     type(prclib_driver_record_t), dimension(:), allocatable :: record
     procedure(prc_get_n_processes), nopass, pointer :: &
          get_n_processes => null ()
     procedure(prc_get_stringptr), nopass, pointer :: &
          get_process_id_ptr => null ()
     procedure(prc_get_stringptr), nopass, pointer :: &
          get_model_name_ptr => null ()
     procedure(prc_get_stringptr), nopass, pointer :: &
          get_md5sum_ptr => null ()
     procedure(prc_get_log), nopass, pointer :: &
          get_openmp_status => null ()
     procedure(prc_get_int), nopass, pointer :: get_n_in  => null ()
     procedure(prc_get_int), nopass, pointer :: get_n_out => null ()
     procedure(prc_get_int), nopass, pointer :: get_n_flv => null ()
     procedure(prc_get_int), nopass, pointer :: get_n_hel => null ()
     procedure(prc_get_int), nopass, pointer :: get_n_col => null ()
     procedure(prc_get_int), nopass, pointer :: get_n_cin => null ()
     procedure(prc_get_int), nopass, pointer :: get_n_cf  => null ()
     procedure(prc_set_int_tab1), nopass, pointer :: &
          set_flv_state_ptr => null ()
     procedure(prc_set_int_tab1), nopass, pointer :: &
          set_hel_state_ptr => null ()
     procedure(prc_set_col_state), nopass, pointer :: &
          set_col_state_ptr => null ()
     procedure(prc_set_color_factors), nopass, pointer :: &
          set_color_factors_ptr => null ()
     procedure(prc_get_fptr), nopass, pointer :: get_fptr => null ()
   contains
     procedure :: write => prclib_driver_write
     procedure :: init => prclib_driver_init
     procedure :: set_md5sum => prclib_driver_set_md5sum
     procedure :: set_record => prclib_driver_set_record
     procedure :: write_interfaces => prclib_driver_write_interfaces
     procedure :: generate_makefile => prclib_driver_generate_makefile
     procedure :: generate_driver_code => prclib_driver_generate_code
     procedure, nopass :: write_module => prclib_driver_write_module
     procedure :: write_lib_md5sum_fun => prclib_driver_write_lib_md5sum_fun
     procedure :: write_get_n_processes_fun
     procedure, nopass :: write_string_to_array_fun
     procedure :: write_get_process_id_fun
     procedure :: write_get_model_name_fun
     procedure :: write_get_md5sum_fun
     procedure :: get_process_id => prclib_driver_get_process_id
     procedure :: get_model_name => prclib_driver_get_model_name
     procedure :: get_md5sum => prclib_driver_get_md5sum
     procedure :: write_get_openmp_status_fun
     procedure :: write_get_int_fun
     procedure :: write_set_int_sub
     procedure :: write_set_col_state_sub
     procedure :: write_set_color_factors_sub
     procedure :: set_flv_state => prclib_driver_set_flv_state
     procedure :: set_hel_state => prclib_driver_set_hel_state
     procedure :: set_col_state => prclib_driver_set_col_state
     procedure :: set_color_factors => prclib_driver_set_color_factors
     procedure :: write_get_fptr_sub
     procedure :: make_source => prclib_driver_make_source
     procedure :: make_compile => prclib_driver_make_compile
     procedure :: make_link => prclib_driver_make_link
     procedure :: clean_library => prclib_driver_clean_library
     procedure :: clean_objects => prclib_driver_clean_objects
     procedure :: clean_source => prclib_driver_clean_source
     procedure :: clean_driver => prclib_driver_clean_driver
     procedure :: clean_makefile => prclib_driver_clean_makefile
     procedure :: clean => prclib_driver_clean
     procedure :: distclean => prclib_driver_distclean
     procedure :: clean_proc => prclib_driver_clean_proc
     procedure :: makefile_exists => prclib_driver_makefile_exists
     procedure :: load => prclib_driver_load
     procedure :: unload => prclib_driver_unload
     procedure (prclib_driver_get_c_funptr), deferred :: get_c_funptr
     procedure :: get_md5sum_makefile => prclib_driver_get_md5sum_makefile
     procedure :: get_md5sum_driver => prclib_driver_get_md5sum_driver
     procedure :: get_md5sum_source => prclib_driver_get_md5sum_source
  end type prclib_driver_t

  type, extends (prclib_driver_t) :: prclib_driver_dynamic_t
     type(dlaccess_t) :: dlaccess
   contains
     procedure :: check_dlerror => prclib_driver_check_dlerror
     procedure :: get_c_funptr => prclib_driver_dynamic_get_c_funptr
  end type prclib_driver_dynamic_t


  abstract interface
     function get_const_string () result (string)
       import
       type(string_t) :: string
     end function get_const_string
  end interface

  abstract interface
     subroutine write_code_file (writer, id)
       import
       class(prc_writer_t), intent(in) :: writer
       type(string_t), intent(in) :: id
     end subroutine write_code_file
  end interface

  abstract interface
     subroutine write_code (writer, unit, id)
       import
       class(prc_writer_t), intent(in) :: writer
       integer, intent(in) :: unit
       type(string_t), intent(in) :: id
     end subroutine write_code
  end interface

  abstract interface
     subroutine write_code_os &
          (writer, unit, id, os_data, verbose, testflag)
       import
       class(prc_writer_t), intent(in) :: writer
       integer, intent(in) :: unit
       type(string_t), intent(in) :: id
       type(os_data_t), intent(in) :: os_data
       logical, intent(in) :: verbose
       logical, intent(in), optional :: testflag
     end subroutine write_code_os
  end interface

  abstract interface
     subroutine write_feature_code (writer, unit, id, feature)
       import
       class(prc_writer_t), intent(in) :: writer
       integer, intent(in) :: unit
       type(string_t), intent(in) :: id, feature
     end subroutine write_feature_code
  end interface

  abstract interface
     subroutine prc_write_wrapper (writer, unit, id, feature)
       import
       class(prc_writer_f_module_t), intent(in) :: writer
       integer, intent(in) :: unit
       type(string_t), intent(in) :: id, feature
     end subroutine prc_write_wrapper
  end interface

  abstract interface
     function prc_get_n_processes () result (n) bind(C)
       import
       integer(c_int) :: n
     end function prc_get_n_processes
  end interface
  abstract interface
     subroutine prc_get_stringptr (i, cptr, len) bind(C)
       import
       integer(c_int), intent(in) :: i
       type(c_ptr), intent(out) :: cptr
       integer(c_int), intent(out) :: len
     end subroutine prc_get_stringptr
  end interface
  abstract interface
     function prc_get_log (pid) result (l) bind(C)
       import
       integer(c_int), intent(in) :: pid
       logical(c_bool) :: l
     end function prc_get_log
  end interface
  abstract interface
     function prc_get_int (pid) result (n) bind(C)
       import
       integer(c_int), intent(in) :: pid
       integer(c_int) :: n
     end function prc_get_int
  end interface
  abstract interface
     subroutine prc_set_int_tab1 (pid, tab, shape) bind(C)
       import
       integer(c_int), intent(in) :: pid
       integer(c_int), dimension(*), intent(out) :: tab
       integer(c_int), dimension(2), intent(in) :: shape
     end subroutine prc_set_int_tab1
  end interface
  abstract interface
     subroutine prc_set_col_state (pid, col_state, ghost_flag, shape) bind(C)
       import
       integer(c_int), intent(in) :: pid
       integer(c_int), dimension(*), intent(out) :: col_state
       logical(c_bool), dimension(*), intent(out) :: ghost_flag
       integer(c_int), dimension(3), intent(in) :: shape
     end subroutine prc_set_col_state
  end interface
  abstract interface
     subroutine prc_set_color_factors &
          (pid, cf_index1, cf_index2, color_factors, shape) bind(C)
       import
       integer(c_int), intent(in) :: pid
       integer(c_int), dimension(*), intent(out) :: cf_index1, cf_index2
       complex(c_default_complex), dimension(*), intent(out) :: color_factors
       integer(c_int), dimension(1), intent(in) :: shape
     end subroutine prc_set_color_factors
  end interface

  abstract interface
     subroutine prc_get_fptr (pid, fid, fptr) bind(C)
       import
       integer(c_int), intent(in) :: pid
       integer(c_int), intent(in) :: fid
       type(c_funptr), intent(out) :: fptr
     end subroutine prc_get_fptr
  end interface

  abstract interface
     subroutine write_driver_code (unit, prefix, id, procname)
       import
       integer, intent(in) :: unit
       type(string_t), intent(in) :: prefix
       type(string_t), intent(in) :: id
       type(string_t), intent(in) :: procname
     end subroutine write_driver_code
  end interface

  abstract interface
     subroutine prclib_unload_hook (libname)
       import
       type(string_t), intent(in) :: libname
     end subroutine prclib_unload_hook

     subroutine prclib_reload_hook (libname)
       import
       type(string_t), intent(in) :: libname
     end subroutine prclib_reload_hook
  end interface
  abstract interface
     function prclib_driver_get_c_funptr (driver, feature) result (c_fptr)
       import
       class(prclib_driver_t), intent(inout) :: driver
       type(string_t), intent(in) :: feature
       type(c_funptr) :: c_fptr
     end function prclib_driver_get_c_funptr
  end interface


  interface
    module subroutine prc_writer_write_use_line (writer, unit, id, feature)
      class(prc_writer_f_module_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t) :: id, feature
    end subroutine prc_writer_write_use_line
    module subroutine prc_writer_init_test (writer)
      class(prc_writer_t), intent(out) :: writer
    end subroutine prc_writer_init_test
    module subroutine prclib_driver_record_write (object, unit)
      class(prclib_driver_record_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine prclib_driver_record_write
    module subroutine prclib_driver_record_write_use_line (record, unit, feature)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
      type(string_t), intent(in) :: feature
    end subroutine prclib_driver_record_write_use_line
    module subroutine prclib_driver_record_write_interface (record, unit, feature)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
      type(string_t), intent(in) :: feature
    end subroutine prclib_driver_record_write_interface
    module subroutine prclib_driver_record_write_interfaces (record, unit)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
    end subroutine prclib_driver_record_write_interfaces
    module subroutine prclib_driver_record_write_wrappers (record, unit)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
    end subroutine prclib_driver_record_write_wrappers
    module subroutine prclib_driver_record_write_makefile_code &
         (record, unit, os_data, verbose, testflag)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      logical, intent(in), optional :: testflag
    end subroutine prclib_driver_record_write_makefile_code
    module subroutine prclib_driver_record_write_source_code (record)
      class(prclib_driver_record_t), intent(in) :: record
    end subroutine prclib_driver_record_write_source_code
    module subroutine prclib_driver_record_before_compile (record)
      class(prclib_driver_record_t), intent(in) :: record
    end subroutine prclib_driver_record_before_compile
    module subroutine prclib_driver_record_after_compile (record)
      class(prclib_driver_record_t), intent(in) :: record
    end subroutine prclib_driver_record_after_compile
    module subroutine prclib_driver_write (object, unit, libpath)
      class(prclib_driver_t), intent(in) :: object
      integer, intent(in) :: unit
      logical, intent(in), optional :: libpath
    end subroutine prclib_driver_write
    module subroutine prclib_driver_init (driver, n_processes)
      class(prclib_driver_t), intent(inout) :: driver
      integer, intent(in) :: n_processes
    end subroutine prclib_driver_init
    module subroutine prclib_driver_set_md5sum (driver, md5sum)
      class(prclib_driver_t), intent(inout) :: driver
      character(32), intent(in) :: md5sum
    end subroutine prclib_driver_set_md5sum
    module subroutine prclib_driver_set_record (driver, i, &
         id, model_name, features, writer)
      class(prclib_driver_t), intent(inout) :: driver
      integer, intent(in) :: i
      type(string_t), intent(in) :: id
      type(string_t), intent(in) :: model_name
      type(string_t), dimension(:), intent(in) :: features
      class(prc_writer_t), intent(in), pointer :: writer
    end subroutine prclib_driver_set_record
    module subroutine prclib_driver_write_interfaces (driver, unit, feature)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: feature
    end subroutine prclib_driver_write_interfaces
    module subroutine prclib_driver_generate_makefile (driver, unit, os_data, verbose, testflag)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      logical, intent(in), optional :: testflag
    end subroutine prclib_driver_generate_makefile
    module subroutine prclib_driver_generate_code (driver, unit)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
    end subroutine prclib_driver_generate_code
    module subroutine prclib_driver_write_module (unit, prefix)
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine prclib_driver_write_module
    module subroutine prclib_driver_write_lib_md5sum_fun (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine prclib_driver_write_lib_md5sum_fun
    module subroutine write_get_n_processes_fun (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_get_n_processes_fun
    module subroutine write_string_to_array_fun (unit, prefix)
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_string_to_array_fun
    module subroutine write_get_process_id_fun (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_get_process_id_fun
    module subroutine write_get_model_name_fun (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_get_model_name_fun
    module subroutine write_get_md5sum_fun (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_get_md5sum_fun
    module subroutine prclib_driver_record_write_md5sum_call (record, unit)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
    end subroutine prclib_driver_record_write_md5sum_call
    module subroutine prc_writer_f_module_write_md5sum_call (writer, unit, id)
      class(prc_writer_f_module_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
    end subroutine prc_writer_f_module_write_md5sum_call
    module subroutine prc_writer_c_lib_write_md5sum_call (writer, unit, id)
      class(prc_writer_c_lib_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
    end subroutine prc_writer_c_lib_write_md5sum_call
    module function prclib_driver_get_process_id (driver, i) result (string)
      type(string_t) :: string
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
    end function prclib_driver_get_process_id
    module function prclib_driver_get_model_name (driver, i) result (string)
      type(string_t) :: string
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
    end function prclib_driver_get_model_name
    module function prclib_driver_get_md5sum (driver, i) result (string)
      type(string_t) :: string
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
    end function prclib_driver_get_md5sum
    module subroutine write_get_openmp_status_fun (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_get_openmp_status_fun
    module subroutine write_get_int_fun (driver, unit, prefix, feature)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
      type(string_t), intent(in) :: feature
    end subroutine write_get_int_fun
    module subroutine write_set_int_sub (driver, unit, prefix, feature)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
      type(string_t), intent(in) :: feature
    end subroutine write_set_int_sub
    module subroutine prclib_driver_record_write_int_sub_call (record, unit, feature)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
      type(string_t), intent(in) :: feature
    end subroutine prclib_driver_record_write_int_sub_call
    module subroutine prc_writer_f_module_write_int_sub_call (writer, unit, id, feature)
      class(prc_writer_f_module_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id, feature
    end subroutine prc_writer_f_module_write_int_sub_call
    module subroutine prc_writer_c_lib_write_int_sub_call (writer, unit, id, feature)
      class(prc_writer_c_lib_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id, feature
    end subroutine prc_writer_c_lib_write_int_sub_call
    module subroutine write_set_col_state_sub (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_set_col_state_sub
    module subroutine prclib_driver_record_write_col_state_call (record, unit)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
    end subroutine prclib_driver_record_write_col_state_call
    module subroutine prc_writer_f_module_write_col_state_call (writer, unit, id)
      class(prc_writer_f_module_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
    end subroutine prc_writer_f_module_write_col_state_call
    module subroutine prc_writer_c_lib_write_col_state_call (writer, unit, id)
      class(prc_writer_c_lib_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
    end subroutine prc_writer_c_lib_write_col_state_call
    module subroutine write_set_color_factors_sub (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_set_color_factors_sub
    module subroutine prclib_driver_record_write_color_factors_call (record, unit)
      class(prclib_driver_record_t), intent(in) :: record
      integer, intent(in) :: unit
    end subroutine prclib_driver_record_write_color_factors_call
    module subroutine prc_writer_f_module_write_color_factors_call (writer, unit, id)
      class(prc_writer_f_module_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
    end subroutine prc_writer_f_module_write_color_factors_call
    module subroutine prc_writer_c_lib_write_color_factors_call (writer, unit, id)
      class(prc_writer_c_lib_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
    end subroutine prc_writer_c_lib_write_color_factors_call
    module subroutine prc_writer_c_lib_write_interface (writer, unit, id, feature)
      class(prc_writer_c_lib_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id, feature
    end subroutine prc_writer_c_lib_write_interface
    module subroutine prclib_driver_set_flv_state (driver, i, flv_state)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
      integer, dimension(:,:), allocatable, intent(out) :: flv_state
    end subroutine prclib_driver_set_flv_state
    module subroutine prclib_driver_set_hel_state (driver, i, hel_state)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
      integer, dimension(:,:), allocatable, intent(out) :: hel_state
    end subroutine prclib_driver_set_hel_state
    module subroutine prclib_driver_set_col_state (driver, i, col_state, ghost_flag)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
      integer, dimension(:,:,:), allocatable, intent(out) :: col_state
      logical, dimension(:,:), allocatable, intent(out) :: ghost_flag
    end subroutine prclib_driver_set_col_state
    module subroutine prclib_driver_set_color_factors (driver, i, color_factors, cf_index)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
      complex(default), dimension(:), allocatable, intent(out) :: color_factors
      integer, dimension(:,:), allocatable, intent(out) :: cf_index
    end subroutine prclib_driver_set_color_factors
    module subroutine write_get_fptr_sub (driver, unit, prefix)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: unit
      type(string_t), intent(in) :: prefix
    end subroutine write_get_fptr_sub
    module function workspace_prefix (workspace) result (prefix)
      type(string_t), intent(in), optional :: workspace
      type(string_t) :: prefix
    end function workspace_prefix
    module function workspace_path (workspace) result (path)
      type(string_t), intent(in), optional :: workspace
      type(string_t) :: path
    end function workspace_path
    module function workspace_cmd (workspace) result (cmd)
      type(string_t), intent(in), optional :: workspace
      type(string_t) :: cmd
    end function workspace_cmd
    module subroutine prclib_driver_make_source (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_make_source
    module subroutine prclib_driver_make_compile (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_make_compile
    module subroutine prclib_driver_make_link (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_make_link
    module subroutine prclib_driver_clean_library (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_clean_library
    module subroutine prclib_driver_clean_objects (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_clean_objects
    module subroutine prclib_driver_clean_source (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_clean_source
    module subroutine prclib_driver_clean_driver (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_clean_driver
    module subroutine prclib_driver_clean_makefile (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_clean_makefile
    module subroutine prclib_driver_clean (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_clean
    module subroutine prclib_driver_distclean (driver, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_distclean
    module subroutine prclib_driver_clean_proc (driver, i, os_data, workspace)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_clean_proc
    module function prclib_driver_makefile_exists (driver, workspace) result (flag)
      class(prclib_driver_t), intent(in) :: driver
      type(string_t), intent(in), optional :: workspace
      logical :: flag
    end function prclib_driver_makefile_exists
    module subroutine prclib_driver_load (driver, os_data, noerror, workspace)
      class(prclib_driver_t), intent(inout) :: driver
      type(os_data_t), intent(in) :: os_data
      logical, intent(in), optional :: noerror
      type(string_t), intent(in), optional :: workspace
    end subroutine prclib_driver_load
    module subroutine prclib_driver_unload (driver)
      class(prclib_driver_t), intent(inout) :: driver
    end subroutine prclib_driver_unload
    module subroutine prclib_driver_check_dlerror (driver)
      class(prclib_driver_dynamic_t), intent(in) :: driver
    end subroutine prclib_driver_check_dlerror
    module function prclib_driver_dynamic_get_c_funptr &
         (driver, feature) result (c_fptr)
      class(prclib_driver_dynamic_t), intent(inout) :: driver
      type(string_t), intent(in) :: feature
      type(c_funptr) :: c_fptr
    end function prclib_driver_dynamic_get_c_funptr
    module function prclib_driver_get_md5sum_makefile &
         (driver, workspace) result (md5sum)
      class(prclib_driver_t), intent(in) :: driver
      type(string_t), intent(in), optional :: workspace
      character(32) :: md5sum
    end function prclib_driver_get_md5sum_makefile
    module function prclib_driver_get_md5sum_driver (driver, workspace) result (md5sum)
      class(prclib_driver_t), intent(in) :: driver
      type(string_t), intent(in), optional :: workspace
      character(32) :: md5sum
    end function prclib_driver_get_md5sum_driver
    module function prclib_driver_get_md5sum_source &
         (driver, i, workspace) result (md5sum)
      class(prclib_driver_t), intent(in) :: driver
      integer, intent(in) :: i
      type(string_t), intent(in), optional :: workspace
      character(32) :: md5sum
    end function prclib_driver_get_md5sum_source
  end interface

contains

  function prc_writer_get_procname (feature) result (name)
    type(string_t) :: name
    type(string_t), intent(in) :: feature
    name = feature
  end function prc_writer_get_procname

  function prc_writer_get_c_procname (writer, id, feature) result (name)
    class(prc_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id, feature
    type(string_t) :: name
    name = id // "_" // feature
  end function prc_writer_get_c_procname

  function prc_writer_get_module_name (id) result (name)
    type(string_t) :: name
    type(string_t), intent(in) :: id
    name = id
  end function prc_writer_get_module_name

  function prclib_driver_record_get_c_procname (record, feature) result (name)
    type(string_t) :: name
    class(prclib_driver_record_t), intent(in) :: record
    type(string_t), intent(in) :: feature
    name = record%writer%get_c_procname (record%id, feature)
  end function prclib_driver_record_get_c_procname

  subroutine dispatch_prclib_driver &
       (driver, basename, modellibs_ldflags)
    class(prclib_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    type(string_t), intent(in), optional :: modellibs_ldflags
    procedure(dispatch_prclib_driver) :: dispatch_prclib_static
    if (allocated (driver))  deallocate (driver)
    call dispatch_prclib_static (driver, basename)
    if (.not. allocated (driver)) then
       allocate (prclib_driver_dynamic_t :: driver)
    end if
    driver%basename = basename
    driver%modellibs_ldflags = modellibs_ldflags
  end subroutine dispatch_prclib_driver


end module prclib_interfaces
