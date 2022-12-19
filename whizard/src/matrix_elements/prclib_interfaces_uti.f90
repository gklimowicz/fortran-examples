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

module prclib_interfaces_uti

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds
  use system_dependencies, only: CC_HAS_QUADMATH, DEFAULT_FC_PRECISION
  use iso_varying_string, string_t => varying_string
  use io_units
  use system_defs, only: TAB
  use os_interface

  use prclib_interfaces

  implicit none
  private

  public :: test_writer_4_t

  public :: prclib_interfaces_1
  public :: prclib_interfaces_2
  public :: prclib_interfaces_3
  public :: prclib_interfaces_4
  public :: prclib_interfaces_5
  public :: prclib_interfaces_6
  public :: prclib_interfaces_7

  type, extends (prc_writer_t) :: test_writer_1_t
   contains
     procedure, nopass :: type_name => test_writer_1_type_name
     procedure :: write_makefile_code => test_writer_1_mk
     procedure :: write_source_code => test_writer_1_src
     procedure :: write_interface => test_writer_1_if
     procedure :: write_md5sum_call => test_writer_1_md5sum
     procedure :: write_int_sub_call => test_writer_1_int_sub
     procedure :: write_col_state_call => test_writer_1_col_state
     procedure :: write_color_factors_call => test_writer_1_col_factors
     procedure :: before_compile => test_writer_1_before_compile
     procedure :: after_compile => test_writer_1_after_compile
  end type test_writer_1_t

  type, extends (prc_writer_f_module_t) :: test_writer_2_t
   contains
     procedure, nopass :: type_name => test_writer_2_type_name
     procedure :: write_makefile_code => test_writer_2_mk
     procedure :: write_source_code => test_writer_2_src
     procedure :: write_interface => test_writer_2_if
     procedure :: write_wrapper => test_writer_2_wr
     procedure :: before_compile => test_writer_2_before_compile
     procedure :: after_compile => test_writer_2_after_compile
  end type test_writer_2_t

  type, extends (prc_writer_f_module_t) :: test_writer_4_t
   contains
     procedure, nopass :: type_name => test_writer_4_type_name
     procedure, nopass :: get_module_name => &
          test_writer_4_get_module_name
     procedure :: write_makefile_code => test_writer_4_mk
     procedure :: write_source_code => test_writer_4_src
     procedure :: write_interface => test_writer_4_if
     procedure :: write_wrapper => test_writer_4_wr
     procedure :: before_compile => test_writer_4_before_compile
     procedure :: after_compile => test_writer_4_after_compile
  end type test_writer_4_t

  type, extends (prc_writer_c_lib_t) :: test_writer_5_t
   contains
     procedure, nopass :: type_name => test_writer_5_type_name
     procedure :: write_makefile_code => test_writer_5_mk
     procedure :: write_source_code => test_writer_5_src
     procedure :: write_interface => test_writer_5_if
     procedure :: before_compile => test_writer_5_before_compile
     procedure :: after_compile => test_writer_5_after_compile
  end type test_writer_5_t

  type, extends (test_writer_5_t) :: test_writer_6_t
   contains
     procedure, nopass :: type_name => test_writer_6_type_name
     procedure :: write_makefile_code => test_writer_6_mk
     procedure :: write_source_code => test_writer_6_src
  end type test_writer_6_t


contains

  subroutine prclib_interfaces_1 (u)
    integer, intent(in) :: u
    class(prclib_driver_t), allocatable :: driver
    character(32), parameter :: md5sum = "prclib_interfaces_1_md5sum      "
    class(prc_writer_t), pointer :: test_writer_1


    write (u, "(A)")  "* Test output: prclib_interfaces_1"
    write (u, "(A)")  "*   Purpose: display the driver object contents"
    write (u, *)
    write (u, "(A)")  "* Create a prclib driver object"
    write (u, "(A)")

    call dispatch_prclib_driver (driver, var_str ("prclib"), var_str (""))
    call driver%init (3)
    call driver%set_md5sum (md5sum)

    allocate (test_writer_1_t :: test_writer_1)

    call driver%set_record (1, var_str ("test1"), var_str ("test_model"), &
         [var_str ("init")], test_writer_1)

    call driver%set_record (2, var_str ("test2"), var_str ("foo_model"), &
         [var_str ("another_proc")], test_writer_1)

    call driver%set_record (3, var_str ("test3"), var_str ("test_model"), &
         [var_str ("init"), var_str ("some_proc")], test_writer_1)

    call driver%write (u)

    deallocate (test_writer_1)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prclib_interfaces_1"
  end subroutine prclib_interfaces_1

  subroutine prclib_interfaces_2 (u)
    integer, intent(in) :: u
    class(prclib_driver_t), allocatable :: driver
    character(32), parameter :: md5sum = "prclib_interfaces_2_md5sum      "
    class(prc_writer_t), pointer :: test_writer_1, test_writer_2

    write (u, "(A)")  "* Test output: prclib_interfaces_2"
    write (u, "(A)")  "*   Purpose: check the generated driver source code"
    write (u, "(A)")
    write (u, "(A)")  "* Create a prclib driver object (2 processes)"
    write (u, "(A)")

    call dispatch_prclib_driver (driver, var_str ("prclib2"), var_str (""))
    call driver%init (2)
    call driver%set_md5sum (md5sum)

    allocate (test_writer_1_t :: test_writer_1)
    allocate (test_writer_2_t :: test_writer_2)

    call driver%set_record (1, var_str ("test1"), var_str ("Test_model"), &
         [var_str ("proc1")], test_writer_1)

    call driver%set_record (2, var_str ("test2"), var_str ("Test_model"), &
         [var_str ("proc1"), var_str ("proc2")], test_writer_2)

    call driver%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write the driver file"
    write (u, "(A)")  "* File contents:"
    write (u, "(A)")

    call driver%generate_driver_code (u)

    deallocate (test_writer_1)
    deallocate (test_writer_2)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prclib_interfaces_2"
  end subroutine prclib_interfaces_2

  subroutine prclib_interfaces_3 (u)
    integer, intent(in) :: u
    class(prclib_driver_t), allocatable :: driver
    type(os_data_t) :: os_data
    character(32), parameter :: md5sum = "prclib_interfaces_3_md5sum      "
    class(prc_writer_t), pointer :: test_writer_1, test_writer_2

    call os_data%init ()
    os_data%fc = "fortran-compiler"
    os_data%whizard_includes = "-I module-dir"
    os_data%fcflags = "-C=all"
    os_data%fcflags_pic = "-PIC"
    os_data%cc = "c-compiler"
    os_data%cflags = "-I include-dir"
    os_data%cflags_pic = "-PIC"
    os_data%whizard_ldflags = ""
    os_data%ldflags = ""
    os_data%whizard_libtool = "my-libtool"
    os_data%latex = "latex -halt-on-error"
    os_data%mpost = "mpost --math=scaled -halt-on-error"
    os_data%dvips = "dvips"
    os_data%ps2pdf = "ps2pdf14"
    os_data%whizard_texpath = ""

    write (u, "(A)")  "* Test output: prclib_interfaces_3"
    write (u, "(A)")  "*   Purpose: check the generated Makefile"
    write (u, *)
    write (u, "(A)")  "* Create a prclib driver object (2 processes)"
    write (u, "(A)")

    call dispatch_prclib_driver (driver, var_str ("prclib3"), var_str (""))
    call driver%init (2)
    call driver%set_md5sum (md5sum)

    allocate (test_writer_1_t :: test_writer_1)
    allocate (test_writer_2_t :: test_writer_2)

    call driver%set_record (1, var_str ("test1"), var_str ("Test_model"), &
         [var_str ("proc1")], test_writer_1)

    call driver%set_record (2, var_str ("test2"), var_str ("Test_model"), &
         [var_str ("proc1"), var_str ("proc2")], test_writer_2)

    call driver%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write Makefile"
    write (u, "(A)")  "* File contents:"
    write (u, "(A)")

    call driver%generate_makefile (u, os_data, verbose = .true.)

    deallocate (test_writer_1)
    deallocate (test_writer_2)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prclib_interfaces_3"
  end subroutine prclib_interfaces_3

  subroutine prclib_interfaces_4 (u)
    integer, intent(in) :: u
    class(prclib_driver_t), allocatable :: driver
    class(prc_writer_t), pointer :: test_writer_4
    type(os_data_t) :: os_data
    integer :: u_file

    integer, dimension(:,:), allocatable :: flv_state
    integer, dimension(:,:), allocatable :: hel_state
    integer, dimension(:,:,:), allocatable :: col_state
    logical, dimension(:,:), allocatable :: ghost_flag
    integer, dimension(:,:), allocatable :: cf_index
    complex(default), dimension(:), allocatable :: color_factors
    character(32), parameter :: md5sum = "prclib_interfaces_4_md5sum      "
    character(32) :: md5sum_file

    type(c_funptr) :: proc1_ptr
    interface
       subroutine proc1_t (n) bind(C)
         import
         integer(c_int), intent(out) :: n
       end subroutine proc1_t
    end interface
    procedure(proc1_t), pointer :: proc1
    integer(c_int) :: n

    write (u, "(A)")  "* Test output: prclib_interfaces_4"
    write (u, "(A)")  "*   Purpose: compile, link, and load process library"
    write (u, "(A)")  "*            with (fake) matrix-element code &
         &as a Fortran module"
    write (u, *)
    write (u, "(A)")  "* Create a prclib driver object (1 process)"
    write (u, "(A)")

    call os_data%init ()

    allocate (test_writer_4_t :: test_writer_4)
    call test_writer_4%init_test ()

    call dispatch_prclib_driver (driver, var_str ("prclib4"), var_str (""))
    call driver%init (1)
    call driver%set_md5sum (md5sum)

    call driver%set_record (1, var_str ("test4"), var_str ("Test_model"), &
         [var_str ("proc1")], test_writer_4)

    call driver%write (u)

    write (u, *)
    write (u, "(A)")  "* Write Makefile"
    u_file = free_unit ()
    open (u_file, file="prclib4.makefile", status="replace", action="write")
    call driver%generate_makefile (u_file, os_data, verbose = .false.)
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Recall MD5 sum from Makefile"
    write (u, "(A)")

    md5sum_file = driver%get_md5sum_makefile ()
    write (u, "(1x,A,A,A)")  "MD5 sum = '", md5sum_file, "'"

    write (u, "(A)")
    write (u, "(A)")  "* Write driver source code"

    u_file = free_unit ()
    open (u_file, file="prclib4.f90", status="replace", action="write")
    call driver%generate_driver_code (u_file)
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Recall MD5 sum from driver source"
    write (u, "(A)")

    md5sum_file = driver%get_md5sum_driver ()
    write (u, "(1x,A,A,A)")  "MD5 sum = '", md5sum_file, "'"

    write (u, "(A)")
    write (u, "(A)")  "* Write matrix-element source code"
    call driver%make_source (os_data)

    write (u, "(A)")
    write (u, "(A)")  "* Recall MD5 sum from matrix-element source"
    write (u, "(A)")

    md5sum_file = driver%get_md5sum_source (1)
    write (u, "(1x,A,A,A)")  "MD5 sum = '", md5sum_file, "'"

    write (u, "(A)")
    write (u, "(A)")  "* Compile source code"
    call driver%make_compile (os_data)

    write (u, "(A)")  "* Link library"
    call driver%make_link (os_data)

    write (u, "(A)")  "* Load library"
    call driver%load (os_data)

    write (u, *)
    call driver%write (u)
    write (u, *)

    if (driver%loaded) then
       write (u, "(A)")  "* Call library functions:"
       write (u, *)
       write (u, "(1x,A,I0)")  "n_processes   = ", driver%get_n_processes ()
       write (u, "(1x,A,A,A)")  "process_id    = '", &
            char (driver%get_process_id (1)), "'"
       write (u, "(1x,A,A,A)")  "model_name    = '", &
            char (driver%get_model_name (1)), "'"
       write (u, "(1x,A,A,A)")  "md5sum (lib)  = '", &
            char (driver%get_md5sum (0)), "'"
       write (u, "(1x,A,A,A)")  "md5sum (proc) = '", &
            char (driver%get_md5sum (1)), "'"
       write (u, "(1x,A,L1)")  "openmp_status = ", driver%get_openmp_status (1)
       write (u, "(1x,A,I0)")  "n_in  = ", driver%get_n_in (1)
       write (u, "(1x,A,I0)")  "n_out = ", driver%get_n_out (1)
       write (u, "(1x,A,I0)")  "n_flv = ", driver%get_n_flv (1)
       write (u, "(1x,A,I0)")  "n_hel = ", driver%get_n_hel (1)
       write (u, "(1x,A,I0)")  "n_col = ", driver%get_n_col (1)
       write (u, "(1x,A,I0)")  "n_cin = ", driver%get_n_cin (1)
       write (u, "(1x,A,I0)")  "n_cf  = ", driver%get_n_cf (1)

       call driver%set_flv_state (1, flv_state)
       write (u, "(1x,A,10(1x,I0))")  "flv_state =", flv_state

       call driver%set_hel_state (1, hel_state)
       write (u, "(1x,A,10(1x,I0))")  "hel_state =", hel_state

       call driver%set_col_state (1, col_state, ghost_flag)
       write (u, "(1x,A,10(1x,I0))")  "col_state =", col_state
       write (u, "(1x,A,10(1x,L1))")  "ghost_flag =", ghost_flag

       call driver%set_color_factors (1, color_factors, cf_index)
       write (u, "(1x,A,10(1x,F5.3))")  "color_factors =", color_factors
       write (u, "(1x,A,10(1x,I0))")  "cf_index =", cf_index

       call driver%get_fptr (1, 1, proc1_ptr)
       call c_f_procpointer (proc1_ptr, proc1)
       if (associated (proc1)) then
          write (u, *)
          call proc1 (n)
          write (u, "(1x,A,I0)")  "proc1(1) = ", n
       end if

    end if

    deallocate (test_writer_4)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prclib_interfaces_4"
  end subroutine prclib_interfaces_4

  subroutine prclib_interfaces_5 (u)
    integer, intent(in) :: u
    class(prclib_driver_t), allocatable :: driver
    class(prc_writer_t), pointer :: test_writer_5
    type(os_data_t) :: os_data
    integer :: u_file

    integer, dimension(:,:), allocatable :: flv_state
    integer, dimension(:,:), allocatable :: hel_state
    integer, dimension(:,:,:), allocatable :: col_state
    logical, dimension(:,:), allocatable :: ghost_flag
    integer, dimension(:,:), allocatable :: cf_index
    complex(default), dimension(:), allocatable :: color_factors
    character(32), parameter :: md5sum = "prclib_interfaces_5_md5sum      "

    type(c_funptr) :: proc1_ptr
    interface
       subroutine proc1_t (n) bind(C)
         import
         integer(c_int), intent(out) :: n
       end subroutine proc1_t
    end interface
    procedure(proc1_t), pointer :: proc1
    integer(c_int) :: n

    write (u, "(A)")  "* Test output: prclib_interfaces_5"
    write (u, "(A)")  "*   Purpose: compile, link, and load process library"
    write (u, "(A)")  "*            with (fake) matrix-element code &
         &as a Fortran bind(C) library"
    write (u, *)
    write (u, "(A)")  "* Create a prclib driver object (1 process)"
    write (u, "(A)")

    call os_data%init ()
    allocate (test_writer_5_t :: test_writer_5)

    call dispatch_prclib_driver (driver, var_str ("prclib5"), var_str (""))
    call driver%init (1)
    call driver%set_md5sum (md5sum)

    call driver%set_record (1, var_str ("test5"), var_str ("Test_model"), &
         [var_str ("proc1")], test_writer_5)

    call driver%write (u)

    write (u, *)
    write (u, "(A)")  "* Write makefile"
    u_file = free_unit ()
    open (u_file, file="prclib5.makefile", status="replace", action="write")
    call driver%generate_makefile (u_file, os_data, verbose = .false.)
    close (u_file)

    write (u, "(A)")  "* Write driver source code"
    u_file = free_unit ()
    open (u_file, file="prclib5.f90", status="replace", action="write")
    call driver%generate_driver_code (u_file)
    close (u_file)

    write (u, "(A)")  "* Write matrix-element source code"
    call driver%make_source (os_data)

    write (u, "(A)")  "* Compile source code"
    call driver%make_compile (os_data)

    write (u, "(A)")  "* Link library"
    call driver%make_link (os_data)

    write (u, "(A)")  "* Load library"
    call driver%load (os_data)

    write (u, *)
    call driver%write (u)
    write (u, *)

    if (driver%loaded) then
       write (u, "(A)")  "* Call library functions:"
       write (u, *)
       write (u, "(1x,A,I0)")  "n_processes   = ", driver%get_n_processes ()
       write (u, "(1x,A,A)")  "process_id    = ", &
            char (driver%get_process_id (1))
       write (u, "(1x,A,A)")  "model_name    = ", &
            char (driver%get_model_name (1))
       write (u, "(1x,A,A)")  "md5sum        = ", &
            char (driver%get_md5sum (1))
       write (u, "(1x,A,L1)")  "openmp_status = ", driver%get_openmp_status (1)
       write (u, "(1x,A,I0)")  "n_in  = ", driver%get_n_in (1)
       write (u, "(1x,A,I0)")  "n_out = ", driver%get_n_out (1)
       write (u, "(1x,A,I0)")  "n_flv = ", driver%get_n_flv (1)
       write (u, "(1x,A,I0)")  "n_hel = ", driver%get_n_hel (1)
       write (u, "(1x,A,I0)")  "n_col = ", driver%get_n_col (1)
       write (u, "(1x,A,I0)")  "n_cin = ", driver%get_n_cin (1)
       write (u, "(1x,A,I0)")  "n_cf  = ", driver%get_n_cf (1)

       call driver%set_flv_state (1, flv_state)
       write (u, "(1x,A,10(1x,I0))")  "flv_state =", flv_state

       call driver%set_hel_state (1, hel_state)
       write (u, "(1x,A,10(1x,I0))")  "hel_state =", hel_state

       call driver%set_col_state (1, col_state, ghost_flag)
       write (u, "(1x,A,10(1x,I0))")  "col_state =", col_state
       write (u, "(1x,A,10(1x,L1))")  "ghost_flag =", ghost_flag

       call driver%set_color_factors (1, color_factors, cf_index)
       write (u, "(1x,A,10(1x,F5.3))")  "color_factors =", color_factors
       write (u, "(1x,A,10(1x,I0))")  "cf_index =", cf_index

       call driver%get_fptr (1, 1, proc1_ptr)
       call c_f_procpointer (proc1_ptr, proc1)
       if (associated (proc1)) then
          write (u, *)
          call proc1 (n)
          write (u, "(1x,A,I0)")  "proc1(1) = ", n
       end if

    end if

    deallocate (test_writer_5)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prclib_interfaces_5"
  end subroutine prclib_interfaces_5

  subroutine prclib_interfaces_6 (u)
    integer, intent(in) :: u
    class(prclib_driver_t), allocatable :: driver
    class(prc_writer_t), pointer :: test_writer_6
    type(os_data_t) :: os_data
    integer :: u_file

    integer, dimension(:,:), allocatable :: flv_state
    integer, dimension(:,:), allocatable :: hel_state
    integer, dimension(:,:,:), allocatable :: col_state
    logical, dimension(:,:), allocatable :: ghost_flag
    integer, dimension(:,:), allocatable :: cf_index
    complex(default), dimension(:), allocatable :: color_factors
    character(32), parameter :: md5sum = "prclib_interfaces_6_md5sum      "

    type(c_funptr) :: proc1_ptr
    interface
       subroutine proc1_t (n) bind(C)
         import
         integer(c_int), intent(out) :: n
       end subroutine proc1_t
    end interface
    procedure(proc1_t), pointer :: proc1
    integer(c_int) :: n

    write (u, "(A)")  "* Test output: prclib_interfaces_6"
    write (u, "(A)")  "*   Purpose: compile, link, and load process library"
    write (u, "(A)")  "*            with (fake) matrix-element code &
         &as a C library"
    write (u, *)
    write (u, "(A)")  "* Create a prclib driver object (1 process)"
    write (u, "(A)")

    call os_data%init ()
    allocate (test_writer_6_t :: test_writer_6)

    call dispatch_prclib_driver (driver, var_str ("prclib6"), var_str (""))
    call driver%init (1)
    call driver%set_md5sum (md5sum)

    call driver%set_record (1, var_str ("test6"), var_str ("Test_model"), &
         [var_str ("proc1")], test_writer_6)

    call driver%write (u)

    write (u, *)
    write (u, "(A)")  "* Write makefile"
    u_file = free_unit ()
    open (u_file, file="prclib6.makefile", status="replace", action="write")
    call driver%generate_makefile (u_file, os_data, verbose = .false.)
    close (u_file)

    write (u, "(A)")  "* Write driver source code"
    u_file = free_unit ()
    open (u_file, file="prclib6.f90", status="replace", action="write")
    call driver%generate_driver_code (u_file)
    close (u_file)

    write (u, "(A)")  "* Write matrix-element source code"
    call driver%make_source (os_data)

    write (u, "(A)")  "* Compile source code"
    call driver%make_compile (os_data)

    write (u, "(A)")  "* Link library"
    call driver%make_link (os_data)

    write (u, "(A)")  "* Load library"
    call driver%load (os_data)

    write (u, *)
    call driver%write (u)
    write (u, *)

    if (driver%loaded) then
       write (u, "(A)")  "* Call library functions:"
       write (u, *)
       write (u, "(1x,A,I0)")  "n_processes   = ", driver%get_n_processes ()
       write (u, "(1x,A,A)")  "process_id    = ", &
            char (driver%get_process_id (1))
       write (u, "(1x,A,A)")  "model_name    = ", &
            char (driver%get_model_name (1))
       write (u, "(1x,A,A)")  "md5sum        = ", &
            char (driver%get_md5sum (1))
       write (u, "(1x,A,L1)")  "openmp_status = ", driver%get_openmp_status (1)
       write (u, "(1x,A,I0)")  "n_in  = ", driver%get_n_in (1)
       write (u, "(1x,A,I0)")  "n_out = ", driver%get_n_out (1)
       write (u, "(1x,A,I0)")  "n_flv = ", driver%get_n_flv (1)
       write (u, "(1x,A,I0)")  "n_hel = ", driver%get_n_hel (1)
       write (u, "(1x,A,I0)")  "n_col = ", driver%get_n_col (1)
       write (u, "(1x,A,I0)")  "n_cin = ", driver%get_n_cin (1)
       write (u, "(1x,A,I0)")  "n_cf  = ", driver%get_n_cf (1)

       call driver%set_flv_state (1, flv_state)
       write (u, "(1x,A,10(1x,I0))")  "flv_state =", flv_state

       call driver%set_hel_state (1, hel_state)
       write (u, "(1x,A,10(1x,I0))")  "hel_state =", hel_state

       call driver%set_col_state (1, col_state, ghost_flag)
       write (u, "(1x,A,10(1x,I0))")  "col_state =", col_state
       write (u, "(1x,A,10(1x,L1))")  "ghost_flag =", ghost_flag

       call driver%set_color_factors (1, color_factors, cf_index)
       write (u, "(1x,A,10(1x,F5.3))")  "color_factors =", color_factors
       write (u, "(1x,A,10(1x,I0))")  "cf_index =", cf_index

       call driver%get_fptr (1, 1, proc1_ptr)
       call c_f_procpointer (proc1_ptr, proc1)
       if (associated (proc1)) then
          write (u, *)
          call proc1 (n)
          write (u, "(1x,A,I0)")  "proc1(1) = ", n
       end if

    end if

    deallocate (test_writer_6)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prclib_interfaces_6"
  end subroutine prclib_interfaces_6

  subroutine prclib_interfaces_7 (u)
    integer, intent(in) :: u
    class(prclib_driver_t), allocatable :: driver
    class(prc_writer_t), pointer :: test_writer_4
    type(os_data_t) :: os_data
    integer :: u_file
    character(32), parameter :: md5sum = "1234567890abcdef1234567890abcdef"

    write (u, "(A)")  "* Test output: prclib_interfaces_7"
    write (u, "(A)")  "*   Purpose: compile and link process library"
    write (u, "(A)")  "*            with (fake) matrix-element code &
         &as a Fortran module"
    write (u, "(A)")  "*            then clean up generated files"
    write (u, *)
    write (u, "(A)")  "* Create a prclib driver object (1 process)"

    allocate (test_writer_4_t :: test_writer_4)

    call os_data%init ()
    call dispatch_prclib_driver (driver, var_str ("prclib7"), var_str (""))
    call driver%init (1)
    call driver%set_md5sum (md5sum)
    call driver%set_record (1, var_str ("test7"), var_str ("Test_model"), &
         [var_str ("proc1")], test_writer_4)

    write (u, "(A)")  "* Write makefile"
    u_file = free_unit ()
    open (u_file, file="prclib7.makefile", status="replace", action="write")
    call driver%generate_makefile (u_file, os_data, verbose = .false.)
    close (u_file)

    write (u, "(A)")  "* Write driver source code"
    u_file = free_unit ()
    open (u_file, file="prclib7.f90", status="replace", action="write")
    call driver%generate_driver_code (u_file)
    close (u_file)

    write (u, "(A)")  "* Write matrix-element source code"
    call driver%make_source (os_data)

    write (u, "(A)")  "* Compile source code"
    call driver%make_compile (os_data)

    write (u, "(A)")  "* Link library"
    call driver%make_link (os_data)


    write (u, "(A)")  "* File check"
    write (u, *)
    call check_file (u, "test7.f90")
    call check_file (u, "tpr_test7.mod")
    call check_file (u, "test7.lo")
    call check_file (u, "prclib7.makefile")
    call check_file (u, "prclib7.f90")
    call check_file (u, "prclib7.lo")
    call check_file (u, "prclib7.la")

    write (u, *)
    write (u, "(A)")  "* Delete library"
    write (u, *)
    call driver%clean_library (os_data)
    call check_file (u, "prclib7.la")

    write (u, *)
    write (u, "(A)")  "* Delete object code"
    write (u, *)
    call driver%clean_objects (os_data)
    call check_file (u, "test7.lo")
    call check_file (u, "tpr_test7.mod")
    call check_file (u, "prclib7.lo")

    write (u, *)
    write (u, "(A)")  "* Delete source code"
    write (u, *)
    call driver%clean_source (os_data)
    call check_file (u, "test7.f90")

    write (u, *)
    write (u, "(A)")  "* Delete driver source code"
    write (u, *)
    call driver%clean_driver (os_data)
    call check_file (u, "prclib7.f90")

    write (u, *)
    write (u, "(A)")  "* Delete makefile"
    write (u, *)
    call driver%clean_makefile (os_data)
    call check_file (u, "prclib7.makefile")

    deallocate (test_writer_4)

    write (u, *)
    write (u, "(A)")  "* Test output end: prclib_interfaces_7"
  end subroutine prclib_interfaces_7


  function test_writer_1_type_name () result (string)
    type(string_t) :: string
    string = "test_1"
  end function test_writer_1_type_name

  subroutine test_writer_1_mk &
       (writer, unit, id, os_data, verbose, testflag)
    class(test_writer_1_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    write (unit, "(5A)")  "# Makefile code for process ", char (id), &
         " goes here."
  end subroutine test_writer_1_mk

  subroutine test_writer_1_src (writer, id)
    class(test_writer_1_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_1_src

  subroutine test_writer_1_if (writer, unit, id, feature)
    class(test_writer_1_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    write (unit, "(2x,9A)")  "! Interface code for ", &
       char (id), "_", char (writer%get_procname (feature)), &
       " goes here."
  end subroutine test_writer_1_if

  subroutine test_writer_1_md5sum (writer, unit, id)
    class(test_writer_1_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(5x,9A)")  "! MD5sum call for ", char (id), " goes here."
  end subroutine test_writer_1_md5sum

  subroutine test_writer_1_int_sub (writer, unit, id, feature)
    class(test_writer_1_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    write (unit, "(5x,9A)")  "! ", char (feature), " call for ", &
         char (id), " goes here."
  end subroutine test_writer_1_int_sub

  subroutine test_writer_1_col_state (writer, unit, id)
    class(test_writer_1_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(5x,9A)")  "! col_state call for ", &
         char (id), " goes here."
  end subroutine test_writer_1_col_state

  subroutine test_writer_1_col_factors (writer, unit, id)
    class(test_writer_1_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    write (unit, "(5x,9A)")  "! color_factors call for ", &
         char (id), " goes here."
  end subroutine test_writer_1_col_factors

  subroutine test_writer_1_before_compile (writer, id)
    class(test_writer_1_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_1_before_compile

  subroutine test_writer_1_after_compile (writer, id)
    class(test_writer_1_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_1_after_compile

  function test_writer_2_type_name () result (string)
    type(string_t) :: string
    string = "test_2"
  end function test_writer_2_type_name

  subroutine test_writer_2_mk &
       (writer, unit, id, os_data, verbose, testflag)
    class(test_writer_2_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    write (unit, "(5A)")  "# Makefile code for process ", char (id), &
         " goes here."
  end subroutine test_writer_2_mk

  subroutine test_writer_2_src (writer, id)
    class(test_writer_2_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_2_src

  subroutine test_writer_2_if (writer, unit, id, feature)
    class(test_writer_2_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    write (unit, "(2x,9A)")  "! Interface code for ", &
       char (writer%get_module_name (id)), "_", &
       char (writer%get_procname (feature)), " goes here."
  end subroutine test_writer_2_if

  subroutine test_writer_2_wr (writer, unit, id, feature)
    class(test_writer_2_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    write (unit, *)
    write (unit, "(9A)")  "! Wrapper code for ", &
       char (writer%get_c_procname (id, feature)), " goes here."
  end subroutine test_writer_2_wr

  subroutine test_writer_2_before_compile (writer, id)
    class(test_writer_2_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_2_before_compile

  subroutine test_writer_2_after_compile (writer, id)
    class(test_writer_2_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_2_after_compile

  function test_writer_4_type_name () result (string)
    type(string_t) :: string
    string = "test_4"
  end function test_writer_4_type_name

  function test_writer_4_get_module_name (id) result (name)
    type(string_t), intent(in) :: id
    type(string_t) :: name
    name = "tpr_" // id
  end function test_writer_4_get_module_name

  subroutine test_writer_4_mk &
       (writer, unit, id, os_data, verbose, testflag)
    class(test_writer_4_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    write (unit, "(5A)")  "SOURCES += ", char (id), ".f90"
    write (unit, "(5A)")  "OBJECTS += ", char (id), ".lo"
    write (unit, "(5A)")  "CLEAN_SOURCES += ", char (id), ".f90"
    write (unit, "(5A)")  "CLEAN_OBJECTS += tpr_", char (id), ".mod"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), ".lo"
    write (unit, "(5A)")  char (id), ".lo: ", char (id), ".f90"
    if (.not. verbose) then
       write (unit, "(5A)")  TAB // '@echo  "  FC       " $@'
    end if
    write (unit, "(5A)")  TAB, "$(LTFCOMPILE) $<"
  end subroutine test_writer_4_mk

  subroutine test_writer_4_src (writer, id)
    class(test_writer_4_t), intent(in) :: writer
    type(string_t), intent(in) :: id
    call write_test_module_file (id, var_str ("proc1"), writer%md5sum)
  end subroutine test_writer_4_src

  subroutine test_writer_4_if (writer, unit, id, feature)
    class(test_writer_4_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    write (unit, "(2x,9A)")  "interface"
    write (unit, "(5x,9A)")  "subroutine ", &
       char (writer%get_c_procname (id, feature)), &
       " (n) bind(C)"
    write (unit, "(7x,9A)")  "import"
    write (unit, "(7x,9A)")  "implicit none"
    write (unit, "(7x,9A)")  "integer(c_int), intent(out) :: n"
    write (unit, "(5x,9A)")  "end subroutine ", &
       char (writer%get_c_procname (id, feature))
    write (unit, "(2x,9A)")  "end interface"
  end subroutine test_writer_4_if

  subroutine test_writer_4_wr (writer, unit, id, feature)
    class(test_writer_4_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    write (unit, *)
    write (unit, "(9A)")  "subroutine ", &
         char (writer%get_c_procname (id, feature)), &
       " (n) bind(C)"
    write (unit, "(2x,9A)")  "use iso_c_binding"
    write (unit, "(2x,9A)")  "use tpr_", char (id), ", only: ", &
         char (writer%get_procname (feature))
    write (unit, "(2x,9A)")  "implicit none"
    write (unit, "(2x,9A)")  "integer(c_int), intent(out) :: n"
    write (unit, "(2x,9A)")  "call ", char (feature), " (n)"
    write (unit, "(9A)")  "end subroutine ", &
       char (writer%get_c_procname (id, feature))
  end subroutine test_writer_4_wr

  subroutine test_writer_4_before_compile (writer, id)
    class(test_writer_4_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_4_before_compile

  subroutine test_writer_4_after_compile (writer, id)
    class(test_writer_4_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_4_after_compile

  subroutine write_test_module_file (basename, feature, md5sum)
    type(string_t), intent(in) :: basename
    type(string_t), intent(in) :: feature
    character(32), intent(in) :: md5sum
    integer :: u
    u = free_unit ()
    open (u, file = char (basename) // ".f90", &
         status = "replace", action = "write")
    write (u, "(A)")  "! (Pseudo) matrix element code file &
         &for WHIZARD self-test"
    write (u, *)
    write (u, "(A)")  "module tpr_" // char (basename)
    write (u, *)
    write (u, "(2x,A)")  "use kinds"
    write (u, "(2x,A)")  "use omega_color, OCF => omega_color_factor"
    write (u, *)
    write (u, "(2x,A)")  "implicit none"
    write (u, "(2x,A)")  "private"
    write (u, *)
    call write_test_me_code_1 (u)
    write (u, *)
    write (u, "(2x,A)")  "public :: " // char (feature)
    write (u, *)
    write (u, "(A)")  "contains"
    write (u, *)
    call write_test_me_code_2 (u, md5sum)
    write (u, *)
    write (u, "(2x,A)")  "subroutine " // char (feature) // " (n)"
    write (u, "(2x,A)")  "  integer, intent(out) :: n"
    write (u, "(2x,A)")  "  n = 42"
    write (u, "(2x,A)")  "end subroutine " // char (feature)
    write (u, *)
    write (u, "(A)")  "end module tpr_" // char (basename)
    close (u)
  end subroutine write_test_module_file

  subroutine write_test_me_code_1 (u)
    integer, intent(in) :: u
    write (u, "(2x,A)")  "public :: md5sum"
    write (u, "(2x,A)")  "public :: openmp_supported"
    write (u, *)
    write (u, "(2x,A)")  "public :: n_in"
    write (u, "(2x,A)")  "public :: n_out"
    write (u, "(2x,A)")  "public :: n_flv"
    write (u, "(2x,A)")  "public :: n_hel"
    write (u, "(2x,A)")  "public :: n_cin"
    write (u, "(2x,A)")  "public :: n_col"
    write (u, "(2x,A)")  "public :: n_cf"
    write (u, *)
    write (u, "(2x,A)")  "public :: flv_state"
    write (u, "(2x,A)")  "public :: hel_state"
    write (u, "(2x,A)")  "public :: col_state"
    write (u, "(2x,A)")  "public :: color_factors"
  end subroutine write_test_me_code_1

  subroutine write_test_me_code_2 (u, md5sum)
    integer, intent(in) :: u
    character(32), intent(in) :: md5sum
    write (u, "(2x,A)")  "pure function md5sum ()"
    write (u, "(2x,A)")  "  character(len=32) :: md5sum"
    write (u, "(2x,A)")  "  md5sum = '" // md5sum // "'"
    write (u, "(2x,A)")  "end function md5sum"
    write (u, *)
    write (u, "(2x,A)")  "pure function openmp_supported () result (status)"
    write (u, "(2x,A)")  "  logical :: status"
    write (u, "(2x,A)")  "  status = .false."
    write (u, "(2x,A)")  "end function openmp_supported"
    write (u, *)
    write (u, "(2x,A)")  "pure function n_in () result (n)"
    write (u, "(2x,A)")  "  integer :: n"
    write (u, "(2x,A)")  "  n = 1"
    write (u, "(2x,A)")  "end function n_in"
    write (u, *)
    write (u, "(2x,A)")  "pure function n_out () result (n)"
    write (u, "(2x,A)")  "  integer :: n"
    write (u, "(2x,A)")  "  n = 2"
    write (u, "(2x,A)")  "end function n_out"
    write (u, *)
    write (u, "(2x,A)")  "pure function n_flv () result (n)"
    write (u, "(2x,A)")  "  integer :: n"
    write (u, "(2x,A)")  "  n = 1"
    write (u, "(2x,A)")  "end function n_flv"
    write (u, *)
    write (u, "(2x,A)")  "pure function n_hel () result (n)"
    write (u, "(2x,A)")  "  integer :: n"
    write (u, "(2x,A)")  "  n = 1"
    write (u, "(2x,A)")  "end function n_hel"
    write (u, *)
    write (u, "(2x,A)")  "pure function n_cin () result (n)"
    write (u, "(2x,A)")  "  integer :: n"
    write (u, "(2x,A)")  "  n = 2"
    write (u, "(2x,A)")  "end function n_cin"
    write (u, *)
    write (u, "(2x,A)")  "pure function n_col () result (n)"
    write (u, "(2x,A)")  "  integer :: n"
    write (u, "(2x,A)")  "  n = 1"
    write (u, "(2x,A)")  "end function n_col"
    write (u, *)
    write (u, "(2x,A)")  "pure function n_cf () result (n)"
    write (u, "(2x,A)")  "  integer :: n"
    write (u, "(2x,A)")  "  n = 1"
    write (u, "(2x,A)")  "end function n_cf"
    write (u, *)
    write (u, "(2x,A)")  "pure subroutine flv_state (a)"
    write (u, "(2x,A)")  "  integer, dimension(:,:), intent(out) :: a"
    write (u, "(2x,A)")  "  a = reshape ([1,2,3], [3,1])"
    write (u, "(2x,A)")  "end subroutine flv_state"
    write (u, *)
    write (u, "(2x,A)")  "pure subroutine hel_state (a)"
    write (u, "(2x,A)")  "  integer, dimension(:,:), intent(out) :: a"
    write (u, "(2x,A)")  "  a = reshape ([0,0,0], [3,1])"
    write (u, "(2x,A)")  "end subroutine hel_state"
    write (u, *)
    write (u, "(2x,A)")  "pure subroutine col_state (a, g)"
    write (u, "(2x,A)")  "  integer, dimension(:,:,:), intent(out) :: a"
    write (u, "(2x,A)")  "  logical, dimension(:,:), intent(out) :: g"
    write (u, "(2x,A)")  "  a = reshape ([0,0, 0,0, 0,0], [2,3,1])"
    write (u, "(2x,A)")  "  g = reshape ([.false., .false., .false.], [3,1])"
    write (u, "(2x,A)")  "end subroutine col_state"
    write (u, *)
    write (u, "(2x,A)")  "pure subroutine color_factors (cf)"
    write (u, "(2x,A)")  "  type(OCF), dimension(:), intent(out) :: cf"
    write (u, "(2x,A)")  "  cf = [ OCF(1,1,+1._default) ]"
    write (u, "(2x,A)")  "end subroutine color_factors"
  end subroutine write_test_me_code_2

  function test_writer_5_type_name () result (string)
    type(string_t) :: string
    string = "test_5"
  end function test_writer_5_type_name

  subroutine test_writer_5_mk &
       (writer, unit, id, os_data, verbose, testflag)
    class(test_writer_5_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    write (unit, "(5A)")  "SOURCES += ", char (id), ".f90"
    write (unit, "(5A)")  "OBJECTS += ", char (id), ".lo"
    write (unit, "(5A)")  char (id), ".lo: ", char (id), ".f90"
    if (.not. verbose) then
       write (unit, "(5A)")  TAB // '@echo  "  FC       " $@'
    end if
    write (unit, "(5A)")  TAB, "$(LTFCOMPILE) $<"
  end subroutine test_writer_5_mk

  subroutine test_writer_5_src (writer, id)
    class(test_writer_5_t), intent(in) :: writer
    type(string_t), intent(in) :: id
    call write_test_f_lib_file (id, var_str ("proc1"))
  end subroutine test_writer_5_src

  subroutine test_writer_5_if (writer, unit, id, feature)
    class(test_writer_5_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    select case (char (feature))
    case ("proc1")
       write (unit, "(2x,9A)")  "interface"
       write (unit, "(5x,9A)")  "subroutine ", &
            char (writer%get_c_procname (id, feature)), &
            " (n) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "implicit none"
       write (unit, "(7x,9A)")  "integer(c_int), intent(out) :: n"
       write (unit, "(5x,9A)")  "end subroutine ", &
            char (writer%get_c_procname (id, feature))
       write (unit, "(2x,9A)")  "end interface"
    case default
       call writer%write_standard_interface (unit, id, feature)
    end select
  end subroutine test_writer_5_if

  subroutine test_writer_5_before_compile (writer, id)
    class(test_writer_5_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_5_before_compile

  subroutine test_writer_5_after_compile (writer, id)
    class(test_writer_5_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine test_writer_5_after_compile

  subroutine write_test_f_lib_file (basename, feature)
    type(string_t), intent(in) :: basename
    type(string_t), intent(in) :: feature
    integer :: u
    u = free_unit ()
    open (u, file = char (basename) // ".f90", &
         status = "replace", action = "write")
    write (u, "(A)")  "! (Pseudo) matrix element code file &
         &for WHIZARD self-test"
    call write_test_me_code_3 (u, char (basename))
    write (u, *)
    write (u, "(A)")  "subroutine " // char (basename) // "_" &
         // char (feature) // " (n) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int), intent(out) :: n"
    write (u, "(A)")  "  n = 42"
    write (u, "(A)")  "end subroutine " // char (basename) // "_" &
         // char (feature)
    close (u)
  end subroutine write_test_f_lib_file

  subroutine write_test_me_code_3 (u, id)
    integer, intent(in) :: u
    character(*), intent(in) :: id
    write (u, "(A)")  "function " // id // "_get_md5sum () &
         &result (cptr) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  type(c_ptr) :: cptr"
    write (u, "(A)")  "  character(c_char), dimension(32), &
         &target, save :: md5sum"
    write (u, "(A)")  "  md5sum = copy (c_char_&
         &'1234567890abcdef1234567890abcdef')"
    write (u, "(A)")  "  cptr = c_loc (md5sum)"
    write (u, "(A)")  "contains"
    write (u, "(A)")  "  function copy (md5sum)"
    write (u, "(A)")  "    character(c_char), dimension(32) :: copy"
    write (u, "(A)")  "    character(c_char), dimension(32), intent(in) :: &
         &md5sum"
    write (u, "(A)")  "    copy = md5sum"
    write (u, "(A)")  "  end function copy"
    write (u, "(A)")  "end function " // id // "_get_md5sum"
    write (u, *)
    write (u, "(A)")  "function " // id // "_openmp_supported () &
         &result (status) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  logical(c_bool) :: status"
    write (u, "(A)")  "  status = .false."
    write (u, "(A)")  "end function " // id // "_openmp_supported"
    write (u, *)
    write (u, "(A)")  "function " // id // "_n_in () result (n) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int) :: n"
    write (u, "(A)")  "  n = 1"
    write (u, "(A)")  "end function " // id // "_n_in"
    write (u, *)
    write (u, "(A)")  "function " // id // "_n_out () result (n) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int) :: n"
    write (u, "(A)")  "  n = 2"
    write (u, "(A)")  "end function " // id // "_n_out"
    write (u, *)
    write (u, "(A)")  "function " // id // "_n_flv () result (n) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int) :: n"
    write (u, "(A)")  "  n = 1"
    write (u, "(A)")  "end function " // id // "_n_flv"
    write (u, *)
    write (u, "(A)")  "function " // id // "_n_hel () result (n) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int) :: n"
    write (u, "(A)")  "  n = 1"
    write (u, "(A)")  "end function " // id // "_n_hel"
    write (u, *)
    write (u, "(A)")  "function " // id // "_n_cin () result (n) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int) :: n"
    write (u, "(A)")  "  n = 2"
    write (u, "(A)")  "end function " // id // "_n_cin"
    write (u, *)
    write (u, "(A)")  "function " // id // "_n_col () result (n) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int) :: n"
    write (u, "(A)")  "  n = 1"
    write (u, "(A)")  "end function " // id // "_n_col"
    write (u, *)
    write (u, "(A)")  "function " // id // "_n_cf () result (n) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int) :: n"
    write (u, "(A)")  "  n = 1"
    write (u, "(A)")  "end function " // id // "_n_cf"
    write (u, *)
    write (u, "(A)")  "subroutine " // id // "_flv_state (flv_state) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int), dimension(*), intent(out) :: flv_state"
    write (u, "(A)")  "  flv_state(1:3) = [1,2,3]"
    write (u, "(A)")  "end subroutine " // id // "_flv_state"
    write (u, *)
    write (u, "(A)")  "subroutine " // id // "_hel_state (hel_state) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int), dimension(*), intent(out) :: hel_state"
    write (u, "(A)")  "  hel_state(1:3) = [0,0,0]"
    write (u, "(A)")  "end subroutine " // id // "_hel_state"
    write (u, *)
    write (u, "(A)")  "subroutine " // id // "_col_state &
         &(col_state, ghost_flag) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int), dimension(*), intent(out) &
         &:: col_state"
    write (u, "(A)")  "  logical(c_bool), dimension(*), intent(out) &
         &:: ghost_flag"
    write (u, "(A)")  "  col_state(1:6) = [0,0, 0,0, 0,0]"
    write (u, "(A)")  "  ghost_flag(1:3) = [.false., .false., .false.]"
    write (u, "(A)")  "end subroutine " // id // "_col_state"
    write (u, *)
    write (u, "(A)")  "subroutine " // id // "_color_factors &
         &(cf_index1, cf_index2, color_factors) bind(C)"
    write (u, "(A)")  "  use iso_c_binding"
    write (u, "(A)")  "  use kinds"
    write (u, "(A)")  "  implicit none"
    write (u, "(A)")  "  integer(c_int), dimension(*), intent(out) :: cf_index1"
    write (u, "(A)")  "  integer(c_int), dimension(*), intent(out) :: cf_index2"
    write (u, "(A)")  "  complex(c_default_complex), dimension(*), &
         &intent(out) :: color_factors"
    write (u, "(A)")  "  cf_index1(1:1) = [1]"
    write (u, "(A)")  "  cf_index2(1:1) = [1]"
    write (u, "(A)")  "  color_factors(1:1) = [1]"
    write (u, "(A)")  "end subroutine " // id // "_color_factors"
  end subroutine write_test_me_code_3

  function test_writer_6_type_name () result (string)
    type(string_t) :: string
    string = "test_6"
  end function test_writer_6_type_name

  subroutine test_writer_6_mk &
       (writer, unit, id, os_data, verbose, testflag)
    class(test_writer_6_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    write (unit, "(5A)")  "SOURCES += ", char (id), ".c"
    write (unit, "(5A)")  "OBJECTS += ", char (id), ".lo"
    write (unit, "(5A)")  char (id), ".lo: ", char (id), ".c"
    if (.not. verbose) then
       write (unit, "(5A)")  TAB // '@echo  "  FC       " $@'
    end if
    write (unit, "(5A)")  TAB, "$(LTCCOMPILE) $<"
  end subroutine test_writer_6_mk

  subroutine test_writer_6_src (writer, id)
    class(test_writer_6_t), intent(in) :: writer
    type(string_t), intent(in) :: id
    call write_test_c_lib_file (id, var_str ("proc1"))
  end subroutine test_writer_6_src

  subroutine write_test_c_lib_file (basename, feature)
    type(string_t), intent(in) :: basename
    type(string_t), intent(in) :: feature
    integer :: u
    u = free_unit ()
    open (u, file = char (basename) // ".c", &
         status = "replace", action = "write")
    write (u, "(A)")  "/* (Pseudo) matrix element code file &
         &for WHIZARD self-test */"
    write (u, "(A)")  "#include <stdbool.h>"
    if (CC_HAS_QUADMATH) then
       write (u, "(A)")  "#include <quadmath.h>"
    end if
    write (u, *)
    call write_test_me_code_4 (u, char (basename))
    write (u, *)
    write (u, "(A)")  "void " // char (basename) // "_" &
         // char (feature) // "(int* n) {"
    write (u, "(A)")  "  *n = 42;"
    write (u, "(A)")  "}"
    close (u)
  end subroutine write_test_c_lib_file

  subroutine write_test_me_code_4 (u, id)
    integer, intent(in) :: u
    character(*), intent(in) :: id
    write (u, "(A)")  "char* " // id // "_get_md5sum() {"
    write (u, "(A)")  "  return ""1234567890abcdef1234567890abcdef"";"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "bool " // id // "_openmp_supported() {"
    write (u, "(A)")  "  return false;"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "int " // id // "_n_in() {"
    write (u, "(A)")  "  return 1;"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "int " // id // "_n_out() {"
    write (u, "(A)")  "  return 2;"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "int " // id // "_n_flv() {"
    write (u, "(A)")  "  return 1;"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "int " // id // "_n_hel() {"
    write (u, "(A)")  "  return 1;"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "int " // id // "_n_cin() {"
    write (u, "(A)")  "  return 2;"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "int " // id // "_n_col() {"
    write (u, "(A)")  "  return 1;"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "int " // id // "_n_cf() {"
    write (u, "(A)")  "  return 1;"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "void " // id // "_flv_state( int (*a)[] ) {"
    write (u, "(A)")  "  static int flv_state[1][3] =  { { 1, 2, 3 } };"
    write (u, "(A)")  "  int j;"
    write (u, "(A)")  "  for (j = 0; j < 3; j++) { (*a)[j] &
         &= flv_state[0][j]; }"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "void " // id // "_hel_state( int (*a)[] ) {"
    write (u, "(A)")  "  static int hel_state[1][3] =  { { 0, 0, 0 } };"
    write (u, "(A)")  "  int j;"
    write (u, "(A)")  "  for (j = 0; j < 3; j++) { (*a)[j] &
         &= hel_state[0][j]; }"
    write (u, "(A)")  "}"
    write (u, *)
    write (u, "(A)")  "void " // id // "_col_state&
         &( int (*a)[], bool (*g)[] ) {"
    write (u, "(A)")  "  static int col_state[1][3][2] = &
         &{ { {0, 0}, {0, 0}, {0, 0} } };"
    write (u, "(A)")  "  static bool ghost_flag[1][3] =  &
         &{ { false, false, false } };"
    write (u, "(A)")  "  int j,k;"
    write (u, "(A)")  "  for (j = 0; j < 3; j++) {"
    write (u, "(A)")  "    for (k = 0; k < 2; k++) {"
    write (u, "(A)")  "       (*a)[j*2+k] = col_state[0][j][k];"
    write (u, "(A)")  "    }"
    write (u, "(A)")  "    (*g)[j] = ghost_flag[0][j];"
    write (u, "(A)")  "  }"
    write (u, "(A)")  "}"
    write (u, *)
    select case (DEFAULT_FC_PRECISION)
    case ("quadruple")
       write (u, "(A)")  "void " // id // "_color_factors&
            &( int (*cf_index1)[], int (*cf_index2)[], &
            &__complex128 (*color_factors)[] ) {"
    case ("extended")
       write (u, "(A)")  "void " // id // "_color_factors&
            &( int (*cf_index1)[], int (*cf_index2)[], &
            &long double _Complex (*color_factors)[] ) {"
    case default
       write (u, "(A)")  "void " // id // "_color_factors&
            &( int (*cf_index1)[], int (*cf_index2)[], &
            &double _Complex (*color_factors)[] ) {"
    end select
    write (u, "(A)")  "  (*color_factors)[0] = 1;"
    write (u, "(A)")  "  (*cf_index1)[0] = 1;"
    write (u, "(A)")  "  (*cf_index2)[0] = 1;"
    write (u, "(A)")  "}"
  end subroutine write_test_me_code_4

  subroutine check_file (u, file)
    integer, intent(in) :: u
    character(*), intent(in) :: file
    logical :: exist
    inquire (file=file, exist=exist)
    write (u, "(2x,A,A,L1)")  file, " = ", exist
  end subroutine check_file


end module prclib_interfaces_uti
