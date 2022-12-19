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

module process_libraries_uti

  use, intrinsic :: iso_c_binding !NODEP!

  use iso_varying_string, string_t => varying_string
  use io_units
  use os_interface
  use particle_specifiers, only: new_prt_spec
  use process_constants
  use prclib_interfaces
  use prc_core_def

  use process_libraries

  use prclib_interfaces_ut, only: test_writer_4_t

  implicit none
  private

  public :: process_libraries_1
  public :: process_libraries_2
  public :: process_libraries_3
  public :: process_libraries_4
  public :: process_libraries_5
  public :: process_libraries_6
  public :: process_libraries_7
  public :: process_libraries_8

  type, extends (prc_core_def_t) :: prcdef_2_t
     integer :: data = 0
     logical :: file = .false.
   contains
     procedure, nopass :: type_string => prcdef_2_type_string
     procedure :: write => prcdef_2_write
     procedure :: read => prcdef_2_read
     procedure, nopass :: get_features => prcdef_2_get_features
     procedure :: generate_code => prcdef_2_generate_code
     procedure :: allocate_driver => prcdef_2_allocate_driver
     procedure :: connect => prcdef_2_connect
  end type prcdef_2_t

  type, extends (process_driver_internal_t) :: prctest_2_t
   contains
     procedure, nopass :: type_name => prctest_2_type_name
     procedure :: fill_constants => prctest_2_fill_constants
  end type prctest_2_t

  type, extends (prc_core_def_t) :: prcdef_5_t
   contains
     procedure, nopass :: type_string => prcdef_5_type_string
     procedure :: init => prcdef_5_init
     procedure :: write => prcdef_5_write
     procedure :: read => prcdef_5_read
     procedure :: allocate_driver => prcdef_5_allocate_driver
     procedure, nopass :: needs_code => prcdef_5_needs_code
     procedure, nopass :: get_features => prcdef_5_get_features
     procedure :: connect => prcdef_5_connect
  end type prcdef_5_t

  type, extends (prc_core_driver_t) :: prctest_5_t
   contains
     procedure, nopass :: type_name => prctest_5_type_name
  end type prctest_5_t

  type, extends (prc_core_def_t) :: prcdef_6_t
   contains
     procedure, nopass :: type_string => prcdef_6_type_string
     procedure :: init => prcdef_6_init
     procedure :: write => prcdef_6_write
     procedure :: read => prcdef_6_read
     procedure :: allocate_driver => prcdef_6_allocate_driver
     procedure, nopass :: needs_code => prcdef_6_needs_code
     procedure, nopass :: get_features => prcdef_6_get_features
     procedure :: connect => prcdef_6_connect
  end type prcdef_6_t

  abstract interface
     subroutine proc1_t (n) bind(C)
       import
       integer(c_int), intent(out) :: n
     end subroutine proc1_t
  end interface

  type, extends (prc_core_driver_t) :: prctest_6_t
     procedure(proc1_t), nopass, pointer :: proc1 => null ()
   contains
     procedure, nopass :: type_name => prctest_6_type_name
  end type prctest_6_t


contains

  subroutine process_libraries_1 (u)
    integer, intent(in) :: u
    type(process_def_list_t) :: process_def_list

    write (u, "(A)")  "* Test output: process_libraries_1"
    write (u, "(A)")  "*   Purpose: Display an empty process definition list"
    write (u, "(A)")

    call process_def_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_libraries_1"
  end subroutine process_libraries_1

  subroutine process_libraries_2 (u)
    integer, intent(in) :: u
    type(prc_template_t), dimension(:), allocatable :: process_core_templates
    type(process_def_list_t) :: process_def_list
    type(process_def_entry_t), pointer :: entry => null ()
    class(prc_core_def_t), allocatable :: test_def
    integer :: scratch_unit

    write (u, "(A)")  "* Test output: process_libraries_2"
    write (u, "(A)")  "* Purpose: Construct a process definition list,"
    write (u, "(A)")  "*          write it to file and reread it"
    write (u, "(A)")  ""
    write (u, "(A)")  "* Construct a process definition list"
    write (u, "(A)")  "*   First process definition: empty"
    write (u, "(A)")  "*   Second process definition: two components"
    write (u, "(A)")  "*     First component: empty"
    write (u, "(A)")  "*     Second component: test data"
    write (u, "(A)")  "*   Third process definition:"
    write (u, "(A)")  "*     Embedded decays and polarization"
    write (u, "(A)")

    allocate (process_core_templates (1))
    allocate (prcdef_2_t :: process_core_templates(1)%core_def)

    allocate (entry)
    call entry%init (var_str ("first"), n_in = 0, n_components = 0)
    call entry%compute_md5sum ()
    call process_def_list%append (entry)

    allocate (entry)
    call entry%init (var_str ("second"), model_name = var_str ("Test"), &
         n_in = 1, n_components = 2)
    allocate (prcdef_2_t :: test_def)
    select type (test_def)
    type is (prcdef_2_t);  test_def%data = 42
    end select
    call entry%import_component (2, n_out = 2, &
         prt_in  = new_prt_spec ([var_str ("a")]), &
         prt_out = new_prt_spec ([var_str ("b"), var_str ("c")]), &
         method  = var_str ("test"), &
         variant = test_def)
    call entry%compute_md5sum ()
    call process_def_list%append (entry)

    allocate (entry)
    call entry%init (var_str ("third"), model_name = var_str ("Test"), &
         n_in = 2, n_components = 1)
    allocate (prcdef_2_t :: test_def)
    call entry%import_component (1, n_out = 3, &
         prt_in  = &
           new_prt_spec ([var_str ("a"), var_str ("b")]), &
         prt_out = &
           [new_prt_spec (var_str ("c")), &
            new_prt_spec (var_str ("d"), .true.), &
            new_prt_spec (var_str ("e"), [var_str ("e_decay")])], &
         method  = var_str ("test"), &
         variant = test_def)
    call entry%compute_md5sum ()
    call process_def_list%append (entry)
    call process_def_list%write (u)

    write (u, "(A)")  ""
    write (u, "(A)")  "* Write the process definition list to (scratch) file"

    scratch_unit = free_unit ()
    open (unit = scratch_unit, status="scratch", action = "readwrite")
    call process_def_list%write (scratch_unit)
    call process_def_list%final ()

    write (u, "(A)")  "* Reread it"
    write (u, "(A)")  ""

    rewind (scratch_unit)
    call process_def_list%read (scratch_unit, process_core_templates)
    close (scratch_unit)

    call process_def_list%write (u)
    call process_def_list%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_libraries_2"
  end subroutine process_libraries_2

  subroutine process_libraries_3 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    type(process_def_entry_t), pointer :: entry
    class(prc_core_driver_t), allocatable :: driver_template

    write (u, "(A)")  "* Test output: process_libraries_3"
    write (u, "(A)")  "* Purpose: Construct a process library object &
         &with entries"
    write (u, "(A)")  ""
    write (u, "(A)")  "* Construct and display a process library object"
    write (u, "(A)")  "*   with 5 entries"
    write (u, "(A)")  "*   associated with 3 matrix element codes"
    write (u, "(A)")  "*   corresponding to 3 process definitions"
    write (u, "(A)")  "*   with 2, 1, 1 components, respectively"
    write (u, "(A)")

    call lib%init (var_str ("testlib"))

    call lib%set_status (STAT_ACTIVE)
    call lib%allocate_entries (5)

    allocate (entry)
    call entry%init (var_str ("test_a"), n_in = 2, n_components = 2)
    allocate (prctest_2_t :: driver_template)
    call lib%init_entry (3, STAT_SOURCE, entry%process_def_t, 2, 2, &
         driver_template)
    call lib%init_entry (4, STAT_COMPILED, entry%process_def_t, 1, 0)
    call lib%append (entry)

    allocate (entry)
    call entry%init (var_str ("test_b"), n_in = 2, n_components = 1)
    call lib%init_entry (2, STAT_CONFIGURED, entry%process_def_t, 1, 1)
    call lib%append (entry)

    allocate (entry)
    call entry%init (var_str ("test_c"), n_in = 2, n_components = 1)
    allocate (prctest_2_t :: driver_template)
    call lib%init_entry (5, STAT_LINKED, entry%process_def_t, 1, 3, &
         driver_template)
    call lib%append (entry)

    call lib%write (u)
    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_libraries_3"
  end subroutine process_libraries_3

  subroutine process_libraries_4 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    type(process_def_entry_t), pointer :: entry
    class(prc_core_def_t), allocatable :: core_def
    type(os_data_t) :: os_data

    write (u, "(A)")  "* Test output: process_libraries_4"
    write (u, "(A)")  "* Purpose: build a process library with an &
         &internal (pseudo) matrix element"
    write (u, "(A)")  "*          No Makefile or code should be generated"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry &
         &(no external code)"
    write (u, "(A)")
    call os_data%init ()
    call lib%init (var_str ("proclibs4"))

    allocate (prcdef_2_t :: core_def)

    allocate (entry)
    call entry%init (var_str ("proclibs4_a"), n_in = 1, n_components = 1)
    call entry%import_component (1, n_out = 2, variant = core_def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure library"
    write (u, "(A)")
    call lib%configure (os_data)

    write (u, "(A)")  "* Compute MD5 sum"
    write (u, "(A)")
    call lib%compute_md5sum ()

    write (u, "(A)")  "* Write makefile (no-op)"
    write (u, "(A)")
    call lib%write_makefile (os_data, force = .true., verbose = .true.)

    write (u, "(A)")  "* Write driver source code (no-op)"
    write (u, "(A)")
    call lib%write_driver (force = .true.)

    write (u, "(A)")  "* Write process source code (no-op)"
    write (u, "(A)")
    call lib%make_source (os_data)

    write (u, "(A)")  "* Compile (no-op)"
    write (u, "(A)")
    call lib%make_compile (os_data)

    write (u, "(A)")  "* Link (no-op)"
    write (u, "(A)")
    call lib%make_link (os_data)

    write (u, "(A)")  "* Load (no-op)"
    write (u, "(A)")
    call lib%load (os_data)

    call lib%write (u)
    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_libraries_4"
  end subroutine process_libraries_4

  subroutine process_libraries_5 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    type(process_def_entry_t), pointer :: entry
    class(prc_core_def_t), allocatable :: core_def
    type(os_data_t) :: os_data

    write (u, "(A)")  "* Test output: process_libraries_5"
    write (u, "(A)")  "* Purpose: build a process library with an &
         &external (pseudo) matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("proclibs5"))
    call os_data%init ()

    allocate (prcdef_5_t :: core_def)
    select type (core_def)
    type is (prcdef_5_t)
       call core_def%init ()
    end select

    allocate (entry)
    call entry%init (var_str ("proclibs5_a"), &
         model_name = var_str ("Test_Model"), &
         n_in = 1, n_components = 1)
    call entry%import_component (1, n_out = 2, &
         prt_in  = new_prt_spec ([var_str ("a")]), &
         prt_out = new_prt_spec ([var_str ("b"), var_str ("c")]), &
         method  = var_str ("test"), &
         variant = core_def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure library"
    write (u, "(A)")
    call lib%configure (os_data)

    write (u, "(A)")  "* Compute MD5 sum"
    write (u, "(A)")
    call lib%compute_md5sum ()

    write (u, "(A)")  "* Write makefile"
    write (u, "(A)")
    call lib%write_makefile (os_data, force = .true., verbose = .false.)

    write (u, "(A)")  "* Write driver source code"
    write (u, "(A)")
    call lib%write_driver (force = .true.)

    write (u, "(A)")  "* Write process source code"
    write (u, "(A)")
    call lib%make_source (os_data)

    write (u, "(A)")  "* Compile"
    write (u, "(A)")
    call lib%make_compile (os_data)

    write (u, "(A)")  "* Link"
    write (u, "(A)")
    call lib%make_link (os_data)

    call lib%write (u, libpath = .false.)

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_libraries_5"
  end subroutine process_libraries_5

  subroutine process_libraries_6 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    type(process_def_entry_t), pointer :: entry
    class(prc_core_def_t), allocatable :: core_def
    type(os_data_t) :: os_data
    type(string_t), dimension(:), allocatable :: name_list
    type(process_constants_t) :: data
    class(prc_core_driver_t), allocatable :: proc_driver
    integer :: i
    integer(c_int) :: n

    write (u, "(A)")  "* Test output: process_libraries_6"
    write (u, "(A)")  "* Purpose: build and load a process library"
    write (u, "(A)")  "*          with an external (pseudo) matrix element"
    write (u, "(A)")  "*          Check single-call linking"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("proclibs6"))
    call os_data%init ()

    allocate (prcdef_6_t :: core_def)
    select type (core_def)
    type is (prcdef_6_t)
       call core_def%init ()
    end select

    allocate (entry)
    call entry%init (var_str ("proclibs6_a"), &
         model_name = var_str ("Test_model"), &
         n_in = 1, n_components = 1)
    call entry%import_component (1, n_out = 2, &
         prt_in  = new_prt_spec ([var_str ("a")]), &
         prt_out = new_prt_spec ([var_str ("b"), var_str ("c")]), &
         method  = var_str ("test"), &
         variant = core_def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure library"
    write (u, "(A)")
    call lib%configure (os_data)

    write (u, "(A)")  "* Write makefile"
    write (u, "(A)")
    call lib%write_makefile (os_data, force = .true., verbose = .false.)

    write (u, "(A)")  "* Write driver source code"
    write (u, "(A)")
    call lib%write_driver (force = .true.)

    write (u, "(A)")  "* Write process source code, compile, link, load"
    write (u, "(A)")
    call lib%load (os_data)

    call lib%write (u, libpath = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,A,A)")  "name                      = '", &
         char (lib%get_name ()), "'"
    write (u, "(1x,A,L1)")  "is active                 = ", &
         lib%is_active ()
    write (u, "(1x,A,I0)")  "n_processes               = ", &
         lib%get_n_processes ()
    write (u, "(1x,A)", advance="no")  "processes                 ="
    call lib%get_process_id_list (name_list)
    do i = 1, size (name_list)
       write (u, "(1x,A)", advance="no")  char (name_list(i))
    end do
    write (u, *)
    write (u, "(1x,A,L1)")  "proclibs6_a is process    = ", &
         lib%contains (var_str ("proclibs6_a"))
    write (u, "(1x,A,I0)")  "proclibs6_a has index     = ", &
         lib%get_entry_index (var_str ("proclibs6_a"))
    write (u, "(1x,A,L1)")  "foobar is process         = ", &
         lib%contains (var_str ("foobar"))
    write (u, "(1x,A,I0)")  "foobar has index          = ", &
         lib%get_entry_index (var_str ("foobar"))
    write (u, "(1x,A,I0)")  "n_in(proclibs6_a)         = ", &
         lib%get_n_in (var_str ("proclibs6_a"))
    write (u, "(1x,A,A)")   "model_name(proclibs6_a)   = ", &
         char (lib%get_model_name (var_str ("proclibs6_a")))
    write (u, "(1x,A)", advance="no")  "components(proclibs6_a)   ="
    call lib%get_component_list (var_str ("proclibs6_a"), name_list)
    do i = 1, size (name_list)
       write (u, "(1x,A)", advance="no")  char (name_list(i))
    end do
    write (u, *)

    write (u, "(A)")
    write (u, "(A)")  "* Constants of proclibs6_a_i1:"
    write (u, "(A)")

    call lib%connect_process (var_str ("proclibs6_a"), 1, data, proc_driver)

    write (u, "(1x,A,A)")  "component ID     = ", char (data%id)
    write (u, "(1x,A,A)")  "model name       = ", char (data%model_name)
    write (u, "(1x,A,A,A)")  "md5sum           = '", data%md5sum, "'"
    write (u, "(1x,A,L1)") "openmp supported = ", data%openmp_supported
    write (u, "(1x,A,I0)") "n_in  = ", data%n_in
    write (u, "(1x,A,I0)") "n_out = ", data%n_out
    write (u, "(1x,A,I0)") "n_flv = ", data%n_flv
    write (u, "(1x,A,I0)") "n_hel = ", data%n_hel
    write (u, "(1x,A,I0)") "n_col = ", data%n_col
    write (u, "(1x,A,I0)") "n_cin = ", data%n_cin
    write (u, "(1x,A,I0)") "n_cf  = ", data%n_cf
    write (u, "(1x,A,10(1x,I0))") "flv state =", data%flv_state
    write (u, "(1x,A,10(1x,I0))") "hel state =", data%hel_state
    write (u, "(1x,A,10(1x,I0))") "col state =", data%col_state
    write (u, "(1x,A,10(1x,L1))") "ghost flag =", data%ghost_flag
    write (u, "(1x,A,10(1x,F5.3))") "color factors =", data%color_factors
    write (u, "(1x,A,10(1x,I0))") "cf index =", data%cf_index

    write (u, "(A)")
    write (u, "(A)")  "* Call feature of proclibs6_a:"
    write (u, "(A)")

    select type (proc_driver)
    type is (prctest_6_t)
       call proc_driver%proc1 (n)
       write (u, "(1x,A,I0)") "proc1 = ", n
    end select

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_libraries_6"
  end subroutine process_libraries_6

  subroutine process_libraries_7 (u)
    integer, intent(in) :: u
    type(prc_template_t), dimension(:), allocatable :: process_core_templates
    type(process_def_entry_t), target :: entry
    class(prc_core_def_t), allocatable :: test_def
    class(prc_core_def_t), pointer :: def

    write (u, "(A)")  "* Test output: process_libraries_7"
    write (u, "(A)")  "* Purpose: Construct a process definition list &
         &and check MD5 sums"
    write (u, "(A)")
    write (u, "(A)")  "* Construct a process definition list"
    write (u, "(A)")  "*   Process: two components"
    write (u, "(A)")

    allocate (process_core_templates (1))
    allocate (prcdef_2_t :: process_core_templates(1)%core_def)

    call entry%init (var_str ("first"), model_name = var_str ("Test"), &
         n_in = 1, n_components = 2)
    allocate (prcdef_2_t :: test_def)
    select type (test_def)
    type is (prcdef_2_t);  test_def%data = 31
    end select
    call entry%import_component (1, n_out = 3, &
         prt_in  = new_prt_spec ([var_str ("a")]), &
         prt_out = new_prt_spec ([var_str ("b"), var_str ("c"), &
                                  var_str ("e")]), &
         method  = var_str ("test"), &
         variant = test_def)
    allocate (prcdef_2_t :: test_def)
    select type (test_def)
    type is (prcdef_2_t);  test_def%data = 42
    end select
    call entry%import_component (2, n_out = 2, &
         prt_in  = new_prt_spec ([var_str ("a")]), &
         prt_out = new_prt_spec ([var_str ("b"), var_str ("c")]), &
         method  = var_str ("test"), &
         variant = test_def)
    call entry%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute MD5 sums"
    write (u, "(A)")

    call entry%compute_md5sum ()
    call entry%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recalculate MD5 sums (should be identical)"
    write (u, "(A)")

    call entry%compute_md5sum ()
    call entry%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Modify a component and recalculate MD5 sums"
    write (u, "(A)")

    def => entry%get_core_def_ptr (2)
    select type (def)
    type is (prcdef_2_t)
       def%data = 54
    end select
    call entry%compute_md5sum ()
    call entry%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Modify the model and recalculate MD5 sums"
    write (u, "(A)")

    call entry%set_model_name (var_str ("foo"))
    call entry%compute_md5sum ()
    call entry%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_libraries_7"
  end subroutine process_libraries_7

  subroutine process_libraries_8 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    type(process_def_entry_t), pointer :: entry
    class(prc_core_def_t), allocatable :: core_def
    type(os_data_t) :: os_data

    write (u, "(A)")  "* Test output: process_libraries_8"
    write (u, "(A)")  "* Purpose: build and load a process library"
    write (u, "(A)")  "*          with an external (pseudo) matrix element"
    write (u, "(A)")  "*          Check status updates"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("proclibs8"))
    call os_data%init ()

    allocate (prcdef_6_t :: core_def)
    select type (core_def)
    type is (prcdef_6_t)
       call core_def%init ()
    end select

    allocate (entry)
    call entry%init (var_str ("proclibs8_a"), &
         model_name = var_str ("Test_model"), &
         n_in = 1, n_components = 1)
    call entry%import_component (1, n_out = 2, &
         prt_in  = new_prt_spec ([var_str ("a")]), &
         prt_out = new_prt_spec ([var_str ("b"), var_str ("c")]), &
         method  = var_str ("test"), &
         variant = core_def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure library"
    write (u, "(A)")

    call lib%configure (os_data)
    call lib%compute_md5sum ()

    call lib%test_transfer_md5sum (1, 1, 1)

    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)

    write (u, "(A)")
    write (u, "(A)")  "* Write makefile"
    write (u, "(A)")
    call lib%write_makefile (os_data, force = .true., verbose = .false.)

    write (u, "(A)")  "* Update status"
    write (u, "(A)")

    call lib%update_status (os_data)
    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)

    write (u, "(A)")
    write (u, "(A)")  "* Write driver source code"
    write (u, "(A)")
    call lib%write_driver (force = .false.)

    write (u, "(A)")  "* Write process source code"
    write (u, "(A)")
    call lib%make_source (os_data)

    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)

    write (u, "(A)")
    write (u, "(A)")  "* Compile and load"
    write (u, "(A)")

    call lib%load (os_data)
    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)

    write (u, "(A)")
    write (u, "(A)")  "* Append process and reconfigure"
    write (u, "(A)")

    allocate (prcdef_6_t :: core_def)
    select type (core_def)
    type is (prcdef_6_t)
       call core_def%init ()
    end select

    allocate (entry)
    call entry%init (var_str ("proclibs8_b"), &
         model_name = var_str ("Test_model"), &
         n_in = 1, n_components = 1)
    call entry%import_component (1, n_out = 2, &
         prt_in  = new_prt_spec ([var_str ("a")]), &
         prt_out = new_prt_spec ([var_str ("b"), var_str ("d")]), &
         method  = var_str ("test"), &
         variant = core_def)
    call lib%append (entry)

    call lib%configure (os_data)
    call lib%compute_md5sum ()
    call lib%test_transfer_md5sum (2, 2, 1)
    call lib%write_makefile (os_data, force = .false., verbose = .false.)
    call lib%write_driver (force = .false.)

    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)
    write (u, "(1x,A,I0)")  "proc2 status = ", lib%get_status (2)

    write (u, "(A)")
    write (u, "(A)")  "* Update status"
    write (u, "(A)")

    call lib%update_status (os_data)
    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)
    write (u, "(1x,A,I0)")  "proc2 status = ", lib%get_status (2)

    write (u, "(A)")
    write (u, "(A)")  "* Write source code"
    write (u, "(A)")

    call lib%make_source (os_data)
    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)
    write (u, "(1x,A,I0)")  "proc2 status = ", lib%get_status (2)

    write (u, "(A)")
    write (u, "(A)")  "* Reset status"
    write (u, "(A)")

    call lib%set_status (STAT_CONFIGURED, entries=.true.)
    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)
    write (u, "(1x,A,I0)")  "proc2 status = ", lib%get_status (2)

    write (u, "(A)")
    write (u, "(A)")  "* Update status"
    write (u, "(A)")

    call lib%update_status (os_data)
    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)
    write (u, "(1x,A,I0)")  "proc2 status = ", lib%get_status (2)

    write (u, "(A)")
    write (u, "(A)")  "* Partial cleanup"
    write (u, "(A)")

    call lib%clean (os_data, distclean = .false.)
    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)
    write (u, "(1x,A,I0)")  "proc2 status = ", lib%get_status (2)


    write (u, "(A)")
    write (u, "(A)")  "* Update status"
    write (u, "(A)")

    call lib%update_status (os_data)
    write (u, "(1x,A,L1)")  "library loaded = ", lib%is_loaded ()
    write (u, "(1x,A,I0)")  "lib status   = ", lib%get_status ()
    write (u, "(1x,A,I0)")  "proc1 status = ", lib%get_status (1)
    write (u, "(1x,A,I0)")  "proc2 status = ", lib%get_status (2)

    write (u, "(A)")
    write (u, "(A)")  "* Complete cleanup"

    call lib%clean (os_data, distclean = .true.)
    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: process_libraries_8"
  end subroutine process_libraries_8


  function prcdef_2_type_string () result (string)
    type(string_t) :: string
    string = "test"
  end function prcdef_2_type_string

  subroutine prcdef_2_write (object, unit)
    class(prcdef_2_t), intent(in) :: object
    integer, intent(in) :: unit
    write (unit, "(3x,A,I0)")  "Test data         = ", object%data
  end subroutine prcdef_2_write

  subroutine prcdef_2_read (object, unit)
    class(prcdef_2_t), intent(out) :: object
    integer, intent(in) :: unit
    character(80) :: buffer
    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    read (buffer, *)  object%data
  end subroutine prcdef_2_read

  subroutine prcdef_2_get_features (features)
    type(string_t), dimension(:), allocatable, intent(out) :: features
    allocate (features (0))
  end subroutine prcdef_2_get_features

  subroutine prcdef_2_generate_code (object, &
       basename, model_name, prt_in, prt_out)
    class(prcdef_2_t), intent(in) :: object
    type(string_t), intent(in) :: basename
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in
    type(string_t), dimension(:), intent(in) :: prt_out
  end subroutine prcdef_2_generate_code

  subroutine prcdef_2_allocate_driver (object, driver, basename)
    class(prcdef_2_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    allocate (prctest_2_t :: driver)
  end subroutine prcdef_2_allocate_driver

  subroutine prcdef_2_connect (def, lib_driver, i, proc_driver)
    class(prcdef_2_t), intent(in) :: def
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
  end subroutine prcdef_2_connect

  function prctest_2_type_name () result (type)
    type(string_t) :: type
    type = "test"
  end function prctest_2_type_name

  subroutine prctest_2_fill_constants (driver, data)
    class(prctest_2_t), intent(in) :: driver
    type(process_constants_t), intent(out) :: data
  end subroutine prctest_2_fill_constants

  function prcdef_5_type_string () result (string)
    type(string_t) :: string
    string = "test_file"
  end function prcdef_5_type_string

  subroutine prcdef_5_init (object)
    class(prcdef_5_t), intent(out) :: object
    allocate (test_writer_4_t :: object%writer)
  end subroutine prcdef_5_init

  subroutine prcdef_5_write (object, unit)
    class(prcdef_5_t), intent(in) :: object
    integer, intent(in) :: unit
  end subroutine prcdef_5_write

  subroutine prcdef_5_read (object, unit)
    class(prcdef_5_t), intent(out) :: object
    integer, intent(in) :: unit
  end subroutine prcdef_5_read

  subroutine prcdef_5_allocate_driver (object, driver, basename)
    class(prcdef_5_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    allocate (prctest_5_t :: driver)
  end subroutine prcdef_5_allocate_driver

  function prcdef_5_needs_code () result (flag)
    logical :: flag
    flag = .true.
  end function prcdef_5_needs_code

  subroutine prcdef_5_get_features (features)
    type(string_t), dimension(:), allocatable, intent(out) :: features
    allocate (features (1))
    features = [ var_str ("proc1") ]
  end subroutine prcdef_5_get_features

  subroutine prcdef_5_connect (def, lib_driver, i, proc_driver)
    class(prcdef_5_t), intent(in) :: def
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
  end subroutine prcdef_5_connect

  function prctest_5_type_name () result (type)
    type(string_t) :: type
    type = "test_file"
  end function prctest_5_type_name

  function prcdef_6_type_string () result (string)
    type(string_t) :: string
    string = "test_file"
  end function prcdef_6_type_string

  subroutine prcdef_6_init (object)
    class(prcdef_6_t), intent(out) :: object
    allocate (test_writer_4_t :: object%writer)
    call object%writer%init_test ()
  end subroutine prcdef_6_init

  subroutine prcdef_6_write (object, unit)
    class(prcdef_6_t), intent(in) :: object
    integer, intent(in) :: unit
  end subroutine prcdef_6_write

  subroutine prcdef_6_read (object, unit)
    class(prcdef_6_t), intent(out) :: object
    integer, intent(in) :: unit
  end subroutine prcdef_6_read

  subroutine prcdef_6_allocate_driver (object, driver, basename)
    class(prcdef_6_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    allocate (prctest_6_t :: driver)
  end subroutine prcdef_6_allocate_driver

  function prcdef_6_needs_code () result (flag)
    logical :: flag
    flag = .true.
  end function prcdef_6_needs_code

  subroutine prcdef_6_get_features (features)
    type(string_t), dimension(:), allocatable, intent(out) :: features
    allocate (features (1))
    features = [ var_str ("proc1") ]
  end subroutine prcdef_6_get_features

  subroutine prcdef_6_connect (def, lib_driver, i, proc_driver)
    class(prcdef_6_t), intent(in) :: def
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
    integer(c_int) :: pid, fid
    type(c_funptr) :: fptr
    select type (proc_driver)
    type is  (prctest_6_t)
       pid = i
       fid = 1
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%proc1)
    end select
  end subroutine prcdef_6_connect

  function prctest_6_type_name () result (type)
    type(string_t) :: type
    type = "test_file"
  end function prctest_6_type_name


end module process_libraries_uti
