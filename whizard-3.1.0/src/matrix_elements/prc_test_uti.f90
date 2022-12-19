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

module prc_test_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use particle_specifiers, only: new_prt_spec
  use process_constants
  use prc_core_def
  use process_libraries

  use prc_test

  implicit none
  private

  public :: prc_test_1
  public :: prc_test_2
  public :: prc_test_3
  public :: prc_test_4

contains

  subroutine prc_test_1 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(string_t) :: model_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(process_constants_t) :: data
    class(prc_core_driver_t), allocatable :: driver
    real(default), dimension(0:3,4) :: p
    integer :: i

    write (u, "(A)")  "* Test output: prc_test_1"
    write (u, "(A)")  "*   Purpose: create a trivial process"
    write (u, "(A)")  "*            build a library and &
         &access the matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call os_data%init ()
    call lib%init (var_str ("prc_test1"))

    model_name = "Test"
    allocate (prt_in (2), prt_out (2))
    prt_in  = [var_str ("s"), var_str ("s")]
    prt_out = [var_str ("s"), var_str ("s")]

    allocate (prc_test_def_t :: def)
    select type (def)
    type is (prc_test_def_t)
       call def%init (model_name, prt_in, prt_out)
    end select
    allocate (entry)
    call entry%init (var_str ("prc_test1_a"), model_name = model_name, &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("test_me"), &
         variant = def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure library"
    write (u, "(A)")
    call lib%configure (os_data)

    write (u, "(A)")  "* Load library"
    write (u, "(A)")
    call lib%load (os_data)

    call lib%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active                 = ", &
         lib%is_active ()
    write (u, "(1x,A,I0)")  "n_processes               = ", &
         lib%get_n_processes ()

    write (u, "(A)")
    write (u, "(A)")  "* Constants of prc_test1_a_i1:"
    write (u, "(A)")

    call lib%connect_process (var_str ("prc_test1_a"), 1, data, driver)

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
    write (u, "(1x,A,10(1x,I2))") "hel state =", data%hel_state(:,1)
    write (u, "(1x,A,10(1x,I0))") "col state =", data%col_state
    write (u, "(1x,A,10(1x,L1))") "ghost flag =", data%ghost_flag
    write (u, "(1x,A,10(1x,F5.3))") "color factors =", data%color_factors
    write (u, "(1x,A,10(1x,I0))") "cf index =", data%cf_index

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics:"
    write (u, "(A)")

    p = reshape ([ &
         1.0_default, 0.0_default, 0.0_default, 1.0_default, &
         1.0_default, 0.0_default, 0.0_default,-1.0_default, &
         1.0_default, 1.0_default, 0.0_default, 0.0_default, &
         1.0_default,-1.0_default, 0.0_default, 0.0_default &
         ], [4,4])
    do i = 1, 4
       write (u, "(2x,A,I0,A,4(1x,F7.4))")  "p", i, " =", p(:,i)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Compute matrix element:"
    write (u, "(A)")

    select type (driver)
    type is (prc_test_t)
       write (u, "(1x,A,1x,E11.4)") "|amp| =", abs (driver%get_amplitude (p))
    end select

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_test_1"

  end subroutine prc_test_1

  subroutine prc_test_2 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    class(prc_core_driver_t), allocatable :: driver
    type(process_constants_t) :: data
    real(default), dimension(0:3,4) :: p

    write (u, "(A)")  "* Test output: prc_test_2"
    write (u, "(A)")  "*   Purpose: create a trivial process"
    write (u, "(A)")  "*            build a library and &
         &access the matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a process library with one entry"

    call prc_test_create_library (var_str ("prc_test2"), lib)
    call lib%connect_process (var_str ("prc_test2"), 1, data, driver)

    p = reshape ([ &
         1.0_default, 0.0_default, 0.0_default, 1.0_default, &
         1.0_default, 0.0_default, 0.0_default,-1.0_default, &
         1.0_default, 1.0_default, 0.0_default, 0.0_default, &
         1.0_default,-1.0_default, 0.0_default, 0.0_default &
         ], [4,4])

    write (u, "(A)")
    write (u, "(A)")  "* Compute matrix element:"
    write (u, "(A)")

    select type (driver)
    type is (prc_test_t)
       write (u, "(1x,A,1x,E11.4)") "|amp| =", abs (driver%get_amplitude (p))
    end select

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_test_2"

  end subroutine prc_test_2

  subroutine prc_test_3 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(string_t) :: model_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(process_constants_t) :: data
    class(prc_core_driver_t), allocatable :: driver
    real(default), dimension(0:3,3) :: p
    integer :: i

    write (u, "(A)")  "* Test output: prc_test_3"
    write (u, "(A)")  "*   Purpose: create a trivial decay process"
    write (u, "(A)")  "*            build a library and &
         &access the matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call os_data%init ()
    call lib%init (var_str ("prc_test3"))

    model_name = "Test"
    allocate (prt_in (1), prt_out (2))
    prt_in  = [var_str ("s")]
    prt_out = [var_str ("f"), var_str ("F")]

    allocate (prc_test_def_t :: def)
    select type (def)
    type is (prc_test_def_t)
       call def%init (model_name, prt_in, prt_out)
    end select
    allocate (entry)
    call entry%init (var_str ("prc_test3_a"), model_name = model_name, &
         n_in = 1, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("test_me"), &
         variant = def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure library"
    write (u, "(A)")
    call lib%configure (os_data)

    write (u, "(A)")  "* Load library"
    write (u, "(A)")
    call lib%load (os_data)

    call lib%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active                 = ", &
         lib%is_active ()
    write (u, "(1x,A,I0)")  "n_processes               = ", &
         lib%get_n_processes ()

    write (u, "(A)")
    write (u, "(A)")  "* Constants of prc_test3_a_i1:"
    write (u, "(A)")

    call lib%connect_process (var_str ("prc_test3_a"), 1, data, driver)

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
    write (u, "(1x,A,10(1x,I2))") "hel state =", data%hel_state(:,1)
    write (u, "(1x,A,10(1x,I2))") "hel state =", data%hel_state(:,2)
    write (u, "(1x,A,10(1x,I0))") "col state =", data%col_state
    write (u, "(1x,A,10(1x,L1))") "ghost flag =", data%ghost_flag
    write (u, "(1x,A,10(1x,F5.3))") "color factors =", data%color_factors
    write (u, "(1x,A,10(1x,I0))") "cf index =", data%cf_index

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics:"
    write (u, "(A)")

    p = reshape ([ &
         125._default, 0.0_default, 0.0_default, 0.0_default, &
         62.5_default, 0.0_default, 0.0_default, 62.5_default, &
         62.5_default, 0.0_default, 0.0_default,-62.5_default &
         ], [4,3])
    do i = 1, 3
       write (u, "(2x,A,I0,A,4(1x,F8.4))")  "p", i, " =", p(:,i)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Compute matrix element:"
    write (u, "(A)")

    select type (driver)
    type is (prc_test_t)
       write (u, "(1x,A,1x,E11.4)") "|amp| =", abs (driver%get_amplitude (p))
    end select

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_test_3"

  end subroutine prc_test_3

  subroutine prc_test_4 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    class(prc_core_driver_t), allocatable :: driver
    type(process_constants_t) :: data
    real(default), dimension(0:3,3) :: p

    write (u, "(A)")  "* Test output: prc_test_4"
    write (u, "(A)")  "*   Purpose: create a trivial decay process"
    write (u, "(A)")  "*            build a library and &
         &access the matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a process library with one entry"

    call prc_test_create_library (var_str ("prc_test4"), lib, &
         scattering=.false., decay=.true.)
    call lib%connect_process (var_str ("prc_test4"), 1, data, driver)

    p = reshape ([ &
         125._default, 0.0_default, 0.0_default, 0.0_default, &
         62.5_default, 0.0_default, 0.0_default, 62.5_default, &
         62.5_default, 0.0_default, 0.0_default,-62.5_default &
         ], [4,3])

    write (u, "(A)")
    write (u, "(A)")  "* Compute matrix element:"
    write (u, "(A)")

    select type (driver)
    type is (prc_test_t)
       write (u, "(1x,A,1x,E11.4)") "|amp| =", abs (driver%get_amplitude (p))
    end select

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_test_4"

  end subroutine prc_test_4


end module prc_test_uti
