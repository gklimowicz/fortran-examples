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

module prc_template_me_uti

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use os_interface
  use particle_specifiers, only: new_prt_spec
  use model_data
  use prc_core_def
  use process_constants
  use process_libraries
  use model_testbed, only: prepare_model, cleanup_model

  use prc_template_me

  implicit none
  private

  public :: prc_template_me_1
  public :: prc_template_me_2

contains

  subroutine prc_template_me_1 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(string_t) :: model_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(process_constants_t) :: data
    class(prc_core_driver_t), allocatable :: driver
    integer, parameter :: cdf = c_default_float
    integer, parameter :: ci = c_int
    real(cdf), dimension(4) :: par
    real(cdf), dimension(0:3,4) :: p
    logical(c_bool) :: flag
    complex(c_default_complex) :: amp
    integer :: i

    write (u, "(A)")  "* Test output: prc_template_me_1"
    write (u, "(A)")  "*   Purpose: create a template matrix element,"
    write (u, "(A)")  "*            normalized to give unit integral,"
    write (u, "(A)")  "*            build a library, link, load, and &
         &access the matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("template_me1"))
    call os_data%init ()

    model_name = "QED"
    model => null ()
    call prepare_model (model, model_name)

    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("e+"), var_str ("e-")]
    prt_out = [var_str ("m+"), var_str ("m-")]

    allocate (template_me_def_t :: def)
    select type (def)
    type is (template_me_def_t)
       call def%init (model, prt_in, prt_out, unity = .false.)
    end select
    allocate (entry)
    call entry%init (var_str ("template_me1_a"), model_name = model_name, &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("template"), &
         variant = def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure library"
    write (u, "(A)")
    call lib%configure (os_data)

    write (u, "(A)")  "* Write makefile"
    write (u, "(A)")
    call lib%write_makefile (os_data, force = .true., verbose = .false.)

    write (u, "(A)")  "* Clean any left-over files"
    write (u, "(A)")
    call lib%clean (os_data, distclean = .false.)

    write (u, "(A)")  "* Write driver"
    write (u, "(A)")
    call lib%write_driver (force = .true.)

    write (u, "(A)")  "* Write process source code, compile, link, load"
    write (u, "(A)")
    call lib%load (os_data)

    call lib%write (u, libpath = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active                 = ", &
         lib%is_active ()
    write (u, "(1x,A,I0)")  "n_processes               = ", &
         lib%get_n_processes ()

    write (u, "(A)")
    write (u, "(A)")  "* Constants of template_me1_a_i1:"
    write (u, "(A)")

    call lib%connect_process (var_str ("template_me1_a"), 1, data, driver)

    write (u, "(1x,A,A)")  "component ID     = ", char (data%id)
    write (u, "(1x,A,A)")  "model name       = ", char (data%model_name)
    write (u, "(1x,A,A,A)")  "md5sum           = '", data%md5sum, "'"
    write (u, "(1x,A,I0)") "n_in  = ", data%n_in
    write (u, "(1x,A,I0)") "n_out = ", data%n_out
    write (u, "(1x,A,I0)") "n_flv = ", data%n_flv
    write (u, "(1x,A,I0)") "n_hel = ", data%n_hel
    write (u, "(1x,A,I0)") "n_col = ", data%n_col
    write (u, "(1x,A,I0)") "n_cin = ", data%n_cin
    write (u, "(1x,A,I0)") "n_cf  = ", data%n_cf
    write (u, "(1x,A,10(1x,I0))") "flv state =", data%flv_state
    write (u, "(1x,A,10(1x,I2))") "hel state =", data%hel_state(:,1)
    do i = 2, 16
       write (u, "(12x,4(1x,I2))")  data%hel_state(:,i)
    end do
    write (u, "(1x,A,10(1x,I0))") "col state =", data%col_state
    write (u, "(1x,A,10(1x,L1))") "ghost flag =", data%ghost_flag
    write (u, "(1x,A,10(1x,F5.3))") "color factors =", data%color_factors
    write (u, "(1x,A,10(1x,I0))") "cf index =", data%cf_index

    write (u, "(A)")
    write (u, "(A)")  "* Set parameters for template_me1_a and initialize:"
    write (u, "(A)")

    par = [0.3_cdf, 0.0_cdf, 0.0_cdf, 0.0_cdf]
    write (u, "(2x,A,F6.4)")  "ee   = ", par(1)
    write (u, "(2x,A,F6.4)")  "me   = ", par(2)
    write (u, "(2x,A,F6.4)")  "mmu  = ", par(3)
    write (u, "(2x,A,F6.4)")  "mtau = ", par(4)

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics:"
    write (u, "(A)")

    p = reshape ([ &
         1.0_cdf, 0.0_cdf, 0.0_cdf, 1.0_cdf, &
         1.0_cdf, 0.0_cdf, 0.0_cdf,-1.0_cdf, &
         1.0_cdf, 1.0_cdf, 0.0_cdf, 0.0_cdf, &
         1.0_cdf,-1.0_cdf, 0.0_cdf, 0.0_cdf &
         ], [4,4])
    do i = 1, 4
       write (u, "(2x,A,I0,A,4(1x,F7.4))")  "p", i, " =", p(:,i)
    end do

    select type (driver)
    type is (template_me_driver_t)
       call driver%init (par, 0)

       call driver%new_event (p)

       write (u, "(A)")
       write (u, "(A)")  "* Compute matrix element:"
       write (u, "(A)")

       call driver%is_allowed (1_ci, 6_ci, 1_ci, flag)
       write (u, "(1x,A,L1)") "is_allowed (1, 6, 1) = ", flag

       call driver%get_amplitude (1_ci, 6_ci, 1_ci, amp)
       write (u, "(1x,A,1x,E11.4)") "|amp (1, 6, 1)| =", abs (amp)
    end select

    call lib%final ()
    call cleanup_model (model)
    deallocate (model)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_template_me_1"

  end subroutine prc_template_me_1

  subroutine prc_template_me_2 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(string_t) :: model_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(process_constants_t) :: data
    class(prc_core_driver_t), allocatable :: driver
    integer, parameter :: cdf = c_default_float
    integer, parameter :: ci = c_int
    real(cdf), dimension(4) :: par
    real(cdf), dimension(0:3,4) :: p
    logical(c_bool) :: flag
    complex(c_default_complex) :: amp
    integer :: i

    write (u, "(A)")  "* Test output: prc_template_me_1"
    write (u, "(A)")  "*   Purpose: create a template matrix element,"
    write (u, "(A)")  "*            being identical to unity,"
    write (u, "(A)")  "*            build a library, link, load, and &
         &access the matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("template_me2"))
    call os_data%init ()

    model_name = "QED"
    model => null ()
    call prepare_model (model, model_name)

    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("e+"), var_str ("e-")]
    prt_out = [var_str ("m+"), var_str ("m-")]

    allocate (template_me_def_t :: def)
    select type (def)
    type is (template_me_def_t)
       call def%init (model, prt_in, prt_out, unity = .true.)
    end select
    allocate (entry)
    call entry%init (var_str ("template_me2_a"), model_name = model_name, &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("template_unity"), &
         variant = def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure library"
    write (u, "(A)")
    call lib%configure (os_data)

    write (u, "(A)")  "* Write makefile"
    write (u, "(A)")
    call lib%write_makefile (os_data, force = .true., verbose = .false.)

    write (u, "(A)")  "* Clean any left-over files"
    write (u, "(A)")
    call lib%clean (os_data, distclean = .false.)

    write (u, "(A)")  "* Write driver"
    write (u, "(A)")
    call lib%write_driver (force = .true.)

    write (u, "(A)")  "* Write process source code, compile, link, load"
    write (u, "(A)")
    call lib%load (os_data)

    call lib%write (u, libpath = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active                 = ", &
         lib%is_active ()
    write (u, "(1x,A,I0)")  "n_processes               = ", &
         lib%get_n_processes ()

    write (u, "(A)")
    write (u, "(A)")  "* Constants of template_me2_a_i1:"
    write (u, "(A)")

    call lib%connect_process (var_str ("template_me2_a"), 1, data, driver)

    write (u, "(1x,A,A)")  "component ID     = ", char (data%id)
    write (u, "(1x,A,A)")  "model name       = ", char (data%model_name)
    write (u, "(1x,A,A,A)")  "md5sum           = '", data%md5sum, "'"
    write (u, "(1x,A,I0)") "n_in  = ", data%n_in
    write (u, "(1x,A,I0)") "n_out = ", data%n_out
    write (u, "(1x,A,I0)") "n_flv = ", data%n_flv
    write (u, "(1x,A,I0)") "n_hel = ", data%n_hel
    write (u, "(1x,A,I0)") "n_col = ", data%n_col
    write (u, "(1x,A,I0)") "n_cin = ", data%n_cin
    write (u, "(1x,A,I0)") "n_cf  = ", data%n_cf
    write (u, "(1x,A,10(1x,I0))") "flv state =", data%flv_state
    write (u, "(1x,A,10(1x,I2))") "hel state =", data%hel_state(:,1)
    do i = 2, 16
       write (u, "(12x,4(1x,I2))")  data%hel_state(:,i)
    end do
    write (u, "(1x,A,10(1x,I0))") "col state =", data%col_state
    write (u, "(1x,A,10(1x,L1))") "ghost flag =", data%ghost_flag
    write (u, "(1x,A,10(1x,F5.3))") "color factors =", data%color_factors
    write (u, "(1x,A,10(1x,I0))") "cf index =", data%cf_index

    write (u, "(A)")
    write (u, "(A)")  "* Set parameters for template_me2_a and initialize:"
    write (u, "(A)")

    par = [0.3_cdf, 0.0_cdf, 0.0_cdf, 0.0_cdf]
    write (u, "(2x,A,F6.4)")  "ee   = ", par(1)
    write (u, "(2x,A,F6.4)")  "me   = ", par(2)
    write (u, "(2x,A,F6.4)")  "mmu  = ", par(3)
    write (u, "(2x,A,F6.4)")  "mtau = ", par(4)

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics:"
    write (u, "(A)")

    p = reshape ([ &
         1.0_cdf, 0.0_cdf, 0.0_cdf, 1.0_cdf, &
         1.0_cdf, 0.0_cdf, 0.0_cdf,-1.0_cdf, &
         1.0_cdf, 1.0_cdf, 0.0_cdf, 0.0_cdf, &
         1.0_cdf,-1.0_cdf, 0.0_cdf, 0.0_cdf &
         ], [4,4])
    do i = 1, 4
       write (u, "(2x,A,I0,A,4(1x,F7.4))")  "p", i, " =", p(:,i)
    end do

    select type (driver)
    type is (template_me_driver_t)
       call driver%init (par, 0)

       call driver%new_event (p)

       write (u, "(A)")
       write (u, "(A)")  "* Compute matrix element:"
       write (u, "(A)")

       call driver%is_allowed (1_ci, 6_ci, 1_ci, flag)
       write (u, "(1x,A,L1)") "is_allowed (1, 6, 1) = ", flag

       call driver%get_amplitude (1_ci, 6_ci, 1_ci, amp)
       write (u, "(1x,A,1x,E11.4)") "|amp (1, 6, 1)| =", abs (amp)
    end select

    call lib%final ()
    call cleanup_model (model)
    deallocate (model)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_template_me_2"

  end subroutine prc_template_me_2


end module prc_template_me_uti
