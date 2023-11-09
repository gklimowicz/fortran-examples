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

module prc_omega_uti

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use io_units
  use file_utils, only: delete_file
  use os_interface
  use sm_qcd
  use lorentz
  use model_data
  use var_base
  use particle_specifiers, only: new_prt_spec
  use prc_core_def
  use process_constants
  use process_libraries
  use prc_core
  use model_testbed, only: prepare_model, cleanup_model

  use prc_omega

  implicit none
  private

  public :: prc_omega_1
  public :: prc_omega_2
  public :: prc_omega_3
  public :: prc_omega_4
  public :: prc_omega_5
  public :: prc_omega_6
  public :: prc_omega_diags_1

contains

  subroutine prc_omega_1 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
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

    write (u, "(A)")  "* Test output: prc_omega_1"
    write (u, "(A)")  "*   Purpose: create a simple process with OMega"
    write (u, "(A)")  "*            build a library, link, load, and &
         &access the matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("omega1"))
    call os_data%init ()

    model_name = "QED"
    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("e+"), var_str ("e-")]
    prt_out = [var_str ("m+"), var_str ("m-")]

    allocate (omega_def_t :: def)
    select type (def)
    type is (omega_def_t)
       call def%init (model_name, prt_in, prt_out, &
            ufo = .false., ovm = .false.)
    end select
    allocate (entry)
    call entry%init (var_str ("omega1_a"), model_name = model_name, &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("omega"), &
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
    write (u, "(A)")  "* Constants of omega1_a_i1:"
    write (u, "(A)")

    call lib%connect_process (var_str ("omega1_a"), 1, data, driver)

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
    do i = 2, 16
       write (u, "(12x,4(1x,I2))")  data%hel_state(:,i)
    end do
    write (u, "(1x,A,10(1x,I0))") "col state =", data%col_state
    write (u, "(1x,A,10(1x,L1))") "ghost flag =", data%ghost_flag
    write (u, "(1x,A,10(1x,F5.3))") "color factors =", data%color_factors
    write (u, "(1x,A,10(1x,I0))") "cf index =", data%cf_index

    write (u, "(A)")
    write (u, "(A)")  "* Set parameters for omega1_a and initialize:"
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
    type is (omega_driver_t)
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

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_omega_1"

  end subroutine prc_omega_1

  subroutine prc_omega_2 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    type(string_t) :: model_name
    class(model_data_t), pointer :: model
    class(vars_t), pointer :: vars
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(string_t) :: restrictions
    type(process_component_def_t), pointer :: config
    type(prc_omega_t) :: prc1, prc2
    type(process_constants_t) :: data
    integer, parameter :: cdf = c_default_float
    integer, parameter :: ci = c_int
    real(cdf), dimension(:), allocatable :: par
    real(cdf), dimension(0:3,4) :: p
    complex(c_default_complex) :: amp
    integer :: i
    logical :: exist

    write (u, "(A)")  "* Test output: prc_omega_2"
    write (u, "(A)")  "*   Purpose: create simple processes with OMega"
    write (u, "(A)")  "*            use the prc_omega wrapper for this"
    write (u, "(A)")  "*            and check OMega options"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with two entries, &
         &different options."
    write (u, "(A)")  "* (1) e- e+ -> e- e+   &
         &(all diagrams, no OpenMP, report progress)"
    write (u, "(A)")  "* (2) e- e+ -> e- e+   &
         &(s-channel only, with OpenMP, report progress to file)"

    call lib%init (var_str ("omega2"))
    call os_data%init ()

    model_name = "QED"
    model => null ()
    call prepare_model (model, model_name, vars)

    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("e-"), var_str ("e+")]
    prt_out = prt_in
    restrictions = "3+4~A"

    allocate (entry)
    call entry%init (var_str ("omega2_a"), &
         model, n_in = 2, n_components = 2)

    call omega_make_process_component (entry, 1, &
         model_name, prt_in, prt_out, &
         report_progress=.true.)
    call omega_make_process_component (entry, 2, &
         model_name, prt_in, prt_out, &
         restrictions=restrictions, openmp_support=.true., &
         extra_options=var_str ("-fusion:progress_file omega2.log"))

    call lib%append (entry)

    write (u, "(A)")
    write (u, "(A)")  "* Remove left-over file"
    write (u, "(A)")

    call delete_file ("omega2.log")
    inquire (file="omega2.log", exist=exist)
    write (u, "(1x,A,L1)")  "omega2.log exists = ", exist

    write (u, "(A)")
    write (u, "(A)")  "* Build and load library"

    call lib%configure (os_data)
    call lib%write_makefile (os_data, force = .true., verbose = .false.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)

    write (u, "(A)")
    write (u, "(A)")  "* Check extra output of OMega"
    write (u, "(A)")

    inquire (file="omega2.log", exist=exist)
    write (u, "(1x,A,L1)")  "omega2.log exists = ", exist

    write (u, "(A)")
    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active                 = ", &
         lib%is_active ()
    write (u, "(1x,A,I0)")  "n_processes               = ", &
         lib%get_n_processes ()

    write (u, "(A)")
    write (u, "(A)")  "* Set parameters for omega2_a and initialize:"
    write (u, "(A)")

    call vars%set_rval (var_str ("ee"), 0.3_default)
    call vars%set_rval (var_str ("me"), 0._default)
    call vars%set_rval (var_str ("mmu"), 0._default)
    call vars%set_rval (var_str ("mtau"), 0._default)
    allocate (par (model%get_n_real ()))
    call model%real_parameters_to_c_array (par)

    write (u, "(2x,A,F6.4)")  "ee   = ", par(1)
    write (u, "(2x,A,F6.4)")  "me   = ", par(2)
    write (u, "(2x,A,F6.4)")  "mmu  = ", par(3)
    write (u, "(2x,A,F6.4)")  "mtau = ", par(4)

    call prc1%set_parameters (model)
    call prc2%set_parameters (model)

    write (u, "(A)")
    write (u, "(A)")  "* Constants of omega2_a_i1:"
    write (u, "(A)")

    entry => lib%get_process_def_ptr (var_str ("omega2_a"))
    config => entry%get_component_def_ptr (1)
    call prc1%init (config%get_core_def_ptr (), &
         lib, var_str ("omega2_a"), 1)
    call prc1%get_constants (data, 1)

    write (u, "(1x,A,A)")  "component ID     = ", &
         char (data%id)
    write (u, "(1x,A,L1)") "openmp supported = ", &
         data%openmp_supported
    write (u, "(1x,A,A,A)") "model name       = '", &
         char (data%model_name), "'"

    write (u, "(A)")
    write (u, "(A)")  "* Constants of omega2_a_i2:"
    write (u, "(A)")

    config => entry%get_component_def_ptr (2)
    call prc2%init (config%get_core_def_ptr (), &
         lib, var_str ("omega2_a"), 2)
    call prc2%get_constants (data, 1)

    write (u, "(1x,A,A)")  "component ID     = ", &
         char (data%id)
    write (u, "(1x,A,L1)") "openmp supported = ", &
         data%openmp_supported
    write (u, "(1x,A,A,A)") "model name       = '", &
         char (data%model_name), "'"

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

    write (u, "(A)")
    write (u, "(A)")  "* Compute matrix element:"
    write (u, "(A)")

    select type (driver => prc1%driver)
    type is (omega_driver_t)
       call driver%new_event (p)
       call driver%get_amplitude (1_ci, 6_ci, 1_ci, amp)
       write (u, "(2x,A,1x,E11.4)") "(1) |amp (1, 6, 1)| =", abs (amp)
    end select

    select type (driver => prc2%driver)
    type is (omega_driver_t)
       call driver%new_event (p)
       call driver%get_amplitude (1_ci, 6_ci, 1_ci, amp)
       write (u, "(2x,A,1x,E11.4)") "(2) |amp (1, 6, 1)| =", abs (amp)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics:"
    write (u, "(A)")

    p = reshape ([ &
         1.0_cdf, 0.0_cdf, 0.0_cdf, 1.0_cdf, &
         1.0_cdf, 0.0_cdf, 0.0_cdf,-1.0_cdf, &
         1.0_cdf, sqrt(0.5_cdf), 0.0_cdf, sqrt(0.5_cdf), &
         1.0_cdf,-sqrt(0.5_cdf), 0.0_cdf,-sqrt(0.5_cdf) &
         ], [4,4])
    do i = 1, 4
       write (u, "(2x,A,I0,A,4(1x,F7.4))")  "p", i, " =", p(:,i)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Compute matrix element:"
    write (u, "(A)")

    select type (driver => prc1%driver)
    type is (omega_driver_t)
       call driver%new_event (p)
       call driver%get_amplitude (1_ci, 6_ci, 1_ci, amp)
       write (u, "(2x,A,1x,E11.4)") "(1) |amp (1, 6, 1)| =", abs (amp)
    end select

    select type (driver => prc2%driver)
    type is (omega_driver_t)
       call driver%new_event (p)
       call driver%get_amplitude (1_ci, 6_ci, 1_ci, amp)
       write (u, "(2x,A,1x,E11.4)") "(2) |amp (1, 6, 1)| =", abs (amp)
    end select

    call lib%final ()
    call cleanup_model (model)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_omega_2"

  end subroutine prc_omega_2

  subroutine prc_omega_3 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    type(string_t) :: model_name
    class(model_data_t), pointer :: model
    class(vars_t), pointer :: vars => null ()
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(process_component_def_t), pointer :: config
    type(prc_omega_t) :: prc1
    type(process_constants_t) :: data
    integer, parameter :: cdf = c_default_float
    real(cdf), dimension(:), allocatable :: par
    real(cdf), dimension(0:3,4) :: p
    type(helicity_selection_t) :: helicity_selection
    integer :: i, h

    write (u, "(A)")  "* Test output: prc_omega_3"
    write (u, "(A)")  "*   Purpose: create simple process with OMega"
    write (u, "(A)")  "*            and check helicity selection"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library."
    write (u, "(A)")  "* (1) e- e+ -> e- e+   (all diagrams, no OpenMP)"

    call lib%init (var_str ("omega3"))
    call os_data%init ()

    model_name = "QED"
    model => null ()
    call prepare_model (model, model_name, vars)

    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("e-"), var_str ("e+")]
    prt_out = prt_in

    allocate (entry)
    call entry%init (var_str ("omega3_a"), &
         model, n_in = 2, n_components = 1)

    call omega_make_process_component (entry, 1, &
         model_name, prt_in, prt_out)
    call lib%append (entry)

    write (u, "(A)")
    write (u, "(A)")  "* Build and load library"

    call lib%configure (os_data)
    call lib%write_makefile (os_data, force = .true., verbose = .false.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)

    write (u, "(A)")
    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active                 = ", &
         lib%is_active ()
    write (u, "(1x,A,I0)")  "n_processes               = ", &
         lib%get_n_processes ()

    write (u, "(A)")
    write (u, "(A)")  "* Set parameters for omega3_a and initialize:"
    write (u, "(A)")

    call vars%set_rval (var_str ("ee"), 0.3_default)
    call vars%set_rval (var_str ("me"), 0._default)
    call vars%set_rval (var_str ("mmu"), 0._default)
    call vars%set_rval (var_str ("mtau"), 0._default)
    allocate (par (model%get_n_real ()))
    call model%real_parameters_to_c_array (par)

    write (u, "(2x,A,F6.4)")  "ee   = ", par(1)
    write (u, "(2x,A,F6.4)")  "me   = ", par(2)
    write (u, "(2x,A,F6.4)")  "mmu  = ", par(3)
    write (u, "(2x,A,F6.4)")  "mtau = ", par(4)

    call prc1%set_parameters (model, helicity_selection=helicity_selection)

    write (u, "(A)")
    write (u, "(A)")  "* Helicity states of omega3_a_i1:"
    write (u, "(A)")

    entry => lib%get_process_def_ptr (var_str ("omega3_a"))
    config => entry%get_component_def_ptr (1)
    call prc1%init (config%get_core_def_ptr (), &
         lib, var_str ("omega3_a"), 1)
    call prc1%get_constants (data, 1)

    do i = 1, data%n_hel
       write (u, "(3x,I2,':',4(1x,I2))") i, data%hel_state(:,i)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Initially allowed helicities:"
    write (u, "(A)")

    write (u, "(4x,16(1x,I2))")  [(h, h = 1, data%n_hel)]
    write (u, "(4x)", advance = "no")
    do h = 1, data%n_hel
       write (u, "(2x,L1)", advance = "no")  prc1%is_allowed (1, 1, h, 1)
    end do
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "* Reset helicity selection (cutoff = 4)"
    write (u, "(A)")

    helicity_selection%active = .true.
    helicity_selection%threshold = 1e10_default
    helicity_selection%cutoff = 4
    call helicity_selection%write (u)

    call prc1%set_parameters (model, helicity_selection=helicity_selection)
    call prc1%reset_helicity_selection ()

    write (u, "(A)")
    write (u, "(A)")  "* Allowed helicities:"
    write (u, "(A)")

    write (u, "(4x,16(1x,I2))")  [(h, h = 1, data%n_hel)]
    write (u, "(4x)", advance = "no")
    do h = 1, data%n_hel
       write (u, "(2x,L1)", advance = "no")  prc1%is_allowed (1, 1, h, 1)
    end do
    write (u, "(A)")

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

    write (u, "(A)")
    write (u, "(A)")  "* Compute scattering matrix 5 times"
    write (u, "(A)")

    write (u, "(4x,16(1x,I2))")  [(h, h = 1, data%n_hel)]

    select type (driver => prc1%driver)
    type is (omega_driver_t)
       do i = 1, 5
          call driver%new_event (p)
          write (u, "(2x,I2)", advance = "no")  i
          do h = 1, data%n_hel
             write (u, "(2x,L1)", advance = "no")  prc1%is_allowed (1, 1, h, 1)
          end do
          write (u, "(A)")
       end do
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Reset helicity selection again"
    write (u, "(A)")

    call prc1%activate_parameters ()

    write (u, "(A)")  "* Allowed helicities:"
    write (u, "(A)")

    write (u, "(4x,16(1x,I2))")  [(h, h = 1, data%n_hel)]
    write (u, "(4x)", advance = "no")
    do h = 1, data%n_hel
       write (u, "(2x,L1)", advance = "no")  prc1%is_allowed (1, 1, h, 1)
    end do
    write (u, "(A)")

    call lib%final ()
    call cleanup_model (model)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_omega_3"

  end subroutine prc_omega_3

  subroutine prc_omega_4 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    type(string_t) :: model_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(process_constants_t) :: data
    class(prc_core_driver_t), allocatable :: driver
    integer, parameter :: cdf = c_default_float
    integer, parameter :: ci = c_int
    real(cdf), dimension(8) :: par
    real(cdf), dimension(0:3,4) :: p
    logical(c_bool) :: flag
    complex(c_default_complex) :: amp
    integer :: i
    real(cdf) :: alpha_s

    write (u, "(A)")  "* Test output: prc_omega_4"
    write (u, "(A)")  "*   Purpose: create a QCD process with OMega"
    write (u, "(A)")  "*            and check alpha_s dependence"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("prc_omega_4_lib"))
    call os_data%init ()

    model_name = "QCD"
    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("u"), var_str ("ubar")]
    prt_out = [var_str ("d"), var_str ("dbar")]

    allocate (omega_def_t :: def)
    select type (def)
    type is (omega_def_t)
       call def%init (model_name, prt_in, prt_out, &
            ufo = .false., ovm = .false.)
    end select
    allocate (entry)
    call entry%init (var_str ("prc_omega_4_p"), model_name = model_name, &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("omega"), &
         variant = def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure and compile process"
    write (u, "(A)")
    call lib%configure (os_data)
    call lib%write_makefile (os_data, force = .true., verbose = .false.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)

    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active = ", lib%is_active ()

    write (u, "(A)")
    write (u, "(A)")  "* Set parameters:"
    write (u, "(A)")

    alpha_s = 0.1178_cdf

    par = [alpha_s, &
         0._cdf, 0._cdf, 0._cdf, 0._cdf, 0._cdf, 173.1_cdf, 1.523_cdf]
    write (u, "(2x,A,F8.4)")  "alpha_s = ", par(1)
    write (u, "(2x,A,F8.4)")  "md      = ", par(2)
    write (u, "(2x,A,F8.4)")  "mu      = ", par(3)
    write (u, "(2x,A,F8.4)")  "ms      = ", par(4)
    write (u, "(2x,A,F8.4)")  "mc      = ", par(5)
    write (u, "(2x,A,F8.4)")  "mb      = ", par(6)
    write (u, "(2x,A,F8.4)")  "mtop    = ", par(7)
    write (u, "(2x,A,F8.4)")  "wtop    = ", par(8)

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics:"
    write (u, "(A)")

    p = reshape ([ &
         100.0_cdf, 0.0_cdf, 0.0_cdf, 100.0_cdf, &
         100.0_cdf, 0.0_cdf, 0.0_cdf,-100.0_cdf, &
         100.0_cdf, 100.0_cdf, 0.0_cdf, 0.0_cdf, &
         100.0_cdf,-100.0_cdf, 0.0_cdf, 0.0_cdf &
         ], [4,4])
    do i = 1, 4
       write (u, "(2x,A,I0,A,4(1x,F7.1))")  "p", i, " =", p(:,i)
    end do

    call lib%connect_process (var_str ("prc_omega_4_p"), 1, data, driver)

    select type (driver)
    type is (omega_driver_t)
       call driver%init (par, 0)

       write (u, "(A)")
       write (u, "(A)")  "* Compute matrix element:"
       write (u, "(A)")

       call driver%new_event (p)

       call driver%is_allowed (1_ci, 6_ci, 1_ci, flag)
       write (u, "(1x,A,L1)") "is_allowed (1, 6, 1) = ", flag

       call driver%get_amplitude (1_ci, 6_ci, 1_ci, amp)
       write (u, "(1x,A,1x,E11.4)") "|amp (1, 6, 1)| =", abs (amp)

       write (u, "(A)")
       write (u, "(A)")  "* Double alpha_s and compute matrix element again:"
       write (u, "(A)")

       call driver%update_alpha_s (2 * alpha_s)
       call driver%new_event (p)

       call driver%is_allowed (1_ci, 6_ci, 1_ci, flag)
       write (u, "(1x,A,L1)") "is_allowed (1, 6, 1) = ", flag

       call driver%get_amplitude (1_ci, 6_ci, 1_ci, amp)
       write (u, "(1x,A,1x,E11.4)") "|amp (1, 6, 1)| =", abs (amp)
    end select

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_omega_4"

  end subroutine prc_omega_4

  subroutine prc_omega_5 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_component_def_t), pointer :: cdef_ptr
    class(prc_core_def_t), pointer :: def_ptr
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(string_t) :: model_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(qcd_t) :: qcd
    class(prc_core_t), allocatable :: core
    class(prc_core_state_t), allocatable :: core_state
    type(vector4_t), dimension(4) :: p
    complex(default) :: amp
    real(default) :: ren_scale
    real(default), allocatable :: alpha_qcd_forced
    integer :: i

    write (u, "(A)")  "* Test output: prc_omega_5"
    write (u, "(A)")  "*   Purpose: create a QCD process with OMega"
    write (u, "(A)")  "*            and check alpha_s dependence"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("prc_omega_5_lib"))
    call os_data%init ()

    model_name = "QCD"
    model => null ()
    call prepare_model (model, model_name)

    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("u"), var_str ("ubar")]
    prt_out = [var_str ("d"), var_str ("dbar")]

    allocate (omega_def_t :: def)
    select type (def)
    type is (omega_def_t)
       call def%init (model_name, prt_in, prt_out, &
            ufo = .false., ovm = .false.)
    end select
    allocate (entry)
    call entry%init (var_str ("prc_omega_5_p"), model_name = model_name, &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("omega"), &
         variant = def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure and compile process"
    write (u, "(A)")
    call lib%configure (os_data)
    call lib%write_makefile (os_data, force = .true., verbose = .false.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)

    write (u, "(A)")  "* Probe library API"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active = ", lib%is_active ()

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics"
    write (u, "(A)")

    p(1) = vector4_moving (100._default, 100._default, 3)
    p(2) = vector4_moving (100._default,-100._default, 3)
    p(3) = vector4_moving (100._default, 100._default, 1)
    p(4) = vector4_moving (100._default,-100._default, 1)
    do i = 1, 4
       call vector4_write (p(i), u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Setup QCD data"
    write (u, "(A)")

    allocate (alpha_qcd_from_scale_t :: qcd%alpha)

    write (u, "(A)")  "* Setup process core"
    write (u, "(A)")

    allocate (prc_omega_t :: core)
    entry => lib%get_process_def_ptr (var_str ("prc_omega_5_p"))
    cdef_ptr => entry%get_component_def_ptr (1)
    def_ptr => cdef_ptr%get_core_def_ptr ()

    select type (core)
    type is (prc_omega_t)
       call core%allocate_workspace (core_state)
       call core%set_parameters (model, qcd = qcd)
       call core%init (def_ptr, lib, var_str ("prc_omega_5_p"), 1)
       call core%write (u)

       write (u, "(A)")
       write (u, "(A)")  "* Compute matrix element"
       write (u, "(A)")

       ren_scale = 100
       write (u, "(1x,A,F4.0)")  "renormalization scale = ", ren_scale

       amp = core%compute_amplitude &
            (1, p, 1, 6, 1, 100._default, ren_scale, alpha_qcd_forced)

       write (u, "(1x,A,1x,E11.4)") "|amp (1, 6, 1)| =", abs (amp)

       write (u, "(A)")
       write (u, "(A)")  "* Modify renormalization scale and &
            &compute matrix element again"
       write (u, "(A)")

       ren_scale = 200
       write (u, "(1x,A,F4.0)")  "renormalization scale = ", ren_scale

       amp = core%compute_amplitude &
            (1, p, 1, 6, 1, 100._default, ren_scale, alpha_qcd_forced)

       write (u, "(1x,A,1x,E11.4)") "|amp (1, 6, 1)| =", abs (amp)

       write (u, "(A)")
       write (u, "(A)")  "* Set alpha(QCD) directly and &
            &compute matrix element again"
       write (u, "(A)")

       allocate (alpha_qcd_forced, source = 0.1_default)
       write (u, "(1x,A,F6.4)")  "alpha_qcd = ", alpha_qcd_forced

       amp = core%compute_amplitude &
            (1, p, 1, 6, 1, 100._default, ren_scale, alpha_qcd_forced)

       write (u, "(1x,A,1x,E11.4)") "|amp (1, 6, 1)| =", abs (amp)

    end select

    call lib%final ()
    call cleanup_model (model)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_omega_5"

  end subroutine prc_omega_5

  subroutine prc_omega_6 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    type(string_t) :: model_name
    class(model_data_t), pointer :: model
    class(vars_t), pointer :: vars
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(string_t) :: restrictions
    type(process_component_def_t), pointer :: config
    type(prc_omega_t) :: prc1, prc2
    type(process_constants_t) :: data
    integer, parameter :: cdf = c_default_float
    integer, parameter :: ci = c_int
    real(cdf), dimension(:), allocatable :: par
    real(cdf), dimension(0:3,4) :: p
    complex(c_default_complex) :: amp
    integer :: i
    logical :: exist

    write (u, "(A)")  "* Test output: prc_omega_6"
    write (u, "(A)")  "*   Purpose: create simple process with OMega / UFO file"
    write (u, "(A)")

    call os_data%init ()

    model_name = "SM"
    model => null ()

    os_data%whizard_modelpath_ufo = "../models/UFO"

    write (u, "(A)")  "* Create process library entry"
    write (u, "(A)")

    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("e-"), var_str ("e+")]
    prt_out = prt_in
    restrictions = "3+4~A"

    allocate (entry)
    call entry%init (var_str ("omega_6_a"), &
         model_name = model_name, n_in = 2, n_components = 1)

    call omega_make_process_component (entry, 1, &
         model_name, prt_in, prt_out, &
         ufo=.true., ufo_path=os_data%whizard_modelpath_ufo, &
         report_progress=.true.)

    call entry%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Build and load library"

    call lib%init (var_str ("omega_6"))
    call lib%append (entry)

    call lib%configure (os_data)
    call lib%write_makefile (os_data, force = .true., verbose = .false.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)

    write (u, "(A)")
    write (u, "(A)")  "* Probe library API:"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active                 = ", &
         lib%is_active ()
    write (u, "(1x,A,I0)")  "n_processes               = ", &
         lib%get_n_processes ()

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_omega_6"

  end subroutine prc_omega_6

  subroutine prc_omega_diags_1 (u)
    integer, intent(in) :: u
    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    type(string_t) :: model_name
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(string_t) :: diags_file, pdf_file, ps_file
    logical :: exist, exist_pdf, exist_ps
    integer :: iostat, u_diags
    character(128) :: buffer

    write (u, "(A)")  "* Test output: prc_omega_diags_1"
    write (u, "(A)")  "*   Purpose: generate Feynman diagrams"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process library with one entry"
    write (u, "(A)")
    call lib%init (var_str ("prc_omega_diags_1_lib"))
    call os_data%init ()

    model_name = "SM"

    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("u"), var_str ("ubar")]
    prt_out = [var_str ("d"), var_str ("dbar")]

    allocate (omega_def_t :: def)
    select type (def)
    type is (omega_def_t)
       call def%init (model_name, prt_in, prt_out, &
            ufo = .false., ovm = .false., &
            diags = .true., diags_color = .true.)
    end select
    allocate (entry)
    call entry%init (var_str ("prc_omega_diags_1_p"), model_name = model_name, &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("omega"), &
         variant = def)
    call lib%append (entry)

    write (u, "(A)")  "* Configure and compile process"
    write (u, "(A)")  "    and generate diagrams"
    write (u, "(A)")
    call lib%configure (os_data)
    call lib%write_makefile &
         (os_data, force = .true., verbose = .false., testflag = .true.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)

    write (u, "(A)")  "* Probe library API"
    write (u, "(A)")

    write (u, "(1x,A,L1)")  "is active = ", lib%is_active ()

    write (u, "(A)")  "* Check produced diagram files"
    write (u, "(A)")

    diags_file = "prc_omega_diags_1_p_i1_diags.tex"
    ps_file  = "prc_omega_diags_1_p_i1_diags.ps"
    pdf_file = "prc_omega_diags_1_p_i1_diags.pdf"
    inquire (file = char (diags_file), exist = exist)
    if (exist) then
       u_diags = free_unit ()
       open (u_diags, file = char (diags_file), action = "read", status = "old")
       iostat = 0
       do while (iostat == 0)
          read (u_diags, "(A)", iostat = iostat)  buffer
          if (iostat == 0)  write (u, "(A)")  trim (buffer)
       end do
       close (u_diags)
    else
       write (u, "(A)")  "[Feynman diagrams LaTeX file is missing]"
    end if
    inquire (file = char (ps_file), exist = exist_ps)
    if (exist_ps) then
       write (u, "(A)")  "[Feynman diagrams Postscript file exists and is nonempty]"
    else
       write (u, "(A)")  "[Feynman diagrams Postscript file is missing/non-regular]"
    end if
    inquire (file = char (pdf_file), exist = exist_pdf)
    if (exist_pdf) then
       write (u, "(A)")  "[Feynman diagrams PDF file exists and is nonempty]"
    else
       write (u, "(A)")  "[Feynman diagrams PDF file is missing/non-regular]"
    end if

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    call lib%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prc_omega_diags_1"

  end subroutine prc_omega_diags_1


end module prc_omega_uti
