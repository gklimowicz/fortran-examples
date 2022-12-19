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

module prc_recola_uti

  use recola_wrapper !NODEP!

  use, intrinsic :: iso_c_binding !NODEP!
  use kinds
  use iso_varying_string, string_t => varying_string

  use constants
  use format_utils, only: write_separator
  use numeric_utils, only: assert_equal
  use os_interface
  use particle_specifiers, only: new_prt_spec
  use prc_core_def
  use process_constants
  use process_libraries
  use prc_core
  use prc_omega

  implicit none
  private

  public :: prc_recola_1
  public :: prc_recola_2

contains

  function get_omega_parameter_array () result (par)
    real(default), dimension(25) :: par
    par = zero

    par(1) = 1.16637d-5 ! gf
    par(2) = 91.153480619182744_default ! mZ
    par(3) = 80.357973609877547_default ! mW
    par(4) = 125._default ! mH
    par(5) = rclwrap_get_alpha_s () ! alpha_s
    par(12) = 173.2_default ! mt
    par(14) = 2.4942663787728243_default ! wZ
    par(15) = 2.0842989982782196_default ! wW
    par(22) = one / sqrt (sqrt (two) * par(1)) ! par%v - Higgs expectation value
    par(23) = par(3) / par(2) ! par%cw
    par(24) = sqrt (one - par(23)**2) ! par%sw
    par(25) = two * par(24) * par(3) / par(22)
  end function get_omega_parameter_array


  subroutine prc_recola_1 (u)
    integer, intent(in) :: u
    real(double) :: p(0:3,1:4)
    real(double) :: sqrts = 500._double
    real(double) :: m_e = 0._double
    real(double) :: m_mu = 0._double
    real(double) :: p_x_out, p_y_out, p_z_out, p_z_in
    integer      :: h_e_p, h_e_m, h_mu_p, h_mu_m, counter
    real(double) :: sqme
    integer :: i
    integer, dimension(:), allocatable :: col_recola, hel_recola
    complex(double) :: amp_recola
    complex(default) :: amp_recola_default

    real(default), parameter :: ee = 0.3 !!! Electromagnetic coupling

    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(os_data_t) :: os_data
    type(process_constants_t) :: data
    class(prc_core_driver_t), allocatable :: driver
    complex(default) :: amp

    integer, dimension(:,:), allocatable :: helicities



    write (u, "(A)") "* Test output: prc_recola_1"
    write (u, "(A)") "* Purpose: Test interface to RECOLA and compare matrix elements with O'Mega"
    write (u, "(A)")

    p_z_in = sqrt ((sqrts / 2)**2 - m_e**2)
    p_z_out = 0._double
    p_y_out = sqrts / 10._default
    p_x_out = sqrt ((sqrts / 2)**2 - p_y_out**2 - p_z_out**2 - m_mu**2)
    p(:,1) = [sqrts / 2,  0._double,  0._double,  p_z_in]
    p(:,2) = [sqrts / 2,  0._double,  0._double, -p_z_in]
    p(:,3) = [sqrts / 2,  p_x_out,  p_y_out,  p_z_out]
    p(:,4) = [sqrts / 2, -p_x_out, -p_y_out, -p_z_out]

    write (u, "(A)") "Use phase-space point: "
    do i = 1, 4
       write (u, "(4(F12.3,1x))") p(:,1)
    end do
    write (u, "(A)")
    call write_separator (u)
    write (u, "(A)")
    write (u, "(A)") "* RECOLA: Evaluate process"
    counter  = 1
    call rclwrap_request_generate_processes ()
    write (u, "(A)") "*  RECOLA: Define process e+ e- -> mu+ mu- at leading order"
    call rclwrap_add_process (counter, var_str ('e+ e- -> mu+ mu-'), var_str ('LO'))
    call rclwrap_define_processes ()
    write (u, "(A)") "* RECOLA: generate process"
    call rclwrap_generate_processes ()
    call rclwrap_compute_process (1, p, 'LO')
    call rclwrap_get_helicity_configurations (1, helicities)
    allocate (hel_recola (4), col_recola (4))
    col_recola = [0,0,0,0]

    write (u, "(A)") "* Setting up Omega to compute the same amplitude"

    call lib%init (var_str ("omega1"))
    allocate (prt_in (2), prt_out (2))
    prt_in = [var_str ("e+"), var_str ("e-")]
    prt_out = [var_str ("mu+"), var_str ("mu-")]

    allocate (omega_def_t :: def)
    select type (def)
    type is (omega_def_t)
       call def%init (var_str ("SM"), prt_in, prt_out, &
            ufo = .false., ovm = .false., cms_scheme = .true.)
    end select

    allocate (entry)
    call entry%init (var_str ("omega1_a"), model_name = var_str ("SM"), &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = 2, &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("omega"), &
         variant = def)
    call lib%append (entry)

    call os_data%init ()
    call lib%configure (os_data)
    call lib%write_makefile (os_data, force = .true., verbose = .false.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)
    call lib%connect_process (var_str ("omega1_a"), 1, data, driver)

    select type (driver)
    type is (omega_driver_t)
       call driver%init (get_omega_parameter_array (), 3)
       call driver%new_event (real(p, kind =  default))
       do i = 1, 6
          call rclwrap_get_amplitude (1, 0, 'LO', col_recola, helicities (:,i), amp_recola)
       end do
       do i = 1, 16
           call rclwrap_get_amplitude (1, 0, 'LO', col_recola, data%hel_state (:,i), amp_recola)
           amp_recola = amp_recola * cmplx (0, -1, double)
           amp_recola_default = amp_recola
           call driver%get_amplitude (1, i, 1, amp)
           write(u,"(A,4(I2),A)") "Helicity: [",data%hel_state (:,i),"]"
           call assert_equal (u, amp, amp_recola_default, rel_smallness = 1.E-7_default)
       end do

    end select

    call rclwrap_reset_recola ()

    write (u, "(A)")
    write (u, "(A)") "* End of test output: prc_recola_1"

  end subroutine prc_recola_1

  subroutine prc_recola_2 (u)
    integer, intent(in) :: u
    real(double) :: p(0:3,1:5)
    real(double) :: sqrts = 700._double
    real(double) :: m_e = 0._double
    real(double) :: m_mu = 0._double
    real(double) :: p_x_out, p_y_out, p_z_out, p_z_in
    real(double) :: sqme
    integer :: i
    integer, dimension(:), allocatable :: col_recola, hel_recola
    integer, dimension(:,:), allocatable :: helicities
    complex(double) :: amp_recola
    complex(default) :: amp_recola_default

    real(default), parameter :: ee = 0.3 !!! Electromagnetic coupling

    type(process_library_t) :: lib
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    type(os_data_t) :: os_data
    type(process_constants_t) :: data
    class(prc_core_driver_t), allocatable :: driver
    complex(default) :: amp
    integer :: n_allowed

    write (u, "(A)") "* Test output: prc_recola_2"
    write (u, "(A)") "* Purpose: Test interface to RECOLA and compare matrix elements with O'Mega for 2->3 process"
    write (u, "(A)")

    p_z_in = sqrt ((sqrts / 2)**2 - m_e**2)
    p(:,1) = [sqrts / 2,  0._double,  0._double,  p_z_in]
    p(:,2) = [sqrts / 2,  0._double,  0._double, -p_z_in]
    p(:,3) = [243.49323116_double, -141.69619338_double, -108.30640321_double,  165.77353656_double]
    p(:,4) = [337.53250628_double,  143.95931207_double,  110.19717026_double, -284.71124482_double]
    p(:,5) = [118.97426257_double, -2.2631186860_double, -1.8907670459_double,  118.93770827_double]

    write (u, "(A)") "Use phase-space point: "
    do i = 1, 5
       write (u, "(4(F12.3,1x))") p(:,1)
    end do
    write (u, "(A)")
    call write_separator (u)
    write (u, "(A)")
    write (u, "(A)") "* RECOLA: Evaluate process"
    call rclwrap_request_generate_processes ()
    write (u, "(A)") "*  RECOLA: Define process e+ e- -> mu+ mu- A at leading order"
    call rclwrap_add_process (2, var_str ('e+ e- -> mu+ mu- A'), var_str ('LO'))
    call rclwrap_define_processes ()
    write (u, "(A)") "* RECOLA: generate process"
    call rclwrap_generate_processes ()
    call rclwrap_compute_process (2, p, 'LO')
    call rclwrap_get_helicity_configurations (2, helicities)

    allocate (hel_recola (5), col_recola (5))
    col_recola = [0,0,0,0,0]


    write (u, "(A)") "* Setting up Omega to compute the same amplitude"

    call lib%init (var_str ("omega2"))
    allocate (prt_in (2), prt_out (3))
    prt_in = [var_str ("e+"), var_str ("e-")]
    prt_out = [var_str ("mu+"), var_str ("mu-"), var_str("A")]

    allocate (omega_def_t :: def)
    select type (def)
    type is (omega_def_t)
       call def%init (var_str ("SM"), prt_in, prt_out, &
            ufo = .false., ovm = .false.)
    end select

    allocate (entry)
    call entry%init (var_str ("omega2_a"), model_name = var_str ("SM"), &
       n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = 3, &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("omega"), &
         variant = def)
    call lib%append (entry)

    call os_data%init ()
    call lib%configure (os_data)
    call lib%write_makefile (os_data, force = .true., verbose = .false.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)
    call lib%connect_process (var_str ("omega2_a"), 1, data, driver)


    select type (driver)
    type is (omega_driver_t)
       call driver%init (get_omega_parameter_array (), 3)
       call driver%new_event (real(p, kind = default))
       do i = 1, 32
           call rclwrap_get_amplitude &
                (2, 0, 'LO', col_recola, data%hel_state (:,i), amp_recola)
           if (data%hel_state(3,i) * data%hel_state(4,i) * &
                data%hel_state(5,i) == -1) then
              amp_recola = amp_recola * cmplx (0, -1, double)
           else
              amp_recola = amp_recola * cmplx (0, 1, double)
           end if
           amp_recola_default = amp_recola
           call driver%get_amplitude (1, i, 1, amp)
           write(u,"(A,5(I2),A)") "Helicity: [", data%hel_state (:,i),"]"
           write(u,"(A,2(F12.7,1x),A,2(F12.7,1x))") "RECOLA:", &
                amp_recola,", O'MEGA:", amp
           call assert_equal &
                (u, amp, amp_recola_default, rel_smallness = 1.E-6_default)
       end do

    end select

    call rclwrap_reset_recola ()

    write (u, "(A)")
    write (u, "(A)") "* End of test output: prc_recola_2"

  end subroutine prc_recola_2

end module prc_recola_uti

