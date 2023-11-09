! compare_lib_recola.f90 --
! compare_lib_recola.f90 -- compare two O'Mega versions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     Christian Speckner <cnspeckn@googlemail.com>
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

module compare_lib_recola
  ! use ieee_arithmetic
  use kinds
  use constants
  use tao_random_numbers
  use recola
  use omega95
  use omega_interface
  use omega_testtools
  use compare_lib
  use parameters_sm_higgs_recola
  implicit none
  private
  public :: check_recola
  public :: set_recola_parameters_from_omega
contains

  subroutine set_recola_parameters_from_omega
    real(kind=double), dimension(size(mass)) :: m
    real(kind=double), dimension(size(width)) :: w
    m = real (mass, kind=double)
    w = real (width, kind=double)
    call set_output_file_rcl ('recola.log')
    call set_print_level_squared_amplitude_rcl (0)
    call set_on_shell_scheme_rcl
    call set_pole_mass_down_rcl      (m( 1))
    call set_pole_mass_up_rcl        (m( 2))
    call set_pole_mass_charm_rcl     (m( 3), w( 3))
    call set_pole_mass_strange_rcl   (m( 4))
    call set_pole_mass_bottom_rcl    (m( 5), w( 5))
    call set_pole_mass_top_rcl       (m( 6), w( 6))
    call set_pole_mass_electron_rcl  (m(11))
    call set_pole_mass_muon_rcl      (m(13), w(13))
    call set_pole_mass_tau_rcl       (m(15), w(15))
    call set_pole_mass_z_rcl         (m(23), w(23))
    call set_pole_mass_w_rcl         (m(24), w(24))
    call set_pole_mass_h_rcl         (m(25), w(25))
    ! call switchoff_coupling3_rcl ('H','mu+','mu-')
    ! call unset_light_muon_rcl
  end subroutine set_recola_parameters_from_omega

  elemental function ieee_is_nan (x) result (yorn)
    logical :: yorn
    real (kind=default), intent(in) :: x
    yorn = (x /= x)
  end function ieee_is_nan

  elemental function ieee_is_nan_double (x) result (yorn)
    logical :: yorn
    real (kind=double), intent(in) :: x
    yorn = (x /= x)
  end function ieee_is_nan_double

  elemental function mass_of_pdg (pdg) result (m)
    integer, intent(in) :: pdg
    real(kind=default) :: m
    m = mass(abs(pdg))
  end function mass_of_pdg

  subroutine check_recola (proc, process, roots, threshold, n, &
                    failures, attempts, seed, abs_threshold)
    type(omega_procedures), intent(in) :: proc
    character(*), intent(in) :: process
    real(kind=default), intent(in) :: roots, threshold
    integer, intent(in) :: n
    integer, intent(out) :: failures, attempts
    integer, intent(in), optional :: seed
    real(kind=default), intent(in), optional :: abs_threshold
    logical :: match, passed, error
    integer :: i
    real(kind=default) :: asq_sum_omega
    real(kind=double) :: asq_sum_recola
    real(kind=default) :: sum_omega, sum_recola
    character(len=80) :: msg
    integer, dimension(:,:), allocatable :: flavors
    real(kind=default), dimension(:,:), allocatable :: p
    real(kind=default), dimension(:), allocatable :: m
    logical, parameter :: UNWEIGHTED = .true.
    logical, parameter :: QUIET = .true.
    real(kind=default) :: weight

    call define_process_rcl (1, trim(process), 'LO')
    call generate_processes_rcl

    call omega_flavor_states (proc, flavors)
    allocate (m(size(flavors,dim=1)))
    allocate (p(0:3,size(flavors,dim=1)))
    m = mass_of_pdg (flavors(:,1))
  
    call beams (roots, m(1), m(2), p(:,1), p(:,2))
    call proc%reset_helicity_selection (-1.0_default, -1)

    failures = 0
    attempts = 0
    sum_omega = 0
    sum_recola = 0

    if (present (seed)) then
       call tao_random_seed (seed)
    end if

    do i = 1, n
       attempts = attempts + 1
       if (size(p,dim=2) > 3) then
          call rambo (roots, m(3:), p(:,3:), weight, UNWEIGHTED)
          call rambo_check (roots, m(3:), p(:,3:), quiet=.true.)
       else
          p(:,3) = p(:,1) + p(:,2)
       end if

       call compute_process_rcl (1, real (p, kind=double), 'LO')
       call get_squared_amplitude_rcl (1, 0, 'LO', asq_sum_recola)

       call omega_squared_matrix_element (proc, p, asq_sum_omega, error)
       if (error) then
          write (*, "(1X,'evt=',I5,': ', A)") i, "O'Mega failed"
          failures = failures + 1
          cycle
       end if

       passed = .true.
       if (ieee_is_nan_double (asq_sum_recola)) then
          write (*, "(1X,'evt=',I5,': ', A)") i, "squared recola amplitude NaN"
          passed = .false.
       end if
       if (ieee_is_nan (asq_sum_omega)) then
          write (*, "(1X,'evt=',I5,': ', A)") i, "squared O'Mega amplitude NaN"
          passed = .false.
       end if
       write (msg, "(1X,'evt=',I5)") i
       call expect (real (asq_sum_recola, kind=default), &
                    asq_sum_omega, trim(msg), passed, &
                    quiet=.true., threshold=threshold, &
                    abs_threshold=abs_threshold)

       if (.not.passed) then
          failures = failures + 1
          cycle
       end if

       sum_omega = sum_omega + asq_sum_omega
       sum_recola = sum_recola + asq_sum_recola

    end do

    call reset_recola_rcl
    deallocate (m, p)

    print *, 'Summed results: '
    print *, 'omega, recola =    ', sum_omega, sum_recola

  end subroutine check_recola

end module compare_lib_recola
