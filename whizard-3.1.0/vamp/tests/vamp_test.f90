! vamp_test.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vamp_test_functions
  use kinds
  use constants, only: PI
  use coordinates
  use vamp, only: vamp_grid, vamp_multi_channel
  use vamp, only: vamp_data_t
  implicit none
  private
  public :: f, j, phi, ihp, w
  public :: lorentzian
  private :: lorentzian_normalized
  real(kind=default), public :: width
contains
  pure function lorentzian_normalized (x, x0, x1, x2, a) result (f)
    real(kind=default), intent(in) :: x, x0, x1, x2, a
    real(kind=default) :: f
    if (x1 <= x .and. x <= x2) then
       f = 1 / ((x - x0)**2 + a**2) &
            * a / (atan2 (x2 - x0, a) - atan2 (x1 - x0, a))
    else
       f = 0
    end if
  end function lorentzian_normalized
  pure function lorentzian (x, x0, x1, x2, r0, a) result (f)
    real(kind=default), dimension(:), intent(in) :: x, x0, x1, x2
    real(kind=default), intent(in) :: r0, a
    real(kind=default) :: f
    real(kind=default) :: r, r1, r2
    integer :: n
    n = size (x)
    if (n > 1) then
       r = sqrt (dot_product (x-x0, x-x0))
       r1 = 0.4_default
       r2 = min (minval (x2-x0), minval (x0-x1))
       if (r1 <= r .and. r <= r2) then
          f = lorentzian_normalized (r, r0, r1, r2, a) * r**(1-n) / surface (n)
       else
          f = 0
       end if
    else
       f = lorentzian_normalized (x(1), x0(1), x1(1), x2(1), a)
    endif
  end function lorentzian
  pure function f (x, data, weights, channel, grids) result (f_x)
    real(kind=default), dimension(:), intent(in) :: x
    class(vamp_data_t), intent(in) :: data
    real(kind=default), dimension(:), intent(in), optional :: weights
    integer, intent(in), optional :: channel
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    real(kind=default) :: f_x
    real(kind=default), dimension(size(x)) :: minus_one, plus_one, zero, w_i, f_i
    integer :: n, i
    n = size(x)
    minus_one = -1
    zero = 0
    plus_one = 1
    w_i = 1
    do i = 1, n
       if (all (abs (x(i+1:)) <= 1)) then
          f_i = lorentzian (x(1:i), zero(1:i), minus_one(1:i), plus_one(1:i), &
                            0.7_default, width) &
               / 2.0_default**(n-i)
       else
          f_i = 0
       end if
    end do
    f_x = dot_product (w_i, f_i) / sum (w_i)
  end function f
  pure function phi (xi, channel) result (x)
    real(kind=default), dimension(:), intent(in) :: xi
    integer, intent(in) :: channel
    real(kind=default), dimension(size(xi)) :: x
    real(kind=default) :: r
    real(kind=default), dimension(0) :: dummy
    integer :: n
    n = size(x)
    if (channel == 1) then
       x = xi
    else if (channel == 2) then
       r = (xi(1) + 1) / 2 * sqrt (2.0_default)
       x(1:2) = spherical_cos_to_cartesian (r, PI * xi(2), dummy)
       x(3:) = xi(3:)
    else if (channel < n) then
       r = (xi(1) + 1) / 2 * sqrt (real (channel, kind=default))
       x(1:channel) = spherical_cos_to_cartesian (r, PI * xi(2), xi(3:channel))
       x(channel+1:) = xi(channel+1:)
    else if (channel == n) then
       r = (xi(1) + 1) / 2 * sqrt (real (channel, kind=default))
       x = spherical_cos_to_cartesian (r, PI * xi(2), xi(3:))
    else
       x = 0
    end if
  end function phi
  pure function ihp (x, channel) result (xi)
    real(kind=default), dimension(:), intent(in) :: x
    integer, intent(in) :: channel
    real(kind=default), dimension(size(x)) :: xi
    real(kind=default) :: r, phi
    integer :: n
    n = size(x)
    if (channel == 1) then
       xi = x
    else if (channel == 2) then
       call cartesian_to_spherical_cos (x(1:2), r, phi)
       xi(1) = 2 * r / sqrt (2.0_default) - 1
       xi(2) = phi / PI
       xi(3:) = x(3:)
    else if (channel < n) then
       call cartesian_to_spherical_cos (x(1:channel), r, phi, xi(3:channel))
       xi(1) = 2 * r / sqrt (real (channel, kind=default)) - 1
       xi(2) = phi / PI
       xi(channel+1:) = x(channel+1:)
    else if (channel == n) then
       call cartesian_to_spherical_cos (x, r, phi, xi(3:))
       xi(1) = 2 * r / sqrt (real (channel, kind=default)) - 1
       xi(2) = phi / PI
    else
       xi = 0
    end if
  end function ihp
  pure function j (x, data, channel) result (j_x)
    real(kind=default), dimension(:), intent(in) :: x
    class(vamp_data_t), intent(in) :: data
    integer, intent(in) :: channel
    real(kind=default) :: j_x
    if (channel == 1) then
       j_x = 1
    else if (channel > 1) then
       j_x = 2 / sqrt (real (channel, kind=default)) !: $1/|\mathrm{d}r/\mathrm{d}\xi_1|$
       j_x = j_x / PI                                !: $1/|\mathrm{d}\phi/\mathrm{d}\xi_2|$
       j_x = j_x * cartesian_to_spherical_cos_j (x(1:channel))
    else
       j_x = 0
    end if
  end function j
  function w (x, data, weights, channel, grids) result (w_x)
    real(kind=default), dimension(:), intent(in) :: x
    class(vamp_data_t), intent(in) :: data
    real(kind=default), dimension(:), intent(in), optional :: weights
    integer, intent(in), optional :: channel
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    real(kind=default) :: w_x
    w_x = vamp_multi_channel (f, data, phi, ihp, j, x, weights, channel, grids)
  end function w
end module vamp_test_functions
module vamp_tests
  use kinds
  use exceptions
  use histograms
  use tao_random_numbers
  use coordinates
  use vamp
  use vamp_test_functions !NODEP!
  implicit none
  private
  ! public :: check_jacobians, check_inverses, check_inverses3
  public :: check_inverses, check_inverses3
  public :: single_channel, multi_channel
  public :: print_results
contains
  subroutine check_inverses (rng, region, weights, samples)
    type(tao_random_state), intent(inout) :: rng
    real(kind=default), dimension(:,:), intent(in) :: region
    real(kind=default), dimension(:), intent(in) :: weights
    integer, intent(in) :: samples
    real(kind=default), dimension(size(region,dim=2)) :: x1, x2, x_dx
    real(kind=default) :: dx, dx_max
    integer :: ch, i
    dx_max = 0
    x_dx = 0
    do ch = 1, size(weights)
       do i = 1, samples
          call tao_random_number (rng, x1)
          x2 = ihp (phi (x1, ch), ch)
          dx = sqrt (dot_product (x1-x2, x1-x2))
          if (dx > dx_max) then
             dx_max = dx
             x_dx = x1
          end if
       end do
       print *, "channel", ch, ": |x-x|=", real(dx), ", @x=", real (x_dx)
    end do
  end subroutine check_inverses
  subroutine check_inverses3 (rng, region, samples)
    type(tao_random_state), intent(inout) :: rng
    real(kind=default), dimension(:,:), intent(in) :: region
    integer, intent(in) :: samples
    real(kind=default), dimension(size(region,dim=2)) :: x1, x2, x_dx, x_dj
    real(kind=default) :: r, phi, jac, caj, dx, dx_max, dj, dj_max
    real(kind=default), dimension(size(x1)-2) :: cos_theta
    integer :: i
    dx_max = 0
    x_dx = 0
    dj_max = 0
    x_dj = 0
    do i = 1, samples
       call tao_random_number (rng, x1)
       call cartesian_to_spherical_cos_2 (x1, r, phi, cos_theta, jac)
       call spherical_cos_to_cartesian_2 (r, phi, cos_theta, x2, caj)
       dx = sqrt (dot_product (x1-x2, x1-x2))
       dj = jac*caj - 1
       if (dx > dx_max) then
          dx_max = dx
          x_dx = x1
       end if
       if (dj > dj_max) then
          dj_max = dj
          x_dj = x1
       end if
    end do
    print *, "channel 3 : j*j-1=", real(dj), ", @x=", real (x_dj)
    print *, "channel 3 : |x-x|=", real(dx), ", @x=", real (x_dx)
  end subroutine check_inverses3
  subroutine single_channel (rng, region, samples, iterations, &
       integral, standard_dev, chi_squared)
    type(tao_random_state), intent(inout) :: rng
    real(kind=default), dimension(:,:), intent(in) :: region
    integer, dimension(:), intent(in) :: samples, iterations
    real(kind=default), intent(out) :: integral, standard_dev, chi_squared
    type(vamp_grid) :: gr
    type(vamp_history), dimension(iterations(1)+iterations(2)) :: history
    call vamp_create_history (history)
    call vamp_create_grid (gr, region, samples(1))
    call vamp_sample_grid (rng, gr, f, NO_DATA, iterations(1), history = history)
    call vamp_discard_integral (gr, samples(2))
    call vamp_sample_grid &
         (rng, gr, f, NO_DATA, iterations(2), &
          integral, standard_dev, chi_squared, &
          history = history(iterations(1)+1:))
    call vamp_write_grid (gr, "vamp_test.grid")
    call vamp_delete_grid (gr)
    call vamp_print_history (history, "single")
    call vamp_delete_history (history)
  end subroutine single_channel
  subroutine multi_channel (rng, region, weights, samples, iterations, powers, &
       integral, standard_dev, chi_squared)
    type(tao_random_state), intent(inout) :: rng
    real(kind=default), dimension(:,:), intent(in) :: region
    real(kind=default), dimension(:), intent(inout) :: weights
    integer, dimension(:), intent(in) :: samples, iterations
    real(kind=default), dimension(:), intent(in) :: powers
    real(kind=default), intent(out) :: integral, standard_dev, chi_squared
    type(vamp_grids) :: grs
    type(vamp_history), dimension(iterations(1)+iterations(2)+size(powers)-1) :: &
         history
    type(vamp_history), dimension(size(history),size(weights)) :: histories
    integer :: it, nit
    nit = size (powers)
    call vamp_create_history (history)
    call vamp_create_history (histories)
    call vamp_create_grids (grs, region, samples(1), weights)
    call vamp_sample_grids (rng, grs, w, NO_DATA, iterations(1) - 1, &
                            history = history, histories = histories)
    call vamp_print_history (history, "multi")
    call vamp_print_history (histories, "multi")
    do it = 1, nit
       call vamp_sample_grids (rng, grs, w, NO_DATA, 1, &
                               history = history(iterations(1)+it-1:), &
                               histories = histories(iterations(1)+it-1:,:))
       call vamp_print_history (history(iterations(1)+it-1:), "multi")
       call vamp_print_history (histories(iterations(1)+it-1:,:), "multi")
       call vamp_refine_weights (grs, powers(it))
    end do
    call vamp_discard_integrals (grs, samples(2))
    call vamp_sample_grids &
         (rng, grs, w, NO_DATA, iterations(2), &
          integral, standard_dev, chi_squared, &
          history = history(iterations(1)+nit:), &
          histories = histories(iterations(1)+nit:,:))
    call vamp_print_history (history(iterations(1)+nit:), "multi")
    call vamp_print_history (histories(iterations(1)+nit:,:), "multi")
    call vamp_write_grids (grs, "vamp_test.grids")
    call vamp_delete_grids (grs)
    call vamp_print_history (history, "multi")
    call vamp_print_history (histories, "multi")
    call vamp_delete_history (history)
    call vamp_delete_history (histories)
  end subroutine multi_channel
  subroutine print_results (prefix, prev_ticks, &
       integral, std_dev, chi2, acceptable, failures)
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: prev_ticks
    real(kind=default), intent(in) :: integral, std_dev, chi2, acceptable
    integer, intent(inout) :: failures
    integer :: ticks, ticks_per_second
    real(kind=default) :: pull
    call system_clock (ticks, ticks_per_second)
    pull = (integral - 1) / std_dev
    print "(1X,A,A,F6.2,A)", prefix, &
         ": time = ", real (ticks - prev_ticks) / ticks_per_second, " secs"
    print *, prefix, ":    int, err, chi2: ", &
         real (integral), real (std_dev), real (chi2)
    if (abs (pull) > acceptable) then
       failures = failures + 1
       print *, prefix, ": inacceptable pull:", real (pull)
    else
       print *, prefix, ":   acceptable pull:", real (pull)
    end if
  end subroutine print_results
end module vamp_tests
program vamp_test
  use kinds
  use tao_random_numbers
  use coordinates
  use vamp
  use vamp_test_functions !NODEP!
  use vamp_tests !NODEP!
  implicit none
  integer :: start_ticks, status
  integer, dimension(2) :: iterations, samples
  real(kind=default), dimension(2,5) :: region
  real(kind=default), dimension(5) :: weight_vector
  real(kind=default), dimension(10) :: powers
  real(kind=default) :: single_integral, single_standard_dev, single_chi_squared
  real(kind=default) :: multi_integral, multi_standard_dev, multi_chi_squared
  type(tao_random_state) :: rng
  real(kind=default), parameter :: ACCEPTABLE = 4
  integer :: failures
  failures = 0
  call tao_random_create (rng, 0)
  call get_environment_variable (name="VAMP_RANDOM_TESTS", status=status)
  if (status == 0) then
     call system_clock (start_ticks)
  else
     start_ticks = 42
  end if
  call tao_random_seed (rng, start_ticks)
  iterations = (/ 4, 3 /)
  samples = (/ 20000, 200000 /)
  region(1,:)  = -1.0
  region(2,:)  =  1.0
  width = 0.0001
  print *, "Starting VAMP 1.0 self test..."
  print *, "serial code"
  call system_clock (start_ticks)
  call single_channel (rng, region, samples, iterations, &
          single_integral, single_standard_dev, single_chi_squared)
  call print_results ("SINGLE", start_ticks, &
          single_integral, single_standard_dev, single_chi_squared, &
          10*ACCEPTABLE, failures)
  weight_vector = 1
  powers = 0.25_default
  call system_clock (start_ticks)
  call multi_channel (rng, region, weight_vector, samples, iterations, &
          powers, multi_integral, multi_standard_dev, multi_chi_squared)
  call print_results ("MULTI", start_ticks, &
           multi_integral, multi_standard_dev, multi_chi_squared, &
          ACCEPTABLE, failures)
  call system_clock (start_ticks)
! call check_jacobians (rng, region, weight_vector, samples(1))
  call check_inverses (rng, region, weight_vector, samples(1))
  call check_inverses3 (rng, region, samples(1))
  if (failures == 0) then
     stop 0
  else if (failures == 1) then
     stop 1
  else
     stop 2
  end if
end program vamp_test
