! vamp_test0.f90 --
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
module vamp_test0_functions
  use kinds
  use vamp, only: vamp_grid, vamp_multi_channel0
  use vamp, only: vamp_data_t
  implicit none
  private
  public :: f, g, phi, w
  public :: create_sample, delete_sample
  private :: f0, psi, g0, f_norm
  real(kind=default), dimension(:), allocatable, private :: c, x_min, x_max
  real(kind=default), dimension(:,:,:), allocatable, public :: x0, gamma
contains
  pure function f0 (x, x_min, x_max, x0, g) result (f_x)
    real(kind=default), intent(in) :: x, x_min, x_max
    real(kind=default), dimension(:), intent(in) :: x0, g
    real(kind=default) :: f_x
    complex(kind=default) :: amp
    real(kind=default) :: norm
    integer :: i, j
    amp = sum (1.0 / cmplx (x - x0, g, kind=default))
    norm = 0
    do i = 1, size (x0)
       norm = norm + f_norm (x_min, x_max, x0(i), g(i), x0(i), g(i))
       do j = i + 1, size (x0)
          norm = norm + 2 * f_norm (x_min, x_max, x0(i), g(i), x0(j), g(j))
       end do
    end do
    f_x = amp * conjg (amp) / norm
  end function f0
  pure function f_norm (x_min, x_max, x0p, gp, x0q, gq) &
       result (norm)
    real(kind=default), intent(in) :: x_min, x_max, x0p, gp, x0q, gq
    real(kind=default) :: norm
    norm = real ((   log (   cmplx (x_max - x0p,   gp, kind=default) &
                           / cmplx (x_min - x0p,   gp, kind=default)) &
                   - log (   cmplx (x_max - x0q, - gq, kind=default) &
                           / cmplx (x_min - x0q, - gq, kind=default))) &
                   / cmplx (x0p - x0q, - gp - gq, kind=default), &
                 kind=default)
  end function f_norm
  pure function f (x, data, weights, channel, grids) result (f_x)
    real(kind=default), dimension(:), intent(in) :: x
    class(vamp_data_t), intent(in) :: data
    real(kind=default), dimension(:), intent(in), optional :: weights
    integer, intent(in), optional :: channel
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    real(kind=default) :: f_x
    real(kind=default) :: fi_x
    integer :: i, j
    f_x = 0.0
    do i = 1, size (c)
       fi_x = 1.0
       do j = 1, size (x)
          if (all (gamma(:,i,j) > 0)) then
             fi_x = fi_x * f0 (x(j), x_min(j), x_max(j), &
                               x0(:,i,j), gamma(:,i,j))
          else
             fi_x = fi_x / (x_max(j) - x_min(j))
          end if
       end do
       f_x = f_x + c(i) * fi_x
    end do
    f_x = f_x / sum (c)
  end function f
  subroutine delete_sample ()
    deallocate (c, x_min, x_max, x0, gamma)
  end subroutine delete_sample
  subroutine create_sample (num_poles, weights, region)
    integer, intent(in) :: num_poles
    real(kind=default), dimension(:), intent(in) :: weights
    real(kind=default), dimension(:,:), intent(in) :: region
    integer :: nd, nc
    nd = size (region, dim=2)
    nc = size (weights)
    allocate (c(nc), x_min(nd), x_max(nd))
    allocate (x0(num_poles,nc,nd), gamma(num_poles,nc,nd))
    x_min = region(1,:)
    x_max = region(2,:)
    c = weights
  end subroutine create_sample
  pure function psi (xi, x_min, x_max, x0, gamma) result (x)
    real(kind=default), intent(in) :: xi, x_min, x_max, x0, gamma
    real(kind=default) :: x
    x = x0 + gamma &
         * tan (((xi - x_min) * atan ((x_max - x0) / gamma) &
                   - (x_max - xi) * atan ((x0 - x_min) / gamma)) &
                 / (x_max - x_min))
  end function psi
  pure function g0 (x, x_min, x_max, x0, gamma) result (g_x)
    real(kind=default), intent(in) :: x, x_min, x_max, x0, gamma
    real(kind=default) :: g_x
    g_x = gamma / (atan ((x_max - x0) / gamma) - atan ((x_min - x0) / gamma)) &
           * (x_max - x_min) / ((x - x0)**2 + gamma**2)
  end function g0
  pure function phi (xi, channel) result (x)
    real(kind=default), dimension(:), intent(in) :: xi
    integer, intent(in) :: channel
    real(kind=default), dimension(size(xi)) :: x
    integer, dimension(size(xi)) :: p
    integer :: j, ch, np, nch, nd, channels
    np = size (x0, dim = 1)
    nch = size (x0, dim = 2)
    nd = size (x0, dim = 3)
    channels = nch * np**nd
    if (channel >= 1 .and. channel <= channels) then
       ch = channel - 1
       do j = 1, size (x)
          p(j) = 1 + modulo (ch, np)
          ch = ch / np
       end do
       ch = ch + 1
       do j = 1, size (xi)
          if (all (gamma(:,ch,j) > 0)) then
             x(j) = psi (xi(j), x_min(j), x_max(j), &
                         x0(p(j),ch,j), gamma(p(j),ch,j))
          else
             x = xi
          end if
       end do
    else if (channel == channels + 1) then
       x = xi
    else
       x = 0
    end if
  end function phi
  pure recursive function g (x, data, channel) result (g_x)
    real(kind=default), dimension(:), intent(in) :: x
    class(vamp_data_t), intent(in) :: data
    integer, intent(in) :: channel
    real(kind=default) :: g_x
    integer, dimension(size(x)) :: p
    integer :: j, ch, np, nch, nd, channels
    np = size (x0, dim = 1)
    nch = size (x0, dim = 2)
    nd = size (x0, dim = 3)
    channels = nch * np**nd
    if (channel >= 1 .and. channel <= channels) then
       ch = channel - 1
       do j = 1, size (x)
          p(j) = 1 + modulo (ch, np)
          ch = ch / np
       end do
       ch = ch + 1
       g_x = 1.0
       do j = 1, size (x)
          if (all (gamma(:,ch,j) > 0)) then
             g_x = g_x * g0 (x(j), x_min(j), x_max(j), &
                             x0(p(j),ch,j), gamma(p(j),ch,j))
          end if
       end do
    else if (channel == channels + 1) then
       g_x = 1.0
    else
       g_x = 0
    end if
  end function g
  function w (x, data, weights, channel, grids) result (w_x)
    real(kind=default), dimension(:), intent(in) :: x
    class(vamp_data_t), intent(in) :: data
    real(kind=default), dimension(:), intent(in), optional :: weights
    integer, intent(in), optional :: channel
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    real(kind=default) :: w_x
    w_x = vamp_multi_channel0 (f, data, phi, g, x, weights, channel)
  end function w
end module vamp_test0_functions
module vamp_tests0
  use kinds
  use exceptions
  use histograms
  use tao_random_numbers
  use vamp_test0_functions !NODEP!
  use vamp
  implicit none
  private
  public :: single_channel, multi_channel
  public :: single_channel_generator, multi_channel_generator
contains
  subroutine single_channel (do_print, region, iterations, samples, rng, &
       acceptable, failures)
    logical, intent(in) :: do_print
    real(kind=default), dimension(:,:), intent(in) :: region
    integer, dimension(:), intent(in) :: iterations, samples
    type(tao_random_state), intent(inout) :: rng
    real(kind=default), intent(in) :: acceptable
    integer, intent(inout) :: failures
    type(vamp_grid) :: gr
    type(vamp_history), dimension(iterations(1)+iterations(2)) :: history
    real(kind=default) :: integral, standard_dev, chi_squared, pull
    call vamp_create_history (history)
    call vamp_create_grid (gr, region, samples(1))
    call vamp_sample_grid (rng, gr, f, NO_DATA, iterations(1), history = history)
    call vamp_discard_integral (gr, samples(2))
    call vamp_sample_grid &
         (rng, gr, f, NO_DATA, iterations(2), &
          integral, standard_dev, chi_squared, &
          history = history(iterations(1)+1:))
    call vamp_write_grid (gr, "vamp_test0.grid")
    call vamp_delete_grid (gr)
    call vamp_print_history (history, "single")
    call vamp_delete_history (history)
    pull = (integral - 1) / standard_dev
    if (do_print) then
       print *, "    int, err, chi2:", integral, standard_dev, chi_squared
    end if
    if (abs (pull) > acceptable) then
       failures = failures + 1
       print *, " unacceptable pull:", pull
    else
       print *, "   acceptable pull:", pull
    end if
  end subroutine single_channel
  subroutine multi_channel (do_print, region, iterations, samples, rng, &
       acceptable, failures)
    logical, intent(in) :: do_print
    real(kind=default), dimension(:,:), intent(in) :: region
    integer, dimension(:), intent(in) :: iterations, samples
    type(tao_random_state), intent(inout) :: rng
    real(kind=default), intent(in) :: acceptable
    type(vamp_grids) :: grs
    integer, intent(inout) :: failures
    real(kind=default), &
         dimension(size(x0,dim=2)*size(x0,dim=1)**size(x0,dim=3)+1) :: &
         weight_vector
    type(vamp_history), dimension(iterations(1)+iterations(2)+4) :: history
    type(vamp_history), dimension(size(history),size(weight_vector)) :: histories
    real(kind=default) :: integral, standard_dev, chi_squared, pull
    integer :: it
    weight_vector = 1.0
    call vamp_create_history (history)
    call vamp_create_history (histories)
    call vamp_create_grids (grs, region, samples(1), weight_vector)
    call vamp_sample_grids (rng, grs, w, NO_DATA, iterations(1) - 1, &
                            history = history, histories = histories)
    do it = 1, 5
       call vamp_sample_grids (rng, grs, w, NO_DATA, 1, &
                               history = history(iterations(1)+it-1:), &
                               histories = histories(iterations(1)+it-1:,:))
       call vamp_refine_weights (grs)
    end do
    call vamp_discard_integrals (grs, samples(2))
    call vamp_sample_grids &
         (rng, grs, w, NO_DATA, iterations(2), &
          integral, standard_dev, chi_squared, &
          history = history(iterations(1)+5:), &
          histories = histories(iterations(1)+5:,:))
    call vamp_write_grids (grs, "vamp_test0.grids")
    call vamp_delete_grids (grs)
    call vamp_print_history (history, "multi")
    call vamp_print_history (histories, "multi")
    call vamp_delete_history (history)
    call vamp_delete_history (histories)
    if (do_print) then
       print *, integral, standard_dev, chi_squared
    end if
    pull = (integral - 1) / standard_dev
    if (abs (pull) > acceptable) then
       failures = failures + 1
       print *, " unacceptable pull:", pull
    else
       print *, "   acceptable pull:", pull
    end if
  end subroutine multi_channel
  subroutine single_channel_generator (do_print, region, iterations, samples, rng)
    logical, intent(in) :: do_print
    real(kind=default), dimension(:,:), intent(in) :: region
    integer, dimension(:), intent(in) :: iterations, samples
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grid) :: gr
    type(vamp_history), dimension(iterations(1)+iterations(2)) :: history
    type(histogram) :: unweighted, reweighted, weighted, weights
    type(exception) :: exc
    real(kind=default) :: weight, integral, standard_dev
    integer :: i
    real(kind=default), dimension(size(region,dim=2)) :: x
    call vamp_create_grid (gr, region, samples(1))
    call vamp_sample_grid (rng, gr, f, NO_DATA, iterations(1), history = history)
    call vamp_discard_integral (gr, samples(2))
    call vamp_warmup_grid &
         (rng, gr, f, NO_DATA, iterations(2), history = history(iterations(1)+1:))
    call vamp_print_history (history, "single")
    call vamp_delete_history (history)
    call create_histogram (unweighted, region(1,1), region(2,1), 100)
    call create_histogram (reweighted, region(1,1), region(2,1), 100)
    call create_histogram (weighted, region(1,1), region(2,1), 100)
    call create_histogram (weights, 0.0_default, 10.0_default, 100)
    ! do i = 1, 1000000
    do i = 1, 100
       call clear_exception (exc)
       call vamp_next_event (x, rng, gr, f, NO_DATA, exc = exc)
       call handle_exception (exc)
       call fill_histogram (unweighted, x(1))
       call fill_histogram (reweighted, x(1), 1.0_default / f (x, NO_DATA))
    end do
    integral = 0.0
    standard_dev = 0.0
    do i = 1, 10000
       call clear_exception (exc)
       call vamp_next_event (x, rng, gr, f, NO_DATA, weight, exc = exc)
       call handle_exception (exc)
       call fill_histogram (weighted, x(1), weight / f (x, NO_DATA))
       call fill_histogram (weights, x(1), weight)
       integral = integral + weight
       standard_dev = standard_dev + weight**2
    end do
    if (do_print) then
       print *, integral / (i-1), sqrt (standard_dev) / (i-1)
       call write_histogram (unweighted, "u_s.d")
       call write_histogram (reweighted, "r_s.d")
       call write_histogram (weighted, "w_s.d")
       call write_histogram (weights, "ws_s.d")
    end if
    call delete_histogram (unweighted)
    call delete_histogram (reweighted)
    call delete_histogram (weighted)
    call delete_histogram (weights)
    call vamp_delete_grid (gr)
  end subroutine single_channel_generator
  subroutine multi_channel_generator (do_print, region, iterations, samples, rng)
    logical, intent(in) :: do_print
    real(kind=default), dimension(:,:), intent(in) :: region
    integer, dimension(:), intent(in) :: iterations, samples
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grids) :: grs
    real(kind=default), &
         dimension(size(x0,dim=2)*size(x0,dim=1)**size(x0,dim=3)+1) :: &
         weight_vector
    type(vamp_history), dimension(iterations(1)+iterations(2)+4) :: history
    type(vamp_history), dimension(size(history),size(weight_vector)) :: histories
    type(histogram) :: unweighted, reweighted, weighted, weights
    type(exception) :: exc
    real(kind=default) :: weight, integral, standard_dev
    real(kind=default), dimension(size(region,dim=2)) :: x
    character(len=5) :: pfx
    integer :: it, i, j
    weight_vector = 1.0
    call vamp_create_history (history)
    call vamp_create_history (histories)
    call vamp_create_grids (grs, region, samples(1), weight_vector)
    call vamp_sample_grids (rng, grs, w, NO_DATA, iterations(1) - 1, &
                            history = history, histories = histories)
    do it = 1, 5
       call vamp_sample_grids (rng, grs, w, NO_DATA, 1, &
                               history = history(iterations(1)+it-1:), &
                               histories = histories(iterations(1)+it-1:,:))
       call vamp_refine_weights (grs)
    end do
    call vamp_discard_integrals (grs, samples(2))
    call vamp_warmup_grids &
         (rng, grs, w, NO_DATA, iterations(2), &
          history = history(iterations(1)+5:), &
          histories = histories(iterations(1)+5:,:))
    call vamp_print_history (history, "multi")
    call vamp_print_history (histories, "multi")
    call vamp_delete_history (history)
    call vamp_delete_history (histories)
    !!! do i = 1, size (grs%grids)
    !!!    do j = 1, size (grs%grids(i)%div)
    !!!       write (pfx, "(I2.2,'/',I2.2)") i, j
    !!!       call dump_division (grs%grids(i)%div(j), pfx)
    !!!    end do
    !!! end do
    call create_histogram (unweighted, region(1,1), region(2,1), 100)
    call create_histogram (reweighted, region(1,1), region(2,1), 100)
    call create_histogram (weighted, region(1,1), region(2,1), 100)
    call create_histogram (weights, 0.0_default, 10.0_default, 100)
    ! do i = 1, 1000000
    do i = 1, 100
       call clear_exception (exc)
       call vamp_next_event (x, rng, grs, f, NO_DATA, phi, exc = exc)
       call handle_exception (exc)
       call fill_histogram (unweighted, x(1))
       call fill_histogram (reweighted, x(1), 1.0_default / f (x, NO_DATA))
    end do
    integral = 0.0
    standard_dev = 0.0
    do i = 1, 10000
       call clear_exception (exc)
       call vamp_next_event (x, rng, grs, f, NO_DATA, phi, weight, exc = exc)
       call handle_exception (exc)
       call fill_histogram (weighted, x(1), weight / f (x, NO_DATA))
       call fill_histogram (weights, x(1), weight)
       integral = integral + weight
       standard_dev = standard_dev + weight**2
    end do
    if (do_print) then
       print *, integral / (i-1), sqrt (standard_dev) / (i-1)
       call write_histogram (unweighted, "u_m.d")
       call write_histogram (reweighted, "r_m.d")
       call write_histogram (weighted, "w_m.d")
       call write_histogram (weights, "ws_m.d")
    end if
    call delete_histogram (unweighted)
    call delete_histogram (reweighted)
    call delete_histogram (weighted)
    call delete_histogram (weights)
    call vamp_delete_grids (grs)
  end subroutine multi_channel_generator
end module vamp_tests0
program vamp_test0
  use kinds
  use tao_random_numbers
  use vamp_test0_functions !NODEP!
  use vamp_tests0 !NODEP!
  implicit none
  logical :: do_print
  integer :: i, j, ticks, ticks_per_second, ticks0, status
  integer, dimension(2) :: iterations, samples
  real(kind=default), dimension(:,:), allocatable :: region
  type(tao_random_state) :: rng
  real(kind=default), parameter :: ACCEPTABLE = 4
  integer :: failures
  do_print = .true.
  print *, "Starting VAMP 1.0 self test..."
  print *, "serial code"
  call tao_random_create (rng, 0)
  call get_environment_variable (name="VAMP_RANDOM_TESTS", status=status)
  if (status == 0) then
     call system_clock (ticks0)
  else
     ticks0 = 42
  end if
  call tao_random_seed (rng, ticks0)
  iterations = (/ 4, 3 /)
  samples = (/ 10000, 50000 /)
  allocate (region(2,2))
  region(1,:)  = -1.0
  region(2,:)  =  2.0
  call create_sample &
       (num_poles = 2, weights = (/ 1.0_default, 2.0_default /), region = region)
  do i = 1, size (x0, dim=2)
     do j = 1, size (x0, dim=3)
        call tao_random_number (rng, x0(:,i,j))
     end do
  end do
  gamma = 0.001
  x0(1,:,:) = 0.2
  x0(2:,:,:) = 0.8
  failures = 0
  call system_clock (ticks0)
  call single_channel (do_print, region, iterations, samples, rng, 10*ACCEPTABLE, failures)
  call system_clock (ticks, ticks_per_second)
  print "(1X,A,F6.2,A)", &
       "time = ", real (ticks - ticks0) / ticks_per_second, " secs"
  call system_clock (ticks0)
  call single_channel_generator &
         (do_print, region, iterations, samples, rng)
  call system_clock (ticks, ticks_per_second)
  print "(1X,A,F6.2,A)", &
       "time = ", real (ticks - ticks0) / ticks_per_second, " secs"
  call system_clock (ticks0)
  call multi_channel (do_print, region, iterations, samples, rng, ACCEPTABLE, failures)
  call system_clock (ticks, ticks_per_second)
  print "(1X,A,F6.2,A)", &
       "time = ", real (ticks - ticks0) / ticks_per_second, " secs"
  call system_clock (ticks0)
  call multi_channel_generator &
         (do_print, region, iterations, samples, rng)
  call system_clock (ticks, ticks_per_second)
  print "(1X,A,F6.2,A)", &
       "time = ", real (ticks - ticks0) / ticks_per_second, " secs"
  call system_clock (ticks0)
  ! call check_jacobians (do_print, region, samples, rng)
  call system_clock (ticks, ticks_per_second)
  print "(1X,A,F6.2,A)", &
         "time = ", real (ticks - ticks0) / ticks_per_second, " secs"

  call delete_sample ()
  deallocate (region)
  if (failures == 0) then
     stop 0
  else if (failures == 1) then
     stop 1
  else
     stop 2
  end if
end program vamp_test0
