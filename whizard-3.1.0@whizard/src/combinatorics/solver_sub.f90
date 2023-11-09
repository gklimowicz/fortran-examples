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

submodule (solver) solver_s

  use numeric_utils
  use diagnostics

  implicit none

contains

  module function solve_secant &
       (func, lower_start, upper_start, success, precision) result (x0)
    class(solver_function_t), intent(in) :: func
    real(default) :: x0
    real(default), intent(in) :: lower_start, upper_start
    real(default), intent(in), optional :: precision
    logical, intent(out) :: success
    real(default) :: desired, x_curr, x_next, f_curr, f_next, x_new
    integer :: n_iter
    desired = DEFAULT_PRECISION; if (present(precision)) desired = precision
    x_curr = lower_start
    x_next = upper_start
    n_iter = 0
    success = .false.
    SEARCH: do
       n_iter = n_iter + 1
       f_curr = real( func%evaluate (x_curr) )
       f_next = real( func%evaluate (x_next) )
       if (abs (f_next) < desired) then
          x0 = x_next
          exit
       end if
       if (n_iter > MAX_TRIES) then
          call msg_warning ("solve: Couldn't find root of function")
          return
       end if
       if (vanishes (f_next - f_curr)) then
          x_next = x_next + (x_next - x_curr) / 10
          cycle
       end if

       x_new = x_next - (x_next - x_curr) / (f_next - f_curr) * f_next
       x_curr = x_next
       x_next = x_new
    end do SEARCH
    if (x0 < lower_start .or. x0 > upper_start) then
       call msg_warning ("solve: The root of the function is not in boundaries")
       return
    end if
    success = .true.
  end function solve_secant

  module function solve_interval &
       (func, lower_start, upper_start, success, precision) result (x0)
    class(solver_function_t), intent(in) :: func
    real(default) :: x0
    real(default), intent(in) :: lower_start, upper_start
    real(default), intent(in), optional :: precision
    logical, intent(out) :: success
    real(default) :: desired
    real(default) :: x_low, x_high, x_half
    real(default) :: f_low, f_high, f_half
    integer :: n_iter
    success = .false.
    desired = DEFAULT_PRECISION; if (present(precision)) desired = precision
    x0 = lower_start
    x_low = lower_start
    x_high = upper_start
    f_low = real( func%evaluate (x_low) )
    f_high = real( func%evaluate (x_high) )
    if (f_low * f_high > 0) return
    if (x_low > x_high) then
       call display_solver_status()
       call msg_fatal ("Interval solver: Upper bound must be &
                       &greater than lower bound")
    end if
    n_iter = 0
    do n_iter = 1, MAX_TRIES
       x_half = (x_high + x_low)/2
       f_half = real( func%evaluate (x_half) )
       if (abs (f_half) <= desired) then
          x0 = x_half
          exit
       end if
       if (f_low * f_half > 0._default) then
          x_low = x_half
          f_low = f_half
       else
          x_high = x_half
          f_high = f_half
       end if
    end do
    if (x0 < lower_start .or. x0 > upper_start) then
       call msg_warning ("Interval solver: The root of the function&
                          & is out of boundaries")
       return
    end if
    success = .true.
  contains
    subroutine display_solver_status ()
       print *, '================='
       print *, 'Status of interval solver: '
       print *, 'initial values: ', lower_start, upper_start
       print *, 'iteration: ', n_iter
       print *, 'x_low: ', x_low, 'f_low: ', f_low
       print *, 'x_high: ', x_high, 'f_high: ', f_high
       print *, 'x_half: ', x_half, 'f_half: ', f_half
    end subroutine display_solver_status
  end function solve_interval

  module function solve_qgaus (integrand, grid) result (integral)
    class(solver_function_t), intent(in) :: integrand
    complex(default) :: integral
    real(default), dimension(:), intent(in) :: grid
    integer :: i, j
    real(default) :: xm, xr
    real(default), dimension(5) :: dx, &
      w = (/ 0.2955242247_default, 0.2692667193_default, &
        0.2190863625_default, 0.1494513491_default, 0.0666713443_default /), &
      x = (/ 0.1488743389_default, 0.4333953941_default, 0.6794095682_default, &
        0.8650633666_default, 0.9739065285_default /)
    integral = 0.0_default
    if ( size(grid) < 2 ) then
      call msg_warning ("solve_qgaus: size of integration grid smaller than 2.")
      return
    end if
    do i=1, size(grid)-1
      xm = 0.5_default * ( grid(i+1) + grid(i) )
      xr = 0.5_default * ( grid(i+1) - grid(i) )
      do j=1, 5
        dx(j) = xr * x(j)
        integral = integral + xr * w(j) * &
          ( integrand%evaluate (xm+dx(j)) + integrand%evaluate (xm-dx(j)) )
      end do
    end do
  end function solve_qgaus


end submodule solver_s

