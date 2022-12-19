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

module solver

  use kinds, only: default
  use constants, only: tiny_10

  implicit none
  private

  public :: solver_function_t
  public :: solve_secant
  public :: solve_interval
  public :: solve_qgaus

  real(default), parameter, public :: DEFAULT_PRECISION = tiny_10
  integer, parameter :: MAX_TRIES = 10000

  type, abstract :: solver_function_t
  contains
    procedure(solver_function_evaluate), deferred :: evaluate
  end type solver_function_t


  abstract interface
     function solver_function_evaluate (solver_f, x) result (f)
       import
       complex(default) :: f
       class(solver_function_t), intent(in) :: solver_f
       real(default), intent(in) :: x
     end function
  end interface


  interface
    module function solve_secant &
         (func, lower_start, upper_start, success, precision) result (x0)
      class(solver_function_t), intent(in) :: func
      real(default) :: x0
      real(default), intent(in) :: lower_start, upper_start
      real(default), intent(in), optional :: precision
      logical, intent(out) :: success
    end function solve_secant
    module function solve_interval &
         (func, lower_start, upper_start, success, precision) result (x0)
      class(solver_function_t), intent(in) :: func
      real(default) :: x0
      real(default), intent(in) :: lower_start, upper_start
      real(default), intent(in), optional :: precision
      logical, intent(out) :: success
    end function solve_interval
    module function solve_qgaus (integrand, grid) result (integral)
      class(solver_function_t), intent(in) :: integrand
      complex(default) :: integral
      real(default), dimension(:), intent(in) :: grid
    end function solve_qgaus
  end interface

end module solver
