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

module solver_uti

  use kinds, only: default
  use constants, only: zero, one, two
  use numeric_utils

  use solver

  implicit none
  private

  public :: solver_1

  type, extends (solver_function_t) :: test_function_1_t
  contains
    procedure :: evaluate => test_func_1
  end type test_function_1_t

  type, extends (solver_function_t) :: test_function_2_t
  contains
    procedure :: evaluate => test_func_2
  end type test_function_2_t

  type, extends (solver_function_t) :: test_function_3_t
  contains
    procedure :: evaluate => test_func_3
  end type test_function_3_t

  type, extends (solver_function_t) :: test_function_4_t
  contains
    procedure :: evaluate => test_func_4
  end type test_function_4_t


contains

  subroutine solver_1 (u)
    integer, intent(in) :: u
    real(default) :: zero_position
    logical :: success
    type(test_function_1_t) :: test_func_1
    type(test_function_2_t) :: test_func_2
    type(test_function_3_t) :: test_func_3
    type(test_function_4_t) :: test_func_4
    write (u, "(A)")  "* Test output: solver_1"
    write (u, "(A)")  "*   Purpose: Solve trivial functions"
    write (u, "(A)")

    zero_position = solve_interval (test_func_1, -one, one, success)
    call assert (u, success, "success")
    call assert_equal (u, zero_position, zero, "test_func_1: zero_position")

    zero_position = solve_interval (test_func_4, two, 10.0_default, success)
    call assert (u, success, "success")
    call assert_equal (u, zero_position, &
         3.5216674011865940283397224_default, &
         "test_func_4: zero_position", rel_smallness=1000*DEFAULT_PRECISION)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: solver_1"
  end subroutine solver_1


  function test_func_1 (solver_f, x) result (f)
    complex(default) :: f
    class(test_function_1_t), intent(in) :: solver_f
    real(default), intent(in) :: x
    f = x
  end function test_func_1

  function test_func_2 (solver_f, x) result (f)
    complex(default) :: f
    class(test_function_2_t), intent(in) :: solver_f
    real(default), intent(in) :: x
    f = x ** 2
  end function test_func_2

  function test_func_3 (solver_f, x) result (f)
    complex(default) :: f
    class(test_function_3_t), intent(in) :: solver_f
    real(default), intent(in) :: x
    f = x ** 3
  end function test_func_3

  function test_func_4 (solver_f, x) result (f)
    complex(default) :: f
    class(test_function_4_t), intent(in) :: solver_f
    real(default), intent(in) :: x
    real(default) :: s, cutoff
    s = 100.0_default
    cutoff = 1.01_default
    if (x < cutoff) then
       f = - (log (s) * log (log (s) / log(cutoff**2)) - log (s / cutoff**2)) - &
         log (one/two)
    else
       f = - (log (s) * log (log (s) / log(x**2)) - log (s / x**2)) - &
            log (one/two)
    end if
  end function test_func_4


end module solver_uti
