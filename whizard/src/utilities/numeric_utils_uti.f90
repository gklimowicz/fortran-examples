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

module numeric_utils_uti

  use kinds, only: default
  use constants, only: one, PI
  use numeric_utils

  implicit none
  private

  public :: numeric_utils_1
  public :: numeric_utils_2

 contains

  subroutine numeric_utils_1 (u)
    integer, intent(in) :: u
    real(default) :: int, intabs, err, errabs
    real(default) :: a, b, c
    a = 0._default
    b = 1._default
    c = Pi

    write (u, "(A)") "* Test output: Numeric utils"
    write (u, "(A)") "*   Purpose: test numeric routines"
    write (u, "(A)")

    write (u, "(A)") "*   41-point Gauss-Kronrod integration"
    write (u, "(A)")
    
    call dqk41 (f_x, a, b, int, intabs, err, errabs)
    write (u, "(1x,A,F9.6)")  " Integral (x,[0,1])       = ", &
         int
    call dqk41 (f_x2, a, b, int, intabs, err, errabs)
    write (u, "(1x,A,F9.6)")  " Integral (x**2,[0,1])    = ", &
         int
    call dqk41 (sinx, a, c, int, intabs, err, errabs)
    write (u, "(1x,A,F9.6)")  " Integral (sin(x),[0,Pi]) = ", &
         int

  contains

    function f_x (x) result (f)
      real(default), intent(in) :: x
      real(default) :: f
      f = x
    end function f_x
    function f_x2 (x) result (f)
      real(default), intent(in) :: x
      real(default) :: f
      f = x**2
    end function f_x2
    function sinx (x) result (f)
      real(default), intent(in) :: x
      real(default) :: f
      f = sin(x)
    end function sinx
  end subroutine numeric_utils_1

  subroutine numeric_utils_2 (u)
    integer, intent(in) :: u
    real(default) :: result, abserr
    real(default), parameter :: epsabs = 0.001_default, &
         epsrel = 0.001_default
    real(default), parameter :: a = 0._default, b = 1._default, &
         z = 0.1_default
    integer, parameter :: limit = 10000

    write (u, "(A)") "* Test output: Numeric utils"
    write (u, "(A)") "*   Purpose: test adaptive Gauss-Kronrod 41"
    write (u, "(A)")

    write (u, "(A)") "*   41-point Gauss-Kronrod integration"
    write (u, "(A)")
    
    call gauss_kronrod (GAUSS_KRONROD_41, f1_x, a, b, limit, result, &
         abserr, epsabs, epsrel)
    write (u, "(1x,A,F9.6)")  " Integral (f1_x,[0,1])     = ", &
         result
    call gauss_kronrod (GAUSS_KRONROD_41, f2_x, a, b, limit, result, &
         abserr, epsabs, epsrel)
    write (u, "(1x,A,F9.6)")  " Integral (f2_x,[0,1])     = ", &
         result

  contains

    function f1_x (x) result (f)
      real(default), intent(in) :: x
      real(default) :: f
      f = log(one - z/(one + (-one + z)*x))
    end function f1_x
    function f2_x (x) result (f)
      real(default), intent(in) :: x
      real(default) :: f
      f = log(one - z/(z + x))
    end function f2_x
  end subroutine numeric_utils_2


end module numeric_utils_uti
