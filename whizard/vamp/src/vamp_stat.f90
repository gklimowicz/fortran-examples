! vamp_stat.f90 --
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
module vamp_stat
  use kinds
  implicit none
  private
  public :: average, standard_deviation, value_spread
  public :: standard_deviation_percent, value_spread_percent
contains
  pure function average (x) result (a)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default) :: a
    integer :: n
    n = size (x)
    if (n == 0) then
       a = 0.0
    else
       a = sum (x) / n
    end if
  end function average
  pure function standard_deviation (x) result (s)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default) :: s
    integer :: n
    n = size (x)
    if (n < 2) then
       s = huge (s)
    else
       s = sqrt (max ((sum (x**2) / n - (average (x))**2) / (n - 1), &
                      0.0_default))
    end if
  end function standard_deviation
  pure function value_spread (x) result (s)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default) :: s
    s = maxval(x) - minval(x)
  end function value_spread
  pure function standard_deviation_percent (x) result (s)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default) :: s
    real(kind=default) :: abs_avg
    abs_avg = abs (average (x))
    if (abs_avg <= tiny (abs_avg)) then
       s = huge (s)
    else
       s = 100.0 * standard_deviation (x) / abs_avg
    end if
  end function standard_deviation_percent
  pure function value_spread_percent (x) result (s)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default) :: s
    real(kind=default) :: abs_avg
    abs_avg = abs (average (x))
    if (abs_avg <= tiny (abs_avg)) then
       s = huge (s)
    else
       s = 100.0 * value_spread (x) / abs_avg
    end if
  end function value_spread_percent
end module vamp_stat
