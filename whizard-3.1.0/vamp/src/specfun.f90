! specfun.f90 --
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
module specfun
  use kinds
! use constants
  implicit none
  private
  public :: gamma
  real(kind=default), public, parameter :: &
       PI = 3.1415926535897932384626433832795028841972_default
contains
  pure function gamma (x) result (g)
    real(kind=default), intent(in) :: x
    real(kind=default) :: g
    integer :: i
    real(kind=default) :: u, f, alpha, b0, b1, b2
    real(kind=default), dimension(0:15), parameter :: &
         c = (/ 3.65738772508338244_default, &
                1.95754345666126827_default, &
                0.33829711382616039_default, &
                0.04208951276557549_default, &
                0.00428765048212909_default, &
                0.00036521216929462_default, &
                0.00002740064222642_default, &
                0.00000181240233365_default, &
                0.00000010965775866_default, &
                0.00000000598718405_default, &
                0.00000000030769081_default, &
                0.00000000001431793_default, &
                0.00000000000065109_default, &
                0.00000000000002596_default, &
                0.00000000000000111_default, &
                0.00000000000000004_default /)

    u = x
    if (u <= 0.0) then
       if (u == int (u)) then
          g = huge (g)
          return
       else
          u = 1 - u
       end if
    endif
    f = 1
    if (u < 3) then
       do i = 1, int (4 - u)
          f = f / u
          u = u + 1
       end do
    else
       do i = 1, int (u - 3)
          u = u - 1
          f = f * u
       end do
    end if
    g = 2*u - 7
    alpha = 2*g
    b1 = 0
    b2 = 0
    do i = 15, 0, -1
       b0 = c(i) + alpha * b1 - b2
       b2 = b1
       b1 = b0
    end do
    g = f * (b0 - g * b2)
    if (x < 0) then
       g = PI / (sin (PI * x) * g)
    end if
  end function gamma
end module specfun
