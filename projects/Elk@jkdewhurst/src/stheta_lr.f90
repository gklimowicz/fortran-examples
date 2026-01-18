
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

elemental real(8) function stheta_lr(x)
implicit none
! arguments
real(8), intent(in) :: x
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
stheta_lr=0.5d0+atan(2.d0*x)/pi
end function

