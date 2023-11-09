
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine k_vwlb(n,rho,grho2,tau)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rho(n),grho2(n)
real(8), intent(inout) :: tau(n)
! local variables
integer i
real(8) t1
! enforce the von Weizsacker lower bound
do i=1,n
  t1=(1.d0/8.d0)*grho2(i)/rho(i)
  if (tau(i).lt.t1) tau(i)=t1
end do
end subroutine

