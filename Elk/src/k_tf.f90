
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine k_tf(n,rho,dtdr)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rho(n)
real(8), intent(out) :: dtdr(n)
! local variables
integer i
do i=1,n
  call k_tf1(rho(i),dtdr(i))
end do
end subroutine

