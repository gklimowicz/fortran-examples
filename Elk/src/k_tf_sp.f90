
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine k_tf_sp(n,rhoup,rhodn,dtdru,dtdrd)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rhoup(n),rhodn(n)
real(8), intent(out) :: dtdru(n),dtdrd(n)
! local variables
integer i
do i=1,n
  call k_tf1(2.d0*rhoup(i),dtdru(i))
end do
do i=1,n
  call k_tf1(2.d0*rhodn(i),dtdrd(i))
end do
end subroutine

