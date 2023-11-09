
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wsplint(n,x,w)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: w(n)
! local variables
integer i
! automatic arrays
real(8) f(9)
! external functions
real(8), external :: splint
if (n.le.9) then
  do i=1,n
    f(:)=0.d0
    f(i)=1.d0
    w(i)=splint(n,x,f)
  end do
  return
end if
do i=1,4
  f(:)=0.d0
  f(i)=1.d0
  w(i)=splint(9,x,f)
end do
f(:)=0.d0
f(5)=1.d0
do i=5,n-4
  w(i)=splint(9,x(i-4),f)
end do
do i=1,4
  f(:)=0.d0
  f(i+5)=1.d0
  w(n-4+i)=splint(9,x(n-8),f)
end do
end subroutine

