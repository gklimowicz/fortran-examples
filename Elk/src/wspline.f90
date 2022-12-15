
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wspline(n,x,wc)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: wc(4,3,n)
! local variables
integer i,j
real(8) f(4),cf(3,4)
if (n.lt.4) then
  write(*,*)
  write(*,'("Error(wspline): n < 4 : ",I8)') n
  write(*,*)
  stop
end if
f(1)=1.d0
f(2:)=0.d0
call spline(4,x,f,cf)
wc(1,:,1)=cf(:,1)
wc(1,:,2)=cf(:,2)
f(1)=0.d0
f(2)=1.d0
call spline(4,x,f,cf)
wc(2,:,1)=cf(:,1)
wc(2,:,2)=cf(:,2)
f(2)=0.d0
f(3)=1.d0
call spline(4,x,f,cf)
wc(3,:,1)=cf(:,1)
wc(3,:,2)=cf(:,2)
f(3)=0.d0
f(4)=1.d0
call spline(4,x,f,cf)
wc(4,:,1)=cf(:,1)
wc(4,:,2)=cf(:,2)
do i=3,n-3
  j=i-1
  f(1)=1.d0
  f(2:)=0.d0
  call spline(4,x(j),f,cf)
  wc(1,:,i)=cf(:,2)
  f(1)=0.d0
  f(2)=1.d0
  call spline(4,x(j),f,cf)
  wc(2,:,i)=cf(:,2)
  f(2)=0.d0
  f(3)=1.d0
  call spline(4,x(j),f,cf)
  wc(3,:,i)=cf(:,2)
  f(3)=0.d0
  f(4)=1.d0
  call spline(4,x(j),f,cf)
  wc(4,:,i)=cf(:,2)
end do
j=n-3
f(1)=1.d0
f(2:)=0.d0
call spline(4,x(j),f,cf)
wc(1,:,n-2)=cf(:,2)
wc(1,:,n-1)=cf(:,3)
wc(1,:,n)=cf(:,4)
f(1)=0.d0
f(2)=1.d0
call spline(4,x(j),f,cf)
wc(2,:,n-2)=cf(:,2)
wc(2,:,n-1)=cf(:,3)
wc(2,:,n)=cf(:,4)
f(2)=0.d0
f(3)=1.d0
call spline(4,x(j),f,cf)
wc(3,:,n-2)=cf(:,2)
wc(3,:,n-1)=cf(:,3)
wc(3,:,n)=cf(:,4)
f(3)=0.d0
f(4)=1.d0
call spline(4,x(j),f,cf)
wc(4,:,n-2)=cf(:,2)
wc(4,:,n-1)=cf(:,3)
wc(4,:,n)=cf(:,4)
end subroutine

