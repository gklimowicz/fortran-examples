
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine splinew(n,wc,f,cf)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wc(4,3,n),f(n)
real(8), intent(out) :: cf(3,n)
! local variables
integer i
real(8) f1,f2,f3,f4
f1=f(1); f2=f(2); f3=f(3); f4=f(4)
cf(1,1)=wc(1,1,1)*f1+wc(2,1,1)*f2+wc(3,1,1)*f3+wc(4,1,1)*f4
cf(2,1)=wc(1,2,1)*f1+wc(2,2,1)*f2+wc(3,2,1)*f3+wc(4,2,1)*f4
cf(3,1)=wc(1,3,1)*f1+wc(2,3,1)*f2+wc(3,3,1)*f3+wc(4,3,1)*f4
cf(1,2)=wc(1,1,2)*f1+wc(2,1,2)*f2+wc(3,1,2)*f3+wc(4,1,2)*f4
cf(2,2)=wc(1,2,2)*f1+wc(2,2,2)*f2+wc(3,2,2)*f3+wc(4,2,2)*f4
cf(3,2)=wc(1,3,2)*f1+wc(2,3,2)*f2+wc(3,3,2)*f3+wc(4,3,2)*f4
!$OMP SIMD LASTPRIVATE(f1,f2,f3,f4)
do i=3,n-2
  f1=f(i-1); f2=f(i); f3=f(i+1); f4=f(i+2)
  cf(1,i)=wc(1,1,i)*f1+wc(2,1,i)*f2+wc(3,1,i)*f3+wc(4,1,i)*f4
  cf(2,i)=wc(1,2,i)*f1+wc(2,2,i)*f2+wc(3,2,i)*f3+wc(4,2,i)*f4
  cf(3,i)=wc(1,3,i)*f1+wc(2,3,i)*f2+wc(3,3,i)*f3+wc(4,3,i)*f4
end do
i=n-1
cf(1,i)=wc(1,1,i)*f1+wc(2,1,i)*f2+wc(3,1,i)*f3+wc(4,1,i)*f4
cf(2,i)=wc(1,2,i)*f1+wc(2,2,i)*f2+wc(3,2,i)*f3+wc(4,2,i)*f4
cf(3,i)=wc(1,3,i)*f1+wc(2,3,i)*f2+wc(3,3,i)*f3+wc(4,3,i)*f4
cf(1,n)=wc(1,1,n)*f1+wc(2,1,n)*f2+wc(3,1,n)*f3+wc(4,1,n)*f4
cf(2,n)=wc(1,2,n)*f1+wc(2,2,n)*f2+wc(3,2,n)*f3+wc(4,2,n)*f4
cf(3,n)=wc(1,3,n)*f1+wc(2,3,n)*f2+wc(3,3,n)*f3+wc(4,3,n)*f4
end subroutine

