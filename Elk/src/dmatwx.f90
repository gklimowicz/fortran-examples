
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dmatwx(n,w,x,dxx,dwx,xn)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: w(n,n),x(n,n)
complex(8), intent(out) :: dxx(n,n),dwx(n,n)
real(8), intent(out) :: xn(n)
! local variables
integer i
! normal bosonic density matrix XX†
call zgemm('N','C',n,n,n,(1.d0,0.d0),x,n,x,n,(0.d0,0.d0),dxx,n)
! store the X-norm
do i=1,n
  xn(i)=dot_product(x(:,i),x(:,i))
end do
! anomalous density matrix -WX†
call zgemm('N','C',n,n,n,(-1.d0,0.d0),w,n,x,n,(0.d0,0.d0),dwx,n)
end subroutine

