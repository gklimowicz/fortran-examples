
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dmatuv(n,ef,e,u,v,dvv,duv,vn)
use modmain
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: ef,e(n)
complex(8), intent(in) :: u(n,n),v(n,n)
complex(8), intent(out) :: dvv(n,n),duv(n,n)
real(8), intent(out) :: vn(n)
! local variables
integer i
! normal fermionic density matrix VV†
call zgemm('N','C',n,n,n,(1.d0,0.d0),v,n,v,n,(0.d0,0.d0),dvv,n)
do i=1,n
! store the V-norm
  vn(i)=dot_product(v(:,i),v(:,i))
! subtract unperturbed density matrix
  if (e(i).le.ef) dvv(i,i)=dvv(i,i)-1.d0
end do
! anomalous density matrix UV†
call zgemm('N','C',n,n,n,(1.d0,0.d0),u,n,v,n,(0.d0,0.d0),duv,n)
end subroutine

