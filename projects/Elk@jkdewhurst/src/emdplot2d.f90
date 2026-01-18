
! Copyright (C) 2014 D. Ernsting, S. Dugdale and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine emdplot2d(emds)
use modmain
use modpw
use modomp
implicit none
! arguments
real(4), intent(in) :: emds(nhkmax,nkpt)
! local variables
integer nh(3),np,ip,n,i,nthd
real(8) vpnl(3),v1(3),t1
! allocatable arrays
real(8), allocatable :: vpl(:,:),vppc(:,:)
real(8), allocatable :: x(:),wx(:),f(:)
! external functions
real(8), external :: rfhkintp
! allocate local arrays
np=np2d(1)*np2d(2)
allocate(vpl(3,np),vppc(2,np))
! generate the 2D plotting points
call plotpt2d(bvec,binv,vpnl,vpl,vppc)
! determine the number of integration points
nh(:)=int(hkmax*sqrt(avec(1,:)**2+avec(2,:)**2+avec(3,:)**2)/pi)+1
n=2*maxval(nh(:)*ngridk(:))
allocate(x(n),wx(n))
do i=1,n
  t1=2.d0*dble(i-1)/dble(n-1)-1.d0
  x(i)=t1*hkmax
end do
! determine the weights for spline integration
call wsplint(n,x,wx)
open(50,file='EMD2D.OUT',form='FORMATTED')
write(50,'(2I6," : grid size")') np2d(:)
! loop over plotting points in the 2D plane
call holdthd(np,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(f,i,v1,t1) &
!$OMP NUM_THREADS(nthd)
allocate(f(n))
!$OMP DO ORDERED
do ip=1,np
! integrate along normal to plane
  do i=1,n
    v1(:)=vpl(:,ip)+x(i)*vpnl(:)
    f(i)=rfhkintp(v1,emds)
  end do
  t1=dot_product(wx(:),f(:))
!$OMP ORDERED
  write(50,'(3G18.10)') vppc(1,ip),vppc(2,ip),t1
!$OMP END ORDERED
end do
!$OMP END DO
deallocate(f)
!$OMP END PARALLEL
call freethd(nthd)
close(50)
deallocate(vpl,vppc,x,wx)
end subroutine

