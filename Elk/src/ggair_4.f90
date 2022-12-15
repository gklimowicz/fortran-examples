
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ggair_4(gvrho,dtdr,dtdgr2,wx,wc,vx,vc)
use modmain
implicit none
! arguments
real(8), intent(in) :: gvrho(ngtot,3)
real(8), intent(in) :: dtdr(ngtot),dtdgr2(ngtot)
real(8), intent(in) :: wx(ngtot),wc(ngtot)
real(8), intent(inout) :: vx(ngtot),vc(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
real(8), allocatable :: rfir(:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(rfir(ngtot),zfft1(ngtot),zfft2(ngtot))
!------------------!
!     exchange     !
!------------------!
vx(:)=vx(:)+wx(:)*dtdr(:)
rfir(:)=wx(:)*dtdgr2(:)
do i=1,3
  zfft1(:)=rfir(:)*gvrho(:,i)
  call zfftifc(3,ngridg,-1,zfft1)
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  vx(:)=vx(:)-2.d0*dble(zfft2(:))
end do
!---------------------!
!     correlation     !
!---------------------!
vc(:)=vc(:)+wc(:)*dtdr(:)
rfir(:)=wc(:)*dtdgr2(:)
do i=1,3
  zfft1(:)=rfir(:)*gvrho(:,i)
  call zfftifc(3,ngridg,-1,zfft1)
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  vc(:)=vc(:)-2.d0*dble(zfft2(:))
end do
deallocate(rfir,zfft1,zfft2)
end subroutine

