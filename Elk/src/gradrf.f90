
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrf(rfmt,rfir,grfmt,grfir)
use modmain
use modomp
implicit none
! arguments
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
real(8), intent(out) :: grfmt(npmtmax,natmtot,3),grfir(ngtot,3)
! local variables
integer is,ias,ld,i
integer ig,ifg,nthd
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
! muffin-tin gradient
ld=npmtmax*natmtot
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call gradrfmt(nrmt(is),nrmti(is),rlmt(:,-1,is),wcrmt(:,:,is),rfmt(:,ias),ld, &
   grfmt(1,ias,1))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! interstitial gradient
allocate(zfft1(ngtot),zfft2(ngtot))
zfft1(:)=rfir(:)
call zfftifc(3,ngridg,-1,zfft1)
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    z1=zfft1(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(z1),dble(z1),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  grfir(:,i)=dble(zfft2(:))
end do
deallocate(zfft1,zfft2)
end subroutine

