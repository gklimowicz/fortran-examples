
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradzf(zfmt,zfir,gzfmt,gzfir)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) :: zfmt(npmtmax,natmtot),zfir(ngtot)
complex(8), intent(out) :: gzfmt(npmtmax,natmtot,3),gzfir(ngtot,3)
! local variables
integer is,ias,ld,i
integer ig,ifg,nthd
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zfft(:)
! muffin-tin gradient
ld=npmtmax*natmtot
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call gradzfmt(nrmt(is),nrmti(is),rlmt(:,-1,is),wcrmt(:,:,is),zfmt(:,ias),ld, &
   gzfmt(1,ias,1))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! interstitial gradient
allocate(zfft(ngtot))
zfft(:)=zfir(:)
call zfftifc(3,ngridg,-1,zfft)
call holdthd(3,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ig,ifg,z1) &
!$OMP NUM_THREADS(nthd)
do i=1,3
  gzfir(:,i)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    z1=zfft(ifg)
    gzfir(ifg,i)=vgc(i,ig)*cmplx(-aimag(z1),dble(z1),8)
  end do
  call zfftifc(3,ngridg,1,gzfir(:,i))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
deallocate(zfft)
end subroutine

