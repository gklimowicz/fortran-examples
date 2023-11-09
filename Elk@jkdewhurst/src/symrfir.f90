
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrfir
! !INTERFACE:
subroutine symrfir(ngdg,ngt,ngv,igf,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngdg : G-vector grid sizes (in,integer(3))
!   ngt  : total number of G-vectors (in,integer)
!   ngv  : number of G-vectors within cut-off (in,integer)
!   igf  : map from G-vector index to FFT array (in,integer(ngv))
!   rfir : real intersitial function (inout,real(ngt))
! !DESCRIPTION:
!   Symmetrises a real scalar interstitial function. The function is first
!   Fourier transformed to $G$-space, and then averaged over each symmetry by
!   rotating the Fourier coefficients and multiplying them by a phase factor
!   corresponding to the symmetry translation.
!
! !REVISION HISTORY:
!   Created July 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngdg(3),ngt,ngv,igf(ngv)
real(8), intent(inout) :: rfir(ngt)
! local variables
logical tv0
integer isym,lspl,ilspl
integer sym(3,3),ig,ifg,jfg
integer i1,i2,i3,j1,j2,j3
real(8) v1,v2,v3,t1
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(zfft1(ngt),zfft2(ngt))
! Fourier transform function to G-space
zfft1(:)=rfir(:)
call zfftifc(3,ngdg,-1,zfft1)
zfft2(:)=0.d0
! loop over crystal symmetries
do isym=1,nsymcrys
! zero translation vector flag
  tv0=tv0symc(isym)
! translation vector in Cartesian coordinates
  if (.not.tv0) then
    v1=vtcsymc(1,isym)
    v2=vtcsymc(2,isym)
    v3=vtcsymc(3,isym)
  end if
! index to lattice symmetry of spatial rotation
  lspl=lsplsymc(isym)
! inverse rotation required for rotation of G-vectors
  ilspl=isymlat(lspl)
  sym(:,:)=symlat(:,:,ilspl)
  do ig=1,ngv
    ifg=igf(ig)
! multiply the transpose of the inverse symmetry matrix with the G-vector
    if (lspl.eq.1) then
      jfg=ifg
    else
      i1=ivg(1,ig); i2=ivg(2,ig); i3=ivg(3,ig)
      j1=sym(1,1)*i1+sym(2,1)*i2+sym(3,1)*i3
      j2=sym(1,2)*i1+sym(2,2)*i2+sym(3,2)*i3
      j3=sym(1,3)*i1+sym(2,3)*i2+sym(3,3)*i3
      jfg=igf(ivgig(j1,j2,j3))
    end if
    if (tv0) then
! zero translation vector
      zfft2(jfg)=zfft2(jfg)+zfft1(ifg)
    else
! complex phase factor for translation
      t1=-(vgc(1,ig)*v1+vgc(2,ig)*v2+vgc(3,ig)*v3)
      z1=cmplx(cos(t1),sin(t1),8)
      zfft2(jfg)=zfft2(jfg)+z1*zfft1(ifg)
    end if
  end do
end do
! Fourier transform to real-space and normalise
call zfftifc(3,ngdg,1,zfft2)
t1=1.d0/dble(nsymcrys)
rfir(:)=t1*dble(zfft2(:))
deallocate(zfft1,zfft2)
end subroutine
!EOC

