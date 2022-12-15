
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrvfir
subroutine symrvfir(tspin,tnc,ngdg,ngt,ngv,igf,ld,rvfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tspin : .true. if spin rotations should be used (in,logical)
!   tnc   : .true. if the vector field is non-collinear, otherwise it is
!           collinear along the z-axis (in,logical)
!   ngdg  : G-vector grid sizes (in,integer(3))
!   ngt   : total number of G-vectors (in,integer)
!   ngv   : number of G-vectors within cut-off (in,integer)
!   igf   : map from G-vector index to FFT array (in,integer(ngv))
!   ld    : leading dimension (in,integer)
!   rvfir : real interstitial vector function (inout,real(ld,*))
! !DESCRIPTION:
!   Symmetrises a real interstitial vector function. See routines {\tt symrvf}
!   and {\tt symrfir} for details.
!
! !REVISION HISTORY:
!   Created July 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tspin,tnc
integer, intent(in) :: ngdg(3),ngt,ngv,igf(ngv)
integer, intent(in) :: ld
real(8), intent(inout) :: rvfir(ld,*)
! local variables
logical tv0
integer nd,isym,lspl,ilspl,lspn
integer sym(3,3),ig,ifg,jfg
integer i1,i2,i3,j1,j2,j3,i
real(8) sc(3,3),v1,v2,v3,t1
complex(8) x1,x2,x3,z1
! allocatable arrays
complex(8), allocatable :: zfft1(:,:),zfft2(:,:)
! dimension of the vector field
if (tnc) then
  nd=3
else
  nd=1
end if
allocate(zfft1(ngt,nd),zfft2(ngt,nd))
! Fourier transform vector function to G-space
do i=1,nd
  zfft1(1:ngt,i)=rvfir(1:ngt,i)
  call zfftifc(3,ngdg,-1,zfft1(:,i))
end do
zfft2(:,:)=0.d0
do isym=1,nsymcrys
! zero translation vector flag
  tv0=tv0symc(isym)
! translation vector in Cartesian coordinates
  if (.not.tv0) then
    v1=vtcsymc(1,isym)
    v2=vtcsymc(2,isym)
    v3=vtcsymc(3,isym)
  end if
! index to spatial rotation lattice symmetry
  lspl=lsplsymc(isym)
! inverse rotation required for rotation of G-vectors
  ilspl=isymlat(lspl)
  sym(:,:)=symlat(:,:,ilspl)
  if (tspin) then
! global spin proper rotation in Cartesian coordinates
    lspn=lspnsymc(isym)
    sc(:,:)=symlatd(lspn)*symlatc(:,:,lspn)
  else
! set spin rotation equal to spatial rotation
    lspn=lspl
    sc(:,:)=symlatc(:,:,lspl)
  end if
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
! translation, spatial rotation and global spin rotation
    if (tv0) then
! zero translation vector
      if (lspn.eq.1) then
! global spin symmetry is the identity
        zfft2(jfg,:)=zfft2(jfg,:)+zfft1(ifg,:)
      else
        if (tnc) then
! non-collinear case
          x1=zfft1(ifg,1); x2=zfft1(ifg,2); x3=zfft1(ifg,3)
          zfft2(jfg,1)=zfft2(jfg,1)+sc(1,1)*x1+sc(1,2)*x2+sc(1,3)*x3
          zfft2(jfg,2)=zfft2(jfg,2)+sc(2,1)*x1+sc(2,2)*x2+sc(2,3)*x3
          zfft2(jfg,3)=zfft2(jfg,3)+sc(3,1)*x1+sc(3,2)*x2+sc(3,3)*x3
        else
! collinear case
          zfft2(jfg,1)=zfft2(jfg,1)+sc(3,3)*zfft1(ifg,1)
        end if
      end if
    else
! complex phase factor for translation
      t1=-(vgc(1,ig)*v1+vgc(2,ig)*v2+vgc(3,ig)*v3)
      z1=cmplx(cos(t1),sin(t1),8)
      if (lspn.eq.1) then
        zfft2(jfg,:)=zfft2(jfg,:)+z1*zfft1(ifg,:)
      else
        if (tnc) then
          x1=zfft1(ifg,1); x2=zfft1(ifg,2); x3=zfft1(ifg,3)
          zfft2(jfg,1)=zfft2(jfg,1)+z1*(sc(1,1)*x1+sc(1,2)*x2+sc(1,3)*x3)
          zfft2(jfg,2)=zfft2(jfg,2)+z1*(sc(2,1)*x1+sc(2,2)*x2+sc(2,3)*x3)
          zfft2(jfg,3)=zfft2(jfg,3)+z1*(sc(3,1)*x1+sc(3,2)*x2+sc(3,3)*x3)
        else
          zfft2(jfg,1)=zfft2(jfg,1)+sc(3,3)*z1*zfft1(ifg,1)
        end if
      end if
    end if
  end do
end do
! Fourier transform to real-space and normalise
t1=1.d0/dble(nsymcrys)
do i=1,nd
  call zfftifc(3,ngdg,1,zfft2(:,i))
  rvfir(1:ngt,i)=t1*dble(zfft2(1:ngt,i))
end do
deallocate(zfft1,zfft2)
end subroutine
!EOC

