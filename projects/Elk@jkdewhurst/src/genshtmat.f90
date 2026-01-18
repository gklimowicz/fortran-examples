
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genshtmat
! !INTERFACE:
subroutine genshtmat
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the forward and backward spherical harmonic transformation (SHT)
!   matrices using the spherical covering set produced by the routine
!   {\tt sphcover}. These matrices are used to transform a function between its
!   $(l,m)$-expansion coefficients and its values at the $(\theta,\phi)$ points
!   on the sphere. Both real and complex SHT matrices are calculated and stored
!   in global arrays.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer itp
real(8) v(3)
! automatic arrays
real(8) tp(2,lmmaxo),vtp(3,lmmaxo),rlm(lmmaxo)
complex(8) ylm(lmmaxo)
!--------------------------------!
!     SHT matrices for lmaxi     !
!--------------------------------!
! allocate real SHT matrices
if (allocated(rbshti)) deallocate(rbshti)
allocate(rbshti(lmmaxi,lmmaxi))
if (allocated(rfshti)) deallocate(rfshti)
allocate(rfshti(lmmaxi,lmmaxi))
! allocate complex SHT matrices
if (allocated(zbshti)) deallocate(zbshti)
allocate(zbshti(lmmaxi,lmmaxi))
if (allocated(zfshti)) deallocate(zfshti)
allocate(zfshti(lmmaxi,lmmaxi))
! allocate single-precision complex copies
if (allocated(cbshti)) deallocate(cbshti)
allocate(cbshti(lmmaxi,lmmaxi))
if (allocated(cfshti)) deallocate(cfshti)
allocate(cfshti(lmmaxi,lmmaxi))
! generate spherical covering set for lmaxi
call sphcover(lmmaxi,tp)
! convert (theta, phi) angles to vectors
do itp=1,lmmaxi
  call sctovec(tp(:,itp),vtp(:,itp))
end do
! rotate the spherical covering set if required
if (trotsht) then
  do itp=1,lmmaxi
    v(:)=vtp(:,itp)
    call r3mv(rotsht,v,vtp(:,itp))
  end do
end if
! generate real and complex spherical harmonics and set the backward SHT arrays
do itp=1,lmmaxi
  call genrlmv(lmaxi,vtp(:,itp),rlm)
  rbshti(itp,1:lmmaxi)=rlm(1:lmmaxi)
  call genylmv(lmaxi,vtp(:,itp),ylm)
  zbshti(itp,1:lmmaxi)=ylm(1:lmmaxi)
end do
! find the forward SHT arrays
! real
rfshti(:,:)=rbshti(:,:)
call rminv(lmmaxi,rfshti)
! complex
zfshti(:,:)=zbshti(:,:)
call zminv(lmmaxi,zfshti)
! make single-precision complex copies
cbshti(:,:)=zbshti(:,:)
cfshti(:,:)=zfshti(:,:)
!--------------------------------!
!     SHT matrices for lmaxo     !
!--------------------------------!
! allocate real SHT matrices
if (allocated(rbshto)) deallocate(rbshto)
allocate(rbshto(lmmaxo,lmmaxo))
if (allocated(rfshto)) deallocate(rfshto)
allocate(rfshto(lmmaxo,lmmaxo))
! allocate complex SHT matrices
if (allocated(zbshto)) deallocate(zbshto)
allocate(zbshto(lmmaxo,lmmaxo))
if (allocated(zfshto)) deallocate(zfshto)
allocate(zfshto(lmmaxo,lmmaxo))
! allocate single-precision complex copies
if (allocated(cbshto)) deallocate(cbshto)
allocate(cbshto(lmmaxo,lmmaxo))
if (allocated(cfshto)) deallocate(cfshto)
allocate(cfshto(lmmaxo,lmmaxo))
! generate spherical covering set
call sphcover(lmmaxo,tp)
! convert (theta, phi) angles to vectors
do itp=1,lmmaxo
  call sctovec(tp(:,itp),vtp(:,itp))
end do
! rotate the spherical covering set if required
if (trotsht) then
  do itp=1,lmmaxo
    v(:)=vtp(:,itp)
    call r3mv(rotsht,v,vtp(:,itp))
  end do
end if
! generate real and complex spherical harmonics and set the backward SHT arrays
do itp=1,lmmaxo
  call genrlmv(lmaxo,vtp(:,itp),rlm)
  rbshto(itp,1:lmmaxo)=rlm(1:lmmaxo)
  call genylmv(lmaxo,vtp(:,itp),ylm)
  zbshto(itp,1:lmmaxo)=ylm(1:lmmaxo)
end do
! find the forward SHT arrays
! real
rfshto(:,:)=rbshto(:,:)
call rminv(lmmaxo,rfshto)
! complex
zfshto(:,:)=zbshto(:,:)
call zminv(lmmaxo,zfshto)
! make single-precision complex copies
cbshto(:,:)=zbshto(:,:)
cfshto(:,:)=zfshto(:,:)
end subroutine
!EOC

