
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: moment
! !INTERFACE:
subroutine moment
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total moments by integrating the
!   magnetisation.
!
! !REVISION HISTORY:
!   Created January 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ias
real(8) t1
! external functions
real(8), external :: rfmtint
if (.not.spinpol) then
  mommt(:,:)=0.d0
  mommttot(:)=0.d0
  momir(:)=0.d0
  momtot(:)=0.d0
  momtotm=0.d0
  return
end if
! find the muffin-tin moments
mommttot(:)=0.d0
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    mommt(idm,ias)=rfmtint(nrmt(is),nrmti(is),wrmt(:,is),magmt(:,ias,idm))
    mommttot(idm)=mommttot(idm)+mommt(idm,ias)
  end do
end do
! find the interstitial moments
do idm=1,ndmag
  t1=dot_product(magir(1:ngtot,idm),cfunir(1:ngtot))
  momir(idm)=t1*omega/dble(ngtot)
end do
momtot(:)=mommttot(:)+momir(:)
! total moment magnitude
if (ncmag) then
  momtotm=sqrt(momtot(1)**2+momtot(2)**2+momtot(3)**2)
else
  momtotm=abs(momtot(1))
end if
! write total moment magnitude to test file
call writetest(450,'total moment magnitude',tol=1.d-3,rv=momtotm)
end subroutine
!EOC

