
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: charge
! !INTERFACE:
subroutine charge
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total charges by integrating the
!   density.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias
real(8) t1
! external functions
real(8), external :: rfmtint
! find the muffin-tin charges
chgmttot=0.d0
do ias=1,natmtot
  is=idxis(ias)
  chgmt(ias)=rfmtint(nrmt(is),nrmti(is),wrmt(:,is),rhomt(:,ias))
  chgmttot=chgmttot+chgmt(ias)
end do
! find the interstitial charge
t1=dot_product(rhoir(1:ngtot),cfunir(1:ngtot))
chgir=t1*omega/dble(ngtot)
! total calculated charge
chgcalc=chgmttot+chgir
! write total calculated charge to test file
call writetest(400,'calculated total charge',tol=1.d-6,rv=chgcalc)
end subroutine
!EOC

