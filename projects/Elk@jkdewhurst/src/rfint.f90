
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfint(rfmt,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias
! external functions
real(8), external :: rfmtint
! interstitial contribution
rfint=dot_product(rfir(:),cfunir(:))
rfint=rfint*omega/dble(ngtot)
! muffin-tin contribution
do ias=1,natmtot
  is=idxis(ias)
  rfint=rfint+rfmtint(nrmt(is),nrmti(is),wrmt(:,is),rfmt(:,ias))
end do
end function

