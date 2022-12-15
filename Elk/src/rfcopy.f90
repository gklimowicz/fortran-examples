
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine rfcopy(rfmt1,rfir1,rfmt2,rfir2)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt1(npmtmax,natmtot),rfir1(ngtot)
real(8), intent(out) :: rfmt2(npmtmax,natmtot),rfir2(ngtot)
! local variables
integer is,ias
do ias=1,natmtot
  is=idxis(ias)
  rfmt2(1:npmt(is),ias)=rfmt1(1:npmt(is),ias)
end do
rfir2(:)=rfir1(:)
end subroutine

