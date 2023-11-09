
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine rfmtlm(lm,nr,nri,rfmt,fr)
use modmain
implicit none
! arguments
integer, intent(in) :: lm,nr,nri
real(8), intent(in) :: rfmt(npmtmax)
real(8), intent(out) :: fr(nrmtmax)
! local variables
integer iro,i0,i1
if (lm.gt.lmmaxi) then
  fr(1:nri)=0.d0
else
  i1=lmmaxi*(nri-1)+lm
  fr(1:nri)=rfmt(lm:i1:lmmaxi)
end if
iro=nri+1
if (lm.gt.lmmaxo) then
  fr(iro:nr)=0.d0
else
  i0=lmmaxi*nri+lm
  i1=lmmaxo*(nr-iro)+i0
  fr(iro:nr)=rfmt(i0:i1:lmmaxo)
end if
end subroutine

