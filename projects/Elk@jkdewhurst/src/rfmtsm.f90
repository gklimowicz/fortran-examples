
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtsm(m,nr,nri,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: m,nr,nri
real(8), intent(inout) :: rfmt(*)
! local variables
integer nro,iro,lm
integer i1,j0,j1
! automatic arrays
real(8) fr(nr)
if (m.le.0) return
nro=nr-nri
iro=nri+1
do lm=1,lmmaxi
  i1=lmmaxi*(nri-1)+lm
  j0=i1+lmmaxi
  j1=lmmaxo*(nr-iro)+j0
  fr(1:nri)=rfmt(lm:i1:lmmaxi)
  fr(iro:nr)=rfmt(j0:j1:lmmaxo)
  call fsmooth(m,nr,fr)
  rfmt(lm:i1:lmmaxi)=fr(1:nri)
  rfmt(j0:j1:lmmaxo)=fr(iro:nr)
end do
do lm=lmmaxi+1,lmmaxo
  j0=lmmaxi*nri+lm
  j1=lmmaxo*(nr-iro)+j0
  fr(iro:nr)=rfmt(j0:j1:lmmaxo)
  call fsmooth(m,nro,fr(iro))
  rfmt(j0:j1:lmmaxo)=fr(iro:nr)
end do
end subroutine

