
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfint0(rf0,rfmt,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rf0
real(8), intent(inout) :: rfmt(npmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias,nr,nri
integer iro,i0,i1
real(8) t1
! external functions
real(8), external :: rfint
t1=rfint(rfmt,rfir)
t1=rf0-t1/omega
rfir(:)=rfir(:)+t1
t1=t1/y00
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  i1=lmmaxi*(nri-1)+1
  rfmt(1:i1:lmmaxi,ias)=rfmt(1:i1:lmmaxi,ias)+t1
  i0=i1+lmmaxi
  i1=lmmaxo*(nr-iro)+i0
  rfmt(i0:i1:lmmaxo,ias)=rfmt(i0:i1:lmmaxo,ias)+t1
end do
end subroutine

