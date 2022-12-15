
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure real(8) function rfmtint(nr,nri,wr,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr),rfmt(*)
! local variables
integer iro,i0,i1
i1=lmmaxi*(nri-1)+1
rfmtint=sum(wr(1:nri)*rfmt(1:i1:lmmaxi))
iro=nri+1
i0=i1+lmmaxi
i1=lmmaxo*(nr-iro)+i0
rfmtint=rfmtint+sum(wr(iro:nr)*rfmt(i0:i1:lmmaxo))
rfmtint=fourpi*y00*rfmtint
end function

