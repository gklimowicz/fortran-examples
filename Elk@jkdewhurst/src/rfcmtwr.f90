
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine rfcmtwr(nr,nri,wr,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
real(8), intent(inout) :: rfmt(*)
! local variables
integer n,ir,i
real(8) t1
i=1
if (lmaxi.eq.1) then
  do ir=1,nri
    rfmt(i:i+3)=(pi*wr(ir))*rfmt(i:i+3)
    i=i+4
  end do
else
  t1=fourpi/dble(lmmaxi)
  n=lmmaxi-1
  do ir=1,nri
    rfmt(i:i+n)=(t1*wr(ir))*rfmt(i:i+n)
    i=i+lmmaxi
  end do
end if
t1=fourpi/dble(lmmaxo)
n=lmmaxo-1
do ir=nri+1,nr
  rfmt(i:i+n)=(t1*wr(ir))*rfmt(i:i+n)
  i=i+lmmaxo
end do
end subroutine

