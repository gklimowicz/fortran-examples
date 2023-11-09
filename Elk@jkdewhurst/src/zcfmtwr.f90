
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine zcfmtwr(nr,nri,wr,zfmt,cfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
complex(8), intent(in) :: zfmt(*)
complex(4), intent(out) :: cfmt(*)
! local variables
integer n,ir,i
i=1
if (lmaxi.eq.1) then
  do ir=1,nri
    cfmt(i:i+3)=cmplx(wr(ir)*zfmt(i:i+3))
    i=i+4
  end do
else
  n=lmmaxi-1
  do ir=1,nri
    cfmt(i:i+n)=cmplx(wr(ir)*zfmt(i:i+n))
    i=i+lmmaxi
  end do
end if
n=lmmaxo-1
do ir=nri+1,nr
  cfmt(i:i+n)=cmplx(wr(ir)*zfmt(i:i+n))
  i=i+lmmaxo
end do
end subroutine

