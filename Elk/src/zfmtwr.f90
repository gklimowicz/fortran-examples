
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine zfmtwr(nr,nri,wr,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
complex(8), intent(inout) :: zfmt(*)
! local variables
integer n,ir,i
i=1
if (lmaxi.eq.1) then
  do ir=1,nri
    zfmt(i:i+3)=wr(ir)*zfmt(i:i+3)
    i=i+4
  end do
else
  n=lmmaxi-1
  do ir=1,nri
    zfmt(i:i+n)=wr(ir)*zfmt(i:i+n)
    i=i+lmmaxi
  end do
end if
n=lmmaxo-1
do ir=nri+1,nr
  zfmt(i:i+n)=wr(ir)*zfmt(i:i+n)
  i=i+lmmaxo
end do
end subroutine

