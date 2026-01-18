
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure complex(8) function rzfmtinp(nr,nri,wr,rfmt,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
real(8), intent(in) :: rfmt(*)
complex(8), intent(in) :: zfmt(*)
! local variables
integer n,ir,i
complex(8) z1,z2
! compute the dot-products for each radial point and integrate over r
z1=0.d0
i=1
if (lmaxi.eq.1) then
  do ir=1,nri
    z1=z1+wr(ir) &
     *(rfmt(i)*zfmt(i) &
      +rfmt(i+1)*zfmt(i+1) &
      +rfmt(i+2)*zfmt(i+2) &
      +rfmt(i+3)*zfmt(i+3))
    i=i+4
  end do
  z1=pi*z1
else
  n=lmmaxi-1
  do ir=1,nri
    z1=z1+wr(ir)*dot_product(rfmt(i:i+n),zfmt(i:i+n))
    i=i+lmmaxi
  end do
  z1=(fourpi/dble(lmmaxi))*z1
end if
z2=0.d0
n=lmmaxo-1
do ir=nri+1,nr
  z2=z2+wr(ir)*dot_product(rfmt(i:i+n),zfmt(i:i+n))
  i=i+lmmaxo
end do
rzfmtinp=z1+(fourpi/dble(lmmaxo))*z2
end function

