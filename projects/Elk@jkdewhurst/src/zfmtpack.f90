
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine zfmtpack(tpack,nr,nri,zfmt1,zfmt2)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(in) :: nr,nri
complex(8), intent(in) :: zfmt1(*)
complex(8), intent(out) :: zfmt2(*)
! local variables
integer ir,i,j,k,n
n=lmmaxi-1
i=1
j=1
if (tpack) then
  do ir=1,nri
    zfmt2(j:j+n)=zfmt1(i:i+n)
    i=i+lmmaxo
    j=j+lmmaxi
  end do
else
  do ir=1,nri
    zfmt2(j:j+n)=zfmt1(i:i+n)
    i=i+lmmaxi
    k=j+lmmaxi
    j=j+lmmaxo
    zfmt2(k:j-1)=0.d0
  end do
end if
n=lmmaxo*(nr-nri)-1
zfmt2(j:j+n)=zfmt1(i:i+n)
end subroutine

