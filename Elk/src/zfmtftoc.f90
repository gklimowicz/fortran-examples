
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine zfmtftoc(nrc,nrci,zfmt,zfcmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nrc,nrci
complex(8), intent(in) :: zfmt(*)
complex(8), intent(out) :: zfcmt(*)
! local variables
integer irc,i,j,m,n
i=1
j=1
m=lmmaxi*lradstp
n=lmmaxi-1
do irc=1,nrci
  zfcmt(i:i+n)=zfmt(j:j+n)
  i=i+lmmaxi
  j=j+m
end do
j=j+(lradstp-1)*(lmmaxo-lmmaxi)
m=lmmaxo*lradstp
n=lmmaxo-1
do irc=nrci+1,nrc
  zfcmt(i:i+n)=zfmt(j:j+n)
  i=i+lmmaxo
  j=j+m
end do
end subroutine

