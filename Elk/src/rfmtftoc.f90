
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine rfmtftoc(nrc,nrci,rfmt,rfcmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nrc,nrci
real(8), intent(in) :: rfmt(*)
real(8), intent(out) :: rfcmt(*)
! local variables
integer irc,i,j,m,n
i=1
j=1
m=lmmaxi*lradstp
n=lmmaxi-1
do irc=1,nrci
  rfcmt(i:i+n)=rfmt(j:j+n)
  i=i+lmmaxi
  j=j+m
end do
j=j+(lradstp-1)*(lmmaxo-lmmaxi)
m=lmmaxo*lradstp
n=lmmaxo-1
do irc=nrci+1,nrc
  rfcmt(i:i+n)=rfmt(j:j+n)
  i=i+lmmaxo
  j=j+m
end do
end subroutine

