
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaaq(is,ngp,ngpq,apwalm,apwalmq,ld,oq)
use modmain
implicit none
! arguments
integer, intent(in) :: is,ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: oq(*)
! local variables
integer io,l,lm,i
! automatic arrays
complex(8) a(lmoapw(is),ngpq),b(lmoapw(is),ngp)
i=0
do l=0,lmaxapw
  do lm=l**2+1,(l+1)**2
    do io=1,apword(l,is)
      i=i+1
      a(i,1:ngpq)=apwalmq(1:ngpq,io,lm)
      b(i,1:ngp)=apwalm(1:ngp,io,lm)
    end do
  end do
end do
call zmctm(lmoapw(is),ngpq,ngp,a,b,ld,oq)
end subroutine

