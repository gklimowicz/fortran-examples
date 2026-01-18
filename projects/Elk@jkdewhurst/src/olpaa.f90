
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaa(tor,is,ngp,apwalm,ld,o)
use modmain
implicit none
! arguments
logical, intent(in) :: tor
integer, intent(in) :: is,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: o(*)
! local variables
integer io,l,lm,i
! automatic arrays
complex(8) a(lmoapw(is),ngp)
i=0
do l=0,lmaxapw
  do lm=l**2+1,(l+1)**2
    do io=1,apword(l,is)
      i=i+1
      a(i,1:ngp)=apwalm(1:ngp,io,lm)
    end do
  end do
end do
call zmctmu(tor,lmoapw(is),ngp,a,a,ld,o)
end subroutine

