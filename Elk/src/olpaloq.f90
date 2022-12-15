
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaloq(is,ias,ngp,ngpq,apwalm,apwalmq,ld,oq)
use modmain
implicit none
! arguments
integer, intent(in) :: is,ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: oq(ld,*)
! local variables
integer ilo,io,l,lm
integer i0,j0,i,j
real(8) t1
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do lm=l**2+1,(l+1)**2
    i=idxlo(lm,ilo,ias)
    i0=ngpq+i
    j0=ngp+i
    do io=1,apword(l,is)
      t1=oalo(io,ilo,ias)
      do i=1,ngpq
        oq(i,j0)=oq(i,j0)+t1*conjg(apwalmq(i,io,lm))
      end do
      do j=1,ngp
        oq(i0,j)=oq(i0,j)+t1*apwalm(j,io,lm)
      end do
    end do
  end do
end do
end subroutine

