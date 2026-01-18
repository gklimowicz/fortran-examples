
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine dolpalo(is,ias,ngp,ngpq,dapwalm,dapwalmq,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: is,ias,ngp,ngpq
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: od(ld,*)
! local variables
integer ilo,io,l,lm
integer i0,j0,i,j
real(8) t1
if (ias.ne.iasph) return
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do lm=l**2+1,(l+1)**2
    i=idxlo(lm,ilo,ias)
    i0=ngpq+i
    j0=ngp+i
    do io=1,apword(l,is)
      t1=oalo(io,ilo,ias)
      do i=1,ngpq
        od(i,j0)=od(i,j0)+t1*conjg(dapwalmq(i,io,lm))
      end do
      do j=1,ngp
        od(i0,j)=od(i0,j)+t1*dapwalm(j,io,lm)
      end do
    end do
  end do
end do
end subroutine

