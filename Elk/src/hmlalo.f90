
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine hmlalo(is,ias,ngp,apwalm,ld,h)
use modmain
implicit none
! arguments
integer, intent(in) :: is,ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: h(ld,*)
! local variables
integer io,ilo,i,j
integer l1,l2,l3
integer lm1,lm2,lm3
complex(8) z1
do ilo=1,nlorb(is)
  l1=lorbl(ilo,is)
  do lm1=l1**2+1,(l1+1)**2
    j=ngp+idxlo(lm1,ilo,ias)
    do l3=0,lmaxapw
      do lm3=l3**2+1,(l3+1)**2
        do io=1,apword(l3,is)
          z1=0.d0
          do l2=0,lmaxo
            if (mod(l1+l2+l3,2).eq.0) then
              do lm2=l2**2+1,(l2+1)**2
                z1=z1+gntyry(lm2,lm3,lm1)*hloa(lm2,io,l3,ilo,ias)
              end do
            end if
          end do
! note that what is actually computed is the Hermitian conjugate of <lo|H|APW>
          if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
            do i=1,ngp
              h(i,j)=h(i,j)+conjg(z1*apwalm(i,io,lm3))
            end do
          end if
        end do
      end do
    end do
  end do
end do
end subroutine

