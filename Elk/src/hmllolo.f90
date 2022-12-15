
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine hmllolo(is,ias,ngp,ld,h)
use modmain
implicit none
! arguments
integer, intent(in) :: is,ias,ngp,ld
complex(8), intent(inout) :: h(ld,*)
! local variables
integer ilo,jlo,i,j
integer l1,l2,l3
integer lm1,lm2,lm3
complex(8) z1
do jlo=1,nlorb(is)
  l3=lorbl(jlo,is)
  do lm3=l3**2+1,(l3+1)**2
    j=ngp+idxlo(lm3,jlo,ias)
    do ilo=1,nlorb(is)
      l1=lorbl(ilo,is)
      do lm1=l1**2+1,(l1+1)**2
        i=ngp+idxlo(lm1,ilo,ias)
        if (i.gt.j) cycle
        z1=0.d0
        do l2=0,lmaxo
          if (mod(l1+l2+l3,2).eq.0) then
            do lm2=l2**2+1,(l2+1)**2
              z1=z1+gntyry(lm2,lm3,lm1)*hlolo(lm2,jlo,ilo,ias)
            end do
          end if
        end do
        h(i,j)=h(i,j)+z1
      end do
    end do
  end do
end do
end subroutine

