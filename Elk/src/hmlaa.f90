
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlaa
! !INTERFACE:
subroutine hmlaa(thr,is,ias,ngp,apwalm,ld,h)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   thr    : .true. if the matrix h is real valued (in,logical)
!   is     : species number (in,integer)
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients (in,complex(ngkmax,apwordmax,lmmaxapw))
!   ld     : leading dimension of h (in,integer)
!   h      : Hamiltonian matrix (inout,complex(*))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: thr
integer, intent(in) :: is,ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: h(*)
! local variables
integer io,jo,i
integer l1,l2,l3
integer lm1,lm2,lm3
real(8) t0
complex(8) z1
! automatic arrays
complex(8) a(lmoapw(is),ngp),b(lmoapw(is),ngp)
t0=0.5d0*rmt(is)**2
i=0
do l1=0,lmaxapw
  do lm1=l1**2+1,(l1+1)**2
    do io=1,apword(l1,is)
      i=i+1
      b(i,:)=0.d0
      do l3=0,lmaxapw
        do lm3=l3**2+1,(l3+1)**2
          do jo=1,apword(l3,is)
            z1=0.d0
! kinetic and potential contribution
            do l2=0,lmaxo
              if (mod(l1+l2+l3,2).eq.0) then
                do lm2=l2**2+1,(l2+1)**2
                  z1=z1+gntyry(lm2,lm3,lm1)*haa(lm2,jo,l3,io,l1,ias)
                end do
              end if
            end do
! kinetic surface contribution
            if (lm1.eq.lm3) then
              z1=z1+t0*apwfr(nrmt(is),1,io,l1,ias)*apwdfr(jo,l1,ias)
            end if
            if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
              b(i,1:ngp)=b(i,1:ngp)+z1*apwalm(1:ngp,jo,lm3)
            end if
          end do
        end do
      end do
      a(i,1:ngp)=apwalm(1:ngp,io,lm1)
    end do
  end do
end do
call zmctmu(thr,lmoapw(is),ngp,a,b,ld,h)
end subroutine
!EOC

