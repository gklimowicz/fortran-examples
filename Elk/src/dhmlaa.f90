
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dhmlaa(is,ias,ngp,ngpq,apwalm,apwalmq,dapwalm,dapwalmq,ld,dh)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: is,ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: dh(*)
! local variables
integer lmo,io,jo,i
integer l1,l2,l3
integer lm1,lm2,lm3
real(8) t0
complex(8) z1
! automatic arrays
complex(8) a1(lmoapw(is),ngpq),b1(lmoapw(is),ngp)
complex(8) a2(lmoapw(is),ngpq),b2(lmoapw(is),ngp)
lmo=lmoapw(is)
t0=0.5d0*rmt(is)**2
i=0
do l1=0,lmaxapw
  do lm1=l1**2+1,(l1+1)**2
    do io=1,apword(l1,is)
      i=i+1
      b1(i,:)=0.d0
      do l3=0,lmaxapw
        do lm3=l3**2+1,(l3+1)**2
          do jo=1,apword(l3,is)
            z1=0.d0
            do l2=0,lmaxo
              if (mod(l1+l2+l3,2).eq.0) then
                do lm2=l2**2+1,(l2+1)**2
                  z1=z1+gntyyy(lm2,lm3,lm1)*dhaa(lm2,jo,l3,io,l1,ias)
                end do
              end if
            end do
            if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
              call zaxpy(ngp,z1,apwalm(:,jo,lm3),1,b1(i,1),lmo)
            end if
          end do
        end do
      end do
      if (ias.eq.iasph) then
        b2(i,:)=0.d0
        do l3=0,lmaxapw
          do lm3=l3**2+1,(l3+1)**2
            do jo=1,apword(l3,is)
              z1=0.d0
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
                call zaxpy(ngp,z1,dapwalm(:,jo,lm3),1,b1(i,1),lmo)
                call zaxpy(ngp,z1,apwalm(:,jo,lm3),1,b2(i,1),lmo)
              end if
            end do
          end do
        end do
        a2(i,1:ngpq)=dapwalmq(1:ngpq,io,lm1)
      end if
      a1(i,1:ngpq)=apwalmq(1:ngpq,io,lm1)
    end do
  end do
end do
call zmctm(lmo,ngpq,ngp,a1,b1,ld,dh)
if (ias.eq.iasph) then
  call zmctm(lmo,ngpq,ngp,a2,b2,ld,dh)
end if
end subroutine

