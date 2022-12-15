
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine dhmlalo(is,ias,ngp,ngpq,apwalm,apwalmq,dapwalm,dapwalmq,ld,dh)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: is,ias,ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: dh(ld,*)
! local variables
integer io,ilo
integer l1,l2,l3
integer lm1,lm2,lm3
integer i0,j0,i,j
complex(8) z1
do ilo=1,nlorb(is)
  l1=lorbl(ilo,is)
  do lm1=l1**2+1,(l1+1)**2
    i=idxlo(lm1,ilo,ias)
    i0=ngpq+i
    j0=ngp+i
    do l3=0,lmaxapw
      do lm3=l3**2+1,(l3+1)**2
        do io=1,apword(l3,is)
          z1=0.d0
          do l2=0,lmaxo
            if (mod(l1+l2+l3,2).eq.0) then
              do lm2=l2**2+1,(l2+1)**2
                z1=z1+gntyyy(lm2,lm3,lm1)*dhloa(lm2,io,l3,ilo,ias)
              end do
            end if
          end do
          if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
            do i=1,ngpq
              dh(i,j0)=dh(i,j0)+conjg(z1*apwalmq(i,io,lm3))
            end do
            do j=1,ngp
              dh(i0,j)=dh(i0,j)+z1*apwalm(j,io,lm3)
            end do
          end if
        end do
      end do
    end do
    if (ias.eq.iasph) then
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
            if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
              do i=1,ngpq
                dh(i,j0)=dh(i,j0)+conjg(z1*dapwalmq(i,io,lm3))
              end do
              do j=1,ngp
                dh(i0,j)=dh(i0,j)+z1*dapwalm(j,io,lm3)
              end do
            end if
          end do
        end do
      end do
    end if
  end do
end do
end subroutine

