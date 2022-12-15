
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dmatls(dmat,xl,xs)
use modmain
implicit none
! arguments
complex(8), intent(in) :: dmat(lmmaxo,nspinor,lmmaxo,nspinor)
real(8), intent(out) :: xl(3),xs(3)
! local variables
integer ispn,lm
! automatic arrays
complex(8) zlflm(lmmaxo,3)
! compute tr(LD)
xl(:)=0.d0
do ispn=1,nspinor
  do lm=1,lmmaxo
    call lopzflm(lmaxo,dmat(:,ispn,lm,ispn),lmmaxo,zlflm)
    xl(:)=xl(:)+dble(zlflm(lm,:))
  end do
end do
! compute tr(sigma D)
xs(:)=0.d0
if (spinpol) then
  do lm=1,lmmaxo
    xs(1)=xs(1)+dble(dmat(lm,2,lm,1)+dmat(lm,1,lm,2))
    xs(2)=xs(2)+dble(-zi*dmat(lm,2,lm,1)+zi*dmat(lm,1,lm,2))
    xs(3)=xs(3)+dble(dmat(lm,1,lm,1)-dmat(lm,2,lm,2))
  end do
! S = 1/2 sigma
  xs(:)=0.5d0*xs(:)
end if
end subroutine

