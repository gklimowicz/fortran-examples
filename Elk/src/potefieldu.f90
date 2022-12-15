
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potefieldu
use modmain
use modulr
implicit none
! local variables
integer ir
real(8) v0,v(3)
! allocatable arrays
real(8), allocatable :: rfft(:)
if (sum(abs(efielduc(:))).lt.epslat) return
allocate(rfft(nqpt))
! constant added to potential so that it is zero at the ultracell center
v(:)=0.5d0*(avecu(:,1)+avecu(:,2)+avecu(:,3))
v0=dot_product(efielduc(:),v(:))
! calculate the potential in real-space
do ir=1,nqpt
  rfft(ir)=v0-dot_product(efielduc(:),vrcu(:,ir))
end do
! Fourier transform to Q-space
call rzfftifc(3,ngridq,-1,rfft,vclq)
deallocate(rfft)
end subroutine

