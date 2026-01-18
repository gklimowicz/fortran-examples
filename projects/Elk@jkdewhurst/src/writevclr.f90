
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writevclr
use modmain
use modulr
implicit none
! local variables
integer i1,i2,i3,ir
! allocatable arrays
real(8), allocatable :: vclr(:)
allocate(vclr(nqpt))
! Fourier transform external Coulomb potential from Q-space to real-space
call rzfftifc(3,ngridq,1,vclr,vclq)
! write the real-space potential to file
open(50,file='VCLR.OUT',form='FORMATTED')
write(50,'(3I6," : grid size")') ngridq
ir=0
do i3=1,ngridq(3)
  do i2=1,ngridq(2)
    do i1=1,ngridq(1)
      ir=ir+1
      write(50,'(3I6,G18.10)') i1,i2,i3,vclr(ir)
    end do
  end do
end do
close(50)
deallocate(vclr)
end subroutine

