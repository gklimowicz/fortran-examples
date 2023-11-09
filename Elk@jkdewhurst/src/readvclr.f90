
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readvclr
use modmain
use modulr
implicit none
! local variables
integer i1,i2,i3,ir
integer ngridq_(3),i1_,i2_,i3_
! allocatable arrays
real(8), allocatable :: vclr(:)
allocate(vclr(nqpt))
! read the real-space external Coulomb potential from file
open(50,file='VCLR.OUT',form='FORMATTED')
read(50,*) ngridq_(:)
if (any(ngridq(:).ne.ngridq_(:))) then
  write(*,*)
  write(*,'("Error(readvclr): differing ngridq")')
  write(*,'(" current  : ",3I6)') ngridq
  write(*,'(" VCLR.OUT : ",3I6)') ngridq_
  write(*,*)
  stop
end if
ir=0
do i3=1,ngridq(3)
  do i2=1,ngridq(2)
    do i1=1,ngridq(1)
      ir=ir+1
      read(50,*) i1_,i2_,i3_,vclr(ir)
      if ((i1.ne.i1_).or.(i2.ne.i2_).or.(i3.ne.i3_)) then
        write(*,*)
        write(*,'("Error(readvclr): differing i1, i2 or i3")')
        write(*,'(" current  : ",3I6)') i1,i2,i3
        write(*,'(" VCLR.OUT : ",3I6)') i1_,i2_,i3_
        write(*,*)
        stop
      end if
    end do
  end do
end do
close(50)
! Fourier transform external Coulomb potential from real-space to Q-space
call rzfftifc(3,ngridq,-1,vclr,vclq)
deallocate(vclr)
end subroutine

