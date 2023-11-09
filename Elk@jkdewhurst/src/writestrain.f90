
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writestrain
use modmain
implicit none
! local variables
integer i,j,k
! initialise universal variables
call init0
! generate the strain tensors
call genstrain
! write the strain tensors to file
open(50,file='STRAIN.OUT',form='FORMATTED',action='WRITE')
do k=1,nstrain
  write(50,*)
  write(50,'("Strain tensor : ",I1)') k
  do j=1,3
    write(50,'(3G18.10)') (strain(i,j,k),i=1,3)
  end do
end do
close(50)
write(*,*)
write(*,'("Info(writestrain)")')
write(*,'(" Strain tensors written to STRAIN.OUT")')
write(*,'(" (the first strain tensor is isotropic expansion)")')
end subroutine

