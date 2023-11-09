
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevalu
use modmain
use modulr
implicit none
! local variables
integer ik0,ik,ist
! write out the valence eigenvalues
open(50,file='EIGVALU.OUT',form='FORMATTED')
write(50,'(I6," : nkpt0")') nkpt0
write(50,'(I6," : nstulr")') nstulr
do ik0=1,nkpt0
! central k-point
  ik=(ik0-1)*nkpa+1
  write(50,*)
  write(50,'(I6,3G18.10," : k-point, vkl")') ik0,vkl(:,ik)
  write(50,'(" (state, eigenvalue and occupancy below)")')
  do ist=1,nstulr
    write(50,'(I6,2G18.10)') ist,evalu(ist,ik0),occulr(ist,ik0)
  end do
end do
close(50)
end subroutine

