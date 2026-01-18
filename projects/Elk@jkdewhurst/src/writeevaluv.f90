
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevaluv
use modmain
use modbog
implicit none
! local variables
integer ik,ist
real(8) e
! write out the fermionic eigenvalues
open(50,file='EIGVALUV.OUT',form='FORMATTED',action='WRITE')
write(50,'(I6," : nkpt")') nkpt
write(50,'(I6," : nstsv")') nstsv
do ik=1,nkpt
  write(50,*)
  write(50,'(I6,3G18.10," : k-point, vkl")') ik,vkl(:,ik)
  write(50,'(" (state, eigenvalue, V-norm below)")')
  do ist=1,nstsv
    e=evaluv(ist,ik)
    if (vnorm(ist,ik).gt.0.5d0) e=-e
    write(50,'(I6,2G18.10)') ist,e,vnorm(ist,ik)
  end do
  write(50,*)
end do
close(50)
end subroutine

