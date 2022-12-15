
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjtot
use modmain
implicit none
! local variables
integer i
! external functions
real(8), external :: rfint
if (.not.tjr) return
! compute the total current in the unit cell
do i=1,3
  jtot(i)=rfint(jrmt(:,:,i),jrir(:,i))
end do
! total current magnitude
jtotm=sqrt(jtot(1)**2+jtot(2)**2+jtot(3)**2)
end subroutine

