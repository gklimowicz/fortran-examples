
! Copyright (C) 2020 Peter Elliott, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readafindt
use modtddft
implicit none
! local variables
integer ios
open(50,file='AFINDT.OUT',form='FORMATTED',action='READ',status='OLD', &
 iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readafindt): error opening AFINDT.OUT")')
  write(*,*)
  stop
end if
read(50,*) afindt(:,:)
close(50)
end subroutine

