
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine checkwrite(twrite)
use modmain
implicit none
! arguments
logical, intent(out) :: twrite
! check for WRITE file
inquire(file='WRITE',exist=twrite)
if (twrite) then
  write(*,*)
  write(*,'("Info(checkwrite): WRITE file exists")')
  open(50,file='WRITE')
  close(50,status='DELETE')
end if
end subroutine

