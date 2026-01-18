
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine becfext(is,ia,ip,fext)
use modmain
implicit none
! arguments
integer, intent(in) :: is,ia,ip
character(*), intent(out) :: fext
write(fext,'("_S",I2.2,"_A",I3.3,"_P",I1,".OUT")') is,ia,ip
end subroutine

