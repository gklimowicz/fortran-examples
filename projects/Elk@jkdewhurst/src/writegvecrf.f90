
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writegvecrf
use modmain
implicit none
! local variables
integer ig
open(50,file='GVECRF.OUT',form='FORMATTED',action='WRITE')
write(50,'(G18.10," : gmaxrf")') gmaxrf
write(50,'(I8," : ngrf; G-vector index, ivg below")') ngrf
do ig=1,ngrf
  write(50,'(I8,3I6)') ig,ivg(1:3,ig)
end do
close(50)
end subroutine

