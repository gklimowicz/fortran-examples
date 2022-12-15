
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writekpa
use modmain
use modulr
implicit none
! local variables
integer ikpa
open(50,file='KAPPA.OUT',form='FORMATTED')
write(50,'(I6," : nkpa; kappa-point, vector in unit cell reciprocal lattice &
 &coordinates below")') nkpa
do ikpa=1,nkpa
  write(50,'(I6,3G18.10)') ikpa,vql(:,ikpa)
end do
close(50)
end subroutine

