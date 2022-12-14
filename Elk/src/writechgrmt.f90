
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writechgrmt
use modmain
use modulr
implicit none
! local variables
integer is,ia,ias,ir
open(50,file='CHGMTRU.OUT',form='FORMATTED')
do ir=1,nqpt
  write(50,*)
  write(50,'("R-point (Cartesian coordinates) :")')
  write(50,'(3G18.10)') vrcu(:,ir)
  do is=1,nspecies
    write(50,'("  species : ",I4," (",A,")")') is,trim(spsymb(is))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,'("   atom ",I4,T30,": ",G18.10)') ia,chgmtru(ias,ir)
    end do
  end do
end do
close(50)
end subroutine

