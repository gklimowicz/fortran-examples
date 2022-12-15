
! Copyright (C) 2021 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine writeevalwx
use modmain
use modphonon
use modbog
implicit none
! local variables
integer iq,i
! write out the bosonic eigenvalues
open(50,file='EIGVALWX.OUT',form='FORMATTED',action='WRITE')
write(50,'(I6," : nqpt")') nqpt
write(50,'(I6," : nbph")') nbph
do iq=1,nqpt
  write(50,*)
  write(50,'(I6,3G18.10," : q-point, vql")') iq,vql(:,iq)
  write(50,'(" (phonon branch, eigenvalue, X-norm below)")')
  do i=1,nbph
    write(50,'(I6,2G18.10)') i,evalwx(i,iq),xnorm(i,iq)
  end do
end do
close(50)
end subroutine

