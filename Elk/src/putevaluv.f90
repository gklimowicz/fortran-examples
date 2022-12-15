
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevaluv(ik,evaluvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: evaluvp(nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,evaluvp
!$OMP CRITICAL(u320)
open(320,file='EVALUV.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
write(320,rec=ik) vkl(:,ik),nstsv,evaluvp
close(320)
!$OMP END CRITICAL(u320)
end subroutine

