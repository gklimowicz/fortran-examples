
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevalsv(fext,ik,evalsvp)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
real(8), intent(in) :: evalsvp(nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,evalsvp
!$OMP CRITICAL(u204)
open(204,file='EVALSV'//trim(fext),form='UNFORMATTED',access='DIRECT',recl=recl)
write(204,rec=ik) vkl(:,ik),nstsv,evalsvp
close(204)
!$OMP END CRITICAL(u204)
end subroutine

