
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine putevalwx(iq,evalwxp)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq
real(8), intent(in) :: evalwxp(nbph)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vql(:,iq),nbph,evalwxp
!$OMP CRITICAL(u330)
open(330,file='EVALWX.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
write(330,rec=iq) vql(:,iq),nbph,evalwxp
close(330)
!$OMP END CRITICAL(u330)
end subroutine

