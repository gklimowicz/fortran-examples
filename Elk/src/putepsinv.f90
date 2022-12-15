
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putepsinv(iq,epsi)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
complex(8), intent(in) :: epsi(ngrf,ngrf,nwrf)
! local variables
integer recl
! determine the record length for EPSINV.OUT
inquire(iolength=recl) vql(:,iq),ngrf,nwrf,epsi
!$OMP CRITICAL(u245)
open(245,file='EPSINV.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
write(245,rec=iq) vql(:,iq),ngrf,nwrf,epsi
close(245)
!$OMP END CRITICAL(u245)
end subroutine

