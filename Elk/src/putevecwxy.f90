
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine putevecwxy(iq,w,x,y)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq
complex(8), intent(in) :: w(nbph,nbph),x(nbph,nbph),y(nbph)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vql(:,iq),nbph,w,x,y
!$OMP CRITICAL(u332)
open(332,file='EVECWXY.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
write(332,rec=iq) vql(:,iq),nbph,w,x,y
close(332)
!$OMP END CRITICAL(u332)
end subroutine

