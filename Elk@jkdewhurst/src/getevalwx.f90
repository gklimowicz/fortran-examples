
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine getevalwx(iq,evalwxp)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq
real(8), intent(out) :: evalwxp(nbph)
! local variables
integer recl,nbph_
real(8) vql_(3),t1
! find the record length
inquire(iolength=recl) vql_,nbph_,evalwxp
!$OMP CRITICAL(u330)
open(330,file='EVALWX.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(330,rec=iq) vql_,nbph_,evalwxp
close(330)
!$OMP END CRITICAL(u330)
t1=abs(vql(1,iq)-vql_(1))+abs(vql(2,iq)-vql_(2))+abs(vql(3,iq)-vql_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevalwx): differing vectors for q-point ",I8)') iq
  write(*,'(" current    : ",3G18.10)') vql(:,iq)
  write(*,'(" EVALWX.OUT : ",3G18.10)') vql_
  write(*,*)
  stop
end if
if (nbph.ne.nbph_) then
  write(*,*)
  write(*,'("Error(getevalwx): differing nbph for q-point ",I8)') iq
  write(*,'(" current    : ",I8)') nbph
  write(*,'(" EVALWX.OUT : ",I8)') nbph_
  write(*,*)
  stop
end if
end subroutine

