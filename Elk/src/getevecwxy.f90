
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine getevecwxy(iq,w,x,y)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq
complex(8), intent(out) :: w(nbph,nbph),x(nbph,nbph),y(nbph)
! local variables
integer recl,nbph_
real(8) vql_(3),t1
! find the record length
inquire(iolength=recl) vql_,nbph_,w,x,y
!$OMP CRITICAL(u332)
open(332,file='EVECWXY.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(332,rec=iq) vql_,nbph_,w,x,y
close(332)
!$OMP END CRITICAL(u332)
t1=abs(vql(1,iq)-vql_(1))+abs(vql(2,iq)-vql_(2))+abs(vql(3,iq)-vql_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevecwxy): differing vectors for q-point ",I8)') iq
  write(*,'(" current     : ",3G18.10)') vql(:,iq)
  write(*,'(" EVECWXY.OUT : ",3G18.10)') vql_
  write(*,*)
  stop
end if
if (nbph.ne.nbph_) then
  write(*,*)
  write(*,'("Error(getevecwxy): differing nbph for q-point ",I8)') iq
  write(*,'(" current     : ",I8)') nbph
  write(*,'(" EVECWXY.OUT : ",I8)') nbph_
  write(*,*)
  stop
end if
end subroutine

