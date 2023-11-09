
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevaluv(ik,evaluvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(out) :: evaluvp(nstsv)
! local variables
integer recl,nstsv_
real(8) vkl_(3),t1
! find the record length
inquire(iolength=recl) vkl_,nstsv_,evaluvp
!$OMP CRITICAL(u320)
open(320,file='EVALUV.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(320,rec=ik) vkl_,nstsv_,evaluvp
close(320)
!$OMP END CRITICAL(u320)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevaluv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVALUV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getevaluv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EVALUV.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
end subroutine

