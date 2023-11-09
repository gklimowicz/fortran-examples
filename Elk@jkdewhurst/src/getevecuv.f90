
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevecuv(ikp,vpl,u,v)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: u(nstsv,nstsv),v(nstsv,nstsv)
! local variables
integer isym,ik
integer recl,nstsv_
real(8) vkl_(3),t1
if (ikp.gt.0) then
  ik=ikp
else
! find the equivalent k-point number and symmetry which rotates vkl to vpl
  call findkpt(vpl,isym,ik)
end if
! find the record length
inquire(iolength=recl) vkl_,nstsv_,u,v
!$OMP CRITICAL(u322)
open(322,file='EVECUV.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(322,rec=ik) vkl_,nstsv_,u,v
close(322)
!$OMP END CRITICAL(u322)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevecuv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVECUV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getevecuv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EVECUV.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
end subroutine

