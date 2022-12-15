
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevalfv(fext,ikp,vpl,evalfv)
use modmain
use modramdisk
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ikp
real(8), intent(in) :: vpl(3)
real(8), intent(out) :: evalfv(nstfv,nspnfv)
! local variables
logical tgs
integer isym,ik
integer recl,nstfv_,nspnfv_
real(8) vkl_(3),t1
character(256) fname
if (ikp.gt.0) then
  ik=ikp
else
! find the k-point number
  call findkpt(vpl,isym,ik)
end if
! construct the filename
fname='EVALFV'//trim(fext)
! find the record length
inquire(iolength=recl) vkl_,nstfv_,nspnfv_,evalfv
!$OMP CRITICAL(u200)
! read from RAM disk if required
if (ramdisk) then
  call getrd(fname,ik,tgs,v1=vkl_,n1=nstfv_,n2=nspnfv_,nrv=nstfv*nspnfv, &
   rva=evalfv)
  if (tgs) goto 10
end if
open(200,file=fname,form='UNFORMATTED',access='DIRECT',recl=recl)
read(200,rec=ik) vkl_,nstfv_,nspnfv_,evalfv
close(200)
10 continue
!$OMP END CRITICAL(u200)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevalfv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVALFV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstfv.ne.nstfv_) then
  write(*,*)
  write(*,'("Error(getevalfv): differing nstfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstfv
  write(*,'(" EVALFV.OUT : ",I8)') nstfv_
  write(*,*)
  stop
end if
if (nspnfv.ne.nspnfv_) then
  write(*,*)
  write(*,'("Error(getevalfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nspnfv
  write(*,'(" EVALFV.OUT : ",I8)') nspnfv_
  write(*,*)
  stop
end if
end subroutine

