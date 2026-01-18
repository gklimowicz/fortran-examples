
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gethdbulr(ik0,hdb)
use modmain
use modulr
use modramdisk
implicit none
! arguments
integer, intent(in) :: ik0
complex(8), intent(out) :: hdb(nstsv,nstsv,2:nkpa)
! local variables
logical tgs
integer nstsv_,nkpa_
integer ik,recl
real(8) vkl_(3),t1
! find the record length
inquire(iolength=recl) vkl_,nstsv_,nkpa_,hdb
!$OMP CRITICAL(u300)
! read from RAM disk if required
if (ramdisk) then
  call getrd('HDBULR.OUT',ik0,tgs,v1=vkl_,n1=nstsv_,n2=nkpa_, &
   nzv=nstsv*nstsv*(nkpa-1),zva=hdb)
  if (tgs) goto 10
end if
open(300,file='HDBULR.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(300,rec=ik0) vkl_,nstsv_,nkpa_,hdb
close(300)
10 continue
!$OMP END CRITICAL(u300)
! central k-point
ik=(ik0-1)*nkpa+1
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(gethdbulr): differing vectors for k-point ",I8)') ik0
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" HDBULR.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(gethdbulr): differing nstsv for k-point ",I8)') ik0
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" HDBULR.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
if (nkpa.ne.nkpa_) then
  write(*,*)
  write(*,'("Error(gethdbulr): differing nkpa for k-point ",I8)') ik0
  write(*,'(" current    : ",I8)') nkpa
  write(*,'(" HDBULR.OUT : ",I8)') nkpa_
  write(*,*)
  stop
end if
end subroutine

