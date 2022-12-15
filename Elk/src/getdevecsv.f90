
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getdevecsv(ik,iq,is,ia,ip,devecsv)
use modmain
use modphonon
use modramdisk
implicit none
! arguments
integer, intent(in) :: ik,iq,is,ia,ip
complex(8), intent(out) :: devecsv(nstsv,nstsv)
! local variables
logical tgs
integer recl,nstsv_
real(8) vkl_(3),t1
character(256) fext,fname
! construct the dynamical matrix file extension
call dynfext(iq,is,ia,ip,fext)
! construct filename
fname=trim(scrpath)//'DEVECSV'//trim(fext)
! find the record length
inquire(iolength=recl) vkl_,nstsv_,devecsv
!$OMP CRITICAL(u226)
! read from RAM disk if required
if (ramdisk) then
  call getrd(fname,ik,tgs,v1=vkl_,n1=nstsv_,nzv=nstsv*nstsv,zva=devecsv)
  if (tgs) goto 10
end if
open(226,file=fname,form='UNFORMATTED',access='DIRECT',recl=recl)
read(226,rec=ik) vkl_,nstsv_,devecsv
close(226)
10 continue
!$OMP END CRITICAL(u226)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getdevecsv): differing vectors for k-point ",I8)') ik
  write(*,'(" current : ",3G18.10)') vkl(:,ik)
  write(*,'(" ",A," : ",3G18.10)') trim(fname),vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getdevecsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current : ",I8)') nstsv
  write(*,'(" ",A," : ",I8)') trim(fname),nstsv_
  write(*,*)
  stop
end if
end subroutine

