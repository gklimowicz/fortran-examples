
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevecsv(fext,ik,evecsv)
use modmain
use modramdisk
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer recl
character(256) fname
! construct the filename
fname=trim(scrpath)//'EVECSV'//trim(fext)
!$OMP CRITICAL(u206)
! write to RAM disk if required
if (ramdisk) then
  call putrd(fname,ik,v1=vkl(:,ik),n1=nstsv,nzv=nstsv*nstsv,zva=evecsv)
end if
! write to disk if required
if (wrtdsk) then
! find the record length
  inquire(iolength=recl) vkl(:,ik),nstsv,evecsv
  open(206,file=fname,form='UNFORMATTED',access='DIRECT',recl=recl)
  write(206,rec=ik) vkl(:,ik),nstsv,evecsv
  close(206)
end if
!$OMP END CRITICAL(u206)
end subroutine

