
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putdevecsv(ik,devecsv)
use modmain
use modphonon
use modramdisk
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: devecsv(nstsv,nstsv)
! local variables
integer recl
character(256) fext,fname
! construct the dynamical matrix file extension
call dynfext(iqph,isph,iaph,ipph,fext)
! construct filename
fname=trim(scrpath)//'DEVECSV'//trim(fext)
!$OMP CRITICAL(u226)
! write to RAM disk if required
if (ramdisk) then
  call putrd(fname,ik,v1=vkl(:,ik),n1=nstsv,nzv=nstsv*nstsv,zva=devecsv)
end if
! write to disk if required
if (wrtdsk) then
! find the record length
  inquire(iolength=recl) vkl(:,ik),nstsv,devecsv
  open(226,file=fname,form='UNFORMATTED',access='DIRECT',recl=recl)
  write(226,rec=ik) vkl(:,ik),nstsv,devecsv
  close(226)
end if
!$OMP END CRITICAL(u226)
end subroutine

