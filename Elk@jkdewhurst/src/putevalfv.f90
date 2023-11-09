
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevalfv(fext,ik,evalfv)
use modmain
use modramdisk
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
real(8), intent(in) :: evalfv(nstfv,nspnfv)
! local variables
integer recl
character(256) fname
! construct the filename
fname='EVALFV'//trim(fext)
!$OMP CRITICAL(u200)
! write to RAM disk if required
if (ramdisk) then
  call putrd(fname,ik,v1=vkl(:,ik),n1=nstfv,n2=nspnfv,nrv=nstfv*nspnfv, &
   rva=evalfv)
end if
! write to disk if required
if (wrtdsk) then
! find the record length
  inquire(iolength=recl) vkl(:,ik),nstfv,nspnfv,evalfv
  open(200,file=fname,form='UNFORMATTED',access='DIRECT',recl=recl)
  write(200,rec=ik) vkl(:,ik),nstfv,nspnfv,evalfv
  close(200)
end if
!$OMP END CRITICAL(u200)
end subroutine

