
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevecfv(fext,ik,evecfv)
use modmain
use modramdisk
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
! local variables
integer recl
character(256) fname
! construct the filename
fname=trim(scrpath)//'EVECFV'//trim(fext)
!$OMP CRITICAL(u202)
! write to RAM disk if required
if (ramdisk) then
  call putrd(fname,ik,v1=vkl(:,ik),n1=nmatmax,n2=nstfv,n3=nspnfv, &
   nzv=nmatmax*nstfv*nspnfv,zva=evecfv)
end if
! write to disk if required
if (wrtdsk) then
! find the record length
  inquire(iolength=recl) vkl(:,ik),nmatmax,nstfv,nspnfv,evecfv
  open(202,file=fname,form='UNFORMATTED',access='DIRECT',recl=recl)
  write(202,rec=ik) vkl(:,ik),nmatmax,nstfv,nspnfv,evecfv
  close(202)
end if
!$OMP END CRITICAL(u202)
end subroutine

