
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putdevecfv(ik,devecfv)
use modmain
use modphonon
use modramdisk
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: devecfv(nmatmax,nstfv,nspnfv)
! local variables
integer recl
character(256) fext,fname
! construct the dynamical matrix file extension
call dynfext(iqph,isph,iaph,ipph,fext)
! construct filename
fname=trim(scrpath)//'DEVECFV'//trim(fext)
!$OMP CRITICAL(u222)
! write to RAM disk if required
if (ramdisk) then
  call putrd(fname,ik,v1=vkl(:,ik),n1=nmatmax,n2=nstfv,n3=nspnfv, &
   nzv=nmatmax*nstfv*nspnfv,zva=devecfv)
end if
! write to disk if required
if (wrtdsk) then
! find the record length
  inquire(iolength=recl) vkl(:,ik),nmatmax,nstfv,nspnfv,devecfv
  open(222,file=fname,form='UNFORMATTED',access='DIRECT',recl=recl)
  write(222,rec=ik) vkl(:,ik),nmatmax,nstfv,nspnfv,devecfv
  close(222)
end if
!$OMP END CRITICAL(u222)
end subroutine

