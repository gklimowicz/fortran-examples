
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putpmat(ik)
use modmain
use modmpi
use modramdisk
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ispn,recl
! automatic arrays
complex(8) pmat(nstsv,nstsv,3)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:)
! get the eigenvectors from file
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
do ispn=1,nspnfv
  call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
! calculate the wavefunctions for all states
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngkmax,nspinor,nstsv))
call genwfsv(.true.,.true.,nstsv,0,ngridg,igfft,ngk(:,ik),igkig(:,:,ik),apwalm,&
 evecfv,evecsv,wfmt,ngkmax,wfgk)
deallocate(evecfv,evecsv,apwalm)
! calculate the momentum matrix elements
call genpmatk(ngk(:,ik),igkig(:,:,ik),vgkc(:,:,:,ik),wfmt,wfgk,pmat)
deallocate(wfmt,wfgk)
! write the matrix elements in the second-variational basis
!$OMP CRITICAL(u230)
! write to RAM disk if required
if (ramdisk) then
  call putrd('PMAT.OUT',ik,v1=vkl(:,ik),n1=nstsv,nzv=nstsv*nstsv*3,zva=pmat)
end if
! write to disk if required
if (wrtdsk) then
! determine the record length
  inquire(iolength=recl) vkl(:,1),nstsv,pmat
  open(230,file='PMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
  write(230,rec=ik) vkl(:,ik),nstsv,pmat
  close(230)
end if
!$OMP END CRITICAL(u230)
end subroutine

