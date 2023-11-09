
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putkmat(tfv,tvclcr,ik,vmt,vir,bmt,bir)
use modmain
use modmpi
use modramdisk
implicit none
! arguments
logical, intent(in) :: tfv,tvclcr
integer, intent(in) :: ik
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
! local variables
integer ist,ispn,recl
! automatic arrays
complex(8) kmat(nstsv,nstsv),a(nstsv,nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
complex(4), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:)
! get the eigenvalues/vectors from file for input reduced k-point
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
do ispn=1,nspnfv
  call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
! calculate the wavefunctions for all states of the input k-point
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngkmax,nspinor,nstsv))
call genwfsv_sp(.false.,.true.,nstsv,0,ngridg,igfft,ngk(:,ik),igkig(:,:,ik), &
 apwalm,evecfv,evecsv,wfmt,ngkmax,wfgk)
deallocate(apwalm,evecfv)
! compute Kohn-Sham potential matrix elements
if (spinpol) then
  call genvbmatk(vmt,vir,bmt,bir,ngk(:,ik),igkig(:,:,ik),wfmt,ngkmax,wfgk,kmat)
else
  call genvmatk(vmt,vir,ngk(:,ik),igkig(:,:,ik),wfmt,ngkmax,wfgk,kmat)
end if
deallocate(wfgk)
! negate the potential matrix elements because we have to subtract them
kmat(:,:)=-kmat(:,:)
! add second-variational eigenvalues along the diagonal
do ist=1,nstsv
  kmat(ist,ist)=kmat(ist,ist)+evalsv(ist,ik)
end do
! add scissor correction if required
if (scissor.ne.0.d0) then
  do ist=1,nstsv
    if (evalsv(ist,ik).gt.efermi) then
      kmat(ist,ist)=kmat(ist,ist)+scissor
    end if
  end do
end if
! add the Coulomb core matrix elements if required
if (tvclcr) call vclcore(wfmt,kmat)
! rotate kinetic matrix elements to first-variational basis if required
if (tfv) then
  call zgemm('N','C',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsv,nstsv,zzero,a, &
   nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,a,nstsv,zzero,kmat, &
   nstsv)
end if
deallocate(evecsv,wfmt)
!$OMP CRITICAL(u220)
! write to RAM disk if required
if (ramdisk) then
  call putrd('KMAT.OUT',ik,v1=vkl(:,ik),n1=nstsv,nzv=nstsv*nstsv,zva=kmat)
end if
! write to disk if required
if (wrtdsk) then
! determine the record length
  inquire(iolength=recl) vkl(:,1),nstsv,kmat
  open(220,file='KMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
  write(220,rec=ik) vkl(:,ik),nstsv,kmat
  close(220)
end if
!$OMP END CRITICAL(u220)
end subroutine

