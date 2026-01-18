
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmat(vmt,vir,vmat)
! generates potential matrix elements for all states and k-points
use modmain
use modmpi
use modomp
implicit none
! arguments
real(8), intent(in) :: vmt(npmtmax,natmtot),vir(ngtot)
complex(8), intent(out) :: vmat(nstsv,nstsv,nkpt)
! local variables
integer ik,ispn
integer is,ias,nrc,nrci
integer n,lp,nthd
! automatic arrays
real(8) rfmt(npcmtmax)
! allocatable arrays
real(8), allocatable :: vmt1(:,:),vir1(:)
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(4), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:)
! allocate local arrays
allocate(vmt1(npcmtmax,natmtot),vir1(ngtot))
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nrc,nrci) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
! convert muffin-tin potential to a coarse radial mesg
  call rfmtftoc(nrc,nrci,vmt(:,ias),rfmt)
! convert potential to spherical coordinates
  call rbsht(nrc,nrci,rfmt,vmt1(:,ias))
! multiply by the radial integration weights
  call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt1(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! multiply interstitial potential by characteristic function
vir1(:)=vir(:)*cfunir(:)
! loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv) &
!$OMP PRIVATE(wfmt,wfgk,ispn) &
!$OMP NUM_THREADS(nthd)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngkmax,nspinor,nstsv))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors from file
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! calculate the wavefunctions for all states of the input k-point
  call genwfsv_sp(.false.,.true.,nstsv,0,ngridg,igfft,ngk(:,ik),igkig(:,:,ik), &
   apwalm,evecfv,evecsv,wfmt,ngkmax,wfgk)
  call genvmatk(vmt1,vir1,ngk(:,ik),igkig(:,:,ik),wfmt,ngkmax,wfgk,vmat(:,:,ik))
end do
!$OMP END DO
deallocate(apwalm,evecfv,evecsv,wfmt,wfgk)
!$OMP END PARALLEL
call freethd(nthd)
! broadcast matrix elements to every process
if (np_mpi.gt.1) then
  n=nstsv*nstsv
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(vmat(:,:,ik),n,mpi_double_complex,lp,mpicom,ierror)
  end do
end if
deallocate(vmt1,vir1)
end subroutine

