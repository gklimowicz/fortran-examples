
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentau
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,ispn,is,ias
integer np,npc,n,nthd
! allocatable arrays
real(8), allocatable :: rfmt(:,:),rfir(:)
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
! set the kinetic energy density to zero
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    taumt(1:npcmt(is),ias,ispn)=0.d0
  end do
end do
tauir(1:ngtc,:)=0.d0
! tau cannot be computed if wavefunctions do not exist
if (iscl.le.0) return
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call gentauk(ik)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
allocate(rfmt(npcmtmax,natmtot))
! convert taumt to spherical harmonics
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    rfmt(1:npcmt(is),1)=taumt(1:npcmt(is),ias,ispn)
    call rfsht(nrcmt(is),nrcmti(is),rfmt,taumt(:,ias,ispn))
  end do
end do
! symmetrise tau
if (spinpol) then
! spin-polarised case: convert to scalar-vector form
  allocate(rvfmt(npcmtmax,natmtot,ndmag))
  allocate(rvfir(ngtc,ndmag),rfir(ngtc))
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    rfmt(1:npc,ias)=taumt(1:npc,ias,1)+taumt(1:npc,ias,2)
    rvfmt(1:npc,ias,1:ndmag-1)=0.d0
    rvfmt(1:npc,ias,ndmag)=taumt(1:npc,ias,1)-taumt(1:npc,ias,2)
  end do
  rfir(1:ngtc)=tauir(1:ngtc,1)+tauir(1:ngtc,2)
  rvfir(1:ngtc,1:ndmag-1)=0.d0
  rvfir(1:ngtc,ndmag)=tauir(1:ngtc,1)-tauir(1:ngtc,2)
  call symrf(nrcmt,nrcmti,npcmt,ngdgc,ngtc,ngvc,igfc,npcmtmax,rfmt,rfir)
  call symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,ngdgc,ngtc,ngvc,igfc,npcmtmax,&
   rvfmt,ngtc,rvfir)
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    taumt(1:npc,ias,1)=0.5d0*(rfmt(1:npc,ias)+rvfmt(1:npc,ias,ndmag))
    taumt(1:npc,ias,2)=0.5d0*(rfmt(1:npc,ias)-rvfmt(1:npc,ias,ndmag))
  end do
  tauir(1:ngtc,1)=0.5d0*(rfir(1:ngtc)+rvfir(1:ngtc,ndmag))
  tauir(1:ngtc,2)=0.5d0*(rfir(1:ngtc)-rvfir(1:ngtc,ndmag))
  deallocate(rvfmt,rvfir,rfir)
else
! spin-unpolarised case
  call symrf(nrcmt,nrcmti,npcmt,ngdgc,ngtc,ngvc,igfc,npmtmax,taumt,tauir)
end if
! convert muffin-tin tau from coarse to fine radial mesh
do ispn=1,nspinor
  call rfmtctof(taumt(:,:,ispn))
end do
! convert interstitial tau from coarse to fine grid
do ispn=1,nspinor
  call rfirctof(tauir(:,ispn),tauir(:,ispn))
end do
! add tau from each process and redistribute
if (np_mpi.gt.1) then
  n=npmtmax*natmtot*nspinor
  call mpi_allreduce(mpi_in_place,taumt,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
  n=ngtot*nspinor
  call mpi_allreduce(mpi_in_place,tauir,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
end if
! generate the core kinetic energy density
call gentaucr
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
! add the core contribution
    taumt(1:np,ias,ispn)=taumt(1:np,ias,ispn)+taucr(1:np,ias,ispn)
! zero tau on the inner part of the muffin-tin
    taumt(1:npmti(is),ias,ispn)=0.d0
  end do
end do
deallocate(rfmt)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

