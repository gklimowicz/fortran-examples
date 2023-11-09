
! Copyright (C) 2018 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjpr
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,is,ias
integer nrc,nrci,npc
integer n,i,nthd
! automatic arrays
integer(8) lock(natmtot)
real(8) rfmt(npcmtmax)
! set the current density to zero
do i=1,3
  do ias=1,natmtot
    is=idxis(ias)
    jrmt(1:npcmt(is),ias,i)=0.d0
  end do
end do
jrir(1:ngtc,1:3)=0.d0
! current density cannot be computed if wavefunctions do not exist
if (iscl.le.0) return
! initialise the OpenMP locks
do ias=1,natmtot
  call omp_init_lock(lock(ias))
end do
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call genjprk(ik,lock)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! destroy the OpenMP locks
do ias=1,natmtot
  call omp_destroy_lock(lock(ias))
end do
! convert muffin-tin current density to spherical harmonics
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nrc,nrci,npc,i) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do i=1,3
    rfmt(1:npc)=jrmt(1:npc,ias,i)
    call rfsht(nrc,nrci,rfmt,jrmt(:,ias,i))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! symmetrise the current density
call symrvf(.false.,.true.,nrcmt,nrcmti,npcmt,ngdgc,ngtc,ngvc,igfc,npmtmax, &
 jrmt,ngtot,jrir)
! convert muffin-tin and interstitial current density from coarse to fine grids
call holdthd(3,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do i=1,3
  call rfmtctof(jrmt(:,:,i))
end do
!$OMP END DO NOWAIT
!$OMP DO
do i=1,3
  call rfirctof(jrir(:,i),jrir(:,i))
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! add current densities from each process and redistribute
if (np_mpi.gt.1) then
  n=npmtmax*natmtot*3
  call mpi_allreduce(mpi_in_place,jrmt,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
  n=ngtot*3
  call mpi_allreduce(mpi_in_place,jrir,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

