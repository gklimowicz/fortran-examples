
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epsinv
use modmain
use modmpi
use modomp
implicit none
! local variables
integer iq,ik,ig,iw
integer n,nthd
! allocatable arrays
integer(8), allocatable :: lock(:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:),epsi(:,:,:)
! allocate local arrays
allocate(vgqc(3,ngrf),gqc(ngrf),gclgq(ngrf))
allocate(jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(epsi(ngrf,ngrf,nwrf))
! initialise the OpenMP locks
allocate(lock(nwrf))
do iw=1,nwrf
  call omp_init_lock(lock(iw))
end do
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(epsinv): ",I6," of ",I6," q-points")') iq,nqpt
! generate the G+q-vectors and related functions
  call gengqf(ngrf,vqc(:,iq),vgqc,gqc,jlgqr,ylmgq,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngrf,gqc,gclgq)
! use the symmetric form
  gclgq(:)=sqrt(gclgq(:))
! zero the response function (stored in epsi)
  epsi(:,:,:)=0.d0
  call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! compute v^1/2 chi0 v^1/2
    call genvchi0(.false.,ik,lock,0.d0,vql(:,iq),gclgq,jlgqr,ylmgq,sfacgq, &
     ngrf,epsi)
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! add epsi from each process and redistribute
  if (np_mpi.gt.1) then
    n=ngrf*ngrf*nwrf
    call mpi_allreduce(mpi_in_place,epsi,n,mpi_double_complex,mpi_sum,mpicom, &
     ierror)
  end if
! negate and add delta(G,G')
  epsi(:,:,:)=-epsi(:,:,:)
  do ig=1,ngrf
    epsi(ig,ig,:)=epsi(ig,ig,:)+1.d0
  end do
!-------------------------------------!
!     invert epsilon over G-space     !
!-------------------------------------!
  call holdthd(nwrf,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do iw=1,nwrf
    call zminv(ngrf,epsi(:,:,iw))
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! write inverse RPA epsilon to EPSINV.OUT
  if (mp_mpi) call putepsinv(iq,epsi)
! end loop over q-points
end do
! destroy the OpenMP locks
do iw=1,nwrf
  call omp_destroy_lock(lock(iw))
end do
deallocate(lock)
deallocate(vgqc,gqc,gclgq,jlgqr)
deallocate(ylmgq,sfacgq,epsi)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

