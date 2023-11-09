
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genevfsv
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,lp,nthd
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
! begin parallel loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) &
!$OMP NUM_THREADS(nthd)
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
! write the eigenvalues/vectors to file
  call putevalfv(filext,ik,evalfv)
  call putevalsv(filext,ik,evalsv(:,ik))
  call putevecfv(filext,ik,evecfv)
  call putevecsv(filext,ik,evecsv)
end do
!$OMP END DO
deallocate(evalfv,evecfv,evecsv)
!$OMP END PARALLEL
call freethd(nthd)
! broadcast eigenvalue array to every MPI process
if (np_mpi.gt.1) then
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(evalsv(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
  end do
end if
end subroutine

