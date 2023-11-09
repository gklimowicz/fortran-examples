
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genhdbulr
use modmain
use modulr
use modmpi
use modomp
implicit none
! local variables
integer ik0,nthd
! loop over original k-points
call holdthd(nkpt0/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik0=1,nkpt0
! distribute among MPI processes
  if (mod(ik0-1,np_mpi).ne.lp_mpi) cycle
! write the long-range Hamiltonian diagonal blocks to file
  call puthdbulr(ik0)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

