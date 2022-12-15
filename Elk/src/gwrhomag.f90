
! Copyright (C) 2018 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwrhomag
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
integer ik,lp,nthd
! compute the GW density matrices and write the natural orbitals and occupation
! numbers to EVECSV.OUT and OCCSV.OUT, respectively
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call gwdmatk(ik)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! broadcast occupation number array to every MPI process
if (np_mpi.gt.1) then
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(occsv(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
  end do
end if
! write the occupation numbers to file
if (mp_mpi) then
  do ik=1,nkpt
    call putoccsv(filext,ik,occsv(:,ik))
  end do
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! determine the density and magnetisation in the usual way
call rhomag
end subroutine

