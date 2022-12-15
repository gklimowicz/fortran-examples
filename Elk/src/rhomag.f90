
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomag
use modmain
use modmpi
use modomp
implicit none
! local variables
integer nthd
! calculate the valence density and magnetisation
call rhomagv
! add the core density and magnetisation to the total
call rhocore
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
! calculate the charges
call charge
! normalise the density
call rhonorm
!$OMP SECTION
! calculate the moments
if (spinpol) call moment
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

