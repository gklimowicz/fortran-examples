
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmldbse
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik2,nthd
call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik2=1,nkptnr
! distribute among MPI processes
  if (mod(ik2-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(hmldbse_)
  write(*,'("Info(hmldbse): ",I6," of ",I6," k-points")') ik2,nkptnr
!$OMP END CRITICAL(hmldbse_)
  call hmldbsek(ik2)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

