
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writevcl1221
! !INTERFACE:
subroutine writevcl1221
! !USES:
use modmain
use modmpi
use modomp
! !DESCRIPTION:
!   Generates Coulomb matrix elements of the type $V(1,2,2,1)$ and outputs them
!   to the file {\tt VCL1221.OUT}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! allocatable arrays
real(8), allocatable :: vcl1221(:,:,:)
integer recl,ik,nthd
! determine record length for vcl1221 and open file
allocate(vcl1221(nstsv,nstsv,nkpt))
inquire(iolength=recl) vcl1221
deallocate(vcl1221)
open(260,file='VCL1221.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vcl1221) &
!$OMP NUM_THREADS(nthd)
allocate(vcl1221(nstsv,nstsv,nkpt))
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(writevcl1221_)
  write(*,'("Info(writevcl1221): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(writevcl1221_)
! calculate Coulomb matrix elements of the type V(1,2,2,1)
  call genvcl1221(ik,vcl1221)
!$OMP CRITICAL(u260)
  write(260,rec=ik) vcl1221
!$OMP END CRITICAL(u260)
end do
!$OMP END DO
deallocate(vcl1221)
!$OMP END PARALLEL
call freethd(nthd)
close(260)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine
!EOC

