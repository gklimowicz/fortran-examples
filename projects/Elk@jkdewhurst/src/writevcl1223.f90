
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writevcl1223
! !INTERFACE:
subroutine writevcl1223
! !USES:
use modmain
use modmpi
use modomp
! !DESCRIPTION:
!   Generates Coulomb matrix elements of the type $V(1,2,2,3)$ and outputs them
!   to the file {\tt VCL1223.OUT}. Also writes the real diagonal of this matrix,
!   $V(1,2,2,1)$, to {\tt VCL1221.OUT}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,recl,nthd
! allocatable arrays
real(8), allocatable :: vcl1221(:,:,:)
complex(8), allocatable :: vcl1223(:,:,:,:)
! determine record length for vcl1221 and open file
allocate(vcl1221(nstsv,nstsv,nkpt))
inquire(iolength=recl) vcl1221
deallocate(vcl1221)
open(260,file='VCL1221.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
! determine record length for vcl1223 and open file
allocate(vcl1223(nstsv,nstsv,nstsv,nkpt))
inquire(iolength=recl) vcl1223
deallocate(vcl1223)
open(262,file='VCL1223.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vcl1221,vcl1223,ist) &
!$OMP NUM_THREADS(nthd)
allocate(vcl1221(nstsv,nstsv,nkpt),vcl1223(nstsv,nstsv,nstsv,nkpt))
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(writevcl1223_)
  write(*,'("Info(writevcl1223): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(writevcl1223_)
! calculate Coulomb matrix elements of the type V(1,2,2,3)
  call genvcl1223(ik,vcl1223)
! make a copy of the diagonal elements V(1,2,2,1)
  do ist=1,nstsv
    vcl1221(ist,:,:)=dble(vcl1223(ist,ist,:,:))
  end do
!$OMP CRITICAL(u260)
  write(260,rec=ik) vcl1221
!$OMP END CRITICAL(u260)
!$OMP CRITICAL(u262)
  write(262,rec=ik) vcl1223
!$OMP END CRITICAL(u262)
end do
!$OMP END DO
deallocate(vcl1221,vcl1223)
!$OMP END PARALLEL
call freethd(nthd)
close(260)
close(262)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine
!EOC

