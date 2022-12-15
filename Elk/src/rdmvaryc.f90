
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmvaryc
! !INTERFACE:
subroutine rdmvaryc
! !USES:
use modmain
use modrdm
use modmpi
! !DESCRIPTION:
!   Calculates new {\tt evecsv} from old by using the derivatives of the total
!   energy w.r.t. {\tt evecsv}. A single step of steepest-descent is made.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik
! allocatable arrays
complex(8), allocatable :: dedc(:,:,:),evecsv(:,:)
! compute the derivative w.r.t. evecsv
allocate(dedc(nstsv,nstsv,nkpt))
call rdmdedc(dedc)
allocate(evecsv(nstsv,nstsv))
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! get the eigenvectors from file
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! calculate new evecsv
  evecsv(:,:)=evecsv(:,:)-taurdmc*dedc(:,:,ik)
! othogonalise evecsv
  call unitary(nstsv,evecsv)
! write new evecsv to file
  call putevecsv(filext,ik,evecsv)
! end loop over k-points
end do
deallocate(dedc,evecsv)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine
!EOC

