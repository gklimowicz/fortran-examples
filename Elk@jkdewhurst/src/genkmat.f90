
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genkmat
! !INTERFACE:
subroutine genkmat(tfv,tvclcr)
! !USES:
use modmain
use modmpi
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   tfv    : .true. if the matrix elements are to be expressed in the
!            first-variational basis; second-variational otherwise (in,logical)
!   tvclvr : .true. if the non-local Coulomb potential from the core states is
!            to be included in the kinetic matrix elements (in,logical)
! !DESCRIPTION:
!   Computes the kinetic matrix elements in the first- or second-variational
!   basis and stores them in the file {\tt KMAT.OUT}. See routine {\tt putkmat}.
!
! !REVISION HISTORY:
!   Created January 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tfv,tvclcr
! local variables
integer ik,nthd
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:),bmt(:,:,:)
allocate(vmt(npcmtmax,natmtot),vir(ngtot))
if (spinpol) allocate(bmt(npcmtmax,natmtot,ndmag))
! generate the Kohn-Sham potential and magnetic field in spherical coordinates
! and multiply by the radial integration weights; also multiply the interstitial
! potential with the characteristic function
call vblocal(vmt,vir,bmt)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
if (mp_mpi) write(*,*)
! loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(genkmat_)
  write(*,'("Info(genkmat): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(genkmat_)
  call putkmat(tfv,tvclcr,ik,vmt,vir,bmt,bsir)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
deallocate(vmt,vir)
if (spinpol) deallocate(bmt)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine
!EOC

