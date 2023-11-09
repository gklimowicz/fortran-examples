
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: polar
! !INTERFACE:
subroutine polar(pvl)
! !USES:
use modmain
use modmpi
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   pvl : polarisation vector modulo $2\pi$ (out,real(8))
! !DESCRIPTION:
!   Calculates the polarisation vector modulo $2\pi$ in lattice coordinates
!   using the formula of R. D. King-Smith and David Vanderbilt [Phys. Rev. B
!   {\bf 47}, 1651(R) (1993)], namely
!   $$ P_l=\sum_{\bf k}{\rm Im}\ln\det\left(\langle
!    u_{i{\bf k}+\Delta{\bf k}_l}|u_{j{\bf k}}\rangle\right),$$
!   where $\Delta{\bf k}_l=(1/n_l){\bf B}_l$ and ${\bf B}_l$ is a reciprocal
!   lattice vector. The number of points $n_l$ is equal to that of the original
!   $k$-point grid in direction of ${\bf B}_l$, multiplied by {\tt nskpolar}.
!   See also the routines {\tt polark} and {\tt bornechg}.
!
! !REVISION HISTORY:
!   Created May 2020 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(out) :: pvl(3)
! local variables
integer ik,l,nthd
real(8) vc(3),vgqc(3),gqc,pl
! allocatable arrays
real(8), allocatable :: jlgqr(:,:)
complex(8), allocatable :: ylmgq(:),sfacgq(:),expqmt(:,:)
! allocate local arrays
allocate(jlgqr(njcmax,nspecies))
allocate(ylmgq(lmmaxo),sfacgq(natmtot))
allocate(expqmt(npcmtmax,natmtot))
maxscl=1
! loop over reciprocal lattice vectors
do l=1,3
! create fine k-point grid in direction l
  ngridk(:)=ngridk0(:)
  ngridk(l)=nkspolar*ngridk(l)
! run one loop of the ground-state calculation
  call gndstate
! difference between adjacent k-vectors in this reciprocal lattice direction
  vc(:)=bvec(:,l)/dble(ngridk(l))
! calculate the phase factor function exp(iq.r)
  call gengqf(1,vc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
  call genexpmt(1,jlgqr,ylmgq,1,sfacgq,expqmt)
  pl=0.d0
! parallel loop over non-reduced k-points
  call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP REDUCTION(+:pl) &
!$OMP NUM_THREADS(nthd)
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    call polark(ik,l,expqmt,pl)
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! add polarisation from each process and redistribute
  if (np_mpi.gt.1) then
    call mpi_allreduce(mpi_in_place,pl,1,mpi_double_precision,mpi_sum,mpicom, &
     ierror)
  end if
  pvl(l)=pl
end do
! restore original input parameters
ngridk(:)=ngridk0(:)
maxscl=maxscl0
deallocate(jlgqr,ylmgq,sfacgq,expqmt)
end subroutine
!EOC

