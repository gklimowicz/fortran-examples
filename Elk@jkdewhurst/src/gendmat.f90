
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmat(tspndg,tlmdg,lmin,lmax,ld,dmat)
use modmain
use modmpi
use modomp
implicit none
! arguments
logical, intent(in) :: tspndg,tlmdg
integer, intent(in) :: lmin,lmax,ld
complex(8), intent(out) :: dmat(ld,nspinor,ld,nspinor,natmtot)
! local variables
integer ik,ispn
integer nst,ist,jst
integer ias,n,nthd
real(8) wo
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: dmatk(:,:,:,:,:)
! zero the density matrix
dmat(:,:,:,:,:)=0.d0
! begin parallel loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(idx,apwalm,evecfv,evecsv,dmatk) &
!$OMP PRIVATE(ispn,nst,ist,jst,ias,wo) &
!$OMP NUM_THREADS(nthd)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
allocate(dmatk(ld,nspinor,ld,nspinor,nstsv))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors from file
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! index to occupied states
  nst=0
  do ist=1,nstsv
    if (abs(occsv(ist,ik)).gt.epsocc) then
      nst=nst+1
      idx(nst)=ist
    end if
  end do
! begin loop over atoms and species
  do ias=1,natmtot
    call gendmatk(tspndg,tlmdg,lmin,lmax,ias,nst,idx,ngk(:,ik),apwalm,evecfv, &
     evecsv,ld,dmatk)
    do ist=1,nst
      jst=idx(ist)
      wo=occsv(jst,ik)*wkpt(ik)
!$OMP CRITICAL(gendmat_)
      dmat(:,:,:,:,ias)=dmat(:,:,:,:,ias)+wo*dmatk(:,:,:,:,ist)
!$OMP END CRITICAL(gendmat_)
    end do
  end do
end do
!$OMP END DO
deallocate(apwalm,evecfv,evecsv,dmatk)
!$OMP END PARALLEL
call freethd(nthd)
! add density matrices from each process and redistribute
if (np_mpi.gt.1) then
  n=((ld*nspinor)**2)*natmtot
  call mpi_allreduce(mpi_in_place,dmat,n,mpi_double_complex,mpi_sum,mpicom, &
   ierror)
end if
! symmetrise the density matrix
call symdmat(lmax,ld,dmat)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

