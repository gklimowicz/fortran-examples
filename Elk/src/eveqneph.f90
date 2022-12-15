
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine eveqneph
use modmain
use modphonon
use modbog
use modmpi
use modomp
implicit none
! local variables
integer iq,ik
integer n,lp,nthd
! allocatable arrays
complex(8), allocatable :: dw(:,:,:),ex(:,:,:),fy(:,:)
complex(8), allocatable :: au(:,:,:),bv(:,:,:)
!------------------------------------!
!     phonon eigenvalue equation     !
!------------------------------------!
allocate(dw(nbph,nbph,nqpt),ex(nbph,nbph,nqpt),fy(nbph,nqpt))
! parallel loop over reduced q-point set
call holdthd(nqpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do iq=1,nqpt
! distribute among MPI processes
  if (mod(iq-1,np_mpi).ne.lp_mpi) cycle
! generate the matrices D and E
  call hmlephde(iq,dw(:,:,iq),ex(:,:,iq))
! zero the vector F
  fy(:,iq)=0.d0
! solve the phononic Bogoliubov equation
  call eveqnwxy(nbph,pwxpsn,dw(:,:,iq),ex(:,:,iq),fy(:,iq),evalwx(:,iq))
end do
!$OMP END DO
!$OMP DO
do iq=1,nqpt
! distribute among MPI processes
  if (mod(iq-1,np_mpi).ne.lp_mpi) cycle
! compute the density matrices
  call dmatwx(nbph,dw(:,:,iq),ex(:,:,iq),dxx(:,:,iq),dwx(:,:,iq),xnorm(:,iq))
! write the eigenvalues, eigenvectors and X-norms to file
  if (tlast) then
    call putevalwx(iq,evalwx(:,iq))
    call putevecwxy(iq,dw(:,:,iq),ex(:,:,iq),fy(:,iq))
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
deallocate(dw,ex,fy)
! broadcast arrays to every MPI process
if (np_mpi.gt.1) then
  n=nbph*nbph
  do iq=1,nqpt
    lp=mod(iq-1,np_mpi)
    call mpi_bcast(evalwx(:,iq),nbph,mpi_double_precision,lp,mpicom,ierror)
    call mpi_bcast(xnorm(:,iq),nbph,mpi_double_precision,lp,mpicom,ierror)
    call mpi_bcast(dxx(:,:,iq),n,mpi_double_complex,lp,mpicom,ierror)
    call mpi_bcast(dwx(:,:,iq),n,mpi_double_complex,lp,mpicom,ierror)
  end do
end if
!--------------------------------------!
!     electron eigenvalue equation     !
!--------------------------------------!
allocate(au(nstsv,nstsv,nkpt),bv(nstsv,nstsv,nkpt))
! parallel loop over reduced k-point set
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! generate the matrix A
  call hmlepha(ik,au(:,:,ik))
! generate the matrix B
  call hmlephb(ik,bv(:,:,ik))
! solve the electronic Bogoliubov equation
  call eveqnuv(nstsv,au(:,:,ik),bv(:,:,ik),evaluv(:,ik))
end do
!$OMP END DO
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! compute the density matrices
  call dmatuv(nstsv,efermi,evalsv(:,ik),au(:,:,ik),bv(:,:,ik),dvv(:,:,ik), &
   duv(:,:,ik),vnorm(:,ik))
! write the eigenvalues and eigenvectors to file
  if (tlast) then
    call putevaluv(ik,evaluv(:,ik))
    call putevecuv(ik,au(:,:,ik),bv(:,:,ik))
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
deallocate(au,bv)
! broadcast arrays to every MPI process
if (np_mpi.gt.1) then
  n=nstsv*nstsv
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(evaluv(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
    call mpi_bcast(vnorm(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
    call mpi_bcast(dvv(:,:,ik),n,mpi_double_complex,lp,mpicom,ierror)
    call mpi_bcast(duv(:,:,ik),n,mpi_double_complex,lp,mpicom,ierror)
  end do
end if
end subroutine

