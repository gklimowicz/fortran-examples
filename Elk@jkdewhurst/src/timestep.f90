
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine timestep
use modmain
use modtddft
use modmpi
use modomp
implicit none
! local variables
integer ik,i,j,nthd
real(8) ca,dt,t1
complex(8) z1,z2
! automatic arrays
real(8) w(nstsv)
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:),bmt(:,:,:)
complex(8), allocatable :: evecsv(:,:),evectv(:,:),evecsvt(:,:)
complex(8), allocatable :: kmat(:,:),pmat(:,:,:)
complex(8), allocatable :: a(:,:),b(:,:),c(:,:)
if (itimes.ge.ntimes) then
  write(*,*)
  write(*,'("Error(timestep): itimes >= ntimes : ",2I8)') itimes,ntimes
  write(*,*)
  stop
end if
allocate(vmt(npcmtmax,natmtot),vir(ngtot))
if (spinpol) allocate(bmt(npcmtmax,natmtot,ndmag))
! generate the Kohn-Sham potential and magnetic field in spherical coordinates
! and multiply by the radial integration weights; also multiply the interstitial
! potential with the characteristic function
call vblocal(vmt,vir,bmt)
! backup existing time-dependent Kohn-Sham eigenvectors if required
call tdbackup
! time step length
dt=times(itimes+1)-times(itimes)
! zero the kinetic energy
engykn=0.d0
! zero the total current
jtot(:)=0.d0
! loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv,evectv,evecsvt) &
!$OMP PRIVATE(kmat,pmat,w,a,b,c) &
!$OMP PRIVATE(i,j,t1,z1,z2) &
!$OMP NUM_THREADS(nthd)
allocate(evecsv(nstsv,nstsv),evectv(nstsv,nstsv),evecsvt(nstsv,nstsv))
allocate(kmat(nstsv,nstsv),pmat(nstsv,nstsv,3))
allocate(a(nstsv,nstsv),b(nstsv,nstsv),c(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! get the kinetic matrix elements from file
  call getkmat(ik,kmat)
! get the momentum matrix elements from file
  call getpmat(vkl(:,ik),pmat)
! generate the Hamiltonian matrix in the ground-state second-variational basis
  call genhmlt(ik,vmt,vir,bmt,bsir,kmat,pmat,evectv)
! diagonalise the Hamiltonian to get third-variational eigenvectors
  if (spinpol.and.(.not.ncmag)) then
! collinear case requires block diagonalisation
    call eveqnzh(nstfv,nstsv,evectv,w)
    i=nstfv+1
    call eveqnzh(nstfv,nstsv,evectv(i,i),w(i))
    do i=1,nstfv
      do j=1,nstfv
        evectv(i,j+nstfv)=0.d0
        evectv(i+nstfv,j)=0.d0
      end do
    end do
  else
! non-collinear or spin-unpolarised: full diagonalisation
    call eveqnzh(nstsv,nstsv,evectv,w)
  end if
! read in ground-state eigenvectors
  call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! convert third-variational eigenvectors to first-variational basis
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,evectv,nstsv,zzero,a, &
   nstsv)
! time propagate instantaneous eigenvectors across one time step
  if (tdphi.eq.0.d0) then
! real time evolution
    do i=1,nstsv
      t1=-w(i)*dt
      z1=cmplx(cos(t1),sin(t1),8)
      b(:,i)=z1*a(:,i)
    end do
  else
! complex time evolution
    z2=cmplx(sin(tdphi),cos(tdphi),8)
    do i=1,nstsv
      t1=-w(i)*dt
      z1=exp(t1*z2)
      b(:,i)=z1*a(:,i)
    end do
  end if
! read in time-dependent Kohn-Sham eigenvectors (first-variational basis)
  call getevecsv(filext,ik,vkl(:,ik),evecsvt)
! apply time evolution operator
  call zgemm('C','N',nstsv,nstsv,nstsv,zone,a,nstsv,evecsvt,nstsv,zzero,c,nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,b,nstsv,c,nstsv,zzero,evecsvt,nstsv)
! orthonormalise the eigenvectors if required
  if (tdphi.ne.0.d0) call unitary(nstsv,evecsvt)
! add to the kinetic energy
  call engyknk(ik,kmat,evecsv,evecsvt)
! add to the total current
  call jtotk(ik,pmat,evecsv,evecsvt)
! write the new eigenvectors to file
  call putevecsv(filext,ik,evecsvt)
end do
!$OMP END DO
deallocate(evecsv,evectv,evecsvt)
deallocate(kmat,pmat,a,b,c)
!$OMP END PARALLEL
call freethd(nthd)
deallocate(vmt,vir)
if (spinpol) deallocate(bmt)
! add the kinetic energy and total current from each process and redistribute
if (np_mpi.gt.1) then
  call mpi_allreduce(mpi_in_place,engykn,1,mpi_double_precision,mpi_sum,mpicom,&
   ierror)
  call mpi_allreduce(mpi_in_place,jtot,3,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
end if
! add the core kinetic energy
engykn=engykn+engykncr
! coupling constant of the external A-field (-1/c)
ca=-1.d0/solsc
! add the diamagnetic current to total
do i=1,3
  jtot(i)=jtot(i)+ca*afieldt(i,itimes)*(chgtot-chgstot(i))
end do
! symmetrise the vector
call symvec(jtot)
! total current magnitude
jtotm=sqrt(jtot(1)**2+jtot(2)**2+jtot(3)**2)
! write the time step to file
if (mp_mpi) call writetimes
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

