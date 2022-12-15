
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdinit
use modmain
use modtddft
use moddftu
use modmpi
use modomp
use modramdisk
implicit none
! local variables
integer ik
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
if ((task.eq.460).or.(task.eq.462)) then
! calculation starts at t=0
  tdt0=.true.
else
! calculation restarts
  tdt0=.false.
end if
! determine the static density and charge and write to file
if (tdt0) call rhostatic
! average force can be non-zero (allow for translation of atomic basis)
tfav00=tfav0
tfav0=.false.
tforce0=tforce
tjr0=tjr
! currents should be calculated with forces
if (tforce) tjr=.true.
! ensure eigenvectors are written to disk
wrtdsk0=wrtdsk
wrtdsk=.true.
! initialise global variables
call init0
call init1
! read the charge density and potentials from file
call readstate
call genvsig
call gencore
call energykncr
call readfermi
call linengy
call genapwlofr
call gensocfr
if (tdt0) then
! generate eigenvalues and eigenvectors only at t=0 (not for the restart) for
! the k-point set reduced with the symmetries which leave A(t) invariant for all
! time steps
  call genevfsv
else
! only read in the second-variational eigenvalues for restarts
  call readevalsv
! read in the static density and charge
  call readrhos
end if
! compute the occupation numbers
call occupy
! DFT+U
if (dftu.ne.0) then
  call gendmatmt
  call genvmatmt
  call vmatmtsc
  if (tmwrite) call genwkpr0
end if
! generate the kinetic matrix elements in the second-variational basis
call genkmat(.false.,.false.)
! write the momentum matrix elements in the second-variational basis
call genpmat
! write the power density to file
if (mp_mpi) call writeafpdt
! copy EVALFV.OUT, EVECFV.OUT, OCCSV.OUT and EVECSV.OUT to _TD.OUT extension
if (tdt0) then
  allocate(evalfv(nstfv,nspnfv),evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    call getevalfv('.OUT',ik,vkl(:,ik),evalfv)
    call putevalfv('_TD.OUT',ik,evalfv)
    call getevecfv('.OUT',ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call putevecfv('_TD.OUT',ik,evecfv)
    call putoccsv('_TD.OUT',ik,occsv(:,ik))
    call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! randomise eigenvectors at t=0 if required
    call rndevsv(rndevt0,evecsv)
    call putevecsv('_TD.OUT',ik,evecsv)
  end do
  deallocate(evalfv,evecfv,evecsv)
end if
! set global file extension
filext='_TD.OUT'
! output the new k-point set to file
if (mp_mpi) call writekpts
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! Ehrenfest dynamics
if ((task.eq.462).or.(task.eq.463)) then
! forces should not be calculated
  tforce=.false.
! enable small amplitude displacements
  tdatpos=.true.
! zero the displacements and velocities
  datposc(:,:,:,:)=0.d0
! generate the gradient of the nucleus and static density Coulomb potential
  call gengvnsmt
end if
if (tdt0) then
! start from t=0
  itimes0=1
else
! restart if required
  call tdrestart
end if
! deallocate the static density if not required
if (.not.tjr) deallocate(rhosmt,rhosir)
! read the forces calculated during the previous TDDFT run
if (tdatpos) call readforcet
! read the atomic displacements and velocities
if (trddatpos) call readdatposc
wrtdsk=wrtdsk0
end subroutine

