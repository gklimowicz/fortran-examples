
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdbackup
use modmain
use modtddft
use modmpi
use modramdisk
implicit none
! local variables
integer ik
! allocatable arrays
complex(8), allocatable :: evecsvt(:,:)
if (ntsbackup.le.0) return
if (mod(itimes-1,ntsbackup).ne.0) return
! ensure eigenvectors are written to disk
wrtdsk0=wrtdsk
wrtdsk=.true.
allocate(evecsvt(nstsv,nstsv))
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! read in time-dependent Kohn-Sham eigenvectors
  call getevecsv(filext,ik,vkl(:,ik),evecsvt)
! write eigenvectors to backup file
  call putevecsv('_TD_BACKUP.OUT',ik,evecsvt)
end do
deallocate(evecsvt)
! write the time step backup file
if (mp_mpi) then
  open(50,file='TIMESTEP_BACKUP.OUT',form='FORMATTED')
  write(50,'(I8,G18.10)') itimes,times(itimes)
  close(50)
end if
wrtdsk=wrtdsk0
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

