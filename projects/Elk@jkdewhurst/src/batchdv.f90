
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine batchdv
use modmain
use modvars
use modmpi
implicit none
! no increment for first task
if (itask.le.1) return
! only increment for ground-state tasks
if (all(task.ne.[0,1,2,3])) return
! increment selected variables
if (mp_mpi) write(*,*)
if (any(dngridk(:).ne.0)) then
  ngridk(:)=ngridk(:)+dngridk(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented ngridk")')
end if
if (dlmaxapw.ne.0) then
  lmaxapw=lmaxapw+dlmaxapw
  if (mp_mpi) write(*,'("Info(batchdv): incremented lmaxapw")')
end if
if (dlmaxo.ne.0) then
  lmaxo=lmaxo+dlmaxo
  if (mp_mpi) write(*,'("Info(batchdv): incremented lmaxo")')
end if
if (any(davec(:,:).ne.0.d0)) then
  avec(:,:)=avec(:,:)+davec(:,:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented avec")')
end if
if (any(datposl(:,:,:).ne.0.d0)) then
  atposl(:,:,:)=atposl(:,:,:)+datposl(:,:,:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented atposl")')
end if
if (drgkmax.ne.0.d0) then
  rgkmax=rgkmax+drgkmax
  if (mp_mpi) write(*,'("Info(batchdv): incremented rgkmax")')
end if
if (dgmaxvr.ne.0.d0) then
  gmaxvr=gmaxvr+dgmaxvr
  if (mp_mpi) write(*,'("Info(batchdv): incremented gmaxvr")')
end if
if (dnempty0.ne.0.d0) then
  nempty0=nempty0+dnempty0
  if (mp_mpi) write(*,'("Info(batchdv): incremented nempty")')
end if
if (dnrmtscf.ne.0.d0) then
  nrmtscf=nrmtscf+dnrmtscf
  if (mp_mpi) write(*,'("Info(batchdv): incremented nrmtscf")')
end if
if (any(dudufix(:).ne.0.d0)) then
  udufix(:)=udufix(:)+dudufix(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented udufix")')
end if
if (any(dmomfix(:).ne.0.d0)) then
  momfix(:)=momfix(:)+dmomfix(:)
  if (mp_mpi) write(*,'("Info(batchdv): incremented momfix")')
end if
end subroutine

