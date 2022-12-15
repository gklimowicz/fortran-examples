
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine atptstep
use modmain
use modtddft
use modmpi
implicit none
! local variables
integer is,ia,ias
real(8) dt,t1
! time step length
dt=times(itimes+1)-times(itimes)
do is=1,nspecies
  t1=dt/spmass(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! add to the atomic velocities
    datposc(:,1,ia,is)=datposc(:,1,ia,is)+t1*forcet(:,ias,itimes)
! add to the atomic displacements
    datposc(:,0,ia,is)=datposc(:,0,ia,is)+datposc(:,1,ia,is)*dt
  end do
end do
! write the atomic displacements and velocities to file
if (mp_mpi) call writedatposc
end subroutine

