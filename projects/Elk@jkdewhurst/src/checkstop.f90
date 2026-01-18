
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine checkstop
use modmain
use modmpi
implicit none
! check for STOP file (only MPI master process)
if (mp_mpi) then
  inquire(file='STOP',exist=tstop)
  if (tstop) then
    write(*,*)
    write(*,'("Info(checkstop): STOP file exists")')
! delete the STOP file
    open(50,file='STOP')
    close(50,status='DELETE')
  end if
end if
! broadcast tstop from master process to all other processes
call mpi_bcast(tstop,1,mpi_logical,0,mpicom,ierror)
end subroutine

